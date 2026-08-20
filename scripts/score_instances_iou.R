#!/usr/bin/env Rscript
.bs_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.bs_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bs <- Find(file.exists, c(
  if (!is.null(.bs_ofile) && length(.bs_ofile) && nzchar(.bs_ofile))
    file.path(dirname(.bs_ofile), "bootstrap.R"),
  if (length(.bs_file)) file.path(dirname(sub("^--file=", "", .bs_file[1])),
                                  "bootstrap.R"),
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs, .bs_ofile, .bs_file)

# #V1 point-set IoU / Coverage / Panoptic-Quality scorer (GitHub issue #63).
#
# The whole benchmark is graded by greedy_match apex-distance 1:1 matching, which
# model-benchmark-results.md's Caveats already flag as blind to point-level IoU
# and crown-shape quality -- the exact failure modes the meta-pipeline's fusion is
# meant to fix. FGI-EMIT and FOR-instanceV2 grade with IoU>=0.5 + Coverage +
# Panoptic Quality; this arm adds that mask-aware suite ALONGSIDE the distance
# matching (it never replaces it), so the repo's deep arms become comparable to
# those external leaderboards and fusion gets a mask-aware objective (#P1, #P3).
#
# It reuses the SAME persisted GPU instance artefacts the #34 crown arm consumes
# (it NEVER runs a container) but keeps the FULL per-point labelling rather than
# reducing to apexes:
#   SegmentAnyTree  -> work/neon/<SITE>/segmentanytree_instances/<plot>_<rung>.laz
#                      via read_instance_points_laz(PredInstance) (0 = non-tree
#                      -> NA, dropped). algo = "segmentanytree".
#   ForestFormer3D  -> work/neon/<SITE>/forestformer3d_instances/<plot>_<rung>.laz
#                      (UserData = block, PointSourceID = per-cylinder instance);
#                      stacked into (block, inst, X, Y, Z) and run through
#                      dedup_blocks() (cross-block apex-cluster merge -- never
#                      launders within-cylinder over-segmentation) to a global_id
#                      per point. read_instance_points_laz alone is insufficient
#                      because the FF3D block id is required. algo = "forestformer3d".
# (TreeisoNet joins once detect_treeisonet_crowns.R persists its per-point treeOff
# labels; today it writes only treeisonet_crown_metrics.csv, so it is skipped.)
#
# Reference (ground-truth) instances are an explicit VORONOI-ON-STEMS PROXY: each
# point of the frozen NORMALIZED clip is assigned to the nearest field stem within
# that stem's crown radius (maxCrownDiameter / 2; FALLBACK_RADIUS where a stem has
# no measured width). This is the best mask obtainable from stem points + a field
# crown width -- it is a proxy for true field crown masks, not a hand-delineated
# truth, and is reported as such.
#
# Shared substrate: standard instance-IoU is computed on one fixed evaluation
# cloud. The persisted model clouds are NOT the frozen normalized cloud
# point-for-point (different Z datum, container voxel downsampling, FF3D cylinder
# overlap), so each model's per-point labels are PROJECTED onto the normalized
# cloud by nearest-neighbour in XY (transfer_labels, <= XFER_TOL m). The substrate
# is restricted to canopy points (Z >= CANOPY_MIN AGL, matching the detector
# min_height) inside the plot core (plot_half by plot type, mirroring score_plot's
# core precision denominator) so buffer trees never count as false commissions.
#
# Scoring (all pure helpers in model_bench_lib.R, unit-tested with no NEON data):
#   point_set_iou -> iou_match(gate 0.5) -> panoptic_quality, plus coverage_table.
# Per (plot x rung x model) score_instance_cell emits one long-form row of
# SUMMABLE accumulators stratified by crown_class; pool_pq pools by SUMMING them
# (sum matched-IoU, sum TP/FP/FN), mirroring pool()'s sum(TP)/sum(n_ref) rule.
#
# Usage:
#   Rscript scripts/score_instances_iou.R SITE=SOAP
#   Rscript scripts/score_instances_iou.R SITES=SOAP,SJER,TEAK CORES=4
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv},
#   the cached field crown widths work/neon/<SITE>/vst/<site>_vst_allyears.rds, the
#   cached frozen normalized clips work/neon/<SITE>/frozen/<SITE>/<plot>/<rung>/
#   clip_normalized.laz, and the persisted instance clouds above.
# Writes: work/neon/<SITE>/instance_iou_pq.csv (one long-form row per
#   site x plot x rung x model; SUMMABLE accumulators + per-cell rates).
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))           # plot_half, greedy_match (via lib)
source(.find("model_bench_lib.R"))     # the #V1 scorer + pooler
source(.find("io_bridge.R"))           # read_instance_points_laz
source(.find("coverage_lib.R"))        # read_arm_cache (apex-Voronoi proxy, #V6)

## ---- args (KEY=VALUE positional, per repo convention) --------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES     <- as.integer(if (is.null(A$CORES)) 4 else A$CORES)
RUNGS     <- if (is.null(A$RUNGS)) "native" else strsplit(A$RUNGS, ",")[[1]]
IOU_GATE  <- as.numeric(if (is.null(A$IOU_GATE)) 0.5 else A$IOU_GATE)
XFER_TOL  <- as.numeric(if (is.null(A$XFER_TOL)) 0.5 else A$XFER_TOL)
MERGE_TOL <- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
CANOPY_MIN <- as.numeric(if (is.null(A$CANOPY_MIN)) 2.0 else A$CANOPY_MIN)
FALLBACK_RADIUS <- as.numeric(if (is.null(A$FALLBACK_RADIUS)) 2.0 else A$FALLBACK_RADIUS)
MINTREES  <- as.integer(if (is.null(A$MINTREES)) 1 else A$MINTREES)
SAT_ID_FIELD <- "PredInstance"
# #V6: apex-Voronoi proxy rows for the cached best-config arms (mask_source =
# "voronoi_apex"), alongside the native per-point masks. APEX_R mirrors the
# fusion scorer's FUSE_R apex-Voronoi radius.
APEX_PROXY <- is.null(A$APEX_PROXY) || A$APEX_PROXY != "0"
APEX_R     <- as.numeric(if (is.null(A$APEX_R)) 4.0 else A$APEX_R)
APEX_CACHE_ARMS <- c("chm_vwf", "lmfauto", "multichm", "ptrees", "ams3d",
                     "li2012", "lidr_lmf_pc", "lidr_li2012", "lasr_lmax_pc",
                     "treeisonet", "segmentanytree", "forestformer3d")

## ---- per-model loaders -> full labelled point cloud data.frame(X, Y, id) --
# id = the model's per-point instance label (unassigned points dropped). NULL =
# artefact absent / schema failure (skip the cell); 0-row = ran-but-empty.
EMPTY_PTS <- data.frame(X = numeric(), Y = numeric(), id = integer())

# Generic loader for any persisted single-id-field instance LAZ (#V6): the
# classical arms, treeiso, and SAT all reduce to read_instance_points_laz +
# drop-unassigned. FF3D stays special (needs the block id for dedup_blocks).
make_laz_loader <- function(id_field) {
  force(id_field)
  function(path) {
    pts <- tryCatch(suppressWarnings(
      read_instance_points_laz(path, id_field = id_field)),
      error = function(e) NULL)
    if (is.null(pts)) return(NULL)
    pts <- pts[!is.na(pts$crown_id), , drop = FALSE]
    if (!nrow(pts)) return(EMPTY_PTS)
    data.frame(X = pts$X, Y = pts$Y, id = as.integer(pts$crown_id))
  }
}
load_sat_points <- make_laz_loader(SAT_ID_FIELD)

load_ff3d_points <- function(path) {
  las <- tryCatch(suppressWarnings(lidR::readLAS(path)), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  dat <- las@data
  if (!all(c("UserData", "PointSourceID") %in% names(dat))) return(NULL)
  if (lidR::is.empty(las)) return(EMPTY_PTS)
  pts <- data.frame(block = as.integer(dat$UserData),
                    inst  = as.integer(dat$PointSourceID),
                    X = dat$X, Y = dat$Y, Z = dat$Z)
  rel <- dedup_blocks(pts, merge_tol = MERGE_TOL)        # -> global_id per point
  if (!nrow(rel)) return(EMPTY_PTS)
  data.frame(X = rel$X, Y = rel$Y, id = as.integer(rel$global_id))
}

MODELS <- list(
  segmentanytree = list(dir = "segmentanytree_instances", load = load_sat_points),
  forestformer3d = list(dir = "forestformer3d_instances", load = load_ff3d_points),
  ptrees         = list(dir = "ptrees_instances",  load = make_laz_loader("treeID")),
  ams3d          = list(dir = "ams3d_instances",   load = make_laz_loader("crown_id")),
  li2012         = list(dir = "li2012_instances",  load = make_laz_loader("treeID")),
  treeiso        = list(dir = "treeiso_instances", load = make_laz_loader("treeiso")))

## ---- field crown diameter per site (cached vst rds; as crown_metrics_3d) ---
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst",
                   paste0(tolower(site), "_vst_allyears.rds"))
  if (!file.exists(rds)) return(data.frame(individualID = character(),
                                           maxCrownDiameter = numeric()))
  dat <- readRDS(rds)
  ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID), c("individualID", "maxCrownDiameter")]
}

instance_path <- function(idir, pid, rung) {
  cand <- file.path(idir, paste0(pid, "_", rung, c(".laz", ".las")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# Frozen normalized clip for the substrate. frozen_clip() callers pass
# out_root = work/neon/<SITE>/frozen and it writes out_root/<SITE>/<plot>/<rung>/.
frozen_norm_path <- function(nd, site, pid, rung)
  file.path(nd, "frozen", site, pid, rung, "clip_normalized.laz")

## ---- per-plot, per-rung IoU/PQ scoring -----------------------------------
# Build the reference partition once per (plot, rung) from the frozen normalized
# canopy substrate inside the plot core, then for every model with a persisted
# cloud project its labels onto that substrate and emit one accumulator row.
run_plot <- function(site, pid, pc, gt, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < MINTREES) return(NULL)
  # per-stem crown radius: measured maxCrownDiameter/2, else the documented fallback
  radius <- ifelse(is.finite(stems$maxCrownDiameter) & stems$maxCrownDiameter > 0,
                   stems$maxCrownDiameter / 2, FALLBACK_RADIUS)
  n_fallback <- sum(!(is.finite(stems$maxCrownDiameter) & stems$maxCrownDiameter > 0))
  # ref id = stem row index; ref_class keyed by that index (as character)
  ref_class <- setNames(as.character(stems$crown_class), as.character(seq_len(nrow(stems))))

  rows <- list()
  for (rung in RUNGS) {
    fp <- frozen_norm_path(nd, site, pid, rung)
    if (!file.exists(fp)) next                          # no substrate -> skip cell
    sub <- tryCatch(suppressWarnings(lidR::readLAS(fp)), error = function(e) NULL)
    if (is.null(sub) || lidR::is.empty(sub)) next
    sd <- sub@data
    keep <- sd$Z >= CANOPY_MIN &
            abs(sd$X - cx) <= ph & abs(sd$Y - cy) <= ph  # canopy points in core
    sx <- sd$X[keep]; sy <- sd$Y[keep]
    if (length(sx) < 1) next
    ref <- assign_points_to_stems(sx, sy, stems$E, stems$N, radius)
    if (!any(!is.na(ref))) next                          # no stem captured a point

    for (mname in names(MODELS)) {
      model <- MODELS[[mname]]
      ipath <- instance_path(file.path(nd, model$dir), pid, rung)
      if (is.na(ipath)) next                             # model not run on this cell
      pp <- tryCatch(model$load(ipath), error = function(e) NULL)
      # Skip ONLY when the loader returned NULL (artifact absent / schema fail).
      # A 0-row pp means the model ran but predicted nothing -> let it fall
      # through so transfer_labels yields all-NA labels and score_instance_cell
      # records it as TP=0 / FN=n_ref rather than silently dropping the cell.
      if (is.null(pp)) next
      pred <- transfer_labels(pp$X, pp$Y, pp$id, sx, sy, tol = XFER_TOL)
      row <- score_instance_cell(pred, ref, ref_class = ref_class, gate = IOU_GATE)
      row$site <- site; row$plot <- pid; row$rung <- rung; row$model <- mname
      row$mask_source <- "native"
      row$n_stems <- nrow(stems); row$n_fallback_radius <- n_fallback
      rows[[length(rows) + 1]] <- row
    }

    # #V6 apex-Voronoi proxy rows: the cached best-config apex sets become
    # masks by nearest-apex assignment within APEX_R on the same substrate --
    # the fuse_detectors.R scoring proxy, symmetric with the Voronoi-on-stems
    # reference. This puts the apex-only detectors (CHM-VWF, multichm, lmfauto,
    # TreeisoNet, the pc twins) on the board, and double-scores the native-mask
    # arms so proxy inflation is measurable per arm. Cells exist only at each
    # arm's cached best rung (read_arm_cache -> NULL elsewhere).
    if (APEX_PROXY) {
      cache_dir <- file.path(nd, "best_treetop_cache")
      for (arm in APEX_CACHE_ARMS) {
        det <- suppressWarnings(read_arm_cache(cache_dir, arm, site, pid, rung))
        if (is.null(det)) next
        pred <- if (nrow(det))
          assign_points_to_stems(sx, sy, det$x, det$y, rep(APEX_R, nrow(det))) else
          rep(NA_integer_, length(sx))
        row <- score_instance_cell(pred, ref, ref_class = ref_class, gate = IOU_GATE)
        row$site <- site; row$plot <- pid; row$rung <- rung; row$model <- arm
        row$mask_source <- "voronoi_apex"
        row$n_stems <- nrow(stems); row$n_fallback_radius <- n_fallback
        rows[[length(rows) + 1]] <- row
      }
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- per-site driver ------------------------------------------------------
run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv")
  pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) {
    cat(sprintf("[%s] no ground truth / centroids -- skipped\n", site)); return(NULL) }
  gt <- read.csv(gtf, stringsAsFactors = FALSE)
  pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  # COALESCE the crown diameter: the CSV already carries maxCrownDiameter; the
  # cached VST RDS may carry a (more authoritative) value. Keep the CSV value as
  # a fallback and overwrite it only where the RDS provides a finite value, so a
  # missing/incomplete RDS no longer silently collapses every reference radius to
  # FALLBACK_RADIUS. RDS present (the committed run) -> numbers unchanged.
  mcd_csv <- if (is.null(gt$maxCrownDiameter)) rep(NA_real_, nrow(gt)) else
    as.numeric(gt$maxCrownDiameter)
  gt <- gt[, setdiff(names(gt), "maxCrownDiameter"), drop = FALSE]
  gt$maxCrownDiameter <- mcd_csv
  fc <- field_crowns(site)
  names(fc)[names(fc) == "maxCrownDiameter"] <- "maxCrownDiameter_rds"
  gt <- merge(gt, fc, by = "individualID", all.x = TRUE)
  rds_ok <- is.finite(gt$maxCrownDiameter_rds)
  gt$maxCrownDiameter[rds_ok] <- gt$maxCrownDiameter_rds[rds_ok]
  gt$maxCrownDiameter_rds <- NULL

  plots <- intersect(unique(gt$plotID), pc$plotID)
  cat(sprintf("[%s] scoring %d plots over rungs {%s}\n",
              site, length(plots), paste(RUNGS, collapse = ",")))
  if (!length(plots)) return(NULL)

  res_list <- mclapply(plots, function(p)
    tryCatch(run_plot(site, p, pc, gt, nd),
             error = function(e) { message("  ", p, " failed: ",
                                            conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res) || !nrow(res)) {
    cat(sprintf("[%s] no instance cells scored\n", site)); return(NULL) }
  # lead with the identifier columns, then the accumulators
  idc <- c("site", "plot", "rung", "model", "mask_source",
           "n_stems", "n_fallback_radius")
  res <- res[, c(idc, setdiff(names(res), idc)), drop = FALSE]
  o <- file.path(nd, "instance_iou_pq.csv")
  write.csv(res, o, row.names = FALSE)
  cat(sprintf("[%s] wrote %d cell rows -> %s\n", site, nrow(res), o))
  res
}

## ---- pooled report --------------------------------------------------------
print_tables <- function(res) {
  cat("\n===== POINT-SET IoU / COVERAGE / PANOPTIC QUALITY (pooled by SUM) =====\n")
  cat(sprintf("(IoU gate %.2f; reference = Voronoi-on-stems proxy)\n", IOU_GATE))
  if (is.null(res$mask_source)) res$mask_source <- "native"
  # Per-rung sections: arms cover different rung subsets (li2012 native-only,
  # the proxy rows one cached rung each), so one table pooled across rungs
  # would mix denominators unequally across arms.
  for (rung in intersect(c("native", "8", "4", "2", "1"), unique(res$rung))) {
    rr <- res[res$rung == rung, , drop = FALSE]
    key <- paste(rr$model, rr$mask_source, sep = "  ")
    cat(sprintf("\n----- rung %s -----\n", rung))
    cat(sprintf("%-30s %5s %6s %6s %6s %6s %6s %6s %6s\n",
                "model  mask", "cells", "nref", "P", "R", "F1", "Cov", "SQ", "PQ"))
    for (k in sort(unique(key))) {
      p <- pool_pq(rr[key == k, , drop = FALSE])
      cat(sprintf("%-30s %5d %6d %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
                  k, p$n_cells, p$n_ref, p$precision, p$recall, p$F1,
                  p$coverage, p$SQ, p$PQ))
    }
  }
  cat("\n--- per crown class: recall@IoU0.5 / Coverage (pooled, native rung) ---\n")
  nat <- res[res$rung == "native", , drop = FALSE]
  key <- paste(nat$model, nat$mask_source, sep = "  ")
  cat(sprintf("%-30s %-12s %6s %6s %6s\n",
              "model  mask", "class", "nref", "recall", "Cov"))
  for (k in sort(unique(key))) {
    p <- pool_pq(nat[key == k, , drop = FALSE])
    for (cl in c(POOL_CLASSES, "understory")) {
      n <- p[[paste0("n_", cl)]]; rec <- p[[paste0("rec_", cl)]]
      cov <- p[[paste0("cov_", cl)]]
      cat(sprintf("%-30s %-12s %6d %6.3f %6.3f\n", k, cl,
                  if (is.null(n)) 0L else n,
                  if (is.null(rec) || is.na(rec)) NA_real_ else rec,
                  if (is.null(cov) || is.na(cov)) NA_real_ else cov))
    }
  }
}

run_main <- function() {
  t0 <- Sys.time()
  all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- do.call(rbind, all_res)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d cell rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) {
    cat("\nrows per model:\n"); print(table(res$model))
    print_tables(res)
  } else {
    cat("\nNo instance cells scored. This needs the persisted GPU instance\n")
    cat("clouds + frozen normalized clips under work/neon/<SITE>/.\n")
  }
}

if (sys.nframe() == 0L) run_main()
