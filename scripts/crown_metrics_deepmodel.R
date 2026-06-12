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

# Deep-model crown-diameter arm for GitHub issue #34 (#M6 / #M8).
#
# Issue #30 crown-scored the three CHM-free 3-D point-instance segmenters
# (Li 2012, ptrees, AMS3D) against NEON field crown diameter; issue #20 added the
# SOAP-only TreeisoNet treeOff arm. The two deep instance-segmentation DETECTION
# arms -- SegmentAnyTree (#M6) and ForestFormer3D (#M8) -- already produce
# per-point instance labels on the SAME frozen clips but were never crown-scored.
# This script closes that gap by reusing the #7/#30 crown-scoring harness: from
# each model's per-point instance labelling it derives a convex-hull crown
# diameter (d_eq / d_caliper) via crown_diameter_table(), reduces each instance to
# its apex, matches apexes to field stems with the SAME greedy_match harness (4 m
# tol + height gate), and writes canonical crown-metrics rows. analyze_crown_
# metrics.R then unions these per-model CSVs with the #7 CHM arms, the #30 3-D
# arms, and the TreeisoNet negative result and pools per algorithm.
#
# This arm consumes PERSISTED GPU instance artefacts (it NEVER runs a container):
#   SegmentAnyTree  -> work/neon/<SITE>/segmentanytree_instances/<plot>_<rung>.laz
#                      (merged per-point LAS with the PredInstance extra dim;
#                      see detect_segmentanytree_sweep.R for the merged-LAS
#                      layout and ID_FIELD). Read full table via the bridge's
#                      read_instance_points_laz(PredInstance) -> crown_diameter_
#                      table + instance_apex -> match. algo = "segmentanytree".
#   ForestFormer3D  -> work/neon/<SITE>/forestformer3d_instances/<plot>_<rung>.laz
#                      (merged per-cylinder labelled LAZ: UserData = block,
#                      PointSourceID = per-cylinder instance id; see
#                      detect_forestformer3d_sweep.R). Stacked into
#                      (block, inst, X, Y, Z) -> ff3d_crown_table() runs
#                      dedup_blocks() (cross-block apex-cluster merge -- never
#                      launders within-cylinder over-segmentation) then
#                      crown_diameter_table + instance_apex keyed by global_id ->
#                      match. algo = "forestformer3d".
#
# Both LAS variants keep ABSOLUTE UTM Z, but the crown diameter (a horizontal X/Y
# convex hull) is invariant to the Z datum; only the apex z feeds the matching
# height gate, where the field stem heights are AGL. The instance LAS is therefore
# converted to AGL at its apex via the cached frozen clip's ground_dtm.tif (the
# same det_to_agl the detection arms use) before matching, so the height gate is
# apples-to-apples. A plot whose DTM or instance LAS is absent is SKIPPED with a
# message -- never fabricated.
#
# Start SOAP-native (the site/rung with existing GPU artefacts); the script is
# TOLERANT -- a requested site/model with no instance artefacts is skipped with a
# clear message and contributes no rows. SJER + TEAK extend automatically once
# their GPU clips exist.
#
# Usage:
#   Rscript scripts/crown_metrics_deepmodel.R SITE=SOAP
#   Rscript scripts/crown_metrics_deepmodel.R SITES=SOAP,SJER,TEAK CORES=4 TOL=4
#
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv},
#   the cached field-crown widths work/neon/<SITE>/vst/<site>_vst_allyears.rds,
#   the cached frozen clips' ground DTM
#   work/neon/<SITE>/frozen/<SITE>/<plot>/<rung>/ground_dtm.tif (reused, never
#   regenerated), and the persisted instance clouds above.
# Writes (NEW per-model files, never overwrites the #7/#30 CSVs):
#   work/neon/<SITE>/segmentanytree_crown_metrics.csv
#   work/neon/<SITE>/forestformer3d_crown_metrics.csv  (one row per matched tree,
#   canonical cols: site, plot, algo, crown_class, individualID, d_eq, d_caliper,
#   area, field_maxCD, field_ninetyCD; algo in {segmentanytree,forestformer3d}).
suppressMessages({
  library(lidR); library(data.table); library(parallel)
})
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))
source(.find("io_bridge.R"))

## ---- args (KEY=VALUE positional, per repo convention) --------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) {
  strsplit(A$SITES, ",")[[1]]
} else if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES <- as.integer(if (is.null(A$CORES)) 4 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL))  4 else A$TOL)
MERGE_TOL <- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
MINTREES  <- 6
# Crown diameter is scored at ONE rung per plot to keep the pooled table one row
# per matched tree (matching #30, which crown-scores native only) -- mixing rungs
# would double-count a stem. Default native; RUNGS=native,8 widens it (each rung
# then contributes its own rows, so only widen for a density-sensitivity study).
RUNGS <- if (is.null(A$RUNGS)) "native" else strsplit(A$RUNGS, ",")[[1]]
# The merged SAT LAS carries the instance label as the PredInstance extra dim
# (merge_pt_ss_is.py: preds_instance_segmentation -> PredInstance, +1, NA->0),
# matching detect_segmentanytree_sweep.R::ID_FIELD.
SAT_ID_FIELD <- "PredInstance"

## ---- per-model instance loaders -> (diam, apex) on the SAME labelling -----
# Each loader takes a persisted instance LAS/LAZ for one (plot, rung) and returns
# list(diam, apex) keyed by the model's instance id, plus the absolute-UTM apex z
# the caller converts to AGL. NULL = artefact absent / unreadable (skip the
# cell); a present-but-empty cloud yields 0-row diam/apex (a legit no-crown cell).

# SegmentAnyTree: full PredInstance-labelled point table -> per-instance diameter
# + apex (the bridge maps the 0 = non-tree label to NA, dropping it).
load_sat <- function(path) {
  pts <- tryCatch(read_instance_points_laz(path, id_field = SAT_ID_FIELD),
                  error = function(e) NULL)
  if (is.null(pts)) return(NULL)                     # missing/schema -> skip
  list(diam = crown_diameter_table(pts, id_col = "crown_id", min_pts = 5),
       apex = instance_apex(pts, id_col = "crown_id"))
}

# ForestFormer3D: merged per-cylinder labelled LAZ (UserData = block,
# PointSourceID = per-cylinder instance id). Stack into (block, inst, X, Y, Z),
# then ff3d_crown_table() runs dedup_blocks (cross-block apex-cluster merge) +
# crown_diameter_table + instance_apex keyed by global_id -- the #M8 design.
load_ff3d <- function(path) {
  las <- tryCatch(lidR::readLAS(path), error = function(e) NULL)
  if (is.null(las)) return(NULL)                     # missing/unreadable -> skip
  dat <- las@data
  if (!all(c("UserData", "PointSourceID") %in% names(dat))) return(NULL)  # schema
  if (lidR::is.empty(las))
    return(list(diam = crown_diameter_table(
                  data.frame(global_id = integer(), X = numeric(),
                             Y = numeric(), Z = numeric()), id_col = "global_id"),
                apex = instance_apex(
                  data.frame(global_id = integer(), X = numeric(),
                             Y = numeric(), Z = numeric()), id_col = "global_id")))
  pts <- data.frame(block = as.integer(dat$UserData),
                    inst  = as.integer(dat$PointSourceID),
                    X = dat$X, Y = dat$Y, Z = dat$Z)
  ff3d_crown_table(pts, merge_tol = MERGE_TOL, min_pts = 5)
}

# Per-model instance-artefact directory + the loader that consumes one cloud.
MODELS <- list(
  segmentanytree = list(dir = "segmentanytree_instances", load = load_sat),
  forestformer3d = list(dir = "forestformer3d_instances", load = load_ff3d))

## ---- field crown diameter per site (join from cached vst rds) ------------
# Identical to crown_metrics_3d.R::field_crowns: dedup apparentindividual to the
# nearest-to-2021 measurement per individualID; keep max/ninetyCrownDiameter.
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst",
                   paste0(tolower(site), "_vst_allyears.rds"))
  dat <- readRDS(rds)
  ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID),
     c("individualID", "maxCrownDiameter", "ninetyCrownDiameter")]
}

## ---- locate the persisted instance cloud for one (plot, rung) ------------
# Discovers <plot>_<rung>.laz|.las under work/neon/<SITE>/<model dir>/ (the
# layout the GPU detection arms persist their merged instance clouds to). Returns
# the path, or NA if no artefact exists for this cell (-> the cell is skipped).
instance_path <- function(idir, pid, rung) {
  cand <- file.path(idir, paste0(pid, "_", rung, c(".laz", ".las")))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}

# Match model_bench_lib.R::frozen_clip(), whose callers pass
# out_root = work/neon/<SITE>/frozen and which writes out_root/<SITE>/<plot>/<rung>.
frozen_dtm_path <- function(nd, site, pid, rung) {
  file.path(nd, "frozen", site, pid, rung, "ground_dtm.tif")
}

## ---- per-plot deep-model crown benchmark ---------------------------------
# For one model and one plot: for every rung with a persisted instance cloud,
# load (diam, apex) from the SAME labelling, convert the apex z to AGL via the
# cached frozen DTM (absolute UTM Z -> height above ground, matching the field
# stem height gate), match to the plot-core stems with the shared
# score_crowns_against_field glue, and return canonical-cols rows. The rung is
# only used to locate the cloud + its DTM; one matched tree per stem per cell.
run_plot_model <- function(site, pid, mname, model, pc, gt, fc, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
  if (nrow(stems) < 1) return(NULL)
  idir <- file.path(nd, model$dir)
  if (!dir.exists(idir)) return(NULL)               # model not run here -> skip

  rows <- list()
  for (rung in RUNGS) {
    ipath <- instance_path(idir, pid, rung)
    if (is.na(ipath)) next                          # no cloud for this cell
    sg <- tryCatch(model$load(ipath), error = function(e) NULL)
    if (is.null(sg)) next                           # unreadable/schema -> skip
    if (!nrow(sg$apex)) next                        # ran-but-empty -> no crowns
    # The cached frozen DTM lives alongside the clip the cloud was inferred on.
    dtm <- frozen_dtm_path(nd, site, pid, rung)
    if (!file.exists(dtm)) next                     # no DTM -> cannot AGL-gate
    apex_agl <- tryCatch(det_to_agl(sg$apex, dtm), error = function(e) NULL)
    if (is.null(apex_agl) || !nrow(apex_agl)) next  # all apexes off-DTM -> skip
    r <- score_crowns_against_field(sg$diam, apex_agl, stems, fc, tol = TOL,
                                    site = site, plot = pid, algo = mname)
    if (nrow(r)) rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- per-site, per-model driver -------------------------------------------
# Mirrors crown_metrics_3d.R::run_site: live & is_tree stems, authoritative rds
# crown-diameter join (drop any pre-existing CD cols first), >=MINTREES plots,
# mclapply over plots, write the NEW per-model crown-metrics CSV. Returns a named
# list of per-model result frames (only models that produced rows).
run_site <- function(site) {
  nd  <- file.path(d, "neon", site)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  gt  <- gt[, setdiff(names(gt),
                      c("maxCrownDiameter", "ninetyCrownDiameter")), drop = FALSE]
  fc  <- field_crowns(site)
  gt  <- merge(gt, fc, by = "individualID", all.x = TRUE)
  gt  <- gt[!is.na(gt$maxCrownDiameter) | !is.na(gt$ninetyCrownDiameter), ]

  counts <- table(gt$plotID)
  keep <- names(counts)[counts >= MINTREES]
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] plots with >=%d stems w/ field CD: %d (%s)\n",
              site, MINTREES, length(keep), paste(keep, collapse = ",")))
  if (!length(keep)) return(NULL)

  out <- list()
  for (mname in names(MODELS)) {
    model <- MODELS[[mname]]
    idir  <- file.path(nd, model$dir)
    if (!dir.exists(idir)) {
      cat(sprintf("[%s] %s: no instance dir (%s) -- skipped\n",
                  site, mname, model$dir))
      next
    }
    res_list <- mclapply(keep, function(p)
      tryCatch(run_plot_model(site, p, mname, model, pc, gt, fc, nd),
               error = function(e) { message("  ", mname, "/", p, " failed: ",
                                              conditionMessage(e)); NULL }),
      mc.cores = CORES, mc.preschedule = FALSE)
    res <- do.call(rbind, Filter(Negate(is.null), res_list))
    if (is.null(res) || !nrow(res)) {
      cat(sprintf("[%s] %s: no crowns matched\n", site, mname)); next }
    o <- file.path(nd, sprintf("%s_crown_metrics.csv", mname))
    write.csv(res, o, row.names = FALSE)
    cat(sprintf("[%s] %s: wrote %d matched-tree rows -> %s\n",
                site, mname, nrow(res), o))
    out[[mname]] <- res
  }
  if (!length(out)) return(NULL)
  out
}

## ---- scoring helpers (mirror crown_metrics_3d.R) -------------------------
# Pooled error stats over matched trees: RMSE/MAE/bias/R2 from summed squared
# errors / n, never a mean of per-plot rates. det = detected, fld = field.
err_stats <- function(det, fld) {
  ok <- is.finite(det) & is.finite(fld)
  det <- det[ok]; fld <- fld[ok]; n <- length(det)
  if (n < 2) return(data.frame(n = n, rmse = NA, mae = NA, bias = NA, r2 = NA))
  e <- det - fld
  ss_res <- sum((fld - det)^2)
  ss_tot <- sum((fld - mean(fld))^2)
  data.frame(n = n,
             rmse = sqrt(mean(e^2)), mae = mean(abs(e)),
             bias = mean(e),
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_)
}

print_tables <- function(res) {
  algos <- unique(res$algo)
  cat("\n===== DEEP-MODEL CROWN-DIAMETER ACCURACY (pooled, all sites) =====\n")
  for (defn in list(c("d_eq", "field_ninetyCD",
                      "equiv-circle d_eq vs ninetyCrownDiameter"),
                    c("d_caliper", "field_maxCD",
                      "max-caliper d_caliper vs maxCrownDiameter"))) {
    cat(sprintf("\n--- %s ---\n", defn[3]))
    cat(sprintf("%-16s %5s %7s %7s %7s %7s\n",
                "algo", "n", "rmse", "mae", "bias", "r2"))
    for (a in algos) {
      s <- err_stats(res[[defn[1]]][res$algo == a], res[[defn[2]]][res$algo == a])
      cat(sprintf("%-16s %5d %7.2f %7.2f %7.2f %7.3f\n",
                  a, s$n, s$rmse, s$mae, s$bias, s$r2))
    }
  }
}

## ---- run ------------------------------------------------------------------
run_main <- function() {
  t0 <- Sys.time()
  all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- do.call(rbind, r)
  }
  res <- do.call(rbind, all_res)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d matched-tree rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) {
    cat("\nrows per algorithm:\n"); print(table(res$algo))
    cat("\nrows per crown class:\n")
    print(table(res$crown_class, useNA = "ifany"))
    print_tables(res)
  } else {
    cat("\nNo deep-model crown rows produced. This needs the persisted GPU\n")
    cat("instance clouds under work/neon/<SITE>/{segmentanytree,",
        "forestformer3d}_instances/.\n", sep = "")
  }
}

if (sys.nframe() == 0L) run_main()
