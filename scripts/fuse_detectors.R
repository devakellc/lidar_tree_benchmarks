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

# #P1 cross-arm fusion arm (GitHub issue #64).
#
# Every detector arm is scored in isolation and pooled side-by-side, but nothing
# combines them -- and the per-class table shows clear complementarity (CHM-VWF /
# multichm carry overstory, SegmentAnyTree carries understory). This is the first
# actual META-PIPELINE detector: it materializes per-cell apexes for each arm on
# the IDENTICAL frozen cells, clusters them across arms (fuse_apexes), and emits
# three operating points -- union (k>=1, recall-max), majority (k>=ceil(N/2),
# precision-max), and a height-layered mode (CHM apices for the overstory, point/
# deep apices for the understory) -- plus the full k=1..N Pareto frontier. Each
# fused set is scored against field stems with BOTH score_plot (apex distance) and
# the #V1 IoU/PQ scorer (a Voronoi-on-apexes mask proxy on the frozen substrate,
# symmetric with the Voronoi-on-stems reference), pooled by SUM and restricted to
# the common cell set, and compared to the best single arm.
#
# Apex materialization (the *_results.csv are scored summaries, not apex
# inventories, so detectors are re-run / read fresh here):
#   chm_vwf       -- detect_lasr on the frozen normalized clip (density res, a=0.10)
#   multichm      -- lidRplugins::multichm on the same clip
#   li2012        -- lidR li2012 on the same clip (NATIVE only; meaningless sparse)
#   segmentanytree-- persisted segmentanytree_instances/<plot>_<rung>.laz
#                    (PredInstance) -> reduce_instances -> det_to_agl(frozen DTM)
#   forestformer3d-- persisted forestformer3d_instances/<plot>_<rung>.laz
#                    (UserData/PointSourceID) -> ff3d_collapse -> agl_guard (native+8)
# TreeisoNet is deferred: it persists no reusable per-point/apex labels (apex-only
# GPU results.csv), so it cannot be re-fused offline; noted in the results doc.
# All apex z are AGL (normalized clip, or DTM-converted), so the fusion height
# gate is apples-to-apples.
#
# Usage:
#   Rscript scripts/fuse_detectors.R SITE=SOAP CORES=1
#   Rscript scripts/fuse_detectors.R SITES=SOAP,SJER,TEAK RUNGS=native CORES=1
# (CORES=1: detect_lasr uses lasR exec, which can drop dense cells under fork.)
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv},
#   the cached field crown widths vst rds, the frozen normalized clips + DTMs, and
#   the persisted deep instance clouds.
# Writes: work/neon/<SITE>/fusion_results.csv (one row per
#   site x plot x rung x config; pool distance cols with pool(), IoU cols by SUM).
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))
source(.find("io_bridge.R"))

## ---- args -----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES     <- as.integer(if (is.null(A$CORES)) 1 else A$CORES)
RUNGS     <- if (is.null(A$RUNGS)) "native" else strsplit(A$RUNGS, ",")[[1]]
MERGE_TOL <- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
Z_TOL     <- as.numeric(if (is.null(A$Z_TOL)) 5.0 else A$Z_TOL)
TOL       <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
VWF_A     <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
OVERSTORY_FRAC  <- as.numeric(if (is.null(A$OVERSTORY_FRAC)) 0.5 else A$OVERSTORY_FRAC)
FUSE_R    <- as.numeric(if (is.null(A$FUSE_R)) 4.0 else A$FUSE_R)      # apex Voronoi radius
FALLBACK_R <- as.numeric(if (is.null(A$FALLBACK_R)) 2.0 else A$FALLBACK_R)
CANOPY_MIN <- as.numeric(if (is.null(A$CANOPY_MIN)) 2.0 else A$CANOPY_MIN)
MINTREES  <- as.integer(if (is.null(A$MINTREES)) 1 else A$MINTREES)
CHM_ARMS  <- c("chm_vwf", "multichm")           # overstory family for the layered mode
SAT_ID_FIELD <- "PredInstance"

field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst", paste0(tolower(site), "_vst_allyears.rds"))
  if (!file.exists(rds)) return(data.frame(individualID = character(),
                                           maxCrownDiameter = numeric()))
  dat <- readRDS(rds); ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4)); ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021); ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID), c("individualID", "maxCrownDiameter")]
}
fz <- function(nd, site, pid, rung, f) file.path(nd, "frozen", site, pid, rung, f)
inst_path <- function(idir, pid, rung) {
  cand <- file.path(idir, paste0(pid, "_", rung, c(".laz", ".las")))
  hit <- cand[file.exists(cand)]; if (length(hit)) hit[1] else NA_character_
}

## ---- per-arm apex materialization for one (plot, rung) -------------------
# Returns a named list of (x,y,z AGL) detection frames, one per arm that ran on
# this cell (NULL/absent arms are dropped). `las` is the frozen normalized clip,
# `clip` its path, `dtm` the frozen DTM (for the deep arms' absolute->AGL).
materialize <- function(las, clip, dtm, frdens, res, nd, site, pid, rung) {
  out <- list()
  out$chm_vwf  <- tryCatch(detect_lasr(clip, res, VWF_A, frdens), error = function(e) NULL)
  out$multichm <- tryCatch({
    tt <- lidR::locate_trees(las, lidRplugins::multichm(res = res, ws = ws_factory(VWF_A)))
    co <- sf::st_coordinates(tt)
    if (!nrow(co)) data.frame(x = numeric(), y = numeric(), z = numeric())
    else data.frame(x = co[, 1], y = co[, 2],
                    z = if (ncol(co) >= 3) co[, 3] else tt$Z) }, error = function(e) NULL)
  if (rung == "native") out$li2012 <- tryCatch({
    if (sum(las$Z >= 2) < 1) data.frame(x = numeric(), y = numeric(), z = numeric())
    else { seg <- lidR::segment_trees(las, lidR::li2012(dt1 = 1.5, dt2 = 2, R = 2, hmin = 2))
      reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z") }
  }, error = function(e) NULL)
  # deep arms: persisted absolute-Z instance clouds -> AGL via the frozen DTM
  if (file.exists(dtm)) {
    sat <- inst_path(file.path(nd, "segmentanytree_instances"), pid, rung)
    if (!is.na(sat)) out$segmentanytree <- tryCatch({
      det <- read_instances_laz(sat, id_field = SAT_ID_FIELD)
      if (is.null(det)) NULL else det_to_agl(det, dtm) }, error = function(e) NULL)
    if (rung %in% c("native", "8")) {
      ff <- inst_path(file.path(nd, "forestformer3d_instances"), pid, rung)
      if (!is.na(ff)) out$forestformer3d <- tryCatch({
        det <- ff3d_collapse(ff, merge_tol = MERGE_TOL)
        if (is.null(det)) NULL else agl_guard(det, dtm) }, error = function(e) NULL)
    }
  }
  out[!vapply(out, is.null, logical(1))]            # drop arms that did not run
}

## ---- per (plot, rung) fusion + dual scoring ------------------------------
run_plot <- function(site, pid, pc, gt, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < MINTREES) return(NULL)
  rstem <- ifelse(is.finite(stems$maxCrownDiameter) & stems$maxCrownDiameter > 0,
                  stems$maxCrownDiameter / 2, FALLBACK_R)
  ref_class <- setNames(as.character(stems$crown_class),
                        as.character(seq_len(nrow(stems))))
  rows <- list()
  for (rung in RUNGS) {
    clip <- fz(nd, site, pid, rung, "clip_normalized.laz")
    if (!file.exists(clip)) next
    las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
    if (is.null(las) || lidR::is.empty(las)) next
    frdens <- tryCatch(jsonlite::read_json(fz(nd, site, pid, rung, "manifest.json"),
                       simplifyVector = TRUE)$frdens, error = function(e) NA_real_)
    if (is.null(frdens) || is.na(frdens)) frdens <- 8
    res <- if (frdens >= 8) 0.25 else 0.5
    dtm <- fz(nd, site, pid, rung, "ground_dtm.tif")

    arms <- materialize(las, clip, dtm, frdens, res, nd, site, pid, rung)
    n_arms <- length(arms)
    if (n_arms < 2) next                            # nothing to fuse

    # IoU substrate (canopy points in core) + Voronoi-on-stems reference (once)
    sd <- las@data
    kp <- sd$Z >= CANOPY_MIN & abs(sd$X - cx) <= ph & abs(sd$Y - cy) <= ph
    sx <- sd$X[kp]; sy <- sd$Y[kp]
    ref <- if (length(sx)) assign_points_to_stems(sx, sy, stems$E, stems$N, rstem) else integer(0)

    stack <- rbindlist(lapply(names(arms), function(a)
      data.table(arm = a, x = arms[[a]]$x, y = arms[[a]]$y, z = arms[[a]]$z)))
    fp <- fusion_points(stack$arm, stack$x, stack$y, stack$z, n_arms = n_arms,
                        merge_tol = MERGE_TOL, z_tol = Z_TOL,
                        chm_arms = CHM_ARMS, overstory_frac = OVERSTORY_FRAC)
    f_all <- fuse_apexes(stack$arm, stack$x, stack$y, stack$z, MERGE_TOL, Z_TOL)

    # configs: each single arm, the 3 fusion modes, and the k=1..N Pareto frontier
    cfgs <- c(lapply(names(arms), function(a) list(name = a, det = arms[[a]])),
              list(list(name = "union",    det = fp$union),
                   list(name = "majority", det = fp$majority),
                   list(name = "layered",  det = fp$layered)),
              lapply(seq_len(n_arms), function(k)
                list(name = sprintf("k%d", k),
                     det = f_all[f_all$votes >= k, c("x", "y", "z"), drop = FALSE])))
    for (cfg in cfgs) {
      det <- cfg$det
      if (is.null(det)) det <- data.frame(x = numeric(), y = numeric(), z = numeric())
      det <- det[, c("x", "y", "z"), drop = FALSE]
      sp <- score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                       core_half = ph)
      pred <- if (length(sx) && nrow(det))
        assign_points_to_stems(sx, sy, det$x, det$y, rep(FUSE_R, nrow(det))) else
        rep(NA_integer_, length(sx))
      ic <- score_instance_cell(pred, ref, ref_class = ref_class)
      rows[[length(rows) + 1]] <- cbind(
        data.frame(site = site, plot = pid, rung = rung, config = cfg$name,
                   n_arms = n_arms, n_apex = nrow(det), stringsAsFactors = FALSE),
        sp,
        data.frame(iou_n_ref = ic$n_ref, iou_TP = ic$TP, iou_FP = ic$FP,
                   iou_FN = ic$FN, iou_sum_iou = ic$sum_iou,
                   iou_sum_maxiou = ic$sum_maxiou))
    }
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows, fill = TRUE)
}

run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv"); pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) {
    cat(sprintf("[%s] no ground truth / centroids -- skipped\n", site)); return(NULL) }
  gt <- read.csv(gtf, stringsAsFactors = FALSE); pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  gt <- gt[, setdiff(names(gt), "maxCrownDiameter"), drop = FALSE]
  gt <- merge(gt, field_crowns(site), by = "individualID", all.x = TRUE)
  if (is.null(gt$maxCrownDiameter)) gt$maxCrownDiameter <- NA_real_
  plots <- intersect(unique(gt$plotID), pc$plotID)
  cat(sprintf("[%s] fusing %d plots over rungs {%s}\n", site, length(plots),
              paste(RUNGS, collapse = ",")))
  if (!length(plots)) return(NULL)
  res_list <- mclapply(plots, function(p)
    tryCatch(run_plot(site, p, pc, gt, nd),
             error = function(e) { message("  ", p, " failed: ", conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- rbindlist(Filter(Negate(is.null), res_list), fill = TRUE)
  if (!nrow(res)) { cat(sprintf("[%s] no cells fused\n", site)); return(NULL) }
  o <- file.path(nd, "fusion_results.csv"); write.csv(res, o, row.names = FALSE)
  cat(sprintf("[%s] wrote %d cell rows -> %s\n", site, nrow(res), o))
  as.data.frame(res)
}

## ---- pooled report (fusion vs best single arm; Pareto) --------------------
ARMS_ALL <- c("chm_vwf", "multichm", "li2012", "segmentanytree", "forestformer3d")
print_report <- function(res) {
  iou_pool <- function(sub) {                       # pooled IoU recall@.5 / cov / PQ
    TP <- sum(sub$iou_TP); FN <- sum(sub$iou_FN); FP <- sum(sub$iou_FP)
    nref <- sum(sub$iou_n_ref); si <- sum(sub$iou_sum_iou); sm <- sum(sub$iou_sum_maxiou)
    SQ <- if (TP > 0) si / TP else NA_real_
    RQ <- if (TP + 0.5 * FP + 0.5 * FN > 0) TP / (TP + 0.5 * FP + 0.5 * FN) else NA_real_
    c(rec = if (nref) TP / nref else NA_real_, cov = if (nref) sm / nref else NA_real_,
      PQ = if (is.finite(SQ) && is.finite(RQ)) SQ * RQ else 0)
  }
  # Equal-set guard: configs are only emitted on cells where their arms ran
  # (e.g. forestformer3d on native+8 only, single arms on their own subsets),
  # while union/majority/k* span all fused cells. Pooling each over its OWN cell
  # set makes the headline fusion-vs-single dF1 (and the k-of-N Pareto) a
  # cross-denominator comparison. Mirror equal_set_guard: restrict the compared
  # configs to the COMMON (site,plot,rung) cell set before pool(). `cells` keys a
  # row to its cell; `common_cells` intersects keys across a config set.
  cell_key <- function(df) do.call(paste, c(df[c("site", "plot", "rung")], sep = "::"))
  common_cells <- function(rr, configs) {
    keys <- lapply(configs, function(nm) {
      s <- rr[rr$config == nm, , drop = FALSE]; if (!nrow(s)) NULL else unique(cell_key(s)) })
    keys <- Filter(Negate(is.null), keys)
    if (!length(keys)) character(0) else Reduce(intersect, keys)
  }
  pool_on <- function(rr, nm, cells) {                # pool config nm over `cells` only
    s <- rr[rr$config == nm, , drop = FALSE]
    if (nrow(s)) s <- s[cell_key(s) %in% cells, , drop = FALSE]
    s
  }
  for (rung in unique(res$rung)) {
    rr <- res[res$rung == rung, , drop = FALSE]
    cat(sprintf("\n===== FUSION @ rung %s (pooled) =====\n", rung))
    cat(sprintf("%-16s %6s %6s %6s %6s | %7s %7s %7s\n",
                "config", "n_ref", "recall", "prec", "F1", "iou_R", "iou_Cov", "iou_PQ"))
    show <- c(intersect(ARMS_ALL, rr$config), "union", "majority", "layered")
    best_arm_f1 <- -1; best_arm <- NA
    for (nm in show) {
      sub <- rr[rr$config == nm, , drop = FALSE]; if (!nrow(sub)) next
      p <- pool(sub); iu <- iou_pool(sub)
      if (nm %in% ARMS_ALL && !is.na(p$F1) && p$F1 > best_arm_f1) { best_arm_f1 <- p$F1; best_arm <- nm }
      cat(sprintf("%-16s %6d %6.3f %6.3f %6.3f | %7.3f %7.3f %7.3f\n",
                  nm, p$n_ref, p$recall, p$precision, p$F1, iu["rec"], iu["cov"], iu["PQ"]))
    }
    # Headline dF1: pool best_arm / union / majority over their COMMON cell set so
    # the comparison is single-denominator (not best_arm_f1 from the table above).
    if (!is.na(best_arm)) {
      cc <- common_cells(rr, c(best_arm, "union", "majority"))
      if (length(cc)) {
        ba <- pool(pool_on(rr, best_arm, cc))
        u  <- pool(pool_on(rr, "union", cc)); mj <- pool(pool_on(rr, "majority", cc))
        cat(sprintf(paste0("best single arm: %s (F1 %.3f over %d common cells); ",
                           "union dF1 %+.3f; majority dF1 %+.3f\n"),
                    best_arm, ba$F1, length(cc), u$F1 - ba$F1, mj$F1 - ba$F1))
      } else {
        cat(sprintf("best single arm: %s -- no cells shared with union/majority\n", best_arm))
      }
    }
    cat("Pareto (k-of-N): ")
    ks <- sort(grep("^k[0-9]+$", unique(rr$config), value = TRUE))
    kc <- common_cells(rr, ks)                         # k* share a denominator too
    for (kk in ks) { p <- pool(pool_on(rr, kk, kc))
      cat(sprintf("%s R=%.3f P=%.3f F1=%.3f  ", kk, p$recall, p$precision, p$F1)) }
    cat("\n")
  }
}

run_main <- function() {
  t0 <- Sys.time(); all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- if (length(all_res)) do.call(rbind, all_res) else NULL
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d cell rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) print_report(res)
}

if (sys.nframe() == 0L) run_main()
