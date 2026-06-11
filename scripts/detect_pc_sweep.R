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

# Point-cloud-detector vs CHM-VWF comparison at NATIVE density (issue #6).
#
# The density-ladder sweep (run_sweep.R) is CHM-VWF only. A 2.5-D canopy height
# model can only see the top surface, so the one place high point density SHOULD
# help is a *point-cloud* detector able to resolve sub-dominant apexes the CHM
# cannot. This script runs three point-cloud detectors plus the CHM-VWF baseline
# on the SAME native-density clip per plot, and scores each against the SAME
# field stems by crown class -- so understory (intermediate + suppressed) recall
# can be compared apples-to-apples.
#
# Detectors (all on the native, normalized clip):
#   chm_vwf       : detect_lasr(file, res, a, frdens) -- CHM-VWF baseline. The
#                   res is density-derived per plot (0.25 m if frdens>=8 else
#                   0.5 m), mirroring run_sweep.R's res_set -- NOT hardcoded.
#   lidr_lmf_pc   : lidR::locate_trees(las, lmf(ws, hmin, shape="circular")) on
#                   the POINT CLOUD (not the CHM).
#   lidr_li2012   : lidR::segment_trees(las, li2012(...)); apex per segment =
#                   the max-Z point of each treeID (x,y,z).
#   lasr_lmax_pc  : lasR point-based local_maximum(ws, min_height) (NOT
#                   local_maximum_raster); returns an sf of points (x,y,z).
#
# All point-cloud windows use the SAME variable-window allometry as the CHM
# path: ws <- ws_factory(A) (slope a, clamped [3,5]); hmin/min_height = 2.
#
# Usage:
#   Rscript scripts/detect_pc_sweep.R [SITES=SJER,SOAP,TEAK] [CORES=6]
#                                     [TOL=4] [A=0.10]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/pc_detect_results.csv  (one row per
#         plot x detector, with overall + per-crown-class recall columns).
suppressMessages({ library(lidR); library(lasR); library(sf)
                   library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("pc_detect_lib.R"))   # shared point-cloud apex extractors (#6/#38)

## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER","SOAP","TEAK") else
           strsplit(A$SITES, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL))  4.0  else A$TOL)
A_VWF <- as.numeric(if (is.null(A$A))    0.10 else A$A)
MINTREES <- 6

## ---- per-detector apex extractors -----------------------------------------
# det_lidr_lmf_pc / det_lidr_li2012 / det_lasr_lmax_pc now live in
# pc_detect_lib.R (shared with the #38 density ladder, detect_pc_ladder.R) so
# the two scripts cannot drift apart. Each returns a base data.frame(x, y, z).

## ---- per-plot worker: run all 4 detectors, score each --------------------
run_plot <- function(pid, gt, pc, ctg, tmpdir) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)                 # tower +/-20, distributed +/-10
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
  if (nrow(stems) < 1) return(NULL)

  prep <- tryCatch(prepare_clip(ctg, cx, cy, NA, tmpdir, core_half = ph),
                   error = function(e) NULL)
  if (is.null(prep)) return(NULL)
  on.exit(unlink(prep$file), add = TRUE)
  frdens <- prep$frdens; pdens <- prep$pdens
  las <- tryCatch(readLAS(prep$file), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)
  ws <- ws_factory(A_VWF)

  # CHM-VWF baseline resolution is density-derived, NOT hardcoded (CLAUDE.md
  # Step 0 / run_sweep.R res_set): finest rung is 0.25 m where first-return
  # density supports it (>=8 pts/m^2), else 0.5 m. The point-cloud detectors use
  # the ws window, not a raster res, so they are unaffected by this.
  res <- if (frdens >= 8) 0.25 else 0.5

  # Run each detector, timing the point-cloud ones (cost is part of the answer).
  dets <- list()
  tim  <- c()
  t0 <- Sys.time()
  dets$chm_vwf <- tryCatch(detect_lasr(prep$file, res, A_VWF, frdens),
                           error = function(e) NULL)
  tim["chm_vwf"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lidr_lmf_pc <- tryCatch(det_lidr_lmf_pc(las, ws),
                               error = function(e) NULL)
  tim["lidr_lmf_pc"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lidr_li2012 <- tryCatch(det_lidr_li2012(las), error = function(e) NULL)
  tim["lidr_li2012"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lasr_lmax_pc <- tryCatch(det_lasr_lmax_pc(prep$file, ws),
                                error = function(e) NULL)
  tim["lasr_lmax_pc"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out <- list()
  for (nm in names(dets)) {
    det <- dets[[nm]]
    if (is.null(det)) next
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                              core_cy = cy, core_half = ph),
                   error = function(e) NULL)
    if (is.null(sc)) next
    sc <- cbind(data.frame(plot = pid, plotType = ci$plotType,
                           detector = nm,
                           frdens = round(frdens, 2),
                           pdens = round(pdens, 2),
                           chm_res = if (nm == "chm_vwf") res else NA_real_,
                           n_apex = nrow(det),
                           secs = round(tim[[nm]], 2)),
                sc)
    out[[length(out) + 1]] <- sc
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

## ---- per-site driver -----------------------------------------------------
run_site <- function(SITE) {
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"),
                  stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),
                  stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)

  counts <- table(gt$plotID)
  keep   <- names(counts)[counts >= MINTREES]
  keep   <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))
  if (!length(keep)) return(NULL)

  tmpdir <- file.path(tempdir(), paste0("pcsweep_", SITE))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)

  t0 <- Sys.time()
  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p, gt, pc, ctg, tmpdir),
                         error = function(e) {
                           message("plot ", p, " failed: ",
                                   conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(results)) { cat(sprintf("[%s] no results\n", SITE)); return(NULL) }

  out <- file.path(nd, "pc_detect_results.csv")
  write.csv(results, out, row.names = FALSE)
  cat(sprintf("[%s] DONE: %d rows (%d plots x detectors) in %.1f min -> %s\n",
              SITE, nrow(results), length(unique(results$plot)), dt, out))
  results$site <- SITE
  results
}

## ---- run all sites, then pool --------------------------------------------
all_list <- lapply(SITES, run_site)
all <- do.call(rbind, Filter(Negate(is.null), all_list))
if (is.null(all) || !nrow(all)) { cat("no results across any site\n"); quit() }

## ---- equal-plot-set guard ------------------------------------------------
## A detector arm can fail or return an empty/unscorable result on a plot (the
## tryCatch -> next in run_plot). If we pooled each detector over whatever plots
## it happened to score, the deltas would compare detectors over NON-identical
## plot sets. Restrict every detector to the COMMON set of (site,plot) where all
## four detectors produced a scored row; drop any plot where any arm failed, and
## log how many and which were dropped.
detectors <- c("chm_vwf","lidr_lmf_pc","lidr_li2012","lasr_lmax_pc")
all$pkey <- paste(all$site, all$plot, sep = "::")
n_det_per_plot <- tapply(all$detector, all$pkey,
                         function(v) length(unique(v)))
common  <- names(n_det_per_plot)[n_det_per_plot == length(detectors)]
dropped <- setdiff(unique(all$pkey), common)
cat(sprintf("\n=== Equal-plot-set guard ===\n"))
cat(sprintf("plots with a scored row from all %d detectors: %d\n",
            length(detectors), length(common)))
if (length(dropped)) {
  cat(sprintf("DROPPED %d plot(s) where >=1 detector arm failed/empty:\n",
              length(dropped)))
  for (pk in dropped) {
    got <- sort(unique(all$detector[all$pkey == pk]))
    miss <- setdiff(detectors, got)
    cat(sprintf("  %-18s  missing: %s\n", pk, paste(miss, collapse = ",")))
  }
} else {
  cat("DROPPED 0 plots: every detector scored every plot.\n")
}
if (!length(common)) { cat("no common-plot results to pool\n"); quit() }
all <- all[all$pkey %in% common, , drop = FALSE]

## assert equal n_plots across detectors before any delta is reported
np <- tapply(all$pkey, all$detector, function(v) length(unique(v)))
cat("n_plots per detector after restriction: ",
    paste(sprintf("%s=%d", names(np), np), collapse = "  "), "\n")
stopifnot(length(unique(as.integer(np))) == 1L)

## pooling (mirror analyze_sweep.R: sum TP / sum n_ref, never mean-of-rates).
classes <- c("dominant","codominant","intermediate","suppressed")
# recover per-plot per-class TP counts for correct pooling
all$tp_core <- round(all$precision * all$n_det)
for (cl in classes) {
  nc <- all[[paste0("n_", cl)]]; rc <- all[[paste0("rec_", cl)]]
  all[[paste0("tp_", cl)]] <- ifelse(nc > 0, round(rc * nc), 0)
}

pool <- function(df) {
  out <- data.frame(
    n_plots   = length(unique(paste(df$site, df$plot))),
    frdens    = round(mean(df$frdens), 2),
    n_ref     = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall    = sum(df$TP) / sum(df$n_ref),
    precision = if (sum(df$n_det) > 0) sum(df$tp_core, na.rm = TRUE) / sum(df$n_det) else NA_real_,
    secs_med  = round(median(df$secs), 2))
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  for (cl in classes) {
    nref <- sum(df[[paste0("n_", cl)]], na.rm = TRUE)
    tp   <- sum(df[[paste0("tp_", cl)]], na.rm = TRUE)
    out[[paste0("rec_", cl)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_", cl)]]   <- nref
  }
  # combined understory = intermediate + suppressed (pooled counts)
  nref_u <- sum(df$n_intermediate, na.rm = TRUE) +
            sum(df$n_suppressed,   na.rm = TRUE)
  tp_u   <- sum(df$tp_intermediate, na.rm = TRUE) +
            sum(df$tp_suppressed,   na.rm = TRUE)
  out$rec_understory <- if (nref_u) tp_u / nref_u else NA_real_
  out$n_understory   <- nref_u
  out
}

cat("\n=== Pooled across all sites/plots (common plot set), by detector ===\n")
pooled <- do.call(rbind, lapply(detectors, function(dt) {
  s <- all[all$detector == dt, ]; if (!nrow(s)) return(NULL)
  cbind(detector = dt, pool(s)) }))
print(pooled[, c("detector","n_plots","n_ref","n_det","recall","precision","F1",
                 "rec_dominant","rec_codominant","rec_intermediate",
                 "rec_suppressed","rec_understory","n_understory","secs_med")],
      row.names = FALSE, digits = 3)

## per-site pooled understory recall (for the cross-site gradient)
cat("\n=== Understory (intermediate+suppressed) recall by site x detector ===\n")
for (S in unique(all$site)) {
  cat(sprintf("-- %s --\n", S))
  ps <- do.call(rbind, lapply(detectors, function(dt) {
    s <- all[all$site == S & all$detector == dt, ]; if (!nrow(s)) return(NULL)
    cbind(detector = dt, pool(s)) }))
  print(ps[, c("detector","recall","rec_understory","n_understory")],
        row.names = FALSE, digits = 3)
}

## deltas vs CHM-VWF baseline
base <- pooled[pooled$detector == "chm_vwf", ]
cat("\n=== Delta vs CHM-VWF (pooled) ===\n")
for (i in seq_len(nrow(pooled))) {
  p <- pooled[i, ]
  cat(sprintf("  %-13s  recall %+0.3f   understory %+0.3f\n", p$detector,
              p$recall - base$recall, p$rec_understory - base$rec_understory))
}

write.csv(pooled, file.path(d, "neon", "pc_detect_pooled.csv"),
          row.names = FALSE)
cat(sprintf("\npooled -> %s\n", file.path(d, "neon", "pc_detect_pooled.csv")))
