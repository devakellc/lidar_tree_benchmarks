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

# multichm treetop arm of the canonical NEON density ladder (#37).
#
# The canonical density-ladder pipeline (run_sweep.R -> density-ladder-sweep-
# results.md) scores CHM-VWF only (lasR local_maximum_raster on a pit_fill CHM).
# This adds a parallel `multichm` arm (Eysn-style multi-layer CHM local maxima,
# lidRplugins::multichm) built on the SAME prepare_clip lasR-based path the
# cached CHM-VWF sweep_results.csv was built on -- identical clip provenance
# (clip -> decimate_points -> normalize_height), so the head-to-head stays
# internally consistent (unlike the frozen-clip pitfree model benchmark).
#
# Per plot x density rung {native, 8, 4, 2, 1 pts/m^2}, with the same
# no-upsampling guard (rung target vs the plot's native all-return density) and
# plot-type-aware core (plot_half: tower +/-20 m, distributed +/-10 m). multichm
# uses a single density-derived res (0.25 m if first-return density >= 8 else
# 0.5 m) and the clamped variable window ws_factory(0.10) -- same discipline as
# detect_lidrplugins_sweep.R::det_multichm; there is no chm_res/vwf_a grid here.
#
# Dependency note: this arm needs only lidR + lidRplugins (prepare_clip uses lidR
# only and we never call lasR's local_maximum_raster(ws = f)), so CRAN lasR is
# fine -- the lasR pre-devel build is NOT required to run the multichm arm. The
# CHM-VWF baseline comes from the cached sweep_results.csv via analyze_multichm_sweep.R.
#
# Usage:
#   Rscript scripts/detect_multichm_sweep.R [SITE=SOAP] [PLOTS=ALL] [CORES=8]
#       [TOL=4] [A=0.10]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/multichm_sweep_results.csv (one row per
#         plot x rung).
suppressMessages({ library(lidR); library(lidRplugins); library(sf)
                   library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))        # prepare_clip, plot_half, score_plot, ws_factory
source(.find("model_bench_lib.R"))  # assert_detection_contract

## ---- tops sf -> (x,y,z) detection contract -------------------------------
# Pull lowercase x,y,z from a locate_trees() sf. multichm geometry is 2-D, so
# the apex z is read from the tops' Z attribute (the same fallback the
# detect_lidrplugins_sweep.R::.tops_to_det uses). assert_detection_contract
# enforces the exact shape score_plot consumes.
.tops_to_det <- function(tt) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.null(tt) || nrow(tt) == 0) { assert_detection_contract(empty); return(empty) }
  co <- sf::st_coordinates(tt)
  z  <- if (ncol(co) >= 3) co[, 3] else if (!is.null(tt$Z)) tt$Z else co[, 2] * 0
  det <- data.frame(x = co[, 1], y = co[, 2], z = z)
  assert_detection_contract(det)
  det
}

# multichm: multi-layer CHM local maxima at a density-derived res; pass our
# clamped variable window (same allometry as CHM-VWF) through to the internal
# lmf for a fair comparison. Returns the x,y,z contract, or NULL on a detector
# crash (the caller skips that cell -- the equal-set guard drops it downstream).
det_multichm_run <- function(las, res = 0.5, a = 0.10) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::multichm(res = res,
                                                               ws = ws_factory(a))),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  .tops_to_det(tt)
}

## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 8 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
A_VWF <- as.numeric(if (is.null(A$A))   0.10 else A$A)
RUNGS    <- c(8, 4, 2, 1)        # native is added per-plot as the top rung
MINTREES <- 6                    # min live trees to sweep a plot (matches run_sweep.R)

run_main <- function() {
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),     stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  counts <- table(gt$plotID)
  keep   <- names(counts)[counts >= MINTREES]
  if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  keep   <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] multichm plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))

  tmpdir <- file.path(tempdir(), "multichm"); dir.create(tmpdir, showWarnings = FALSE)

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)            # tower +/-20 m, distributed +/-10 m
    stems <- gt[gt$plotID == pid &
                abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    out <- list(); native_pdens <- NA_real_
    # native first (rung = NA) -> captures native density; then decimated rungs.
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(prepare_clip(ctg, cx, cy, rung, tmpdir, core_half = ph),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      # no-upsampling guard: rung TARGET vs the plot's NATIVE all-return density.
      else if (is.na(native_pdens) || rung >= native_pdens) { unlink(prep$file); next }
      las <- tryCatch(readLAS(prep$file), error = function(e) NULL)
      if (is.null(las) || is.empty(las)) { unlink(prep$file); next }
      res <- if (frdens >= 8) 0.25 else 0.5      # density-derived, like CHM-VWF
      det <- det_multichm_run(las, res = res, a = A_VWF)
      if (is.null(det)) { unlink(prep$file); next }  # detector crash -> skip cell
      sc <- score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                       core_half = ph)
      sc <- cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                             rung = ifelse(is.na(rung), "native", as.character(rung)),
                             pdens = round(pdens, 2), frdens = round(frdens, 2),
                             chm_res = res, n_apex = nrow(det)), sc)
      out[[length(out) + 1]] <- sc
      unlink(prep$file)
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }

  t0 <- Sys.time()
  res_list <- mclapply(keep, function(p) tryCatch(run_plot(p), error = function(e) {
                message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(results) || !nrow(results)) {
    cat("no multichm results produced; check input LiDAR, ground truth, filters\n")
    return(invisible())
  }
  # core true positives for the precision denominator (poolers recover them here).
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "multichm_sweep_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] multichm DONE: %d rows in %.1f min -> %s\n", SITE, nrow(results),
              dt, file.path(nd, "multichm_sweep_results.csv")))

  cat("\noverall recall/precision by rung (pooled; res density-derived, a=0.10):\n")
  agg <- do.call(rbind, lapply(c("native", "8", "4", "2", "1"), function(rl) {
    ss <- results[results$rung == rl, ]
    if (!nrow(ss)) return(NULL)
    recall    <- sum(ss$TP) / sum(ss$n_ref)
    precision <- if (sum(ss$n_det) > 0) sum(ss$tp_core, na.rm = TRUE) /
      sum(ss$n_det) else NA_real_
    f1 <- if (!is.na(precision) && (recall + precision) > 0)
      2 * recall * precision / (recall + precision) else NA_real_
    data.frame(rung = rl, n_plots = length(unique(ss$plot)),
               n_ref = sum(ss$n_ref), n_det = sum(ss$n_det),
               recall = recall, precision = precision, F1 = f1)
  }))
  print(agg, row.names = FALSE, digits = 2)
}

if (sys.nframe() == 0L) run_main()
