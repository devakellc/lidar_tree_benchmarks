#!/usr/bin/env Rscript
# Driver for the NEON SOAP density-ladder parameter sweep.
# Usage:
#   Rscript scripts/run_sweep.R [PLOTS=SOAP_031,SOAP_048|ALL] [OUT=results.csv]
#                               [CORES=8] [TOL=2.5]
# Reads work/neon/{ground_truth_stems.csv,plot_centroids.csv} and the LiDAR
# catalog in work/neon/lidar/. Writes a long-form metrics CSV: one row per
# (plot x density-rung x chm_res x vwf_slope).
suppressMessages({ library(lidR); library(parallel) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
source(file.path("scripts", "sweep_lib.R"))

## ---- args ----------------------------------------------------------------
args <- strsplit(commandArgs(TRUE), "=")
A <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE)) "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
nd    <- file.path(d, "neon", SITE)
OUT   <- if (is.null(A$OUT))   file.path(nd, "sweep_results.csv") else A$OUT
CORES <- as.integer(if (is.null(A$CORES)) 8 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)

## ---- data ----------------------------------------------------------------
gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"))
pc  <- read.csv(file.path(nd, "plot_centroids.csv"))
gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                  recursive = TRUE, full.names = TRUE)
ctg <- readLAScatalog(laz, progress = FALSE)

## plots to run: those with >= MINTREES live trees, intersected with PLOTS arg
MINTREES <- 6
counts <- table(gt$plotID)
keep   <- names(counts)[counts >= MINTREES]
if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
keep   <- intersect(keep, pc$plotID)
cat(sprintf("plots to sweep: %d (%s)\n", length(keep), paste(keep, collapse=",")))

RUNGS <- c(8, 4, 2, 1)          # native is added per-plot as the top rung
A_SET <- c(0.05, 0.10, 0.15)    # VWF slope

tmpdir <- file.path(tempdir(), "sweep"); dir.create(tmpdir, showWarnings = FALSE)

## ---- per-plot worker -----------------------------------------------------
run_plot <- function(pid) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)               # tower +/-20 m, distributed +/-10 m
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
    # no-upsampling guard: compare the rung TARGET to the plot's NATIVE density
    # (not the post-decimation density, which would false-skip on undershoot).
    else if (is.na(native_pdens) || rung >= native_pdens) { unlink(prep$file); next }
    res_set <- if (frdens >= 8) c(0.25, 0.5, 1.0) else c(0.5, 1.0)
    for (res in res_set) for (a in A_SET) {
      det <- tryCatch(detect_lasr(prep$file, res, a, frdens),  # frdens gates smoothing
                      error = function(e) NULL)
      if (is.null(det)) next
      sc <- score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                       core_half = ph)
      sc <- cbind(data.frame(plot = pid, plotType = ci$plotType,
                             rung = ifelse(is.na(rung), "native", as.character(rung)),
                             pdens = round(pdens, 2), frdens = round(frdens, 2),
                             chm_res = res, vwf_a = a),
                  sc)
      out[[length(out) + 1]] <- sc
    }
    unlink(prep$file)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

t0 <- Sys.time()
res_list <- mclapply(keep, function(p) tryCatch(run_plot(p), error=function(e) {
              message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
              mc.cores = CORES, mc.preschedule = FALSE)
results <- do.call(rbind, Filter(Negate(is.null), res_list))
dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(results, OUT, row.names = FALSE)
if (is.null(results)) results <- data.frame()
cat(sprintf("\nDONE: %d rows in %.1f min -> %s\n", nrow(results), dt, OUT))
if (nrow(results)) {
  cat("\noverall recall/precision by rung (mean over plots, res=0.5, a=0.10):\n")
  s <- results[results$chm_res == 0.5 & results$vwf_a == 0.10, ]
  agg <- aggregate(cbind(recall, precision, F1) ~ rung, data = s, mean, na.rm = TRUE)
  print(agg, row.names = FALSE)
}
