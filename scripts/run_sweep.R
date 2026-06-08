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

# Driver for the NEON SOAP density-ladder parameter sweep.
# Usage:
#   Rscript scripts/run_sweep.R [PLOTS=SOAP_031,SOAP_048|ALL] [OUT=results.csv]
#                               [CORES=8] [TOL=2.5] [MEAS_YEAR=2021]
# Reads work/neon/{ground_truth_stems.csv,plot_centroids.csv} and the LiDAR
# catalog in work/neon/lidar/. Writes a long-form metrics CSV: one row per
# (plot x density-rung x chm_res x vwf_slope).
#
# MEAS_YEAR (optional): when set (e.g. MEAS_YEAR=2021) the ground truth is
# additionally restricted to stems whose nearest apparentindividual measurement
# falls in that exact calendar year (the gt `meas_year` column). This drives the
# temporal-sensitivity cut (issue #5): re-score using ONLY stems measured in the
# LiDAR-acquisition year vs. the default +/-4 yr nearest-measurement baseline.
# When unset, behaviour is IDENTICAL to the baseline.
#
# OUT defaulting: when MEAS_YEAR is set AND OUT is omitted, OUT defaults to a
# distinct sweep_results_<YEAR>.csv under the site dir, so the exact-year subset
# never overwrites the canonical +/-4 yr baseline sweep_results.csv. When
# MEAS_YEAR is unset, OUT defaults to the canonical sweep_results.csv (byte-for-
# byte the prior behaviour). Pass OUT= explicitly to override either default.
suppressMessages({ library(lidR); library(parallel) })
options(lidR.progress = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))

## ---- args ----------------------------------------------------------------
args <- strsplit(commandArgs(TRUE), "=")
A <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE)) "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
nd    <- file.path(d, "neon", SITE)
CORES <- as.integer(if (is.null(A$CORES)) 8 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
MEAS_YEAR <- if (is.null(A$MEAS_YEAR)) NA_integer_ else as.integer(A$MEAS_YEAR)
## OUT path: when MEAS_YEAR is set but OUT is omitted, default to a distinct
## sweep_results_<YEAR>.csv under the site dir so the exact-year subset never
## clobbers the canonical +/-4 yr baseline sweep_results.csv. When MEAS_YEAR is
## unset, the default stays byte-identical to the baseline (sweep_results.csv).
OUT   <- if (!is.null(A$OUT)) {
  A$OUT
} else if (is.na(MEAS_YEAR)) {
  file.path(nd, "sweep_results.csv")
} else {
  file.path(nd, sprintf("sweep_results_%d.csv", MEAS_YEAR))
}

## ---- data ----------------------------------------------------------------
gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"))
pc  <- read.csv(file.path(nd, "plot_centroids.csv"))
gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
## temporal-sensitivity cut (issue #5): keep ONLY exact-measurement-year stems.
if (!is.na(MEAS_YEAR)) {
  n_before <- nrow(gt)
  gt <- gt[!is.na(gt$meas_year) & gt$meas_year == MEAS_YEAR, ]
  cat(sprintf("MEAS_YEAR filter = %d : kept %d of %d live-tree stems\n",
              MEAS_YEAR, nrow(gt), n_before))
} else {
  cat("MEAS_YEAR filter = none (+/-4 yr nearest-measurement baseline)\n")
}
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

if (is.null(results) || !nrow(results))
  stop("no sweep results produced; check input LiDAR, ground truth, and plot filters",
       call. = FALSE)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
write.csv(results, OUT, row.names = FALSE)
cat(sprintf("\nDONE: %d rows in %.1f min -> %s\n", nrow(results), dt, OUT))

s <- results[results$chm_res == 0.5 & results$vwf_a == 0.10, ]
if (nrow(s)) {
  cat("\noverall recall/precision by rung (pooled, res=0.5, a=0.10):\n")
  s$tp_core <- round(s$precision * s$n_det)
  agg <- do.call(rbind, lapply(c("native", "8", "4", "2", "1"), function(rl) {
    ss <- s[s$rung == rl, ]
    if (!nrow(ss)) return(NULL)
    recall <- sum(ss$TP) / sum(ss$n_ref)
    precision <- if (sum(ss$n_det) > 0) sum(ss$tp_core, na.rm = TRUE) /
      sum(ss$n_det) else NA_real_
    f1 <- if (!is.na(precision) && (recall + precision) > 0)
      2 * recall * precision / (recall + precision) else NA_real_
    data.frame(rung = rl, n_plots = length(unique(ss$plot)),
               n_ref = sum(ss$n_ref), n_det = sum(ss$n_det),
               recall = recall, precision = precision, F1 = f1)
  }))
  print(agg, row.names = FALSE)
}
