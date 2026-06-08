#!/usr/bin/env Rscript
# TreeisoNet (#M7) density-ladder arm. Runs the headless TreeisoNet apex driver
# (gpu/run_treeisonet.py) on the NORMALIZED frozen clip per plot x density rung,
# SERIALLY (single GPU -- no mclapply contention), scoring against field stems
# with the existing harness. Apex-only (treeLoc -> postPeakExtraction ->
# local-canopy-max z-snap); the treeOff crown variant is issue #20. `conf` is a
# fixed zero-shot threshold, calibrated once (NOT per plot). Normalized Z is
# already height-above-ground, so no DTM transform is needed (see the GPU-arm
# plan's implementation findings).
#
# Usage:
#   Rscript scripts/detect_treeisonet_sweep.R [SITE=SOAP] [PLOTS=ALL]
#       [CONF=0.22] [VOXEL=0] [TOL=4]
# Requires the venv + weights from gpu/setup_treeisonet_env.sh + gpu/mirror_weights.sh.
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/treeisonet_results.csv (one row per
#         plot x rung).
suppressMessages({ library(lidR); library(data.table) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
.script_path <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) && nzchar(ofile)) return(ofile)
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit)) return(sub("^--file=", "", hit[1]))
  NA_character_
}
.sp <- .script_path()
.ROOT <- if (!is.na(.sp)) normalizePath(file.path(dirname(.sp), ".."),
                                        mustWork = FALSE) else getwd()
.find <- function(rel) Find(file.exists, c(file.path(.ROOT, "scripts", rel),
                                           file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))
source(.find("model_runner.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CONF  <- if (is.null(A$CONF))  "0.22" else A$CONF
VOXEL <- if (is.null(A$VOXEL)) "0" else A$VOXEL
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
RUNGS <- c(8, 4, 2, 1); MINTREES <- 6
VENV  <- file.path(.ROOT, "gpu/.venv/bin/python")
DRV   <- file.path(.ROOT, "gpu/run_treeisonet.py")
LOC   <- file.path(.ROOT, "gpu/store/treeaibox/als_treeloc.pth")
CFG   <- Sys.glob(file.path(.ROOT, "gpu/store/treeaibox/*reclamation*treeloc*.json"))[1]

run_main <- function() {
  stopifnot(file.exists(VENV), file.exists(DRV), file.exists(LOC), !is.na(CFG))
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
  cat(sprintf("[%s] treeisonet plots: %d (conf=%s voxel=%s)\n",
              SITE, length(keep), CONF, VOXEL))

  out <- list()
  for (pid in keep) {                       # SERIAL -- one GPU, no contention
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) next
    native_pdens <- NA_real_; ncell <- 0L
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                   out_root = file.path(nd, "frozen")),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) next
      ocsv <- file.path(tempdir(), sprintf("ti_%s_%s.csv", pid,
                        ifelse(is.na(rung), "native", rung)))
      det <- run_python_arm(VENV, DRV, prep$normalized, ocsv,
                            extra = c(LOC, CFG, VOXEL, CONF), timeout = 900)
      if (is.null(det)) next                # GPU crash -> skip cell (guard drops)
      sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                core_cy = cy, core_half = ph),
                     error = function(e) NULL)
      if (is.null(sc)) next
      out[[length(out) + 1]] <- cbind(data.frame(site = SITE, plot = pid,
        plotType = ci$plotType, detector = "treeisonet",
        rung = ifelse(is.na(rung), "native", as.character(rung)),
        pdens = round(pdens, 2), frdens = round(frdens, 2), n_apex = nrow(det)), sc)
      ncell <- ncell + 1L
    }
    cat(sprintf("  %s: %d cells\n", pid, ncell))
  }
  results <- do.call(rbind, out)
  if (is.null(results) || !nrow(results)) { cat("no treeisonet results\n"); return(invisible()) }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "treeisonet_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] treeisonet DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "treeisonet_results.csv")))
}

if (sys.nframe() == 0L) run_main()
