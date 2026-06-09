#!/usr/bin/env Rscript
# SegmentAnyTree (#M6) density-ladder arm. Runs the upstream SAT instance
# segmenter zero-shot (Docker, sm_120 rebuild) on the RAW-WITH-GROUND frozen clip
# per plot x density rung, SERIALLY (single GPU -- no mclapply contention),
# reducing per-tree instances to apex detections and scoring against field stems
# with the existing harness.
#
# Per-cell pipeline:
#   rawground.laz --laz_to_ply--> clip.ply
#     --run_docker_arm(sat-sm120-test, run_segmentanytree.py)--> merged.las
#     --read_instances_laz(PredInstance)--> reduce_instances (max-Z per id, drop
#       0/NA) --det_to_agl(ground_dtm.tif)--> apex(x,y,z AGL) --score_plot.
# SAT self-localizes XY but keeps ABSOLUTE UTM Z, so det_to_agl is REQUIRED here
# (unlike the normalized-clip TreeisoNet arm, which needs no DTM transform).
#
# Usage:
#   Rscript scripts/detect_segmentanytree_sweep.R [SITE=SOAP] [PLOTS=ALL]
#       [TOL=4] [IMAGE=sat-sm120-test] [TIMEOUT=1800]
# Requires the sm_120 image:  bash gpu/segmentanytree-sm120/build.sh
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/segmentanytree_results.csv (one row per
#         plot x rung).
suppressMessages({ library(lidR); library(data.table) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
bs <- Find(file.exists, c(
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs)
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))
source(.find("model_runner.R")); source(.find("io_bridge.R"))

args    <- strsplit(commandArgs(TRUE), "=")
A       <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE    <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS   <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
TOL     <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
IMAGE   <- if (is.null(A$IMAGE)) "sat-sm120-test" else A$IMAGE
TIMEOUT <- as.numeric(if (is.null(A$TIMEOUT)) 1800 else A$TIMEOUT)
RUNGS    <- c(8, 4, 2, 1); MINTREES <- 6
DRIVER   <- file.path(.ROOT, "gpu", "run_segmentanytree.py")
# The merged LAS carries the instance label as the PredInstance extra dim
# (merge_pt_ss_is.py: preds_instance_segmentation -> PredInstance, +1, NA->0).
ID_FIELD <- "PredInstance"

run_main <- function() {
  stopifnot(file.exists(DRIVER))
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
  cat(sprintf("[%s] segmentanytree plots: %d (image=%s tol=%s)\n",
              SITE, length(keep), IMAGE, TOL))

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
      tag  <- ifelse(is.na(rung), "native", as.character(rung))
      ply  <- file.path(tempdir(), sprintf("sat_%s_%s.ply", pid, tag))
      olas <- file.path(tempdir(), sprintf("sat_%s_%s.las", pid, tag))
      # Absolute-UTM double PLY (lossless); SAT does its own XY localization.
      tryCatch(laz_to_ply(prep$rawground, ply), error = function(e) NULL)
      if (!file.exists(ply)) next
      det_abs <- run_docker_arm(IMAGE, ply, olas,
                   cmd     = c("python3", DRIVER),
                   mounts  = dirname(DRIVER),
                   # SAT's clustering Pool spawns workers; --ipc=host + a real
                   # shm avoid multiprocessing stalls under the default 64M shm.
                   extra_docker = c("--shm-size=8g", "--ipc=host"),
                   reader  = function(p) read_instances_laz(p, id_field = ID_FIELD),
                   gpus    = "all", timeout = TIMEOUT,
                   label   = sprintf("%s/%s", pid, tag))
      if (is.null(det_abs)) next            # container crash/schema -> skip cell
      det <- det_to_agl(det_abs, prep$dtm)  # absolute Z -> height above ground
      sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                core_cy = cy, core_half = ph),
                     error = function(e) NULL)
      if (is.null(sc)) next
      out[[length(out) + 1]] <- cbind(data.frame(site = SITE, plot = pid,
        plotType = ci$plotType, detector = "segmentanytree", rung = tag,
        pdens = round(pdens, 2), frdens = round(frdens, 2), n_apex = nrow(det)), sc)
      ncell <- ncell + 1L
    }
    cat(sprintf("  %s: %d cells\n", pid, ncell))
  }
  results <- do.call(rbind, out)
  if (is.null(results) || !nrow(results)) {
    cat("no segmentanytree results\n"); return(invisible())
  }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "segmentanytree_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] segmentanytree DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "segmentanytree_results.csv")))
}

if (sys.nframe() == 0L) run_main()
