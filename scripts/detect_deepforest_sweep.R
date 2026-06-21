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

# #X1 DeepForest RGB arm (GitHub issue #68).
#
# Every LiDAR arm degrades as pts/m2 drops; an RGB detector does NOT, so its flat
# accuracy-vs-density curve is the density-INVARIANT anchor that tells the router
# (#P2) which rung each LiDAR arm stops beating optical, and gives #P1 its first
# optical fusion member. This runs DeepForest's NEON-pretrained crown model
# (predict_tile) on the downloaded RGB tiles (neon_download_aop.R, DP3.30010),
# reduces each crown box to its centroid, samples an apex Z from a CHM built on
# demand from the matched frozen normalized clip, and scores against field stems
# by crown_class with the SAME score_plot harness. detector = "deepforest",
# rung = "rgb" (RGB has no density ladder -- that is the point).
#
# Usage:
#   Rscript scripts/detect_deepforest_sweep.R SITE=SOAP
# Env: PYTHON=~/miniconda3/envs/deepforest/bin/python (deepforest 2.x, CPU ok).
# Reads work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv,rgb/*.tif} +
#   the frozen normalized clips. Writes work/neon/<SITE>/deepforest_results.csv
#   and caches per-tile boxes in deepforest_boxes/.
suppressMessages({ library(lidR); library(terra); library(data.table) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))

args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE <- if (!is.null(A$SITE)) A$SITE else "SOAP"
TOL  <- as.numeric(if (is.null(A$TOL)) 4 else A$TOL)
SCORE_THRESH <- as.numeric(if (is.null(A$SCORE_THRESH)) 0.0 else A$SCORE_THRESH)
CHM_RES <- as.numeric(if (is.null(A$CHM_RES)) 0.5 else A$CHM_RES)
PATCH <- as.integer(if (is.null(A$PATCH)) 400 else A$PATCH)
PYTHON <- path.expand(if (!is.null(A$PYTHON)) A$PYTHON else "~/miniconda3/envs/deepforest/bin/python")
RUNNER <- Find(file.exists, c("gpu/run_deepforest.py", file.path(getwd(), "gpu/run_deepforest.py")))
nd <- file.path(d, "neon", SITE)

## ---- run DeepForest on every RGB tile (cached) -> all crown centroids -----
all_boxes <- function() {
  rdir <- file.path(nd, "rgb")
  tifs <- list.files(rdir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
  if (!length(tifs)) { cat("no RGB tiles under", rdir, "\n"); return(NULL) }
  bdir <- file.path(nd, "deepforest_boxes"); dir.create(bdir, showWarnings = FALSE, recursive = TRUE)
  rows <- list()
  for (tif in tifs) {
    ocsv <- file.path(bdir, paste0(tools::file_path_sans_ext(basename(tif)), ".csv"))
    if (!file.exists(ocsv)) {
      st <- tryCatch(system2(PYTHON, c(shQuote(RUNNER), shQuote(tif), shQuote(ocsv),
                     as.character(PATCH)), stdout = FALSE, stderr = FALSE), error = function(e) 1L)
      if (!identical(as.integer(st), 0L) || !file.exists(ocsv)) {
        cat("deepforest failed on", basename(tif), "\n"); next }
    }
    b <- tryCatch(read.csv(ocsv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(b) && nrow(b)) rows[[length(rows) + 1]] <- b
  }
  if (!length(rows)) return(NULL)
  b <- rbindlist(rows, fill = TRUE)
  b[b$score >= SCORE_THRESH, , drop = FALSE]
}

## ---- per-plot CHM (from frozen normalized clip) to give boxes an apex z ----
plot_chm <- function(pid) {
  clip <- file.path(nd, "frozen", SITE, pid, "native", "clip_normalized.laz")
  if (!file.exists(clip)) return(NULL)
  las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)
  tryCatch(suppressWarnings(lidR::rasterize_canopy(las, res = CHM_RES, algorithm = lidR::p2r())),
           error = function(e) NULL)
}

run_main <- function() {
  if (is.null(RUNNER) || !file.exists(PYTHON)) {
    cat(sprintf("PYTHON=%s exists=%s RUNNER=%s\n", PYTHON, file.exists(PYTHON), RUNNER)) }
  gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  boxes <- all_boxes(); if (is.null(boxes)) { cat("no DeepForest boxes\n"); return(invisible()) }
  cat(sprintf("[%s] DeepForest crown boxes (all tiles, score>=%.2f): %d\n",
              SITE, SCORE_THRESH, nrow(boxes)))
  keep <- intersect(unique(gt$plotID), pc$plotID)
  rows <- list()
  for (pid in keep) {
    ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
    if (!nrow(stems)) next
    bp <- boxes[abs(boxes$x - cx) <= ph + TOL & abs(boxes$y - cy) <= ph + TOL, , drop = FALSE]
    chm <- plot_chm(pid)
    z <- if (!is.null(chm) && nrow(bp))
      as.numeric(terra::extract(chm, cbind(bp$x, bp$y))[, 1]) else rep(NA_real_, nrow(bp))
    z[!is.finite(z)] <- 2.0                              # off-CHM crowns: floor at min_height
    det <- data.frame(x = bp$x, y = bp$y, z = z)
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                              core_half = ph), error = function(e) NULL)
    if (is.null(sc)) next
    rows[[length(rows) + 1]] <- cbind(
      data.frame(site = SITE, plot = pid, plotType = ci$plotType, detector = "deepforest",
                 rung = "rgb", frdens = NA_real_, n_apex = nrow(det), stringsAsFactors = FALSE), sc)
  }
  res <- rbindlist(rows, fill = TRUE)
  if (!nrow(res)) { cat("no rows scored\n"); return(invisible()) }
  o <- file.path(nd, "deepforest_results.csv"); write.csv(res, o, row.names = FALSE)
  p <- pool(res)
  cat(sprintf("\n[%s] DeepForest (rgb): n_ref=%d recall=%.3f prec=%.3f F1=%.3f rec_und=%.3f\n",
              SITE, p$n_ref, p$recall, p$precision, p$F1, p$rec_understory))
  cat(sprintf("wrote %d plot rows -> %s\n", nrow(res), o))
}

if (sys.nframe() == 0L) run_main()
