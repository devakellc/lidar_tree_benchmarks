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

# #X3 Detectree2 RGB arm (GitHub issue #76).
#
# A second OPTICAL detector for ensemble diversity: Detectree2 is a Detectron2
# Mask R-CNN crown segmenter, architecturally different from #X1's DeepForest
# (RetinaNet boxes). Its agreement with DeepForest is a confidence signal, and
# its POLYGONS let the optical modality contribute crown WIDTH (d_eq), not just
# detection. Runs the vendored runner (gpu/run_detectree2.py) on a per-plot RGB
# CROP (a 1 km2 mosaic is ~625 Mask R-CNN sub-tiles -- too slow; plot crops are
# fast). Detectron2 built CPU-only against torch 2.12 (no Blackwell CUDA). The
# pretrained weights are tropical-trained, so CA-conifer transfer is the tested
# unknown this arm measures.
#
# Per plot: crop the covering RGB tile -> runner -> crown polygon centroids +
# d_eq -> filter to core -> apex Z from the frozen-clip CHM -> score_plot
# (detector="detectree2", rung="rgb"). Crown d_eq distribution is reported vs
# field maxCrownDiameter on matched stems.
#
# Usage: Rscript scripts/detect_detectree2_sweep.R SITE=SOAP \
#          MODEL=~/.detectree2_models/250312_flexi.pth
# Env: PYTHON=~/miniconda3/envs/detectree2/bin/python. Reads rgb/ + frozen clips;
# writes work/neon/<SITE>/detectree2_results.csv.
suppressMessages({ library(lidR); library(terra); library(data.table) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))

args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE <- if (!is.null(A$SITE)) A$SITE else "SOAP"
PLOTS <- if (!is.null(A$PLOTS)) strsplit(A$PLOTS, ",")[[1]] else NULL
TOL  <- as.numeric(if (is.null(A$TOL)) 4 else A$TOL)
CHM_RES <- as.numeric(if (is.null(A$CHM_RES)) 0.5 else A$CHM_RES)
PYTHON <- path.expand(if (!is.null(A$PYTHON)) A$PYTHON else "~/miniconda3/envs/detectree2/bin/python")
MODEL  <- path.expand(if (!is.null(A$MODEL)) A$MODEL else "~/.detectree2_models/250312_flexi.pth")
RUNNER <- Find(file.exists, c("gpu/run_detectree2.py", file.path(getwd(), "gpu/run_detectree2.py")))
nd <- file.path(d, "neon", SITE)

rgb_tiles <- list.files(file.path(nd, "rgb"), pattern = "image\\.tif$", recursive = TRUE, full.names = TRUE)
cover_tile <- function(cx, cy) {
  for (t in rgb_tiles) { e <- as.vector(terra::ext(terra::rast(t)))
    if (cx >= e[1] && cx <= e[2] && cy >= e[3] && cy <= e[4]) return(t) }
  NA_character_
}
plot_chm <- function(pid) {
  clip <- file.path(nd, "frozen", SITE, pid, "native", "clip_normalized.laz")
  if (!file.exists(clip)) return(NULL)
  las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)
  tryCatch(suppressWarnings(lidR::rasterize_canopy(las, res = CHM_RES, algorithm = lidR::p2r())),
           error = function(e) NULL)
}

run_main <- function() {
  if (!length(rgb_tiles)) { cat("no RGB tiles; run neon_download_aop.R\n"); return(invisible()) }
  if (is.null(RUNNER) || !file.exists(PYTHON) || !file.exists(MODEL))
    cat(sprintf("NOTE PYTHON=%s(%s) MODEL=%s(%s) RUNNER=%s\n", PYTHON, file.exists(PYTHON),
                MODEL, file.exists(MODEL), RUNNER))
  gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  keep <- intersect(unique(gt$plotID), pc$plotID); if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  cdir <- file.path(nd, "detectree2_boxes"); dir.create(cdir, showWarnings = FALSE, recursive = TRUE)
  rows <- list(); deq <- numeric(0); fcd <- numeric(0)
  for (pid in keep) {
    ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
    if (!nrow(stems)) next
    tif <- cover_tile(cx, cy); if (is.na(tif)) next
    ocsv <- file.path(cdir, paste0(pid, ".csv"))
    if (!file.exists(ocsv)) {
      cropf <- tempfile(fileext = ".tif"); h <- ph + 20
      ok <- tryCatch({ r <- terra::crop(terra::rast(tif), terra::ext(cx - h, cx + h, cy - h, cy + h))
        terra::writeRaster(r, cropf, overwrite = TRUE); TRUE }, error = function(e) FALSE)
      if (!ok) next
      st <- tryCatch(system2(PYTHON, c(shQuote(RUNNER), shQuote(cropf), shQuote(ocsv),
              shQuote(MODEL), shQuote(file.path(tempdir(), paste0("dt2_", pid)))),
              stdout = FALSE, stderr = FALSE), error = function(e) 1L)
      unlink(cropf)
      if (!identical(as.integer(st), 0L) || !file.exists(ocsv)) { cat(sprintf("[%s] %s: runner failed\n", SITE, pid)); next }
    }
    b <- tryCatch(read.csv(ocsv, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(b) || !nrow(b)) { cat(sprintf("[%s] %s: 0 crowns\n", SITE, pid)); next }
    bp <- b[abs(b$x - cx) <= ph + TOL & abs(b$y - cy) <= ph + TOL, , drop = FALSE]
    chm <- plot_chm(pid)
    z <- if (!is.null(chm) && nrow(bp)) as.numeric(terra::extract(chm, cbind(bp$x, bp$y))[, 1]) else rep(NA_real_, nrow(bp))
    z[!is.finite(z)] <- 2.0
    det <- data.frame(x = bp$x, y = bp$y, z = z)
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy, core_half = ph),
                   error = function(e) NULL)
    if (is.null(sc)) next
    rows[[length(rows) + 1]] <- cbind(
      data.frame(site = SITE, plot = pid, plotType = ci$plotType, detector = "detectree2",
                 rung = "rgb", frdens = NA_real_, n_apex = nrow(det), stringsAsFactors = FALSE), sc)
    deq <- c(deq, bp$d_eq[is.finite(bp$d_eq)])
    cat(sprintf("[%s] %s: %d crowns in core\n", SITE, pid, nrow(bp)))
  }
  res <- rbindlist(rows, fill = TRUE)
  if (!nrow(res)) { cat("no rows scored\n"); return(invisible()) }
  o <- file.path(nd, "detectree2_results.csv"); write.csv(res, o, row.names = FALSE)
  p <- pool(res)
  cat(sprintf("\n[%s] Detectree2 (rgb): n_ref=%d recall=%.3f prec=%.3f F1=%.3f rec_und=%.3f\n",
              SITE, p$n_ref, p$recall, p$precision, p$F1, p$rec_understory))
  if (length(deq)) cat(sprintf("crown d_eq: median %.1f m, IQR [%.1f, %.1f] (n=%d polygons)\n",
              median(deq), quantile(deq, .25), quantile(deq, .75), length(deq)))
  cat(sprintf("wrote %d plot rows -> %s\n", nrow(res), o))
}

if (sys.nframe() == 0L) run_main()
