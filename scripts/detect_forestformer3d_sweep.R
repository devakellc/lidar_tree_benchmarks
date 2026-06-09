#!/usr/bin/env Rscript
# ForestFormer3D (#M8) density-ladder arm. Tiles each plot core into 16 m-radius
# cylinders, runs FF3D zero-shot ONCE per plot over all cylinders in the
# ff3d-sm120 container (RAW-WITH-GROUND frozen clip), cross-block dedups the
# instances, reduces to apex detections, and scores against field stems. SOAP,
# native + 8 only (the heaviest arm). Serial -- one GPU.
#
#   rawground.laz --clip_circle(16) x N--> in_dir/cyl_*.laz
#     --run_docker_arm(ff3d-sm120, ff3d_entry.sh)--> merged.laz
#     --ff3d_collapse(dedup_blocks + reduce)--> apex(x,y,z UTM)
#     --agl_guard(ground_dtm.tif)--> apex(z AGL) --score_plot.
#
# Usage:
#   Rscript scripts/detect_forestformer3d_sweep.R [SITE=SOAP] [PLOTS=ALL]
#     [SPACING=24] [MERGE_TOL=2.0] [TOL=4] [IMAGE=ff3d-sm120]
#     [REPO=<abs>] [CKPT=<abs>] [TIMEOUT=3600]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/forestformer3d_results.csv (row per plot x rung).
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

args     <- strsplit(commandArgs(TRUE), "=")
A        <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE     <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS    <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
SPACING  <- as.numeric(if (is.null(A$SPACING)) 24 else A$SPACING)
MERGE_TOL<- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
TOL      <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
IMAGE    <- if (is.null(A$IMAGE)) "ff3d-sm120" else A$IMAGE
TIMEOUT  <- as.numeric(if (is.null(A$TIMEOUT)) 3600 else A$TIMEOUT)
RADIUS   <- 16; RUNGS <- c(8); MINTREES <- 6
REPO  <- if (is.null(A$REPO)) file.path(.ROOT, "gpu/store/forestformer3d/ForestFormer3D") else A$REPO
CKPT  <- if (is.null(A$CKPT)) file.path(REPO, "work_dirs/clean_forestformer/epoch_3000_fix.pth") else A$CKPT
ENTRY <- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_entry.sh")
PATCH <- file.path(REPO, "ff3d_repo.patch")          # the #27 patch, shipped in the repo
DRIVER<- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_arm.py")

# Cylinder centers on a square grid covering [-ph, ph]^2 about the plot center,
# spacing SPACING; each cylinder processes radius RADIUS (overlap = 2*RADIUS-SPACING).
cyl_centers <- function(cx, cy, ph, spacing) {
  k <- max(1L, ceiling((2 * ph) / spacing) + 1L)
  off <- seq(-ph, ph, length.out = k)
  g <- expand.grid(dx = off, dy = off)
  data.frame(cx = cx + g$dx, cy = cy + g$dy)
}

run_main <- function() {
  stopifnot(file.exists(ENTRY), file.exists(DRIVER), file.exists(CKPT), file.exists(REPO))
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),     stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  counts <- table(gt$plotID); keep <- names(counts)[counts >= MINTREES]
  if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] forestformer3d plots: %d (image=%s spacing=%g merge_tol=%g)\n",
              SITE, length(keep), IMAGE, SPACING, MERGE_TOL))

  out <- list()
  for (pid in keep) {                       # SERIAL -- one GPU
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
      tag <- ifelse(is.na(rung), "native", as.character(rung))
      # tile the raw clip into cylinders
      raw <- tryCatch(lidR::readLAS(prep$rawground), error = function(e) NULL)
      if (is.null(raw) || lidR::is.empty(raw)) next
      in_dir <- file.path(tempdir(), sprintf("ff3d_%s_%s", pid, tag))
      unlink(in_dir, recursive = TRUE); dir.create(in_dir, recursive = TRUE)
      cc <- cyl_centers(cx, cy, ph, SPACING); n_cyl <- 0L
      for (i in seq_len(nrow(cc))) {
        cyl <- lidR::clip_circle(raw, cc$cx[i], cc$cy[i], RADIUS)
        if (lidR::is.empty(cyl) || lidR::npoints(cyl) < 50) next
        lidR::writeLAS(cyl, file.path(in_dir, sprintf("cyl_%03d.laz", n_cyl)))
        n_cyl <- n_cyl + 1L
      }
      if (n_cyl == 0L) next
      out_laz <- file.path(tempdir(), sprintf("ff3d_%s_%s.laz", pid, tag))
      # run_docker_arm passes no env, so REPO/PATCH/DRIVER ride in `extra` as
      # positional args (entry.sh reads $3..$6); all are identity-mounted.
      det_abs <- run_docker_arm(IMAGE, in_dir, out_laz,
                   cmd    = c("bash", ENTRY),
                   extra  = c(CKPT, REPO, PATCH, DRIVER),
                   mounts = c(REPO, dirname(CKPT), dirname(ENTRY)),
                   reader = function(p) ff3d_collapse(p, merge_tol = MERGE_TOL),
                   gpus = "all", timeout = TIMEOUT,
                   label = sprintf("%s/%s", pid, tag))
      if (is.null(det_abs)) next            # container crash/schema -> skip cell
      det <- agl_guard(det_abs, prep$dtm)
      if (is.null(det)) next                # wholesale off-DTM (frame bug) -> skip
      sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                core_cy = cy, core_half = ph),
                     error = function(e) NULL)
      if (is.null(sc)) next
      out[[length(out) + 1]] <- cbind(data.frame(site = SITE, plot = pid,
        plotType = ci$plotType, detector = "forestformer3d", rung = tag,
        pdens = round(pdens, 2), frdens = round(frdens, 2),
        n_cyl = n_cyl, n_apex = nrow(det)), sc)
      ncell <- ncell + 1L
    }
    cat(sprintf("  %s: %d cells\n", pid, ncell))
  }
  results <- do.call(rbind, out)
  if (is.null(results) || !nrow(results)) {
    cat("no forestformer3d results\n"); return(invisible())
  }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "forestformer3d_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] forestformer3d DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "forestformer3d_results.csv")))
}

if (sys.nframe() == 0L) run_main()
