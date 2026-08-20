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

# #P3 seed->refine arm: CHM-VWF tops as prompts to SAM2Point (GitHub issue #70).
#
# The strongest overstory detector (CHM-VWF) drives a zero-shot promptable 3-D
# segmenter (SAM2Point, Apache-2.0) for crown masks, decoupling "where is the
# tree" from "what is its crown". SAM2Point voxelizes the clip, slices the voxel
# grid into a "video", and runs SAM2's video predictor from each 3-D apex prompt;
# each returned mask is a crown instance. Runs in the gpu/sam2point-sm120 Docker
# image (torch 2.7 / cu128, Blackwell sm_120 -- VERIFIED working).
#
# COST: each prompt is a full 3-axis SAM2 video segmentation (~10 s/prompt on a
# 5090), so a plot's worth of CHM-VWF apexes is ~10-35 min and a full 3-site
# ladder is many GPU-hours. This driver therefore caps prompts (MAXPROMPTS,
# core+tol apexes ranked by height) and accepts a PLOTS subset, for a bounded
# proof-of-concept; a full sweep is future GPU-time.
#
# Per (plot, rung): detect_lasr CHM-VWF apexes on the NORMALIZED clip (z = AGL) ->
# core+tol subset capped to MAXPROMPTS -> prompts.csv -> the container runner
# (run_sam2point_arm.py) voxelizes + segments each prompt -> per-point `sam2point`
# labels (AGL, so NO det_to_agl) -> reduce to apex/crowns -> score_plot vs stems.
# Reports the SEEDED arm next to the BARE CHM-VWF tops (the seeds) so the deep
# refiner's added value is isolated.
#
# Usage:
#   Rscript scripts/detect_sam2point_sweep.R SITE=SOAP PLOTS=SOAP_031,SOAP_021 \
#     MAXPROMPTS=40 IMAGE=sam2point-sm120:test
# Reads the frozen normalized clips; writes work/neon/<SITE>/sam2point_results.csv
# + persists sam2point_instances/<plot>_<rung>.laz.
suppressMessages({ library(lidR); library(data.table) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R")); source(.find("io_bridge.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (!is.null(A$SITE)) A$SITE else "SOAP"
PLOTS <- if (!is.null(A$PLOTS)) strsplit(A$PLOTS, ",")[[1]] else NULL
RUNG  <- if (is.null(A$RUNG)) "native" else A$RUNG
TOL   <- as.numeric(if (is.null(A$TOL)) 4 else A$TOL)
VWF_A <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
MAXPROMPTS <- as.integer(if (is.null(A$MAXPROMPTS)) 40 else A$MAXPROMPTS)
VOXEL <- as.numeric(if (is.null(A$VOXEL)) 0.02 else A$VOXEL)
IMAGE <- if (!is.null(A$IMAGE)) A$IMAGE else "sam2point-sm120:test"
RUNNER <- normalizePath(Find(file.exists, c("gpu/sam2point-sm120/run_sam2point_arm.py",
            file.path(getwd(), "gpu/sam2point-sm120/run_sam2point_arm.py"))), mustWork = FALSE)
nd <- file.path(d, "neon", SITE)
fz <- function(pid, f) file.path(nd, "frozen", SITE, pid, RUNG, f)

gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
keep <- intersect(unique(gt$plotID), pc$plotID)
if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
idir <- file.path(nd, "sam2point_instances"); dir.create(idir, showWarnings = FALSE, recursive = TRUE)

score_one <- function(stems, det, cx, cy, ph, label) {
  sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                            core_half = ph), error = function(e) NULL)
  if (is.null(sc)) return(NULL)
  cbind(data.frame(site = SITE, plot = NA, rung = RUNG, detector = label,
                   n_apex = nrow(det), stringsAsFactors = FALSE), sc)
}

rows <- list()
for (pid in keep) {
  ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  clip <- fz(pid, "clip_normalized.laz"); if (!file.exists(clip) || !nrow(stems)) next
  # Density-derived res/dens (CLAUDE.md: CHM resolution + the <8 pts/m^2 pre-LM
  # smoothing branch must be functions of measured first-return density, never
  # hardcoded). Read frdens from the frozen clip's manifest via the SAME helper
  # frozen_clip() writes (frozen_dir owns the layout), then derive res exactly as
  # CHM-VWF/multichm_seed_tops does and pass frdens as dens so detect_lasr's
  # smoothing branch fires on sparse rungs.
  mf <- file.path(frozen_dir(file.path(nd, "frozen"), SITE, pid, RUNG), "manifest.json")
  frdens <- tryCatch(as.numeric(jsonlite::read_json(mf, simplifyVector = TRUE)$frdens),
                     error = function(e) NA_real_)
  res <- if (!is.na(frdens) && frdens >= 8) 0.25 else 0.5  # density-derived, as CHM-VWF
  det <- tryCatch(detect_lasr(clip, res, VWF_A, frdens), error = function(e) NULL)
  if (is.null(det) || !nrow(det)) next
  # core+tol apexes, capped to MAXPROMPTS by descending height (overstory seeds)
  inreg <- abs(det$x - cx) <= ph + TOL & abs(det$y - cy) <= ph + TOL
  dseed <- det[inreg, , drop = FALSE]; dseed <- dseed[order(-dseed$z), , drop = FALSE]
  if (nrow(dseed) > MAXPROMPTS) dseed <- dseed[seq_len(MAXPROMPTS), , drop = FALSE]
  pf <- tempfile(fileext = ".csv")
  write.table(dseed[, c("x", "y", "z")], pf, sep = ",", row.names = FALSE, col.names = FALSE)
  out <- file.path(idir, sprintf("%s_%s.laz", pid, RUNG))
  cmd <- c("run", "--rm", "--gpus", "all",
           "-v", paste0(RUNNER, ":/workspace/run_arm.py:ro"),
           "-v", paste0(normalizePath(clip), ":/data/in.laz:ro"),
           "-v", paste0(pf, ":/data/prompts.csv:ro"),
           "-v", paste0(idir, ":/out"), IMAGE,
           "python3", "/workspace/run_arm.py", "--input", "/data/in.laz",
           "--prompts", "/data/prompts.csv", "--output", paste0("/out/", basename(out)),
           "--voxel_size", as.character(VOXEL))
  # shQuote the docker args (matches run_docker_arm in model_runner.R): system2
  # pastes into a shell, so an unquoted space/metachar in IMAGE or a
  # CLAUDE_JOB_DIR-derived path would break sh or inject.
  st <- tryCatch(system2("docker", shQuote(cmd), stdout = FALSE, stderr = FALSE), error = function(e) 1L)
  unlink(pf)
  # SEEDED SAM2Point: read per-point labels (AGL), reduce to apex
  if (identical(as.integer(st), 0L) && file.exists(out)) {
    sdet <- tryCatch(read_instances_laz(out, id_field = "sam2point"), error = function(e) NULL)
    if (!is.null(sdet) && nrow(sdet)) {
      r <- score_one(stems, sdet, cx, cy, ph, "sam2point_seeded"); if (!is.null(r)) { r$plot <- pid; rows[[length(rows)+1]] <- r }
    }
  }
  # BARE CHM-VWF tops (the seeds) -- the baseline the refiner must beat
  r2 <- score_one(stems, dseed[, c("x", "y", "z")], cx, cy, ph, "chm_vwf_seeds")
  if (!is.null(r2)) { r2$plot <- pid; rows[[length(rows)+1]] <- r2 }
  cat(sprintf("[%s] %s: scored (%d seeds)\n", SITE, pid, nrow(dseed)))
}
res <- rbindlist(rows, fill = TRUE)
if (nrow(res)) {
  o <- file.path(nd, "sam2point_results.csv"); write.csv(res, o, row.names = FALSE)
  cat(sprintf("\nwrote %d rows -> %s\n", nrow(res), o))
  for (dt in unique(res$detector)) { p <- pool(res[res$detector == dt, , drop = FALSE])
    cat(sprintf("  %-18s n_ref=%d recall=%.3f prec=%.3f F1=%.3f\n", dt, p$n_ref, p$recall, p$precision, p$F1)) }
} else cat("no rows scored\n")
