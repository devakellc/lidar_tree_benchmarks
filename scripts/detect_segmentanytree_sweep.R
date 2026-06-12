#!/usr/bin/env Rscript
# SegmentAnyTree (#M6) density-ladder arm. Runs the upstream SAT instance
# segmenter zero-shot (Docker, sm_120 rebuild) on the RAW-WITH-GROUND frozen clip
# per plot x density rung, reducing per-tree instances to apex detections and
# scoring against field stems with the existing harness. CORES controls how many
# independent SAT containers may run concurrently on the same GPU; default 1,
# CORES=2 is the tested RTX 5090 throughput setting.
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
#       [TOL=4] [IMAGE=sat-sm120-test] [TIMEOUT=1800] [CORES=1] [BATCH=0]
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
BATCH   <- !is.null(A$BATCH) && A$BATCH != "0"
CORES   <- max(1L, as.integer(if (is.null(A$CORES)) 1L else A$CORES))
RUNGS    <- c(8, 4, 2, 1); MINTREES <- 6
DRIVER   <- file.path(.ROOT, "gpu", "run_segmentanytree.py")
BATCH_DRIVER <- file.path(.ROOT, "gpu", "run_segmentanytree_batch.py")
# The merged LAS carries the instance label as the PredInstance extra dim
# (merge_pt_ss_is.py: preds_instance_segmentation -> PredInstance, +1, NA->0).
ID_FIELD <- "PredInstance"

.find_las <- function(dir, stem) {
  direct <- file.path(dir, paste0(stem, c(".las", ".laz")))
  hit <- direct[file.exists(direct)]
  if (length(hit)) return(hit[1])
  esc <- gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", stem)
  hit <- list.files(dir, pattern = paste0("^", esc, "(_out)?\\.la[sz]$"),
                    full.names = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

run_sat_batch <- function(image, input_dir, output_dir, timeout, label = NULL) {
  if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  in_abs  <- normalizePath(input_dir, mustWork = TRUE)
  out_abs <- normalizePath(output_dir, mustWork = FALSE)
  mdirs <- unique(c(in_abs, out_abs, normalizePath(dirname(BATCH_DRIVER), mustWork = TRUE)))
  vol <- as.vector(rbind("-v", paste0(mdirs, ":", mdirs)))
  args <- c("run", "--rm", "--gpus", "all", "--shm-size=8g", "--ipc=host",
            vol, image, "python3", BATCH_DRIVER, in_abs, out_abs)
  out <- tryCatch(suppressWarnings(system2("docker", shQuote(args), stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")
  if (!is.null(st) && st != 0) {
    if (!is.null(label)) message(sprintf("  %s: SAT batch failed (status %s)",
                                         label, st))
    if (length(out)) message(paste(tail(out, 5), collapse = "\n"))
    return(FALSE)
  }
  TRUE
}

run_main <- function() {
  stopifnot(file.exists(DRIVER))
  if (BATCH) stopifnot(file.exists(BATCH_DRIVER))
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
  cat(sprintf("[%s] segmentanytree plots: %d (image=%s tol=%s cores=%d batch=%s)\n",
              SITE, length(keep), IMAGE, TOL, CORES, BATCH))

  result_file <- file.path(nd, "segmentanytree_results.csv")
  results <- if (file.exists(result_file)) {
    read.csv(result_file, stringsAsFactors = FALSE)
  } else data.frame()
  if (nrow(results) && !"tp_core" %in% names(results))
    results$tp_core <- round(results$precision * results$n_det)
  cell_key <- function(site, plot, rung) paste(site, plot, rung, sep = "\r")
  done <- if (nrow(results)) cell_key(results$site, results$plot, results$rung)
          else character()
  write_results <- function() {
    if (!nrow(results)) return(invisible())
    keep_rows <- !duplicated(cell_key(results$site, results$plot, results$rung),
                             fromLast = TRUE)
    write.csv(results[keep_rows, , drop = FALSE], result_file, row.names = FALSE)
  }

  # Durable per-point instance dir the #34 crown-diameter arm
  # (crown_metrics_deepmodel.R) consumes: the merged PredInstance-labelled LAS for
  # each successful (plot, rung) is persisted here as <plot>_<rung>.laz so the
  # crown arm can re-derive crown_diameter_table + instance_apex from the SAME
  # labelling without re-running the container. The merged LAS already carries the
  # PredInstance extra dim (verified by read_instances_laz above), so a plain copy
  # preserves the schema the crown arm reads.
  inst_dir <- file.path(nd, "segmentanytree_instances")
  dir.create(inst_dir, recursive = TRUE, showWarnings = FALSE)
  persist_instances <- function(src, tag) {
    if (is.null(src) || is.na(src) || !file.exists(src)) return(invisible())
    dst <- file.path(inst_dir, sprintf("%s_%s.laz", pid, tag))
    tryCatch(file.copy(src, dst, overwrite = TRUE), error = function(e) FALSE)
    invisible()
  }

  for (pid in keep) {                       # per plot; cells may run in parallel
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) next
    native_pdens <- NA_real_; ncell <- 0L; nskip <- 0L
    cells <- list()
    batch_in  <- file.path(tempdir(), sprintf("sat_%s_batch_in", pid))
    batch_out <- file.path(tempdir(), sprintf("sat_%s_batch_out", pid))
    unlink(c(batch_in, batch_out), recursive = TRUE, force = TRUE)
    dir.create(batch_in, recursive = TRUE, showWarnings = FALSE)
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                   out_root = file.path(nd, "frozen")),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) next
      tag  <- ifelse(is.na(rung), "native", as.character(rung))
      key <- cell_key(SITE, pid, tag)
      if (key %in% done) { nskip <- nskip + 1L; next }
      stem <- sprintf("sat_%s_%s", pid, tag)
      ply  <- file.path(batch_in, paste0(stem, ".ply"))
      # Absolute-UTM double PLY (lossless); SAT does its own XY localization.
      tryCatch(laz_to_ply(prep$rawground, ply), error = function(e) NULL)
      if (!file.exists(ply)) next
      cells[[length(cells) + 1L]] <- list(tag = tag, key = key, prep = prep,
        pdens = pdens, frdens = frdens, stem = stem, ply = ply)
    }

    run_one_cell <- function(cell) {
        olas <- file.path(tempdir(), paste0(cell$stem, "_single.las"))
        det_abs <- run_docker_arm(IMAGE, cell$ply, olas,
                     cmd     = c("python3", DRIVER),
                     mounts  = dirname(DRIVER),
                     # SAT's clustering Pool spawns workers; --ipc=host + a real
                     # shm avoid multiprocessing stalls under the default 64M shm.
                     extra_docker = c("--shm-size=8g", "--ipc=host"),
                     reader  = function(p) read_instances_laz(p, id_field = ID_FIELD),
                     gpus    = "all", timeout = TIMEOUT,
                     label   = sprintf("%s/%s", pid, cell$tag))
        # Persist the merged instance LAS for the #34 crown-diameter arm before the
        # tempdir is reused (run_docker_arm wrote it to `olas`, or, if it relocated
        # it, .find_las recovers it under the same tempdir).
        if (!is.null(det_abs)) {
          src <- if (file.exists(olas)) olas else .find_las(tempdir(), cell$stem)
          persist_instances(src, cell$tag)
        }
        list(cell = cell, det_abs = det_abs)
    }

    cell_results <- list()
    if (length(cells) && BATCH) {
      batch_ok <- run_sat_batch(IMAGE, batch_in, batch_out,
                                timeout = TIMEOUT * length(cells), label = pid)
      cell_results <- lapply(cells, function(cell) {
        det_abs <- NULL
        if (batch_ok) {
          las <- .find_las(batch_out, cell$stem)
          if (!is.na(las)) {
            det_abs <- read_instances_laz(las, id_field = ID_FIELD)
            # Persist the batch-merged LAS for the #34 crown arm (same schema).
            if (!is.null(det_abs)) persist_instances(las, cell$tag)
          }
        }
        if (!is.null(det_abs)) list(cell = cell, det_abs = det_abs)
        else run_one_cell(cell)
      })
    } else if (length(cells) && CORES > 1L) {
      cell_results <- parallel::mclapply(cells, run_one_cell,
        mc.cores = min(CORES, length(cells)), mc.preschedule = FALSE)
    } else if (length(cells)) {
      cell_results <- lapply(cells, run_one_cell)
    }

    for (res in cell_results) {
      cell <- res$cell
      det_abs <- res$det_abs
      if (is.null(det_abs)) next            # container crash/schema -> skip cell
      det <- det_to_agl(det_abs, cell$prep$dtm)  # absolute Z -> height above ground
      sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                core_cy = cy, core_half = ph),
                     error = function(e) NULL)
      if (is.null(sc)) next
      row <- cbind(data.frame(site = SITE, plot = pid,
        plotType = ci$plotType, detector = "segmentanytree", rung = cell$tag,
        pdens = round(cell$pdens, 2), frdens = round(cell$frdens, 2),
        n_apex = nrow(det)), sc)
      row$tp_core <- round(row$precision * row$n_det)
      results <- rbind(results, row)
      done <- c(done, cell$key)
      write_results()
      ncell <- ncell + 1L
    }
    cat(sprintf("  %s: %d new cells, %d skipped\n", pid, ncell, nskip))
  }
  if (is.null(results) || !nrow(results)) {
    cat("no segmentanytree results\n"); return(invisible())
  }
  results$tp_core <- round(results$precision * results$n_det)
  write_results()
  cat(sprintf("[%s] segmentanytree DONE: %d rows -> %s\n", SITE, nrow(results),
              result_file))
}

if (sys.nframe() == 0L) run_main()
