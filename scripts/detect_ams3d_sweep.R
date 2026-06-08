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

# AMS3D (adaptive mean shift, crownsegmentr) arm of the NEON model benchmark.
# Walking skeleton (#B1): runs segment_tree_crowns on the normalized clip per
# plot x density rung, collapses each crown_id to its max-Z apex (x,y,z), and
# scores against field stems with the existing sweep_lib harness.
#
# Usage:
#   Rscript scripts/detect_ams3d_sweep.R [SITE=SOAP] [PLOTS=ALL] [CORES=6]
#       [TOL=4] [CD_RATIO=0.4] [CL_RATIO=0.8]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/ams3d_results.csv (one row per
#         plot x rung) + ams3d_ledger.csv (recorded zero-shot knobs).
suppressMessages({ library(lidR); library(sf); library(data.table)
                   library(parallel); library(crownsegmentr) })
options(lidR.progress = FALSE)
d <- .job_dir()
# Locate sweep_lib.R robustly across run contexts: `Rscript scripts/foo.R`
# (cwd = repo root), `Rscript tests/run_tests.R` and `testthat::test_file`
# (cwd = tests/testthat during tests). Use the first candidate that exists.
.swp <- Find(file.exists, c(file.path("scripts", "sweep_lib.R"),
                            file.path("..", "..", "scripts", "sweep_lib.R"),
                            file.path(getwd(), "scripts", "sweep_lib.R")))
if (is.null(.swp)) stop("detect_ams3d_sweep.R: cannot locate scripts/sweep_lib.R")
source(.swp)
.mbl <- Find(file.exists, c(file.path("scripts", "model_bench_lib.R"),
                            file.path("..", "..", "scripts", "model_bench_lib.R"),
                            file.path(getwd(), "scripts", "model_bench_lib.R")))
if (is.null(.mbl)) stop("detect_ams3d_sweep.R: cannot locate scripts/model_bench_lib.R")
source(.mbl)

## ---- AMS3D apex extractor ------------------------------------------------
# las: a NORMALIZED lidR::LAS (ground at 0). Returns data.frame(x,y,z) of
# per-crown apexes (max-Z point of each crown_id). cd_ratio/cl_ratio are the
# crown-diameter/length-to-tree-height allometric ratios (zero-shot knobs).
det_ams3d <- function(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.empty(las) || npoints(las) < 5) return(empty)
  seg <- tryCatch(
    crownsegmentr::segment_tree_crowns(
      point_cloud = las,
      crown_diameter_to_tree_height = cd_ratio,
      crown_length_to_tree_height   = cl_ratio,
      segment_crowns_only_above     = min_above,
      ground_height                 = NULL,        # input already normalized
      crown_id_column_name          = "crown_id"),
    error = function(e) NULL)
  if (is.null(seg)) return(NULL)  # crash -> NULL; driver SKIPs cell. 0-row = ran-but-empty = legit recall=0.
  det <- reduce_instances(seg@data, id_col = "crown_id", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
CD    <- as.numeric(if (is.null(A$CD_RATIO)) 0.4 else A$CD_RATIO)
CL    <- as.numeric(if (is.null(A$CL_RATIO)) 0.8 else A$CL_RATIO)
RUNGS <- c(8, 4, 2, 1)
MINTREES <- 6

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
  cat(sprintf("[%s] AMS3D plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    out <- list(); native_pdens <- NA_real_
    for (rung in c(NA, RUNGS)) {
      prep <- tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph,
                                   out_root = file.path(nd, "frozen")),
                       error = function(e) NULL)
      if (is.null(prep)) next
      pdens <- prep$pdens; frdens <- prep$frdens
      if (is.na(rung)) native_pdens <- pdens
      else if (is.na(native_pdens) || rung >= native_pdens) next
      las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
      if (!is.null(las) && !is.empty(las)) {
        det <- det_ams3d(las, cd_ratio = CD, cl_ratio = CL, min_above = 2)
        if (!is.null(det)) {                           # NULL = segmenter crashed -> skip (equal-set guard drops it)
          sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                    core_cy = cy, core_half = ph),
                         error = function(e) NULL)
          if (!is.null(sc)) {
            sc <- cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                                   detector = "ams3d",
                                   rung = ifelse(is.na(rung), "native", as.character(rung)),
                                   pdens = round(pdens, 2), frdens = round(frdens, 2),
                                   n_apex = nrow(det)), sc)
            out[[length(out) + 1]] <- sc
          }
        }
      }
    }
    if (!length(out)) return(NULL)
    do.call(rbind, out)
  }

  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p), error = function(e) {
                  message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(results) || !nrow(results)) { cat("no AMS3D results\n"); return(invisible()) }
  write.csv(results, file.path(nd, "ams3d_results.csv"), row.names = FALSE)
  write.csv(data.frame(knob = c("crown_diameter_to_tree_height",
                                "crown_length_to_tree_height",
                                "segment_crowns_only_above"),
                       value = c(CD, CL, 2),
                       source = "literature default (Ferraz 2016 allometry)"),
            file.path(nd, "ams3d_ledger.csv"), row.names = FALSE)
  cat(sprintf("[%s] AMS3D DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "ams3d_results.csv")))
}

if (sys.nframe() == 0L) run_main()
