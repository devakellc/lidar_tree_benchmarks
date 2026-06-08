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

# Native-only Li 2012 arm (#R10) of the NEON model benchmark. Runs lidR's
# li2012 point-cloud segmenter on the SAME native frozen clip per plot,
# collapses per-point treeID through the bridge's reduce_instances(), and
# scores against field stems with the existing harness. Native-only by design:
# li2012 is the dense-input sub-canopy test; decimated rungs are meaningless for
# a point segmenter, and CHM-VWF/ptrees/ams3d already carry the full ladder.
#
# Usage:  Rscript scripts/detect_li2012_native.R [SITE=SOAP] [PLOTS=ALL]
#             [CORES=6] [TOL=4]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/li2012_results.csv (one row per plot,
#         detector "li2012", rung "native").
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- .job_dir()
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))

# li2012 on a normalized clip -> per-point treeID -> reduce_instances apex set.
# Guard the no-canopy case (parity with ptrees): too few returns above hmin
# means no trees anyway, so return a 0-row frame without entering the segmenter.
det_li2012 <- function(las, hmin = 2, dt1 = 1.5, dt2 = 2, R = 2) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (sum(las$Z >= hmin) < 1) { assert_detection_contract(empty); return(empty) }
  seg <- tryCatch(lidR::segment_trees(las,
                    lidR::li2012(dt1 = dt1, dt2 = dt2, R = R, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg)) return(NULL)               # crash -> skip (guard drops cell)
  if (!"treeID" %in% names(seg@data)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
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
  cat(sprintf("[%s] li2012 plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))

  run_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]
    cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) return(NULL)
    prep <- tryCatch(frozen_clip(ctg, SITE, pid, NA, cx, cy, ph,
                                 out_root = file.path(nd, "frozen")),
                     error = function(e) NULL)
    if (is.null(prep)) return(NULL)
    las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
    if (is.null(las) || is.empty(las)) return(NULL)
    det <- det_li2012(las, hmin = 2)
    if (is.null(det)) return(NULL)
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                              core_cy = cy, core_half = ph),
                   error = function(e) NULL)
    if (is.null(sc)) return(NULL)
    cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                     detector = "li2012", rung = "native",
                     pdens = round(prep$pdens, 2), frdens = round(prep$frdens, 2),
                     n_apex = nrow(det)), sc)
  }

  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p), error = function(e) {
                  message("plot ", p, " failed: ", conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(results) || !nrow(results)) { cat("no li2012 results\n"); return(invisible()) }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "li2012_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] li2012 DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "li2012_results.csv")))
}

if (sys.nframe() == 0L) run_main()
