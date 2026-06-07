#!/usr/bin/env Rscript
# lidRplugins competitor arm (#C9) of the NEON model benchmark.
# Runs three classical LiDAR-native detectors (lmfauto + multichm via
# locate_trees; ptrees via segment_trees + reduce_instances) plus the CHM-VWF
# baseline on the SAME frozen normalized clip per plot x density rung, scoring
# each against field stems with the existing harness.
# lmfauto/multichm return treetops directly -> st_coordinates + .tops_to_det;
# ptrees returns per-point treeID -> reduce_instances. All -> assert_detection_contract.
#
# Usage:
#   Rscript scripts/detect_lidrplugins_sweep.R [SITE=SOAP] [PLOTS=ALL]
#       [CORES=6] [TOL=4] [A=0.10]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/lidrplugins_results.csv (one row per
#         plot x rung x detector).
suppressMessages({ library(lidR); library(lidRplugins); library(sf)
                   library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
# Locate the two libs robustly across run contexts (Rscript from root vs
# source()-d under testthat where cwd = tests/testthat).
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))

## ---- tops sf -> (x,y,z) detection contract -------------------------------
# Pull lowercase x,y,z from a locate_trees() sf. Prefer 3D geometry; fall back
# to the Z attribute column if the geometry is 2D (multichm). assert_detection_contract
# enforces the exact shape score_plot consumes.
.tops_to_det <- function(tt) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.null(tt) || nrow(tt) == 0) { assert_detection_contract(empty); return(empty) }
  co <- sf::st_coordinates(tt)
  z  <- if (ncol(co) >= 3) co[, 3] else if (!is.null(tt$Z)) tt$Z else co[, 2] * 0
  det <- data.frame(x = co[, 1], y = co[, 2], z = z)
  assert_detection_contract(det)
  det
}

# lmfauto: parameter-free variable-window LMF (plot mode for plot-scale clips).
det_lmfauto <- function(las, hmin = 2) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::lmfauto(plot = TRUE, hmin = hmin)),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)            # crash -> skip (equal-set guard drops it)
  .tops_to_det(tt)
}

# multichm: multi-layer CHM; pass our clamped variable window (same allometry as
# CHM-VWF) through ... to the internal lmf for a fair comparison.
det_multichm <- function(las, res = 0.5, a = 0.10) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::multichm(res = res,
                                                               ws = ws_factory(a))),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  .tops_to_det(tt)
}

# ptrees: point-based PTrees (Vega 2014). Its locate_trees path is broken under
# lidR 4.3.2 ("object not a >= 2-column array"), but segment_trees works and
# returns per-point treeID -- an instance segmentation -- which collapses through
# the bridge's reduce_instances() exactly like the AMS3D crown_id arm.
det_ptrees <- function(las, hmin = 2, k = c(30, 15)) {
  seg <- tryCatch(lidR::segment_trees(las, lidRplugins::ptrees(k = k, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg)) return(NULL)            # crash -> skip (equal-set guard drops it)
  if (!"treeID" %in% names(seg@data)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
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
A_VWF <- as.numeric(if (is.null(A$A))   0.10 else A$A)
RUNGS <- c(8, 4, 2, 1)
MINTREES <- 6
ARMS  <- c("lmfauto", "multichm", "ptrees", "chm_vwf")

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
  cat(sprintf("[%s] lidRplugins plots: %d (%s)\n", SITE, length(keep),
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
      if (is.null(las) || is.empty(las)) next
      res <- if (frdens >= 8) 0.25 else 0.5     # density-derived, like CHM-VWF
      dets <- list(
        lmfauto  = det_lmfauto(las, hmin = 2),
        multichm = det_multichm(las, res = res, a = A_VWF),
        ptrees   = det_ptrees(las, hmin = 2),
        chm_vwf  = tryCatch(detect_lasr(prep$normalized, res, A_VWF, frdens),
                            error = function(e) NULL))
      for (nm in names(dets)) {
        det <- dets[[nm]]
        if (is.null(det)) next                  # crash -> skip this detector/cell
        sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                                  core_cy = cy, core_half = ph),
                       error = function(e) NULL)
        if (is.null(sc)) next
        sc <- cbind(data.frame(site = SITE, plot = pid, plotType = ci$plotType,
                               detector = nm,
                               rung = ifelse(is.na(rung), "native", as.character(rung)),
                               pdens = round(pdens, 2), frdens = round(frdens, 2),
                               chm_res = if (nm %in% c("multichm","chm_vwf")) res else NA_real_,
                               n_apex = nrow(det)), sc)
        out[[length(out) + 1]] <- sc
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
  if (is.null(results) || !nrow(results)) { cat("no lidRplugins results\n"); return(invisible()) }
  results$tp_core <- round(results$precision * results$n_det)
  write.csv(results, file.path(nd, "lidrplugins_results.csv"), row.names = FALSE)
  cat(sprintf("[%s] lidRplugins DONE: %d rows -> %s\n", SITE, nrow(results),
              file.path(nd, "lidrplugins_results.csv")))

  guarded <- equal_set_guard(results, arms = ARMS)
  if (length(attr(guarded, "dropped")))
    cat(sprintf("equal-set guard dropped %d (plot,rung) cells\n",
                length(attr(guarded, "dropped"))))
  cat("\n=== Pooled recall/precision/F1 by detector (common cells) ===\n")
  pooled <- do.call(rbind, lapply(ARMS, function(a) {
    s <- guarded[guarded$detector == a, ]; if (!nrow(s)) return(NULL)
    cbind(detector = a, pool(s)) }))
  print(pooled[, c("detector","n_plots","n_ref","recall","precision","F1",
                   "rec_dominant","rec_understory")], row.names = FALSE, digits = 3)
}

if (sys.nframe() == 0L) run_main()
