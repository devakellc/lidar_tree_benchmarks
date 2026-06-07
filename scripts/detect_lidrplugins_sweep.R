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
