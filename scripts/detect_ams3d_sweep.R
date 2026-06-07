#!/usr/bin/env Rscript
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
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
# Locate scripts/ dir robustly: works both as `Rscript scripts/foo.R` (project
# root cwd) and when source()-d from tests/testthat/ (two levels up).
local({
  cf <- commandArgs(FALSE)
  ff <- sub("^--file=", "", grep("^--file=", cf, value = TRUE)[1])
  if (!is.na(ff) && nzchar(ff)) {
    .src_dir <<- dirname(normalizePath(ff, mustWork = FALSE))
  } else {
    # sourced — search candidate directories for sweep_lib.R
    # candidates cover both cwd=repo-root (Rscript) and cwd=tests/testthat (sourced under testthat)
    candidates <- c("scripts",
                    file.path("..", "..", "scripts"),
                    file.path(getwd(), "scripts"))
    found <- Filter(function(cand) file.exists(file.path(cand, "sweep_lib.R")),
                    candidates)
    .src_dir <<- if (length(found)) found[[1L]] else "scripts"
  }
})
source(file.path(.src_dir, "sweep_lib.R"))

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
  if (is.null(seg)) return(empty)
  if (!"crown_id" %in% names(seg@data)) return(empty)
  dt <- as.data.table(seg@data)[!is.na(crown_id), .(X, Y, Z, crown_id)]
  if (!nrow(dt)) return(empty)
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)),
           by = crown_id]
  data.frame(x = ap$x, y = ap$y, z = ap$z)
}
