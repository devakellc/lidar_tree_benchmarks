#!/usr/bin/env Rscript
# #I4 I/O bridge. instances_to_det: labeled point table/LAS (instance-id field;
# 0 OR NA = unassigned) -> apex (x,y,z) via the bridge's reduce_instances
# (max-Z per id). det_to_agl: apex absolute elevation -> height above ground via
# the frozen clip's ground_dtm.tif, so z matches score_plot's height gate; it
# drops apexes that fall off the DTM and records the count in attr "n_dropped".
# (LAZ->PLY for the SAT/FF3D Docker arms is added in their issue, not here.)
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))
suppressMessages(library(terra))

instances_to_det <- function(pts, id_field = "pred_itc",
                             x = "X", y = "Y", z = "Z") {
  pts <- as.data.frame(pts)
  if (!id_field %in% names(pts))
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  ids <- pts[[id_field]]
  pts[[id_field]][!is.na(ids) & ids == 0] <- NA      # 0 or NA = unassigned -> dropped
  det <- reduce_instances(pts, id_col = id_field, x = x, y = y, z = z)
  assert_detection_contract(det)
  det
}

read_instances_laz <- function(path, id_field = "pred_itc") {
  las <- lidR::readLAS(path)
  if (is.null(las) || lidR::is.empty(las))
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  instances_to_det(las@data, id_field = id_field)
}

# Subtract ground elevation at each apex's (x,y). Apexes off the DTM are dropped;
# attr(.,"n_dropped") carries the count so callers can refuse a wholesale drop.
det_to_agl <- function(det, dtm_path) {
  if (!nrow(det)) { attr(det, "n_dropped") <- 0L; return(det) }
  g <- terra::extract(terra::rast(dtm_path), cbind(det$x, det$y))[, 1]
  ok <- !is.na(g)
  out <- det[ok, , drop = FALSE]; out$z <- out$z - g[ok]
  rownames(out) <- NULL
  attr(out, "n_dropped") <- as.integer(sum(!ok))
  out
}
