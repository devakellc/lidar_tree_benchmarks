#!/usr/bin/env Rscript
# Shared bridge for the NEON model benchmark (#B2). Sourced by every model arm.
# Provides: reduce_instances (universal predictions->detections collapse),
# crown_diameter_table (optional side metric), seed_for + frozen_clip
# (deterministic two-variant clip provider), pool + equal_set_guard (canonical
# pooling), assert_detection_contract (scorer-contract conformance harness).
suppressMessages({ library(data.table); library(lidR); library(terra)
                   library(jsonlite) })

## ---- universal reducer: labelled points -> data.frame(x,y,z) -------------
# pts: data.frame/data.table with coordinate columns (default X,Y,Z) and an
# instance-id column (default "crown_id"; NA = unassigned, dropped). Apex per
# instance = its max-Z point. Returns a base data.frame with exactly x,y,z.
reduce_instances <- function(pts, id_col = "crown_id",
                             x = "X", y = "Y", z = "Z") {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  dt <- as.data.table(pts)
  if (!nrow(dt) || !id_col %in% names(dt)) return(empty)
  dt <- dt[!is.na(dt[[id_col]]), c(x, y, z, id_col), with = FALSE]
  if (!nrow(dt)) return(empty)
  setnames(dt, c(x, y, z, id_col), c("X", "Y", "Z", "ID"))
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)), by = ID]
  data.frame(x = ap$x, y = ap$y, z = ap$z)
}

## ---- optional crown-diameter side metric (NOT in the detection contract) -
# Per instance: d_eq = 2*sqrt(area/pi) from the 2-D convex hull area;
# d_caliper = max pairwise distance among the instance's points. Both NA when
# the instance has fewer than `min_pts` points (a diameter over a few points is
# noise). Returns data.frame(id, n_pts, d_eq, d_caliper).
crown_diameter_table <- function(pts, id_col = "crown_id", min_pts = 5) {
  dt <- as.data.table(pts)
  dt <- dt[!is.na(dt[[id_col]])]
  if (!nrow(dt)) return(data.frame(id = integer(), n_pts = integer(),
                                   d_eq = numeric(), d_caliper = numeric()))
  setnames(dt, id_col, "ID")
  one <- function(s) {
    n <- nrow(s)
    if (n < min_pts) return(list(n_pts = n, d_eq = NA_real_, d_caliper = NA_real_))
    h  <- grDevices::chull(s$X, s$Y)
    hx <- s$X[h]; hy <- s$Y[h]
    area <- abs(sum(hx * c(hy[-1], hy[1]) - c(hx[-1], hx[1]) * hy)) / 2
    cal  <- max(dist(cbind(s$X, s$Y)))
    list(n_pts = n, d_eq = 2 * sqrt(area / pi), d_caliper = cal)
  }
  res <- dt[, one(.SD), by = ID, .SDcols = c("X", "Y")]
  data.frame(id = res$ID, n_pts = res$n_pts, d_eq = res$d_eq,
             d_caliper = res$d_caliper)
}
