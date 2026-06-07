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
