#!/usr/bin/env Rscript
# Shared point-cloud-detector library for the issue-#6 native comparison
# (detect_pc_sweep.R) and the issue-#38 density ladder (detect_pc_ladder.R).
#
# Holds the three POINT-CLOUD apex extractors (lidR lmf-on-points, lidR Li 2012
# segmentation, lasR point local_maximum) plus the no-upsampling rung guard. The
# CHM-VWF baseline arm is detect_lasr() from sweep_lib.R; it is not duplicated
# here. Every extractor returns a base data.frame with exactly columns x, y, z
# (the score_plot / assert_detection_contract contract) and a 0-row frame -- not
# NULL -- on empty/no-canopy input, so a detector that finds nothing is scored as
# zero recall rather than silently dropped.
suppressMessages({ library(lidR); library(lasR); library(sf)
                   library(data.table) })

.pc_empty <- function() data.frame(x = numeric(), y = numeric(), z = numeric())

## ---- which decimated rungs are runnable without upsampling ---------------
# Mirrors run_sweep.R's no-upsampling guard (`rung >= native_pdens -> next`):
# decimating to a density at or above the plot's native all-return density would
# either upsample (impossible) or alias the rung onto native, so such rungs are
# dropped. Returns the subset of `decim_rungs` strictly below native_pdens; an
# unknown (NA) native density drops all decimated rungs.
pc_rungs_for <- function(native_pdens, decim_rungs = 8) {
  if (length(native_pdens) != 1L || is.na(native_pdens)) return(numeric(0))
  decim_rungs[decim_rungs < native_pdens]
}

## ---- (a) lidR lmf on the POINT CLOUD -------------------------------------
det_lidr_lmf_pc <- function(las, ws) {
  if (lidR::is.empty(las)) return(.pc_empty())
  tt <- lidR::locate_trees(las, lidR::lmf(ws = ws, hmin = 2,
                                          shape = "circular"))
  if (is.null(tt) || nrow(tt) == 0) return(.pc_empty())
  xy <- sf::st_coordinates(tt)
  data.frame(x = xy[, 1], y = xy[, 2], z = xy[, 3])
}

## ---- (b) lidR Li 2012 point-cloud segmentation; apex = max-Z per treeID ---
det_lidr_li2012 <- function(las) {
  if (lidR::is.empty(las)) return(.pc_empty())
  seg <- lidR::segment_trees(las, lidR::li2012(dt1 = 1.5, dt2 = 2, R = 2,
                                               Zu = 15, hmin = 2,
                                               speed_up = 10))
  dt  <- seg@data[!is.na(seg@data$treeID),
                  c("X", "Y", "Z", "treeID"), with = FALSE]
  if (!nrow(dt)) return(.pc_empty())
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)),
           by = treeID]
  data.frame(x = ap$x, y = ap$y, z = ap$z)
}

## ---- (c) lasR POINT-BASED local_maximum (not local_maximum_raster) -------
det_lasr_lmax_pc <- function(las_file, ws) {
  ans <- lasR::exec(lasR::local_maximum(ws, min_height = 2), on = las_file,
                    with = list(progress = FALSE, ncores = 1L))
  g <- if (inherits(ans, "sf")) ans else ans$local_maximum
  if (is.null(g) || nrow(g) == 0) return(.pc_empty())
  co <- sf::st_coordinates(g)
  data.frame(x = co[, "X"], y = co[, "Y"], z = co[, "Z"])
}
