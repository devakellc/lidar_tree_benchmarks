#!/usr/bin/env Rscript
# Library of functions for the NEON SOAP density-ladder parameter sweep.
# Sourced by run_sweep.R. Implements the docs/treetop-detection-approach.md pipeline
# (decimate -> normalize -> pit-free-ish CHM -> variable-window LM) per plot,
# per density rung, per parameter set, and scores against field stems by crown
# class.  Engine: lasR pre-devel (native ws-as-function VWF) for detection;
# lidR for clip/decimate/normalize/segmentation.
suppressMessages({ library(lidR); library(lasR); library(terra); library(sf) })

NEON_EPSG <- 32611          # SOAP = UTM 11N
PLOT_HALF <- 20             # default half-extent; overridden per plot type below
BUF       <- 25             # LiDAR clip buffer beyond the plot core (edge crowns)

# NEON woody-veg mapped extent differs by plot type: tower base plots map the
# full 40x40 m (+/-20 m), distributed plots only the 20x20 m core (+/-10 m).
# Scoring a 40x40 box over a distributed plot would count the unmapped ring as
# false commission, so the core must track the plot type. (Verified empirically:
# tower stems reach +/-20.9 m, distributed stems +/-11 m.)
plot_half <- function(plotType) ifelse(plotType == "tower", 20, 10)

## ---- variable-window allometry (Popescu & Wynne), clamped to [lo, hi] -----
ws_factory <- function(a, lo = 3, hi = 5) {
  force(a); force(lo); force(hi)
  function(h) pmin(pmax(a * h + lo, lo), hi)   # window clamped within [lo, hi]
}

## ---- global nearest-distance greedy 1:1 matcher --------------------------
# a = field stems (ax,ay), b = detections (bx,by). Builds all stem-detection
# pairs within `tol`, sorts by ascending distance, and assigns each greedily so
# the closest pairs match first (a dominant tree's own apex is claimed before a
# nearby understory stem can steal it). Returns, for each stem, the index of its
# matched detection (0 = unmatched). This is the standard ITD matching rule and
# avoids the spatial-order bias of per-stem nearest matching.
# az = stem heights, bz = detection apex heights (optional). When supplied, a
# pair is valid only if the detected apex is height-consistent with the stem:
# bz in [0.5*az, az + tol_z_up].  This stops a short understory stem from
# "stealing" the local-maximum that belongs to a tall neighbour.
greedy_match <- function(ax, ay, bx, by, tol, az = NULL, bz = NULL,
                         tol_z_up = 8) {
  m <- integer(length(ax))
  if (!length(ax) || !length(bx)) return(m)
  # candidate pairs within tol (vectorised over the typically small sets)
  pr <- do.call(rbind, lapply(seq_along(ax), function(i) {
    dd <- sqrt((bx - ax[i])^2 + (by - ay[i])^2)
    j  <- which(dd <= tol)
    if (length(j) && !is.null(az) && !is.null(bz) && !is.na(az[i])) {
      j <- j[bz[j] >= 0.5 * az[i] & bz[j] <= az[i] + tol_z_up]
    }
    if (!length(j)) return(NULL)
    cbind(stem = i, det = j, dist = dd[j])
  }))
  if (is.null(pr)) return(m)
  pr <- pr[order(pr[, "dist"]), , drop = FALSE]
  su <- logical(length(ax)); du <- logical(length(bx))
  for (k in seq_len(nrow(pr))) {
    i <- pr[k, "stem"]; j <- pr[k, "det"]
    if (!su[i] && !du[j]) { su[i] <- TRUE; du[j] <- TRUE; m[i] <- j }
  }
  m
}

## ---- one detection run on a prepared (decimated+normalized) LAS file -----
# Returns a data.frame of treetops (x,y,z) using lasR VWF on a pit-filled CHM.
detect_lasr <- function(las_file, res, a, dens, smooth_below = 8) {
  ws  <- ws_factory(a)
  del <- lasR::triangulate(filter = lasR::keep_first())
  chm <- lasR::rasterize(res, del)
  pf  <- lasR::pit_fill(chm)
  eo <- list(progress = FALSE, ncores = 1L)
  if (dens < smooth_below) {
    sm  <- lasR::focal(pf, size = 3, fun = "mean")
    lm  <- lasR::local_maximum_raster(sm, ws, min_height = 2)
    ans <- lasR::exec(del + chm + pf + sm + lm, on = las_file, with = eo)
  } else {
    lm  <- lasR::local_maximum_raster(pf, ws, min_height = 2)
    ans <- lasR::exec(del + chm + pf + lm, on = las_file, with = eo)
  }
  tt <- ans$local_maximum
  if (is.null(tt) || nrow(tt) == 0) return(data.frame(x=numeric(), y=numeric(), z=numeric()))
  xy <- sf::st_coordinates(tt)
  data.frame(x = xy[,1], y = xy[,2], z = xy[,3])
}

## ---- prepare a plot clip at a target density -----------------------------
# Clips tiles to plot AOI, decimates to `rung` (NA = native, no decimation),
# normalizes height using existing ground class, writes a temp laz. Returns
# list(file, pdens, frdens): pdens = all-return pts/m^2 (the homogenize unit, used
# for the no-upsampling guard); frdens = first-return pts/m^2 (the QL/pulse unit,
# used to gate CHM resolution & the density-tiered smoothing).
prepare_clip <- function(ctg, cx, cy, rung, tmpdir, core_half = PLOT_HALF) {
  half <- core_half + BUF
  las  <- clip_rectangle(ctg, cx - half, cy - half, cx + half, cy + half)
  if (is.empty(las) || npoints(las) < 100) return(NULL)
  if (!is.na(rung)) {
    las <- decimate_points(las, homogenize(density = rung, res = 5))
  }
  if (sum(las$Classification == 2L) < 10) return(NULL)   # need ground for DTM
  nrm <- normalize_height(las, tin(), na.rm = TRUE)
  nrm <- filter_poi(nrm, Z >= -1, Z < 80)                # drop junk
  area  <- (2 * half)^2
  pdens  <- npoints(nrm) / area                          # point density (all returns)
  frdens <- sum(nrm$ReturnNumber == 1L) / area           # first-return (pulse) density
  f <- tempfile(tmpdir = tmpdir, fileext = ".laz")
  writeLAS(nrm, f)
  list(file = f, pdens = pdens, frdens = frdens)
}

## ---- score detections vs field stems within a plot core ------------------
# stems: data.frame with E,N,crown_class,height (live trees, plot core).
# det:   data.frame x,y,z. tol_xy in m. Returns a one-row metrics data.frame
# plus per-crown-class recall.
score_plot <- function(stems, det, tol_xy = 4.0, core_cx, core_cy,
                       core_half = PLOT_HALF) {
  # Recall region: detections within the core expanded by tol_xy, so a stem at
  # the core boundary can match its apex even if the apex sits just outside.
  in_reg  <- abs(det$x - core_cx) <= core_half + tol_xy &
             abs(det$y - core_cy) <= core_half + tol_xy
  detr    <- det[in_reg, , drop = FALSE]
  is_core <- abs(detr$x - core_cx) <= core_half & abs(detr$y - core_cy) <= core_half
  m <- greedy_match(stems$E, stems$N, detr$x, detr$y, tol_xy,
                    az = stems$height, bz = detr$z)
  matched <- m > 0
  TP    <- sum(matched); n_ref <- nrow(stems)
  n_det <- sum(is_core)                       # precision denominator = core dets
  # core detections that matched a stem = true positives for precision
  tp_core <- sum(is_core[m[matched]])
  recall    <- if (n_ref)  TP / n_ref      else NA_real_
  precision <- if (n_det)  tp_core / n_det else NA_real_
  f1 <- if (!is.na(recall) && !is.na(precision) && (recall+precision) > 0)
          2*recall*precision/(recall+precision) else NA_real_
  # height RMSE on matched trees (field height vs detected apex)
  hz <- NA_real_
  if (TP > 0) {
    fh <- stems$height[matched]; dh <- detr$z[m[matched]]
    ok <- !is.na(fh) & !is.na(dh)
    if (any(ok)) hz <- sqrt(mean((dh[ok] - fh[ok])^2))
  }
  # per-crown-class recall
  cc_rec <- tapply(matched, stems$crown_class, function(v) mean(v))
  base <- data.frame(n_ref = n_ref, n_det = n_det, TP = TP,
                     recall = recall, precision = precision, F1 = f1,
                     height_rmse = hz)
  for (cls in c("dominant","codominant","intermediate","suppressed")) {
    base[[paste0("rec_", cls)]] <- if (cls %in% names(cc_rec)) cc_rec[[cls]] else NA_real_
    base[[paste0("n_",  cls)]]  <- sum(stems$crown_class == cls, na.rm = TRUE)
  }
  # height-band stratification (CHM-relevant: a surface model sees tall trees
  # regardless of social class). Bands: short <8 m, mid 8-15 m, tall >=15 m.
  hb     <- cut(stems$height, c(-Inf, 8, 15, Inf), labels = c("short","mid","tall"))
  hb_rec <- tapply(matched, hb, mean)
  for (b in c("short","mid","tall")) {
    base[[paste0("rec_h_", b)]] <- if (b %in% names(hb_rec)) hb_rec[[b]] else NA_real_
    base[[paste0("n_h_",   b)]] <- sum(hb == b, na.rm = TRUE)
  }
  base
}
