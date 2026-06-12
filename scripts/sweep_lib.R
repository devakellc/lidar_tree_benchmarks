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

## ---- multichm treetop seeds for the crown segmenters (issue #31) ----------
# The crown benchmark (crown_metrics_sweep.R) normally seeds every segmenter
# from one CHM-VWF lmf top set. The model benchmark shows `multichm` is the best
# classical *detector* on SOAP (Eysn-style multi-layer CHM local maxima), so this
# helper produces an ALTERNATIVE seed set the lidR segmenters (dalponte2016,
# silva2016) can consume in place of the lmf tops, to test whether better tops
# improve crown-width RMSE (detection and segmentation are decoupled).
#
# multichm runs on the POINT CLOUD `las` (mirrors detect_multichm_sweep.R::
# det_multichm_run) at a density-derived res and the SAME clamped variable window
# ws_factory(a) as the lmf seeds, so the only thing that changes vs the control
# is the detector. multichm geometry is 2-D, so the apex height each segmenter
# needs is read from the SAME pit-free `chm` at the tops' (x, y) (the segmenters
# grow on that CHM, so a CHM-consistent Z is the right seed height; the sf Z, if
# present, is used only as a fallback when a top lands off the CHM extent). NOTE:
# this is NOT the lasR seed path -- lasR region_growing takes a seed *stage*, not
# an external point set, and there is no matching lasR local_maximum that emits
# the multichm tops, so the lasR arm cannot be re-seeded from multichm and stays
# lmf-seeded as the control (see crown_metrics_sweep.R).
#
# Returns an sf POINT object with a `treeID` column and a Z coordinate (the shape
# dalponte2016/silva2016 accept as `treetops`), or NULL if multichm yields no
# tops. Off-CHM tops (NA height with no usable sf Z) are dropped so a segmenter
# never receives a non-finite seed height.
multichm_seed_tops <- function(las, chm, a = 0.10, frdens = NA_real_) {
  res <- if (!is.na(frdens) && frdens >= 8) 0.25 else 0.5  # density-derived, as CHM-VWF
  tt  <- tryCatch(
    lidR::locate_trees(las, lidRplugins::multichm(res = res, ws = ws_factory(a))),
    error = function(e) NULL)
  if (is.null(tt) || nrow(tt) == 0) return(NULL)
  co <- sf::st_coordinates(tt)
  x  <- co[, 1]; y <- co[, 2]
  z_sf <- if (ncol(co) >= 3) co[, 3] else if (!is.null(tt$Z)) tt$Z else rep(NA_real_, length(x))
  z_chm <- tryCatch(as.numeric(terra::extract(chm, cbind(x, y))[, 1]),
                    error = function(e) rep(NA_real_, length(x)))
  z <- ifelse(is.finite(z_chm), z_chm, z_sf)   # CHM-consistent height; sf Z fallback
  keep <- is.finite(x) & is.finite(y) & is.finite(z)
  if (!any(keep)) return(NULL)
  sf::st_as_sf(data.frame(treeID = seq_len(sum(keep)),
                          X = x[keep], Y = y[keep], Z = z[keep]),
               coords = c("X", "Y", "Z"), crs = sf::st_crs(tt))
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

## ---- th_cr-style stop rule for the random-walker crown labels -------------
# The Grady (2006) random walker on a CHM is a *label-everything* partitioner:
# every canopy pixel is argmax-assigned to some seed, so crowns tile the whole
# canopy mask with no "stop growing" criterion (issue #7's honest negative).
# This applies a per-crown height cutoff analogous to lidR dalponte2016's `th_cr`
# (default 0.55): a pixel assigned to crown k is dropped to background when its
# CHM height falls below `th_cr * seed_height_k` — i.e. the crown is truncated at
# a fraction of its own seed's apex height, so it can no longer creep down into
# the inter-crown gaps that inflate diameter.
#
# Pure relabel: takes the argmax label vector `out_lab` (0 = background, k = crown
# k, CHM-cell order), the per-pixel CHM heights `chm_vals` (same length/order),
# and a numeric vector `lab_height` mapping label index k -> its seed apex height
# (lab_height[k]); returns a relabelled integer vector of the same length with
# below-cutoff pixels set to 0. Labels with no seed height (NA, or index beyond
# `lab_height`) are left untouched (no cutoff applied), and pixels already at
# background stay background.
apply_thcr_cutoff <- function(out_lab, chm_vals, lab_height, th_cr = 0.55) {
  stopifnot(length(out_lab) == length(chm_vals))
  out <- as.integer(out_lab)
  fg  <- which(out > 0L)
  if (!length(fg)) return(out)
  k   <- out[fg]
  hk  <- ifelse(k <= length(lab_height), lab_height[k], NA_real_)
  cut <- th_cr * hk
  # drop only where a finite cutoff exists AND the pixel sits below it
  drop <- is.finite(cut) & is.finite(chm_vals[fg]) & chm_vals[fg] < cut
  out[fg[drop]] <- 0L
  out
}

## ---- marker-controlled (seeded) watershed on a CHM (issue #32) ------------
# The marker-FREE watershed in #7 (lidR/EBImage, th_tree=2) finds its own basins
# from CHM minima and systematically over-grows crowns (pooled bias +2.31 m on
# d_eq). The approach doc (sec 6) instead recommends a MARKER-controlled
# watershed: each detected treetop is a basin marker, so exactly one crown grows
# per seed -- no spurious extra basins, and a crown is keyed by the treeID of its
# seed (not matched back by containment).
#
# seeds_to_marker_raster() rasterizes a seed point set into an integer label
# image aligned to `chm`: cell of seed i (terra::cellFromXY) carries label
# treeID[i]; all other cells are 0. When two seeds fall in one cell the first
# seed wins (mirrors the random_walker() seed-collision rule). NA-cell seeds
# (off the CHM extent) are dropped. Returns an integer vector in CHM cell order
# (length == ncell(chm)); 0 == no marker.
seeds_to_marker_raster <- function(chm, seed_xy, treeID) {
  stopifnot(nrow(seed_xy) == length(treeID))
  n   <- terra::ncell(chm)
  lab <- integer(n)                       # 0 = no marker
  if (!length(treeID)) return(lab)
  cells <- terra::cellFromXY(chm, seed_xy)
  for (i in seq_along(cells)) {
    ci <- cells[i]
    if (is.na(ci)) next                   # seed off the CHM extent
    if (lab[ci] == 0L) lab[ci] <- as.integer(treeID[i])   # first seed wins
  }
  lab
}

## Priority-flood seeded watershed. imager::watershed(img, seeds) does exactly
## this priority-flood from a labelled seed image, but imager is not installed in
## this environment (verified at runtime in crown_metrics_sweep.R), so we run the
## flood directly on the CHM. Each marker cell starts a basin; the canopy mask is
## flooded in order of DESCENDING CHM height (a frontier always grows downhill
## from ridges, the watershed convention), and each unlabelled canopy pixel is
## claimed by the FIRST basin to reach it. A canopy pixel that no marker can
## reach (a seed-less island, disconnected from every marker over the mask) stays
## background -- like random_walker(), "floodable" does not mean "force-labelled".
##
## Inputs: `chm` SpatRaster, `markers` integer vector in CHM cell order
## (0 = none, k = label, e.g. from seeds_to_marker_raster), `hmin` canopy
## threshold (m). Returns an integer label vector in CHM cell order (0 =
## background) -- caller writes it onto a copy of `chm` to polygonize. Pure
## (no I/O); the 4-neighbour frontier uses base R only.
priority_flood_watershed <- function(chm, markers, hmin = 2) {
  nr <- nrow(chm); nc <- ncol(chm); n <- nr * nc
  stopifnot(length(markers) == n)
  vals  <- terra::values(chm, mat = FALSE)
  canopy <- is.finite(vals) & vals >= hmin
  lab <- integer(n)
  seeds <- which(markers != 0L & canopy)     # only markers on the canopy mask
  if (!length(seeds)) return(lab)
  lab[seeds] <- markers[seeds]
  # frontier = unlabelled canopy neighbours of the labelled set, processed
  # highest-CHM-first so basins meet along the canopy valleys between crowns.
  # rowcol arithmetic on the linear cell index (terra: cell = (row-1)*nc + col).
  row_of <- function(cell) ((cell - 1L) %/% nc) + 1L
  col_of <- function(cell) ((cell - 1L) %%  nc) + 1L
  neighbours <- function(cell) {
    rr <- row_of(cell); cc <- col_of(cell)
    out <- integer(0)
    if (rr > 1L)  out <- c(out, cell - nc)   # up
    if (rr < nr)  out <- c(out, cell + nc)   # down
    if (cc > 1L)  out <- c(out, cell - 1L)   # left
    if (cc < nc)  out <- c(out, cell + 1L)   # right
    out
  }
  # initial frontier: unlabelled canopy neighbours of any seed cell
  front <- unique(unlist(lapply(seeds, neighbours)))
  front <- front[canopy[front] & lab[front] == 0L]
  enq   <- logical(n); enq[front] <- TRUE     # in-frontier guard (no dup work)
  while (length(front)) {
    # pop the highest-CHM frontier cell (descending flood from the ridges)
    j  <- which.max(vals[front])
    cell <- front[j]; front <- front[-j]; enq[cell] <- FALSE
    if (lab[cell] != 0L) next                 # already claimed (defensive)
    # claim from the highest labelled neighbour (steepest-ascent basin)
    nb  <- neighbours(cell)
    lnb <- nb[lab[nb] != 0L]
    if (length(lnb)) {
      lab[cell] <- lab[lnb[which.max(vals[lnb])]]
      # enqueue this cell's still-unlabelled canopy neighbours
      add <- nb[canopy[nb] & lab[nb] == 0L & !enq[nb]]
      if (length(add)) { front <- c(front, add); enq[add] <- TRUE }
    }
  }
  lab
}

## ---- crown-diameter error stats (pooled SSE, never a mean of per-plot) -----
# Pooled error stats over matched trees: RMSE/MAE/bias/R2 computed from the
# SUMMED squared errors over all rows (n_total = sum of matched trees), NEVER a
# mean of per-plot rates -- the crown-benchmark analogue of the detection
# sum(TP)/sum(n_ref) rule. det = detected diameter, fld = field diameter; rows
# with a non-finite det or fld are dropped before pooling. Mirrors
# crown_metrics_sweep.R::err_stats / analyze_crown_metrics.R::err_stats (kept as
# the single canonical definition so the per-rung pooler below cannot drift from
# the script's own scoring).
crown_err_stats <- function(det, fld) {
  ok <- is.finite(det) & is.finite(fld); det <- det[ok]; fld <- fld[ok]
  n <- length(det)
  if (n < 2) return(data.frame(n = n, rmse = NA_real_, mae = NA_real_,
                               bias = NA_real_, r2 = NA_real_))
  e <- det - fld; ss_res <- sum((fld - det)^2); ss_tot <- sum((fld - mean(fld))^2)
  data.frame(n = n, rmse = sqrt(mean(e^2)), mae = mean(abs(e)), bias = mean(e),
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_)
}

## ---- per-(algo, rung) crown-diameter pooler (issue #33) -------------------
# Crown-diameter density ladder: pool the matched-tree crown rows WITHIN each
# (algo, rung) cell by the summed-squared-error rule (crown_err_stats above), so
# a rung's RMSE/bias is over every matched tree at that rung -- never a mean of
# per-plot RMSEs (a small plot would otherwise dominate, exactly as the detection
# sweep forbids mean-of-rates). One pooled row per (algo, rung, definition).
#
# BACKWARD COMPATIBILITY: a crown_metrics_results.csv written before this issue
# has no `rung` column; such rows are the native-density #7 run, so a missing /
# NA / "" rung defaults to "native". This lets the analyzer pool an old CSV (one
# native rung) and a new ladder CSV identically.
#
# `defs` is a list of c(det_col, fld_col, label) triples (default the two
# canonical diameter pairs). Returns a long data.frame:
#   algo, rung, definition, n, rmse, mae, bias, r2
# ordered by (definition, rung-order native->sparse, algo). Pure: no I/O.
RUNG_LEVELS <- c("native", "8", "4", "2", "1")
normalize_rung <- function(rung) {
  r <- as.character(rung)
  r[is.na(r) | !nzchar(r)] <- "native"
  r
}
pool_crown_by_rung <- function(df, defs = list(
                                 c("d_eq", "field_ninetyCD",
                                   "d_eq vs ninetyCrownDiameter"),
                                 c("d_caliper", "field_maxCD",
                                   "d_caliper vs maxCrownDiameter"))) {
  if (is.null(df) || !nrow(df)) {
    return(data.frame(algo = character(), rung = character(),
                      definition = character(), n = integer(),
                      rmse = numeric(), mae = numeric(), bias = numeric(),
                      r2 = numeric(), stringsAsFactors = FALSE))
  }
  rung <- if (is.null(df$rung)) rep("native", nrow(df)) else normalize_rung(df$rung)
  df$rung <- rung
  out <- list()
  for (defn in defs) {
    det_col <- defn[1]; fld_col <- defn[2]; lbl <- defn[3]
    for (a in sort(unique(df$algo))) for (rg in unique(df$rung)) {
      sel <- df$algo == a & df$rung == rg
      if (!any(sel)) next
      s <- crown_err_stats(df[[det_col]][sel], df[[fld_col]][sel])
      out[[length(out) + 1]] <- data.frame(
        algo = a, rung = rg, definition = lbl,
        n = s$n, rmse = s$rmse, mae = s$mae, bias = s$bias, r2 = s$r2,
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  # native first, then the sparser rungs in known order, then any unknown rung.
  rk <- match(res$rung, RUNG_LEVELS); rk[is.na(rk)] <- length(RUNG_LEVELS) + 1L
  res[order(res$definition, rk, res$algo), , drop = FALSE]
}
