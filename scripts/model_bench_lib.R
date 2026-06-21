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

## ---- cross-block apex-cluster dedup (#M8) --------------------------------
# Stacked per-cylinder labelled points (block, inst, X, Y, Z; inst 0/NA =
# unassigned) -> the SAME table relabelled with a globally-consistent integer
# `global_id`. Per (block, inst) apex (max-Z); union-find over apexes restricted
# to pairs in DIFFERENT blocks within horizontal `merge_tol` (so the model's own
# within-cylinder over-segmentation is NEVER laundered). The driver then calls
# reduce_instances(id_col = "global_id"). Returns a 0-row frame (with global_id)
# when nothing is assigned.
dedup_blocks <- function(pts, merge_tol = 2.0, block = "block", id = "inst",
                         x = "X", y = "Y", z = "Z") {
  empty <- data.frame(block = integer(), inst = integer(), X = numeric(),
                      Y = numeric(), Z = numeric(), global_id = integer())
  dt <- as.data.table(pts)
  if (!nrow(dt) || !all(c(block, id, x, y, z) %in% names(dt))) return(empty)
  ids <- dt[[id]]
  dt <- dt[!is.na(ids) & ids != 0, c(block, id, x, y, z), with = FALSE]
  if (!nrow(dt)) return(empty)
  setnames(dt, c(block, id, x, y, z), c("block", "inst", "X", "Y", "Z"))
  dt[, key := .GRP, by = .(block, inst)]
  ap <- dt[, .(blk = block[1L], cx = X[which.max(Z)], cy = Y[which.max(Z)]),
           by = key][order(key)]
  n <- nrow(ap)
  parent <- seq_len(n)
  find <- function(i) { r <- i; while (parent[r] != r) r <- parent[r]
                        while (parent[i] != r) { nx <- parent[i]; parent[i] <<- r; i <- nx }
                        r }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (ap$blk[i] != ap$blk[j] &&
        (ap$cx[i] - ap$cx[j])^2 + (ap$cy[i] - ap$cy[j])^2 <= merge_tol^2) {
      ri <- find(i); rj <- find(j); if (ri != rj) parent[rj] <- ri
    }
  roots <- vapply(seq_len(n), find, integer(1))
  ap[, global_id := match(roots, sort(unique(roots)))]
  dt <- merge(dt, ap[, .(key, global_id)], by = "key")
  dt[, key := NULL]
  as.data.frame(dt)
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

## ---- per-instance apex keyed by instance id ------------------------------
# Like reduce_instances, but KEEPS the instance id so the apex can be joined
# back to crown_diameter_table (keyed by id) after matching. Apex per instance =
# its max-Z point. Returns data.frame(id, x, y, z); 0-row frame when empty.
instance_apex <- function(pts, id_col = "crown_id",
                          x = "X", y = "Y", z = "Z") {
  empty <- data.frame(id = integer(), x = numeric(), y = numeric(),
                      z = numeric())
  dt <- as.data.table(pts)
  if (!nrow(dt) || !id_col %in% names(dt)) return(empty)
  dt <- dt[!is.na(dt[[id_col]]), c(id_col, x, y, z), with = FALSE]
  if (!nrow(dt)) return(empty)
  setnames(dt, c(id_col, x, y, z), c("ID", "X", "Y", "Z"))
  ap <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)), by = ID]
  data.frame(id = ap$ID, x = ap$x, y = ap$y, z = ap$z)
}

## ---- FF3D dedup -> crown-diameter glue (#M8 / #34) -----------------------
# Stacked per-cylinder labelled points (block, inst, X, Y, Z) -> the cross-block
# deduped table (dedup_blocks, global_id) -> the SAME two crown artefacts the
# 3-D arm derives from a single labelling, keyed by `global_id`:
#   $diam : crown_diameter_table(id_col = "global_id", min_pts)
#   $apex : instance_apex(id_col = "global_id")
# dedup_blocks merges cross-block apex duplicates ONLY (never launders the
# model's within-cylinder over-segmentation), so the diameters are per-tree
# (one row per global_id), not per-cylinder. Returns 0-row $diam/$apex when
# nothing is assigned. The driver feeds $diam + $apex straight into
# score_crowns_against_field; the test asserts per-global_id diameters.
ff3d_crown_table <- function(pts, merge_tol = 2.0, min_pts = 5) {
  rel <- dedup_blocks(pts, merge_tol = merge_tol)
  list(diam = crown_diameter_table(rel, id_col = "global_id", min_pts = min_pts),
       apex = instance_apex(rel, id_col = "global_id"))
}

## ---- crown-diameter scoring glue (#30) -----------------------------------
# Pure helper shared by the 3-D crown-diameter arms (Li 2012 / ptrees / AMS3D)
# and reusable by any point-instance segmenter. Given:
#   diam_table : crown_diameter_table() output (id, n_pts, d_eq, d_caliper)
#   apex       : instance_apex() output (id, x, y, z) -- the matching geometry
#   stems      : field stems (individualID, E, N, height, crown_class) already
#                restricted to the plot core
#   field_cd   : per-stem field crown diameter (individualID, maxCrownDiameter,
#                ninetyCrownDiameter)
#   tol        : matching radius (m)
# Matches each instance apex to a field stem with greedy_match (global
# nearest-distance 1:1 + the SAME height-consistency gate as the detection
# arms), joins each matched stem's diameter table row and field crown diameter,
# and returns the canonical crown-metrics rows. d_eq / d_caliper carry through
# verbatim (NA when the instance fell below crown_diameter_table's min_pts);
# area is the equivalent-circle area implied by d_eq (NA when d_eq is NA). algo
# and site/plot are stamped by the caller-supplied values. Returns a 0-row
# canonical frame when nothing matches. greedy_match is provided by sweep_lib.R
# (sourced alongside this lib in every arm).
CROWN_COLS <- c("site", "plot", "algo", "crown_class", "individualID",
                "d_eq", "d_caliper", "area", "field_maxCD", "field_ninetyCD")
score_crowns_against_field <- function(diam_table, apex, stems, field_cd,
                                       tol = 4, site = NA_character_,
                                       plot = NA_character_, algo = NA_character_) {
  empty <- data.frame(site = character(), plot = character(), algo = character(),
                      crown_class = character(), individualID = character(),
                      d_eq = numeric(), d_caliper = numeric(), area = numeric(),
                      field_maxCD = numeric(), field_ninetyCD = numeric(),
                      stringsAsFactors = FALSE)
  if (is.null(apex) || !nrow(apex) || is.null(stems) || !nrow(stems))
    return(empty)
  # match field stems -> instance apexes (position tol + height gate); m[i] is
  # the apex-row index matched to stem i (0 = unmatched), mirroring the CHM arm.
  m <- greedy_match(stems$E, stems$N, apex$x, apex$y, tol,
                    az = stems$height, bz = apex$z)
  matched <- which(m > 0)
  if (!length(matched)) return(empty)
  dt <- if (!is.null(diam_table) && nrow(diam_table))
    as.data.table(diam_table) else NULL
  fc <- as.data.table(field_cd)
  rows <- lapply(matched, function(si) {
    a_id <- apex$id[m[si]]
    de <- NA_real_; dc <- NA_real_
    if (!is.null(dt)) { r <- dt[dt$id == a_id, ]
      if (nrow(r)) { de <- r$d_eq[1]; dc <- r$d_caliper[1] } }
    iid <- stems$individualID[si]
    fr  <- fc[fc$individualID == iid, ]
    data.frame(site = site, plot = plot, algo = algo,
               crown_class = stems$crown_class[si], individualID = iid,
               d_eq = de, d_caliper = dc,
               area = if (is.finite(de)) pi * (de / 2)^2 else NA_real_,
               field_maxCD = if (nrow(fr)) fr$maxCrownDiameter[1] else NA_real_,
               field_ninetyCD = if (nrow(fr)) fr$ninetyCrownDiameter[1] else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[, CROWN_COLS, drop = FALSE]
}

## ---- deterministic seed from a (site, plot, rung) key --------------------
# Maps the key string to a stable non-negative 31-bit integer (FNV-1a), so
# decimation is reproducible across arms and runs without external deps.
# All key characters are ASCII (0-127), so XOR only touches the low 7 bits;
# h is kept as a double with %% 2^32 wrapping to avoid R's signed int32 limits.
seed_for <- function(site, plot, rung) {
  key <- paste(site, plot, ifelse(is.na(rung), "native", rung), sep = "|")
  h <- 2166136261                                  # FNV offset basis (uint32)
  for (b in utf8ToInt(key)) {
    h_low <- as.integer(h %% 128)
    h <- ((h - h_low + bitwXor(h_low, as.integer(b))) * 16777619) %% 2^32
  }
  as.integer(h %% 2^31)
}

## ---- canonical on-disk path for a frozen (site, plot, rung) cell ----------
# SINGLE source of truth for the frozen-clip directory layout, shared by the
# writer (frozen_clip) and every reader (e.g. crown_metrics_sweep's
# frdens_by_rung). out_root is whatever the caller passed to frozen_clip (the
# detection ladder passes d/neon/<SITE>/frozen, so the realized layout is
# out_root/<SITE>/<plot>/<rung> -- the site segment intentionally appears under
# out_root). rung = NA (or "native") collapses to the "native" segment, exactly
# as frozen_clip writes it. Keeping both sides on this one helper means a reader
# path can never silently drift from the writer path again (the bug that left
# the issue #33 density-robustness PNGs un-emitted).
frozen_dir <- function(out_root, site, plot, rung) {
  seg <- if (length(rung) != 1L || is.na(rung) || identical(rung, "native"))
    "native" else as.character(rung)
  file.path(out_root, site, plot, seg)
}

## ---- frozen clip provider: one decimated set -> two variants + DTM -------
# Decimates the plot clip ONCE (seeded by site/plot/rung), then derives both the
# raw-with-ground and the normalized variant from that SAME decimated set, plus
# a TIN DTM and a JSON manifest. Caches under out_root; reuses on re-call.
# rung = NA means native (no decimation). Returns a list with file paths +
# densities, or NULL if the clip is unusable.
frozen_clip <- function(ctg, site, plot, rung, cx, cy, core_half, out_root,
                        buffer = 25) {
  rdir <- frozen_dir(out_root, site, plot, rung)
  fp <- list(rawground   = file.path(rdir, "clip_rawground.laz"),
             normalized  = file.path(rdir, "clip_normalized.laz"),
             dtm         = file.path(rdir, "ground_dtm.tif"),
             manifest    = file.path(rdir, "manifest.json"))
  if (file.exists(fp$manifest)) {                    # cached -> reuse verbatim
    mf <- jsonlite::read_json(fp$manifest, simplifyVector = TRUE)
    return(c(fp, list(pdens = mf$pdens, frdens = mf$frdens, seed = mf$seed)))
  }
  dir.create(rdir, showWarnings = FALSE, recursive = TRUE)
  half <- core_half + buffer
  las  <- lidR::clip_rectangle(ctg, cx - half, cy - half, cx + half, cy + half)
  if (lidR::is.empty(las) || lidR::npoints(las) < 100) return(NULL)
  seed <- seed_for(site, plot, rung)
  if (!is.na(rung)) { set.seed(seed)
    las <- lidR::decimate_points(las, lidR::homogenize(density = rung, res = 5)) }
  if (sum(las$Classification == 2L) < 10) return(NULL)   # need ground for DTM
  dtm <- lidR::rasterize_terrain(las, res = 1, algorithm = lidR::tin())
  nrm <- lidR::normalize_height(las, lidR::tin(), na.rm = TRUE)
  nrm <- lidR::filter_poi(nrm, Z >= -1, Z < 80)
  area   <- (2 * half)^2
  pdens  <- lidR::npoints(nrm) / area
  frdens <- sum(nrm$ReturnNumber == 1L) / area
  lidR::writeLAS(las, fp$rawground)
  lidR::writeLAS(nrm, fp$normalized)
  terra::writeRaster(dtm, fp$dtm, overwrite = TRUE)
  jsonlite::write_json(list(site = site, plot = plot,
                            rung = ifelse(is.na(rung), "native", rung),
                            seed = seed, n_raw = lidR::npoints(las),
                            n_norm = lidR::npoints(nrm),
                            pdens = round(pdens, 3), frdens = round(frdens, 3),
                            buffer = buffer, core_half = core_half),
                       fp$manifest, auto_unbox = TRUE, pretty = TRUE)
  c(fp, list(pdens = pdens, frdens = frdens, seed = seed))
}

## ---- canonical pooler: sum counts, never average rates -------------------
# df: long-form scored rows (one per site x plot x rung x detector subset).
# Pools to a single row: recall = sum(TP)/sum(n_ref), precision =
# sum(tp_core)/sum(n_det), F1 from the pooled rates; per-class recall recovered
# from round(rec_<cls> * n_<cls>); understory = intermediate + suppressed. When
# the height-band columns (n_h_short/mid/tall) are present it pools those too,
# by the same summed-count rule; absent, no rec_h_*/n_h_* columns are emitted.
POOL_CLASSES <- c("dominant", "codominant", "intermediate", "suppressed")
pool <- function(df, classes = POOL_CLASSES) {
  if (is.null(df$tp_core)) df$tp_core <- round(df$precision * df$n_det)
  out <- data.frame(
    n_plots = length(unique(paste(df$site, df$plot, df$rung, sep = "::"))),
    n_ref = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall = sum(df$TP) / sum(df$n_ref),
    precision = if (sum(df$n_det) > 0) sum(df$tp_core, na.rm = TRUE) / sum(df$n_det) else NA_real_)
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  for (cl in classes) {
    nref <- sum(df[[paste0("n_", cl)]], na.rm = TRUE)
    tp   <- sum(ifelse(df[[paste0("n_", cl)]] > 0,
                       round(df[[paste0("rec_", cl)]] * df[[paste0("n_", cl)]]), 0),
                na.rm = TRUE)
    out[[paste0("rec_", cl)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_", cl)]]   <- nref
  }
  nref_u <- sum(df$n_intermediate, na.rm = TRUE) + sum(df$n_suppressed, na.rm = TRUE)
  tp_u <- sum(ifelse(df$n_intermediate > 0,
                     round(df$rec_intermediate * df$n_intermediate), 0), na.rm = TRUE) +
          sum(ifelse(df$n_suppressed > 0,
                     round(df$rec_suppressed * df$n_suppressed), 0), na.rm = TRUE)
  out$rec_understory <- if (nref_u) tp_u / nref_u else NA_real_
  out$n_understory   <- nref_u
  hbands <- c("short", "mid", "tall")
  if (all(paste0("n_h_", hbands) %in% names(df))) for (b in hbands) {
    nref <- sum(df[[paste0("n_h_", b)]], na.rm = TRUE)
    tp   <- sum(ifelse(df[[paste0("n_h_", b)]] > 0,
                       round(df[[paste0("rec_h_", b)]] * df[[paste0("n_h_", b)]]), 0),
                na.rm = TRUE)
    out[[paste0("rec_h_", b)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_h_", b)]]   <- nref
  }
  out
}

## ---- per-(plot,rung) equal-set guard -------------------------------------
# Keep only (site,plot,rung) cells scored by EVERY arm; drop the rest so a rung's
# cross-arm comparison uses an identical plot population. Returns the filtered
# df with the dropped cell keys in attr(.,"dropped").
equal_set_guard <- function(df, arms, key_cols = c("site", "plot", "rung")) {
  k <- do.call(paste, c(df[key_cols], sep = "::"))
  df$.k <- k
  have_all <- tapply(df$detector, k, function(v) all(arms %in% v))
  common <- names(have_all)[have_all]
  dropped <- setdiff(unique(k), common)
  out <- df[df$.k %in% common, setdiff(names(df), ".k"), drop = FALSE]
  attr(out, "dropped") <- dropped
  out
}

## ---- scorer-contract conformance harness ---------------------------------
# Asserts `det` satisfies what score_plot consumes: a base data.frame (not sf),
# with exactly numeric columns x, y, z (lowercase). Empty must be a 0-row frame,
# never NULL. Returns TRUE invisibly or stops with a specific message.
assert_detection_contract <- function(det) {
  if (is.null(det)) stop("detection table is NULL; return a 0-row data.frame instead")
  if (!is.data.frame(det)) stop("detection table must be a base data.frame")
  if (inherits(det, "sf")) stop("detection table must not be an sf object")
  if (!identical(names(det), c("x", "y", "z")))
    stop("detection table must have exactly columns x, y, z (lowercase)")
  if (!all(vapply(det, is.numeric, logical(1))))
    stop("detection columns x, y, z must be numeric")
  invisible(TRUE)
}

## ==========================================================================
## #V1 point-set IoU / Coverage / Panoptic-Quality scorer ===================
## ==========================================================================
# The whole benchmark is graded by greedy_match apex-distance 1:1 matching,
# which model-benchmark-results.md flags as blind to point-level IoU and crown
# shape. These pure helpers add the mask-aware suite the external leaderboards
# (FGI-EMIT, FOR-instanceV2) use -- point-set IoU>=0.5 matching, Coverage, and
# Panoptic Quality -- alongside that distance matching (NOT replacing it).
#
# All five helpers below are pure set-algebra / geometry over a SHARED point
# substrate: `pred` and `ref` are equal-length integer label vectors (one entry
# per substrate point; NA = unassigned), so they are unit-testable with no NEON
# data. The driver (score_instances_iou.R) builds those vectors -- ref = the
# Voronoi-on-stems partition of the frozen normalized cloud (assign_points_to_
# stems), pred = the model's per-point labels projected onto that SAME cloud
# (transfer_labels) -- then accumulates per (plot, rung, crown_class) and pools
# by SUMMING the accumulators (sum matched-IoU, sum TP/FP/FN), mirroring the
# sum(TP)/sum(n_ref) rule in pool().

## ---- pairwise point-set IoU ----------------------------------------------
# pred, ref: equal-length integer label vectors over the shared substrate (NA =
# unassigned, dropped). Returns one row per (pred instance, ref instance) pair
# that shares >=1 point: inter = shared points; pred_size/ref_size = each label's
# TOTAL point count over the whole substrate (so the union counts points the
# other partition left unassigned); union = pred_size + ref_size - inter;
# iou = inter / union. Pairs with zero intersection have IoU 0 and are omitted
# (they never match and never raise Coverage). 0-row frame when nothing overlaps.
point_set_iou <- function(pred, ref) {
  stopifnot(length(pred) == length(ref))
  empty <- data.frame(pred_id = integer(), ref_id = integer(), inter = integer(),
                      pred_size = integer(), ref_size = integer(),
                      union = integer(), iou = numeric())
  pred <- as.integer(pred); ref <- as.integer(ref)
  if (!length(pred)) return(empty)
  psz <- table(pred[!is.na(pred)])                  # total size of each pred label
  rsz <- table(ref[!is.na(ref)])                    # total size of each ref label
  both <- !is.na(pred) & !is.na(ref)
  if (!any(both)) return(empty)
  agg <- as.data.table(list(p = pred[both], r = ref[both]))[, .(inter = .N),
                                                            by = .(p, r)]
  agg[, pred_size := as.integer(psz[as.character(p)])]
  agg[, ref_size  := as.integer(rsz[as.character(r)])]
  agg[, union := pred_size + ref_size - inter]
  agg[, iou := inter / union]
  data.frame(pred_id = agg$p, ref_id = agg$r, inter = agg$inter,
             pred_size = agg$pred_size, ref_size = agg$ref_size,
             union = agg$union, iou = agg$iou)
}

## ---- greedy max-IoU 1:1 matcher with an IoU>=gate gate --------------------
# Sibling of greedy_match (sweep_lib.R), but ranks candidate (pred, ref) pairs by
# DESCENDING IoU instead of ascending distance: the best-IoU pair is claimed
# first, so a dominant crown's prediction is matched before a worse overlap can
# steal either side. 1:1 (each pred and each ref used once); only pairs with
# IoU >= gate (default 0.5, the FGI-EMIT/FOR-instance threshold) are accepted.
# iou_df: point_set_iou() output. Returns the matched pairs (pred_id, ref_id, iou)
# = the true positives; 0-row frame when none clear the gate.
iou_match <- function(iou_df, gate = 0.5) {
  empty <- data.frame(pred_id = integer(), ref_id = integer(), iou = numeric())
  if (is.null(iou_df) || !nrow(iou_df)) return(empty)
  d <- as.data.frame(iou_df)[, c("pred_id", "ref_id", "iou"), drop = FALSE]
  d <- d[d$iou >= gate, , drop = FALSE]
  if (!nrow(d)) return(empty)
  d <- d[order(-d$iou), , drop = FALSE]
  seen_p <- integer(0); seen_r <- integer(0); keep <- logical(nrow(d))
  for (k in seq_len(nrow(d))) {
    if (!(d$pred_id[k] %in% seen_p) && !(d$ref_id[k] %in% seen_r)) {
      keep[k] <- TRUE
      seen_p <- c(seen_p, d$pred_id[k]); seen_r <- c(seen_r, d$ref_id[k])
    }
  }
  out <- d[keep, , drop = FALSE]; rownames(out) <- NULL
  out
}

## ---- Coverage: per-reference max IoU over predictions ---------------------
# For every reference instance, the maximum IoU against ANY predicted instance
# (0 when no prediction overlaps it). Coverage = mean(max_iou) over references is
# the ungated companion to recall@0.5: how well each tree's points are captured
# by its best-matching prediction, regardless of the 0.5 cutoff. Returns one row
# per reference id (zero-overlap refs included with max_iou 0); 0-row frame when
# there are no references.
coverage_table <- function(pred, ref) {
  empty <- data.frame(ref_id = integer(), max_iou = numeric())
  ref <- as.integer(ref)
  ref_ids <- sort(unique(ref[!is.na(ref)]))
  if (!length(ref_ids)) return(empty)
  io <- point_set_iou(pred, ref)
  mx <- if (nrow(io)) tapply(io$iou, io$ref_id, max) else numeric(0)
  out <- data.frame(ref_id = ref_ids,
                    max_iou = as.numeric(mx[as.character(ref_ids)]))
  out$max_iou[is.na(out$max_iou)] <- 0                # refs with no overlap -> 0
  out
}

## ---- Panoptic Quality (PQ = SQ * RQ) -------------------------------------
# matched: iou_match() output (the TP pairs + their IoUs). n_pred / n_ref: the
# number of predicted / reference instances on the substrate. Accumulators (TP,
# FP, FN, sum_iou) are returned RAW so the driver pools them by summing, then
# recomputes the rates from the pooled sums:
#   SQ (segmentation quality) = mean IoU over matched pairs = sum_iou / TP
#   RQ (recognition quality)  = TP / (TP + 0.5 FP + 0.5 FN)  (== F1 over matches)
#   PQ                        = SQ * RQ
# SQ is NA when there are no matches (mean over nothing); PQ is 0 then (RQ 0).
panoptic_quality <- function(matched, n_pred, n_ref) {
  TP <- if (is.null(matched)) 0L else nrow(matched)
  sum_iou <- if (TP > 0) sum(matched$iou) else 0
  FP <- n_pred - TP; FN <- n_ref - TP
  SQ <- if (TP > 0) sum_iou / TP else NA_real_
  denom <- TP + 0.5 * FP + 0.5 * FN
  RQ <- if (denom > 0) TP / denom else NA_real_
  PQ <- if (TP > 0 && is.finite(SQ) && is.finite(RQ)) SQ * RQ else 0
  precision <- if (n_pred > 0) TP / n_pred else NA_real_
  recall    <- if (n_ref > 0)  TP / n_ref  else NA_real_
  F1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  data.frame(TP = TP, FP = FP, FN = FN, sum_iou = sum_iou,
             n_pred = n_pred, n_ref = n_ref,
             SQ = SQ, RQ = RQ, PQ = PQ,
             precision = precision, recall = recall, F1 = F1)
}

## ---- reference instances: Voronoi-on-stems within a per-stem crown radius -
# Builds the ground-truth point partition the IoU suite scores against. Each
# substrate point (px, py) is assigned to the NEAREST field stem, but only if it
# falls within THAT stem's crown radius (radius = field maxCrownDiameter / 2,
# per stem); otherwise it is left unassigned (NA). This is an explicit
# Voronoi-on-stems PROXY for true field crown masks -- it is the best reference
# obtainable from stem points + a measured crown width, and is documented as such
# in results/instance-iou-pq-results.md. radius may be a scalar or a per-stem
# vector (recycled to length(sx)). Ties go to the first stem (which.min). Returns
# an integer vector (stem index 1..n, or NA) of length(px).
assign_points_to_stems <- function(px, py, sx, sy, radius) {
  n <- length(px); out <- rep(NA_integer_, n); ns <- length(sx)
  if (!n || !ns) return(out)
  radius <- rep_len(radius, ns)
  best_d2 <- rep(Inf, n); best_s <- rep(NA_integer_, n)
  for (s in seq_len(ns)) {
    d2  <- (px - sx[s])^2 + (py - sy[s])^2
    upd <- d2 < best_d2                              # strict -> first stem wins ties
    best_d2[upd] <- d2[upd]; best_s[upd] <- s
  }
  r <- radius[best_s]                                 # NA for unmatched/NA-radius stems
  ok <- !is.na(best_s) & is.finite(r) & best_d2 <= r^2
  out[ok] <- best_s[ok]
  out
}

## ---- dependency-free nearest-neighbour label transfer --------------------
# Projects a model's per-point instance labels onto the reference substrate so
# both partitions live on the SAME points (the standard instance-IoU evaluation
# uses one fixed eval cloud). The persisted model clouds are NOT the frozen
# normalized cloud point-for-point (different Z datum, container voxel
# downsampling, FF3D cylinder overlap), so each substrate point (dx, dy) takes
# the label of the nearest SOURCE point (sx, sy, sid) within `tol` metres in XY;
# beyond tol it stays NA (the prediction did not reach that point). Implemented
# as a grid hash (cell = tol; each source serves its own cell + 8 neighbours, so
# any source within tol of a destination is in the searched 3x3 block) to avoid a
# new spatial-index dependency (RANN/FNN are not installed). Deterministic.
# Returns an integer vector of length(dx).
transfer_labels <- function(sx, sy, sid, dx, dy, tol = 0.5) {
  n <- length(dx); out <- rep(NA_integer_, n)
  if (!n || !length(sx)) return(out)
  cs  <- tol
  src <- data.table(sx = sx, sy = sy, sid = as.integer(sid),
                    cx = floor(sx / cs), cy = floor(sy / cs))
  ns  <- nrow(src)
  off <- expand.grid(ox = -1:1, oy = -1:1)           # 3x3 serving offsets
  srv <- src[rep(seq_len(ns), 9L)]
  srv[, kx := cx + rep(off$ox, each = ns)]
  srv[, ky := cy + rep(off$oy, each = ns)]
  dst <- data.table(di = seq_len(n), dx = dx, dy = dy,
                    kx = floor(dx / cs), ky = floor(dy / cs))
  j <- merge(dst, srv, by = c("kx", "ky"), allow.cartesian = TRUE)
  if (!nrow(j)) return(out)
  j[, d2 := (dx - sx)^2 + (dy - sy)^2]
  j <- j[d2 <= tol^2]
  if (!nrow(j)) return(out)
  setorder(j, di, d2)
  best <- j[, .(sid = sid[1L]), by = di]             # nearest source per dest
  out[best$di] <- best$sid
  out
}

## ---- one (plot x rung x model) accumulator row ---------------------------
# Scores ONE predicted partition against the reference partition on a shared
# substrate and returns a single long-form row of SUMMABLE accumulators (the
# rate columns are per-cell convenience; pool_pq recomputes them from the sums).
# pred, ref: integer label vectors over the substrate (NA = unassigned).
# ref_class: named character vector keyed by ref id (as character) giving each
# reference instance's crown_class (a ref id absent from it -> NA class: counted
# in the cell totals but in no per-class bucket). n_ref counts reference
# instances that captured >=1 substrate point (the natural denominator for a
# point-set metric -- a stem with no canopy points has no mask to detect).
# Per class cl: n_ (ref instances), tp_ (matched), fn_, sumiou_ (matched-IoU sum
# = the class SQ numerator), sumcov_ (max-IoU sum = the class Coverage numerator).
score_instance_cell <- function(pred, ref, ref_class, gate = 0.5,
                               classes = POOL_CLASSES) {
  io  <- point_set_iou(pred, ref)
  mt  <- iou_match(io, gate = gate)
  cov <- coverage_table(pred, ref)
  n_pred <- length(unique(pred[!is.na(pred)]))
  n_ref  <- nrow(cov)
  pq <- panoptic_quality(mt, n_pred = n_pred, n_ref = n_ref)
  sum_maxiou <- if (nrow(cov)) sum(cov$max_iou) else 0
  out <- data.frame(
    n_pred = n_pred, n_ref = n_ref, TP = pq$TP, FP = pq$FP, FN = pq$FN,
    sum_iou = pq$sum_iou, sum_maxiou = sum_maxiou,
    SQ = pq$SQ, RQ = pq$RQ, PQ = pq$PQ, precision = pq$precision,
    recall = pq$recall, F1 = pq$F1,
    coverage = if (n_ref > 0) sum_maxiou / n_ref else NA_real_)
  rc <- function(id) if (length(ref_class)) ref_class[as.character(id)] else
    rep(NA_character_, length(id))
  matched_cls <- if (nrow(mt))  rc(mt$ref_id)  else character(0)
  cov_cls     <- if (nrow(cov)) rc(cov$ref_id) else character(0)
  for (cl in classes) {
    n_c  <- sum(cov_cls == cl, na.rm = TRUE)
    tp_c <- sum(matched_cls == cl, na.rm = TRUE)
    out[[paste0("n_",  cl)]]     <- n_c
    out[[paste0("tp_", cl)]]     <- tp_c
    out[[paste0("fn_", cl)]]     <- n_c - tp_c
    out[[paste0("sumiou_", cl)]] <- if (nrow(mt))
      sum(mt$iou[matched_cls == cl], na.rm = TRUE) else 0
    out[[paste0("sumcov_", cl)]] <- if (nrow(cov))
      sum(cov$max_iou[cov_cls == cl], na.rm = TRUE) else 0
  }
  out
}

## ---- canonical PQ pooler: SUM accumulators, recompute rates ---------------
# The IoU-suite sibling of pool(): pools long-form score_instance_cell rows to a
# single row by SUMMING the panoptic accumulators (TP/FP/FN, sum_iou, sum_maxiou,
# per-class counts) -- never averaging the per-cell rates, the same rule that
# makes pool() recall sum(TP)/sum(n_ref). Pooled rates are then recomputed from
# the sums: SQ = sum_iou/TP, RQ = TP/(TP+0.5FP+0.5FN), PQ = SQ*RQ, precision =
# TP/n_pred, recall = TP/n_ref, Coverage = sum_maxiou/n_ref; per class rec_ =
# tp/n, SQ_ = sumiou/tp, cov_ = sumcov/n. understory = intermediate + suppressed.
pool_pq <- function(df, classes = POOL_CLASSES) {
  s <- function(col) if (is.null(df[[col]])) 0 else sum(df[[col]], na.rm = TRUE)
  TP <- s("TP"); FP <- s("FP"); FN <- s("FN")
  sum_iou <- s("sum_iou"); sum_maxiou <- s("sum_maxiou")
  n_pred <- s("n_pred"); n_ref <- s("n_ref")
  SQ <- if (TP > 0) sum_iou / TP else NA_real_
  denom <- TP + 0.5 * FP + 0.5 * FN
  RQ <- if (denom > 0) TP / denom else NA_real_
  precision <- if (n_pred > 0) TP / n_pred else NA_real_
  recall    <- if (n_ref  > 0) TP / n_ref  else NA_real_
  F1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  out <- data.frame(
    n_cells = nrow(df), n_pred = n_pred, n_ref = n_ref,
    TP = TP, FP = FP, FN = FN, sum_iou = sum_iou, sum_maxiou = sum_maxiou,
    SQ = SQ, RQ = RQ,
    PQ = if (TP > 0 && is.finite(SQ) && is.finite(RQ)) SQ * RQ else 0,
    precision = precision, recall = recall, F1 = F1,
    coverage = if (n_ref > 0) sum_maxiou / n_ref else NA_real_)
  for (cl in classes) {
    n_c  <- s(paste0("n_", cl)); tp_c <- s(paste0("tp_", cl))
    out[[paste0("n_",   cl)]] <- n_c
    out[[paste0("rec_", cl)]] <- if (n_c  > 0) tp_c / n_c else NA_real_
    out[[paste0("SQ_",  cl)]] <- if (tp_c > 0) s(paste0("sumiou_", cl)) / tp_c else NA_real_
    out[[paste0("cov_", cl)]] <- if (n_c  > 0) s(paste0("sumcov_", cl)) / n_c else NA_real_
  }
  n_u  <- s("n_intermediate") + s("n_suppressed")
  tp_u <- s("tp_intermediate") + s("tp_suppressed")
  out$n_understory   <- n_u
  out$rec_understory <- if (n_u  > 0) tp_u / n_u else NA_real_
  out$SQ_understory  <- if (tp_u > 0)
    (s("sumiou_intermediate") + s("sumiou_suppressed")) / tp_u else NA_real_
  out$cov_understory <- if (n_u  > 0)
    (s("sumcov_intermediate") + s("sumcov_suppressed")) / n_u else NA_real_
  out
}

## ==========================================================================
## #P1 cross-arm detector fusion =============================================
## ==========================================================================
# Every arm is scored in isolation today; nothing combines them. The per-class
# table shows complementarity to exploit (CHM-VWF/multichm carry overstory,
# SegmentAnyTree carries understory), so a consensus/NMS layer should trade recall
# against precision. These helpers cluster per-cell apexes across arms and emit
# the union / majority / height-layered operating points the fusion driver scores.

## ---- cross-arm apex clustering (single-linkage, height-gated) -------------
# Stacked per-cell apexes from N detector arms (arm, x, y, z) -> fused trees.
# Union-find over apex pairs from DIFFERENT arms within horizontal merge_tol AND
# |dz| <= z_tol (the dedup_blocks union-find + a fusion height gate, so an
# overstory CHM apex and an understory point apex at the same x,y stay SEPARATE
# trees). Same-arm apexes never merge directly (each arm's detections are distinct
# trees), though single-linkage can transitively chain through a cross-arm bridge.
# Per cluster: representative = the max-Z member; votes = number of DISTINCT arms;
# arms = sorted comma string. Returns data.frame(cluster,x,y,z,votes,arms); 0-row
# when empty.
fuse_apexes <- function(arm, x, y, z, merge_tol = 2.0, z_tol = 5.0) {
  empty <- data.frame(cluster = integer(), x = numeric(), y = numeric(),
                      z = numeric(), votes = integer(), arms = character(),
                      stringsAsFactors = FALSE)
  n <- length(x)
  if (!n) return(empty)
  arm <- as.character(arm)
  parent <- seq_len(n)
  find <- function(i) { r <- i; while (parent[r] != r) r <- parent[r]
                        while (parent[i] != r) { nx <- parent[i]; parent[i] <<- r; i <- nx }
                        r }
  if (n > 1L) for (i in seq_len(n - 1L)) for (j in (i + 1L):n)
    if (arm[i] != arm[j] &&
        (x[i] - x[j])^2 + (y[i] - y[j])^2 <= merge_tol^2 &&
        abs(z[i] - z[j]) <= z_tol) {
      ri <- find(i); rj <- find(j); if (ri != rj) parent[rj] <- ri
    }
  roots <- vapply(seq_len(n), find, integer(1))
  cl <- match(roots, sort(unique(roots)))
  dt <- as.data.table(list(cluster = cl, arm = arm, x = x, y = y, z = z))
  agg <- dt[, { k <- which.max(z)
                .(x = x[k], y = y[k], z = z[k],
                  votes = length(unique(arm)),
                  arms = paste(sort(unique(arm)), collapse = ",")) },
            by = cluster][order(cluster)]
  as.data.frame(agg)
}

## ---- fusion operating points: union / majority / layered -----------------
# Three operating points over the stacked per-cell apexes:
#   union    -- every fused cluster (k>=1): recall-max.
#   majority -- clusters seen by >= ceil(n_arms/2) arms: precision-max.
#   layered  -- height-stratified: chm-family apexes in the overstory band
#               (z >= overstory_frac * max z) + point/deep apexes in the
#               understory band (z < that split), then fused -- trusting the CHM
#               for the overstory and the point/deep arms for the understory.
# n_arms = the number of arms actually scored in this cell (the majority
# denominator). Returns list(union, majority, layered), each data.frame(x,y,z,votes).
fusion_points <- function(arm, x, y, z, n_arms, merge_tol = 2.0, z_tol = 5.0,
                          chm_arms = c("chm_vwf", "multichm"),
                          overstory_frac = 0.5) {
  cols <- function(f) f[, c("x", "y", "z", "votes"), drop = FALSE]
  f_all <- fuse_apexes(arm, x, y, z, merge_tol, z_tol)
  union <- cols(f_all)
  majority <- cols(f_all[f_all$votes >= ceiling(n_arms / 2), , drop = FALSE])
  layered <- if (!length(x)) cols(f_all) else {
    H <- overstory_frac * max(z)
    is_chm <- arm %in% chm_arms
    keep <- (is_chm & z >= H) | (!is_chm & z < H)
    cols(fuse_apexes(arm[keep], x[keep], y[keep], z[keep], merge_tol, z_tol))
  }
  list(union = union, majority = majority, layered = layered)
}

## ==========================================================================
## #P4 per-detection confidence calibration =================================
## ==========================================================================
# Weighted consensus in #P1 needs each detection to carry a TRUSTWORTHY
# probability, but raw scores are not comparable across arms (a CNN's 0.6 != a
# watershed prominence of 0.6). These pure helpers turn a per-arm raw score +
# TP/FP labels into a reliability diagram, an Expected Calibration Error, and a
# monotone isotonic calibrator that maps raw scores to empirical precision -- so
# a calibrated 0.7 means ~70% precision for EVERY arm and fuse_detectors() can
# compare them. precision_at_recall measures the ensemble payoff: ranking the
# pooled multi-arm detections by calibrated (comparable) score vs raw score.

## ---- reliability table (equal-width bins over [0,1]) ----------------------
# prob: predicted probability in [0,1]; label: 0/1 (TP=1). Splits [0,1] into
# `bins` equal-width bins; per bin reports n, mean predicted prob, and observed
# accuracy (TP rate). Empty bins are kept (n=0, NaN stats) so the table is a
# stable length. The companion to expected_calibration_error.
reliability_table <- function(prob, label, bins = 10) {
  prob <- as.numeric(prob); label <- as.numeric(label)
  br <- seq(0, 1, length.out = bins + 1L)
  idx <- findInterval(pmin(pmax(prob, 0), 1), br, rightmost.closed = TRUE,
                      all.inside = TRUE)
  out <- data.frame(bin = seq_len(bins), n = 0L,
                    mean_prob = NA_real_, obs_acc = NA_real_)
  if (length(prob)) for (b in seq_len(bins)) {
    sel <- idx == b
    if (any(sel)) { out$n[b] <- sum(sel)
      out$mean_prob[b] <- mean(prob[sel]); out$obs_acc[b] <- mean(label[sel]) }
  }
  out
}

## ---- Expected Calibration Error ------------------------------------------
# The n-weighted mean gap between predicted confidence and observed accuracy
# across the reliability bins: ECE = sum_b (n_b/N) * |obs_acc_b - mean_prob_b|.
# 0 = perfectly calibrated. NA when there are no observations.
expected_calibration_error <- function(prob, label, bins = 10) {
  rt <- reliability_table(prob, label, bins = bins)
  rt <- rt[rt$n > 0, , drop = FALSE]
  N <- sum(rt$n)
  if (!N) return(NA_real_)
  sum(rt$n / N * abs(rt$obs_acc - rt$mean_prob))
}

## ---- isotonic (monotone) post-hoc calibrator ------------------------------
# Fits a non-decreasing step mapping raw score -> P(TP) via stats::isoreg (PAVA),
# then returns a function that interpolates it to arbitrary new scores (constant
# extrapolation past the trained range). Monotone, so it NEVER changes a single
# arm's detection RANKING (precision-recall is rank-based) -- its value is making
# scores COMPARABLE across arms. With no positives it calibrates everything to 0
# (and to 1 with no negatives). NA scores map to NA.
isotonic_calibrate <- function(score, label) {
  score <- as.numeric(score); label <- as.numeric(label)
  ok <- is.finite(score) & is.finite(label)
  s <- score[ok]; y <- label[ok]
  if (!length(s)) return(function(z) rep(NA_real_, length(z)))
  if (all(y == y[1])) { const <- y[1]; return(function(z) ifelse(is.na(z), NA_real_, const)) }
  o <- order(s); ss <- s[o]
  fit <- stats::isoreg(ss, y[o])
  yf  <- fit$yf                                   # fitted, in ss order, non-decreasing
  # collapse ties in ss to the last fitted value, then linear-interpolate
  ux <- !duplicated(ss, fromLast = TRUE)
  kx <- ss[ux]; ky <- yf[ux]
  function(z) {
    z <- as.numeric(z)
    out <- approx(kx, ky, xout = z, rule = 2, ties = "ordered")$y
    out[is.na(z)] <- NA_real_
    pmin(pmax(out, 0), 1)
  }
}

## ---- precision at a fixed recall over a scored detection set --------------
# Rank detections by DESCENDING score, walk the threshold down, and return the
# precision at the first point where cumulative recall (TP_kept / total_TP)
# reaches target_recall. This is the ensemble calibration payoff: comparing the
# pooled multi-arm ranking under raw vs calibrated scores. NA if no positives.
# Ties in score are broken stably (label order preserved) -- callers compare
# the SAME detections under two scorings, so the comparison is fair.
precision_at_recall <- function(score, label, target_recall = 0.8) {
  score <- as.numeric(score); label <- as.numeric(label)
  P <- sum(label == 1)
  if (!P) return(NA_real_)
  o <- order(-score)
  y <- label[o]
  tp <- cumsum(y == 1); kept <- seq_along(y)
  recall <- tp / P; precision <- tp / kept
  hit <- which(recall >= target_recall)
  if (!length(hit)) return(NA_real_)
  precision[hit[1]]
}
