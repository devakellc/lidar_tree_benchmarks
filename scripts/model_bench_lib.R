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

## ---- frozen clip provider: one decimated set -> two variants + DTM -------
# Decimates the plot clip ONCE (seeded by site/plot/rung), then derives both the
# raw-with-ground and the normalized variant from that SAME decimated set, plus
# a TIN DTM and a JSON manifest. Caches under out_root; reuses on re-call.
# rung = NA means native (no decimation). Returns a list with file paths +
# densities, or NULL if the clip is unusable.
frozen_clip <- function(ctg, site, plot, rung, cx, cy, core_half, out_root,
                        buffer = 25) {
  rdir <- file.path(out_root, site, plot, ifelse(is.na(rung), "native", rung))
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
# from round(rec_<cls> * n_<cls>); understory = intermediate + suppressed.
POOL_CLASSES <- c("dominant", "codominant", "intermediate", "suppressed")
pool <- function(df, classes = POOL_CLASSES) {
  out <- data.frame(
    n_plots = length(unique(paste(df$site, df$plot, df$rung))),
    n_ref = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall = sum(df$TP) / sum(df$n_ref),
    precision = sum(df$tp_core, na.rm = TRUE) / sum(df$n_det))
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
  out
}

## ---- per-(plot,rung) equal-set guard -------------------------------------
# Keep only (site,plot,rung) cells scored by EVERY arm; drop the rest so a rung's
# cross-arm comparison uses an identical plot population. Returns the filtered
# df with the dropped cell keys in attr(.,"dropped").
equal_set_guard <- function(df, arms, key_cols = c("site", "plot", "rung")) {
  k <- do.call(paste, c(df[key_cols], sep = "::"))
  df$.k <- k
  n_arms <- tapply(df$detector, k, function(v) length(unique(v)))
  common <- names(n_arms)[n_arms == length(arms)]
  dropped <- setdiff(unique(k), common)
  out <- df[df$.k %in% common, setdiff(names(df), ".k"), drop = FALSE]
  attr(out, "dropped") <- dropped
  out
}
