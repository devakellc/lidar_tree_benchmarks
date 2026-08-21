#!/usr/bin/env Rscript
# Pure helpers for the #V5 coverage-gap crediting study (coverage_gap.R):
# modality-family mapping, best_treetop_cache readers, and optical-box ->
# detection conversion. No I/O beyond the explicit readers; unit-tested with
# synthetic fixtures in tests/testthat/test-coverage-gap.R.

## ---- modality families -----------------------------------------------------
# Crediting independence is defined at the FAMILY level, not the arm level:
# errors correlate strongly within a family (all CHM local-maxima share the
# same surface artefacts; the RGB arms share the same imagery), so "co-detected
# by >= min_fam families" is the meaningful independence test, and a target
# arm's own family never testifies for it.
FAMILY_MAP <- c(
  chm_vwf = "chm", lmfauto = "chm", multichm = "chm",
  ptrees = "pc", ams3d = "pc", li2012 = "pc", lidr_li2012 = "pc",
  lidr_lmf_pc = "pc", lasr_lmax_pc = "pc",
  treeisonet = "deep", segmentanytree = "deep", forestformer3d = "deep",
  deepforest = "rgb", detectree2 = "rgb")

arm_family <- function(arm) unname(FAMILY_MAP[as.character(arm)])

## ---- best_treetop_cache reader ---------------------------------------------
# Reads one arm's cached (x, y, z-AGL) detections for a (site, plot, rung) cell
# as written by export_best_treetops_geojson.R: either the plain
# <arm>__<site>__<plot>__<rung>.csv or the param-suffixed
# <arm>__<site>__<plot>__<rung>__<params...>.csv (chm_vwf / treeisonet / the
# GPU arms). The arm name anchors the whole basename, so "li2012" can never
# swallow a "lidr_li2012" file. NULL = no cache / unreadable / wrong schema.
read_arm_cache <- function(dir, arm, site, plot, rung, params = NULL) {
  exact <- file.path(dir, sprintf("%s__%s__%s__%s.csv", arm, site, plot, rung))
  hit <- if (file.exists(exact)) exact else NA_character_
  # exact parameter pinning (selection manifest chm_res/vwf_a): never let a
  # stale variant shadow the selected leaderboard cell.
  if (is.na(hit) && !is.null(params) && length(params)) {
    pinned <- file.path(dir, sprintf("%s__%s__%s__%s__%s.csv", arm, site, plot,
                                     rung, paste(params, collapse = "__")))
    if (file.exists(pinned)) hit <- pinned
    else if (length(Sys.glob(file.path(dir, sprintf("%s__%s__%s__%s__*.csv",
                                                    arm, site, plot, rung)))))
      # only a real shadowing risk: the cell IS cached, just not at the pinned
      # parameters. An uncached cell (no variants at all) is silently skipped.
      warning(sprintf(
        "read_arm_cache: pinned variant %s missing; falling back to glob",
        basename(pinned)), call. = FALSE)
  }
  if (is.na(hit)) {
    g <- Sys.glob(file.path(dir, sprintf("%s__%s__%s__%s__*.csv",
                                         arm, site, plot, rung)))
    if (length(g) > 1)
      warning(sprintf(
        "read_arm_cache: %d parameter variants for %s__%s__%s__%s; using %s",
        length(g), arm, site, plot, rung, basename(g[1])), call. = FALSE)
    hit <- if (length(g)) g[1] else NA_character_
  }
  if (is.na(hit)) return(NULL)
  det <- tryCatch(read.csv(hit, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(det) || !all(c("x", "y", "z") %in% names(det))) return(NULL)
  data.frame(x = as.numeric(det$x), y = as.numeric(det$y), z = as.numeric(det$z))
}

## ---- per-cell crediting orchestration ---------------------------------------
# Counts how many of a target arm's ISOLATED core FPs (fps: x, y, isolated from
# fp_points) are co-detected by witness FPs (wit: x, y, fam, isolated pooled
# over the other arms). Only isolated witnesses testify (a near-FP is over-seg
# of a mapped tree, not evidence of an unmapped one) and the target arm's own
# family is struck from the witness pool before co_detect_credit runs.
credit_isolated <- function(fps, wit, target_fam, r = 2.0, min_fam = 2) {
  # Prefer the hardened credit_eligible flag (isolated AND not near ANY mapped
  # stem) on both sides when present; fall back to the #V4 isolated split.
  et <- if (!is.null(fps$credit_eligible)) fps$credit_eligible else fps$isolated
  iso <- fps[et, , drop = FALSE]
  if (!nrow(iso)) return(0L)
  ew <- if (!is.null(wit$credit_eligible)) wit$credit_eligible else wit$isolated
  w <- wit[ew & wit$fam != target_fam, , drop = FALSE]
  cr <- co_detect_credit(iso$x, iso$y, w$x, w$y, w$fam, r = r, min_fam = min_fam)
  if (!any(cr)) return(0L)
  # One credit per probable tree, not per duplicate detection: credited FPs
  # within r of an already-counted credit are over-segmentation of the SAME
  # unmapped tree and must stay in the precision denominator. Greedy
  # suppression in deterministic (x, y) order.
  cx <- iso$x[cr]; cy <- iso$y[cr]
  o <- order(cx, cy)
  n <- 0L; kx <- numeric(0); ky <- numeric(0)
  for (i in o) {
    if (!length(kx) || all(sqrt((kx - cx[i])^2 + (ky - cy[i])^2) > r)) {
      n <- n + 1L; kx <- c(kx, cx[i]); ky <- c(ky, cy[i])
    }
  }
  n
}

## ---- best-configuration selection manifest -----------------------------------
# Reads export_best_treetops_geojson.R's best_treetop_selection.csv and returns
# the (method, rung) rows selected for `site` -- the authoritative "each arm at
# its best tested configuration" manifest, so stale extra cache files can never
# add cells. NULL when the manifest is absent/unreadable or has no rows for the
# site (callers fall back to cached_rungs discovery).
read_selection <- function(path, site) {
  if (!file.exists(path)) return(NULL)
  s <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(s) || !all(c("site", "method", "rung") %in% names(s))) return(NULL)
  keep <- intersect(c("site", "method", "rung", "chm_res", "vwf_a"), names(s))
  s <- s[s$site == site, keep, drop = FALSE]
  if (!nrow(s)) return(NULL)
  s$rung <- as.character(s$rung)
  s
}

# Enumerate the rungs an arm has cached for one plot (each arm is cached at
# its best tested rung, so this is usually one row). Returns data.frame(arm,
# rung); 0 rows when the arm never cached this plot. rung is the 4th "__"
# token with any param suffix and the .csv extension stripped.
cached_rungs <- function(dir, arm, site, plot) {
  none <- data.frame(arm = character(), rung = character(),
                     stringsAsFactors = FALSE)
  g <- basename(Sys.glob(file.path(dir, sprintf("%s__%s__%s__*.csv",
                                                arm, site, plot))))
  if (!length(g)) return(none)
  rung <- vapply(strsplit(sub("\\.csv$", "", g), "__", fixed = TRUE),
                 function(p) if (length(p) >= 4) p[4] else NA_character_,
                 character(1))
  rung <- unique(rung[!is.na(rung)])
  if (!length(rung)) return(none)
  data.frame(arm = arm, rung = rung, stringsAsFactors = FALSE)
}

## ---- optical boxes -> apex detections ---------------------------------------
# Converts an RGB crown-box table (x, y = box centre in UTM; optional score) to
# the (x, y, z) detection contract, mirroring detect_deepforest_sweep.R: apex z
# is read from the plot's native frozen CHM at the box centre and floored at
# z_floor (the detector min_height) where the box falls off the CHM. chm = NULL
# floors everything (callers should treat that as a degraded fallback).
boxes_to_dets <- function(boxes, chm, score_min = 0, z_floor = 2.0) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.null(boxes) || !nrow(boxes)) return(empty)
  keep <- is.finite(boxes$x) & is.finite(boxes$y)
  if (!is.null(boxes$score)) keep <- keep & boxes$score >= score_min
  b <- boxes[keep, , drop = FALSE]
  if (!nrow(b)) return(empty)
  z <- if (!is.null(chm))
    as.numeric(terra::extract(chm, cbind(b$x, b$y))[, 1]) else rep(NA_real_, nrow(b))
  z[!is.finite(z)] <- z_floor
  data.frame(x = b$x, y = b$y, z = z)
}
