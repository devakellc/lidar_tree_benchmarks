# Shared calibration/validation split helpers for the density-ladder sweep.
#
# These define the *plot universe* and the deterministic, stratified
# calibration/validation assignment used by BOTH the CHM-VWF cal/val study
# (calval_split.R, issue #3) and the multichm cal/val addendum (calval_multichm.R,
# issue #39). They live here as the SINGLE SOURCE OF TRUTH for the split so the
# two studies can never drift: the multichm head-to-head must be scored on the
# exact same held-out plots as the CHM-VWF #3 test for the comparison to be fair.
#
# Nothing here reads LiDAR or pools metrics -- it is pure data-frame logic over a
# cached long-form sweep CSV. Pooling lives in each caller (calval_split.R uses
# its own count-summing pool(); calval_multichm.R uses model_bench_lib::pool).
#
# Sourced, not run. No side effects.

# Density rungs (top to bottom) and NEON crown classes, in canonical order.
rung_lab <- c("native", "8", "4", "2", "1")
classes  <- c("dominant", "codominant", "intermediate", "suppressed")

## ---- per-plot stratification covariates ----------------------------------
# Each plot is one row of a sweep parameter combo; crown-class counts and
# plotType are identical across combos for a plot (they describe the ground
# truth), so one row per plot suffices and captures EVERY plot even if a plot is
# missing a specific (chm_res, vwf_a) slice in the grid. `r` is a cached
# long-form sweep data.frame carrying plot, plotType, and n_<class> columns.
plot_table <- function(r) {
  ref <- r[!duplicated(r$plot), ]
  os <- ref$n_dominant + ref$n_codominant             # overstory stems
  us <- ref$n_intermediate + ref$n_suppressed         # understory stems
  data.frame(
    plot     = ref$plot,
    plotType = ref$plotType,
    n_over   = os, n_under = us,
    # understory-present iff understory stems > 0 AND >= a nontrivial share
    # (us >= os/4), else overstory. (At SJER understory is near-absent.)
    crown_mix = ifelse(us > 0 & us >= os / 4, "understory", "overstory"),
    stringsAsFactors = FALSE)
}

## ---- deterministic stratified calib/valid assignment ---------------------
# Strata = plotType x crown_mix. Within each stratum a deterministic shuffle
# (set.seed) then take ceil(FRAC * n) to calibration. A singleton stratum goes to
# calibration (so validation never silently borrows tuning data from an
# unrepresented stratum -- a caveat for tiny sites). Returns `pt` with $stratum
# and $split ("calib"/"valid") columns added.
assign_split <- function(pt, seed, frac) {
  set.seed(seed)
  pt$stratum <- paste(pt$plotType, pt$crown_mix, sep = "/")
  pt$split <- NA_character_
  for (st in unique(pt$stratum)) {
    idx <- which(pt$stratum == st)
    idx <- idx[sample.int(length(idx))]          # deterministic shuffle
    n_cal <- max(1L, ceiling(frac * length(idx)))
    if (length(idx) == 1L) n_cal <- 1L            # singleton -> calibration
    pt$split[idx[seq_len(n_cal)]] <- "calib"
    if (n_cal < length(idx)) pt$split[idx[(n_cal + 1):length(idx)]] <- "valid"
  }
  pt
}
