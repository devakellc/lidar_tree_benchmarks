#!/usr/bin/env Rscript
# #P2 routing helpers (pure data.frame logic; no I/O, no heavy deps), so the
# per-cell detector-selection policy is unit-testable with no NEON data. The
# driver route_detectors.R sources this plus model_bench_lib.R (for pool /
# equal_set_guard) and lidR (for the deploy-time structure features).
#
# Vocabulary: a "cell" is one (site, plot, rung). `long` is the per-cell x per-arm
# scored table (one row per arm per cell, carrying F1 + the pool() count columns).
# A "policy" assigns one arm to each cell; select_policy_rows materializes the
# chosen rows so pool() can score the routed selection by the canonical summed-
# count rule (never by averaging per-cell F1).

CELL_KEYS <- c("site", "plot", "rung")

## ---- per-cell argmax-metric arm (the oracle / label) ---------------------
# For each cell, the arm with the maximum `value` (default F1). Ties resolve to
# the first arm in row order (deterministic). Returns data.frame(site,plot,rung,
# arm) -- one row per cell. A cell whose arms are all NA on `value` is dropped.
oracle_pick <- function(long, value = "F1", keys = CELL_KEYS) {
  long <- as.data.frame(long)
  k <- do.call(paste, c(long[keys], sep = "\r"))
  v <- long[[value]]
  ord <- order(k, -ifelse(is.finite(v), v, -Inf), seq_along(v))   # stable, NA last
  best <- long[ord, , drop = FALSE][!duplicated(k[ord]), , drop = FALSE]
  best <- best[is.finite(best[[value]]), , drop = FALSE]          # drop all-NA cells
  out <- best[, c(keys, "arm"), drop = FALSE]
  rownames(out) <- NULL
  out
}

## ---- materialize a routing policy's chosen rows --------------------------
# Inner-join `long` to `picks` (cell keys + arm) so the result is exactly the
# rows the policy selects -- one per cell that has a matching (cell, arm) row in
# long. A pick naming an arm absent from that cell contributes no row (the cell
# is simply unscored, never fabricated). Order follows `picks`.
select_policy_rows <- function(long, picks, keys = CELL_KEYS) {
  long <- as.data.frame(long); picks <- as.data.frame(picks)
  lk <- do.call(paste, c(long[c(keys, "arm")], sep = "\r"))
  pk <- do.call(paste, c(picks[c(keys, "arm")], sep = "\r"))
  rows <- long[match(pk, lk), , drop = FALSE]
  rows <- rows[!is.na(rows[[keys[1]]]), , drop = FALSE]            # drop unmatched picks
  rownames(rows) <- NULL
  rows
}

## ---- leave-one-plot-out fold ids -----------------------------------------
# Honest CV for a per-cell router: every cell of a plot is held out together (a
# plot's rungs are not independent). Returns an integer fold id per row, one
# distinct id per unique plot. Used by route_detectors.R to predict each plot's
# cells from a model trained on the others.
lopo_folds <- function(plot_vec) as.integer(factor(plot_vec))
