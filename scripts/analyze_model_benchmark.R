#!/usr/bin/env Rscript
# Cross-model density-ladder synthesis (#R10). Unifies every arm scored on the
# shared frozen SOAP clips (AMS3D, lmfauto, multichm, ptrees, CHM-VWF, and the
# native-only Li 2012) into per-class + per-height-band recall tables, density-
# robustness curves, and head-to-head deltas vs CHM-VWF. Pools with the canonical
# pool() (sum counts, never average rates) over equal_set_guard common cells.
#
# Usage:  Rscript scripts/analyze_model_benchmark.R [SITE=SOAP]
# Reads:  $CLAUDE_JOB_DIR/neon/<SITE>/{ams3d,lidrplugins,li2012}_results.csv
# Writes: summary CSVs + figs/ + a generated markdown table fragment under the
#         same dir; the narrative results/model-benchmark-results.md is authored.
suppressMessages({ library(data.table) })
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))

RUNG_LEVELS <- c("native", "8", "4", "2", "1")

# Union arm data.frames with differing schemas; recompute tp_core uniformly so
# pooling is identical across arms regardless of which CSV carried the column.
harmonize_union <- function(dfs) {
  u <- as.data.frame(data.table::rbindlist(dfs, fill = TRUE, use.names = TRUE))
  u$tp_core <- round(u$precision * u$n_det)
  u$rung <- factor(u$rung, levels = RUNG_LEVELS)
  u
}

# Restrict to `arms` and `rungs`, apply equal_set_guard per (site,plot,rung), then
# pool() each (detector,rung). Carries frdens (mean over the cell's plots) as the
# density-curve x. Returns a long table; attr "dropped" = guard's dropped cells.
pool_arms <- function(u, arms, rungs = RUNG_LEVELS) {
  sub <- u[u$detector %in% arms & as.character(u$rung) %in% rungs, ]
  sub$rung <- as.character(sub$rung)
  g <- equal_set_guard(sub, arms = arms)
  out <- list()
  for (a in arms) for (rl in rungs) {
    s <- g[g$detector == a & g$rung == rl, ]
    if (!nrow(s)) next
    out[[length(out) + 1]] <- cbind(
      detector = a, rung = rl, frdens = round(mean(s$frdens), 2), pool(s))
  }
  res <- do.call(rbind, out)
  if (!is.null(res)) res$rung <- factor(res$rung, levels = rungs)
  attr(res, "dropped") <- attr(g, "dropped")
  res
}

# Per rung, each non-baseline arm's recall/F1/understory minus the baseline arm's.
deltas_vs_baseline <- function(pooled, baseline = "chm_vwf") {
  base <- pooled[pooled$detector == baseline, ]
  others <- pooled[pooled$detector != baseline, ]
  m <- merge(others, base[, c("rung", "recall", "F1", "rec_understory")],
             by = "rung", suffixes = c("", "_base"))
  data.frame(detector = m$detector, rung = m$rung,
             d_recall = m$recall - m$recall_base,
             d_F1 = m$F1 - m$F1_base,
             d_understory = m$rec_understory - m$rec_understory_base)
}
