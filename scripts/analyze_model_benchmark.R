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

LADDER_ARMS <- c("ams3d", "lmfauto", "multichm", "ptrees", "chm_vwf")
NATIVE_ARMS <- c("chm_vwf", "ptrees", "ams3d", "li2012")

# Render a long pooled table to a GitHub-markdown block (selected columns).
# Format numeric columns BEFORE apply() (which would coerce the whole frame to a
# character matrix and silently drop the rounding); integer-valued columns
# (counts) stay integers, rates get `digits` decimals.
.md_table <- function(df, cols, digits = 2) {
  d2 <- df[, cols, drop = FALSE]
  for (j in seq_along(d2)) if (is.numeric(d2[[j]])) {
    col <- d2[[j]]
    d2[[j]] <- if (all(col == round(col), na.rm = TRUE)) as.character(col)
               else formatC(col, format = "f", digits = digits)
  }
  hdr  <- paste0("| ", paste(cols, collapse = " | "), " |")
  sep  <- paste0("|", paste(rep(" --- ", length(cols)), collapse = "|"), "|")
  body <- apply(d2, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(hdr, sep, body), collapse = "\n")
}

run_main <- function() {
  args <- strsplit(commandArgs(TRUE), "=")
  A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
  SITE <- if (is.null(A$SITE)) "SOAP" else A$SITE
  d    <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
  nd   <- file.path(d, "neon", SITE)

  arm_files <- c(ams3d = "ams3d_results.csv", lidrplugins = "lidrplugins_results.csv",
                 li2012 = "li2012_results.csv")
  dfs <- lapply(file.path(nd, arm_files), function(f)
    if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL)
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) stop("no arm result CSVs found under ", nd)
  u <- harmonize_union(dfs)
  cat("arms present:", paste(sort(unique(u$detector)), collapse = ", "), "\n")

  ladder <- pool_arms(u, arms = intersect(LADDER_ARMS, unique(u$detector)))
  native <- pool_arms(u, arms = intersect(NATIVE_ARMS, unique(u$detector)),
                      rungs = "native")
  dl <- deltas_vs_baseline(ladder, baseline = "chm_vwf")

  for (nm in c("ladder", "native", "dl"))
    write.csv(get(nm), file.path(nd, paste0("model_bench_", nm, ".csv")),
              row.names = FALSE)
  if (length(attr(ladder, "dropped")))
    cat(sprintf("ladder guard dropped %d (plot,rung) cells\n",
                length(attr(ladder, "dropped"))))

  ## figures: recall and understory recall vs first-return density, per arm
  fig <- file.path(nd, "figs"); dir.create(fig, showWarnings = FALSE)
  arms <- intersect(LADDER_ARMS, unique(ladder$detector))
  pal  <- setNames(c("#1b9e77","#d95f02","#7570b3","#e7298a","#666666")[seq_along(arms)], arms)
  draw <- function(file, ycol, ylab, main) {
    png(file.path(fig, file), 1000, 750, res = 130); on.exit(dev.off())
    par(mar = c(4.2, 4.2, 2.5, 1)); first <- TRUE
    for (a in arms) {
      s <- ladder[ladder$detector == a, ]; s <- s[order(s$frdens), ]
      if (!nrow(s)) next
      if (first) { plot(s$frdens, s[[ycol]], type = "o", pch = 19, ylim = c(0, 1),
                        log = "x", xlab = "first-return density (pulses/m^2)",
                        ylab = ylab, main = main, col = pal[a], lwd = 2); first <- FALSE
      } else lines(s$frdens, s[[ycol]], type = "o", pch = 19, col = pal[a], lwd = 2)
    }
    legend("topleft", arms, col = pal[arms], lwd = 2, pch = 19, bty = "n", cex = 0.85)
  }
  draw("model_recall_vs_density.png", "recall", "recall (overall)",
       "SOAP cross-model recall vs density")
  draw("model_understory_vs_density.png", "rec_understory", "understory recall",
       "SOAP understory recall vs density")

  ## generated markdown fragment (tables only; narrative is authored separately)
  cols_l <- c("detector","rung","frdens","n_plots","n_ref","recall","precision",
              "F1","rec_dominant","rec_codominant","rec_understory")
  cols_h <- c("detector","rung","rec_h_tall","n_h_tall","rec_h_mid","n_h_mid",
              "rec_h_short","n_h_short")
  cols_n <- c("detector","n_plots","n_ref","recall","precision","F1",
              "rec_understory","n_understory")
  cols_d <- c("detector","rung","d_recall","d_F1","d_understory")
  frag <- c("<!-- generated by scripts/analyze_model_benchmark.R; do not edit by hand -->",
            "", "#### Density ladder, per crown class (pooled, equal-set)", "",
            .md_table(ladder[order(ladder$detector, ladder$rung), ], cols_l), "",
            "#### Density ladder, per height band", "",
            .md_table(ladder[order(ladder$detector, ladder$rung), ], cols_h), "",
            "#### Native point-segmenter head-to-head", "",
            .md_table(native[order(native$detector), ], cols_n), "",
            "#### Head-to-head deltas vs CHM-VWF (recall/F1/understory)", "",
            .md_table(dl[order(dl$detector, dl$rung), ], cols_d))
  writeLines(frag, file.path(nd, "model_bench_tables.md"))
  cat(sprintf("synthesis -> %s {ladder,native,dl}.csv, figs/, model_bench_tables.md\n", nd))
}

if (sys.nframe() == 0L) run_main()
