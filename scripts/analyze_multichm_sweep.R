#!/usr/bin/env Rscript
.bs_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.bs_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bs <- Find(file.exists, c(
  if (!is.null(.bs_ofile) && length(.bs_ofile) && nzchar(.bs_ofile))
    file.path(dirname(.bs_ofile), "bootstrap.R"),
  if (length(.bs_file)) file.path(dirname(sub("^--file=", "", .bs_file[1])),
                                  "bootstrap.R"),
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs, .bs_ofile, .bs_file)

# Pool the multichm density-ladder arm (#37) and put it head-to-head with the
# cached CHM-VWF baseline (sweep_results.csv) on the SAME prepare_clip lasR path.
#
# The comparison is apples-to-apples by construction: for each cached CHM-VWF
# (plot, rung) we select the row whose chm_res equals that row's OWN density-
# derived value (0.25 if first-return density >= 8 else 0.5) and vwf_a == 0.10 --
# exactly the discipline multichm follows (it picks res from its own measured
# frdens, fixed a = 0.10). Both arms are then restricted to the COMMON (plot,
# rung) set (a 2-arm equal-set guard) before pooling, and every rate is pooled
# with the canonical model_bench_lib::pool (sum TP / sum n_ref; per-class TP
# recovered as round(rec_<cls> * n_<cls>); understory = intermediate +
# suppressed). Deltas are differences of POOLED rates, never a mean of per-row
# deltas.
#
# Usage:  Rscript scripts/analyze_multichm_sweep.R [SITES=SJER,SOAP,TEAK]
# Writes: work/neon/<SITE>/multichm_summary_by_rung.csv,
#         work/neon/<SITE>/figs/multichm_vs_chmvwf.png,
#         work/neon/multichm_vs_chmvwf.csv (cross-site roll-up),
#         work/neon/multichm_addendum.md (results-doc fragment).
suppressMessages({ library(data.table) })
d <- .job_dir()
source(.find("model_bench_lib.R"))   # pool (canonical sum-counts pooler)

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
RUNGS <- c("native", "8", "4", "2", "1")

## res-rule: the density-derived CHM resolution both arms use, from a row's OWN
## measured first-return density. multichm picks this internally; for the cached
## CHM-VWF sweep we pick the matching row with this same rule (and vwf_a = 0.10).
res_rule <- function(frdens) ifelse(frdens >= 8, 0.25, 0.5)

## Load one site's two arms, restricted to the common (plot, rung) set. Returns
## list(mc, cv) -- both with `site` + `tp_core`, the columns pool() consumes --
## or NULL if either CSV is missing.
load_site <- function(site) {
  nd <- file.path(d, "neon", site)
  fmc <- file.path(nd, "multichm_sweep_results.csv")
  fcv <- file.path(nd, "sweep_results.csv")
  if (!file.exists(fmc) || !file.exists(fcv)) {
    cat(sprintf("[%s] missing CSV (multichm:%s chmvwf:%s) -- skipped\n",
                site, file.exists(fmc), file.exists(fcv)))
    return(NULL)
  }
  mc <- read.csv(fmc, stringsAsFactors = FALSE)
  cv <- read.csv(fcv, stringsAsFactors = FALSE)
  mc$site <- site
  cv$site <- site
  # CHM-VWF comparison cell: same res-rule (from the CHM-VWF row's own frdens)
  # and a = 0.10, so one row per (plot, rung) on the same discipline as multichm.
  cv <- cv[cv$vwf_a == 0.10 & abs(cv$chm_res - res_rule(cv$frdens)) < 1e-9, ]
  if (is.null(mc$tp_core)) mc$tp_core <- round(mc$precision * mc$n_det)
  cv$tp_core <- round(cv$precision * cv$n_det)
  # 2-arm equal-set guard: keep only (plot, rung) cells scored by BOTH arms.
  kmc <- paste(mc$plot, mc$rung, sep = "::")
  kcv <- paste(cv$plot, cv$rung, sep = "::")
  common <- intersect(kmc, kcv)
  dropped <- length(union(kmc, kcv)) - length(common)
  if (dropped) cat(sprintf("[%s] equal-set: %d common (plot,rung); %d dropped\n",
                           site, length(common), dropped))
  list(mc = mc[kmc %in% common, ], cv = cv[kcv %in% common, ])
}

## Pool an arm per rung (canonical pooler); returns one row per present rung with
## the pooled headline + crown-class + height-band rates and mean frdens.
pool_by_rung <- function(df) {
  do.call(rbind, lapply(RUNGS, function(rl) {
    s <- df[df$rung == rl, ]
    if (!nrow(s)) return(NULL)
    p <- pool(s)
    cbind(data.frame(rung = rl, frdens = round(mean(s$frdens), 2)), p)
  }))
}

## ---- per-site pooling + head-to-head -------------------------------------
sites_loaded <- character(0)
xsite <- list()           # cross-site roll-up rows

for (site in SITES) {
  L <- load_site(site)
  if (is.null(L)) next
  sites_loaded <- c(sites_loaded, site)
  pm <- pool_by_rung(L$mc)    # multichm pooled by rung
  pc <- pool_by_rung(L$cv)    # CHM-VWF  pooled by rung

  cat(sprintf("\n=== %s : multichm pooled by rung (res density-derived, a=0.10) ===\n", site))
  print(pm[, c("rung", "frdens", "n_plots", "n_ref", "n_det", "recall",
               "precision", "F1", "rec_dominant", "rec_understory",
               "rec_h_tall", "rec_h_mid", "rec_h_short")],
        row.names = FALSE, digits = 2)

  # head-to-head Δ (multichm − CHM-VWF) per rung, on the common cells.
  rl_common <- intersect(pm$rung, pc$rung)
  delta <- do.call(rbind, lapply(rl_common, function(rl) {
    a <- pm[pm$rung == rl, ]; b <- pc[pc$rung == rl, ]
    data.frame(rung = rl, frdens = a$frdens, n_plots = a$n_plots,
               recall_mc = a$recall, recall_cv = b$recall,
               d_recall = a$recall - b$recall,
               prec_mc = a$precision, prec_cv = b$precision,
               d_prec = a$precision - b$precision,
               F1_mc = a$F1, F1_cv = b$F1, d_F1 = a$F1 - b$F1)
  }))
  delta <- delta[match(intersect(RUNGS, delta$rung), delta$rung), ]
  cat(sprintf("\n=== %s : head-to-head Δ (multichm − CHM-VWF), common cells ===\n", site))
  print(delta[, c("rung", "frdens", "n_plots", "recall_mc", "recall_cv", "d_recall",
                  "F1_mc", "F1_cv", "d_F1")], row.names = FALSE, digits = 2)

  # persist per-site pooled table (both arms, tidy long form)
  pm$arm <- "multichm"; pc$arm <- "chm_vwf"
  by_rung <- rbind(pm, pc)
  by_rung <- by_rung[, c("arm", setdiff(names(by_rung), "arm"))]
  nd <- file.path(d, "neon", site)
  write.csv(by_rung, file.path(nd, "multichm_summary_by_rung.csv"), row.names = FALSE)

  # figure: multichm vs CHM-VWF recall & F1 across measured first-return density.
  fig <- file.path(nd, "figs"); dir.create(fig, showWarnings = FALSE)
  ord <- order(delta$frdens)
  png(file.path(fig, "multichm_vs_chmvwf.png"), 1000, 750, res = 130)
  par(mar = c(4.2, 4.2, 2.5, 1))
  plot(delta$frdens[ord], delta$recall_mc[ord], type = "o", pch = 19, lwd = 2,
       col = "#1b9e77", ylim = c(0, 1), log = "x",
       xlab = "measured first-return density (pulses/m^2)", ylab = "score",
       main = sprintf("%s: multichm vs CHM-VWF (same lasR clip path)", site))
  lines(delta$frdens[ord], delta$recall_cv[ord], type = "o", pch = 1, lwd = 2, col = "#1b9e77")
  lines(delta$frdens[ord], delta$F1_mc[ord], type = "o", pch = 19, lwd = 2, col = "#d95f02")
  lines(delta$frdens[ord], delta$F1_cv[ord], type = "o", pch = 1, lwd = 2, col = "#d95f02")
  legend("topleft", bty = "n", cex = 0.85,
         legend = c("recall multichm", "recall CHM-VWF", "F1 multichm", "F1 CHM-VWF"),
         col = c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02"),
         pch = c(19, 1, 19, 1), lwd = 2)
  dev.off()

  # cross-site roll-up rows (per rung, both arms + delta)
  for (rl in delta$rung)
    xsite[[length(xsite) + 1]] <- cbind(site = site, delta[delta$rung == rl, ])
}

if (!length(sites_loaded)) stop("no sites with both arms found", call. = FALSE)

## ---- cross-site roll-up + markdown fragment ------------------------------
xs <- do.call(rbind, xsite)
write.csv(xs, file.path(d, "neon", "multichm_vs_chmvwf.csv"), row.names = FALSE)

cat("\n=== Cross-site roll-up (per site × rung; multichm vs CHM-VWF) ===\n")
print(xs[, c("site", "rung", "frdens", "recall_mc", "recall_cv", "d_recall",
             "F1_mc", "F1_cv", "d_F1")], row.names = FALSE, digits = 2)

# pooled-over-rungs per site (one Δ number per site, for the addendum TL;DR).
# overstory = dominant + codominant, recovered from the pooled per-class rates
# (pooled rec_<cls> = tp_<cls>/n_<cls>, so tp_<cls> = rec_<cls> * n_<cls>).
overstory <- function(p) {
  n <- p$n_dominant + p$n_codominant
  if (!n) return(NA_real_)
  (p$rec_dominant * p$n_dominant + p$rec_codominant * p$n_codominant) / n
}
cat("\n=== Per-site pooled over all common cells (multichm vs CHM-VWF) ===\n")
site_pool <- do.call(rbind, lapply(sites_loaded, function(site) {
  L <- load_site(site)
  pm <- pool(L$mc); pc <- pool(L$cv)
  data.frame(site = site, n_plots = pm$n_plots, n_ref = pm$n_ref,
             recall_mc = pm$recall, recall_cv = pc$recall, d_recall = pm$recall - pc$recall,
             prec_mc = pm$precision, prec_cv = pc$precision,
             over_mc = overstory(pm), over_cv = overstory(pc),
             und_mc = pm$rec_understory, und_cv = pc$rec_understory,
             F1_mc = pm$F1, F1_cv = pc$F1, d_F1 = pm$F1 - pc$F1)
}))
print(site_pool[, c("site", "n_plots", "n_ref", "recall_mc", "recall_cv",
                    "prec_mc", "prec_cv", "over_mc", "over_cv", "und_mc",
                    "und_cv", "F1_mc", "F1_cv", "d_F1")],
      row.names = FALSE, digits = 2)

## markdown fragment for the results-doc addendum (real numbers; rumdl 80-char).
fmt <- function(x) ifelse(is.na(x), "—", sprintf("%.2f", x))
md <- c(
  "<!-- generated by scripts/analyze_multichm_sweep.R -->",
  "",
  "Per-site pooled over all common (plot, rung) cells:",
  "",
  "| Site | plots | n_ref | recall (mc / vwf) | F1 (mc / vwf) | ΔF1 |",
  "|---|---|---|---|---|---|")
for (i in seq_len(nrow(site_pool))) {
  r <- site_pool[i, ]
  md <- c(md, sprintf("| %s | %d | %d | %s / %s | %s / %s | %+.2f |",
                      r$site, r$n_plots, r$n_ref, fmt(r$recall_mc), fmt(r$recall_cv),
                      fmt(r$F1_mc), fmt(r$F1_cv), r$d_F1))
}
md <- c(md, "", "Per site × density rung (recall / F1, multichm vs CHM-VWF):", "",
        "| Site | rung | frdens | recall mc | recall vwf | ΔF1 |",
        "|---|---|---|---|---|---|")
for (i in seq_len(nrow(xs))) {
  r <- xs[i, ]
  md <- c(md, sprintf("| %s | %s | %.1f | %s | %s | %+.2f |",
                      r$site, r$rung, r$frdens, fmt(r$recall_mc), fmt(r$recall_cv), r$d_F1))
}
writeLines(md, file.path(d, "neon", "multichm_addendum.md"))

cat(sprintf("\nwrote: per-site multichm_summary_by_rung.csv + figs/multichm_vs_chmvwf.png\n  %s\n  %s\n",
            file.path(d, "neon", "multichm_vs_chmvwf.csv"),
            file.path(d, "neon", "multichm_addendum.md")))
