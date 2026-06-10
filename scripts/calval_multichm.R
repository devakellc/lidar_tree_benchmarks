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

# Out-of-sample (calibration/validation) scrutiny of the multichm treetop arm
# (issue #39). The density-ladder §8 head-to-head (analyze_multichm_sweep.R ->
# density-ladder-sweep-results.md) found multichm beats CHM-VWF on F1 at the
# closed-canopy conifer sites (SOAP +0.05, TEAK +0.10 pooled), flat at open SJER
# -- but that comparison pools over ALL plots (in-sample). This script asks
# whether multichm's advantage survives on HELD-OUT plots, using the EXACT same
# stratified calib/valid split as the CHM-VWF #3 study (calval_split.R) so the
# two are directly comparable.
#
# multichm has NO tuning grid (single density-derived res 0.25/0.5 m by frdens,
# fixed window ws_factory(0.10)), so its "held-out" metric is simply its pooled
# rate on the validation plots -- there is no parameter to overfit. The split
# therefore tests robustness of the head-to-head to plot subsampling: on a
# stratified half of the plots that CHM-VWF's #3 test held out, is multichm still
# ahead? We score multichm against TWO CHM-VWF baselines on the same held-out
# plots:
#   * chm_vwf_matched -- density-derived res (0.25/0.5 by frdens) + a=0.10. This
#     is the §8 baseline the "multichm wins" claim was made against; the PRIMARY
#     comparison. Neither arm is tuned -> a clean apples-to-apples subsample test.
#   * chm_vwf_tuned   -- the per-rung best (chm_res, vwf_a) tuned on the
#     CALIBRATION plots (exactly as calval_split.R selects), applied to the
#     validation plots. This is the adversarial test: does multichm still win
#     when CHM-VWF is allowed to tune on calibration data and multichm is not?
#
# Every comparison is PAIRED + complete-case: at each rung only the validation
# plots scored by BOTH arms are pooled (a plot missing one arm's cell is dropped
# and counted), so a rung's Δ is never contaminated by an unequal plot
# population. Pooling is the canonical model_bench_lib::pool (sum TP / sum n_ref;
# precision = sum tp_core / sum n_det; per-class TP = round(rec*n); understory =
# intermediate + suppressed) -- the SAME pooler the in-sample §8 used, so
# in-sample and held-out numbers are produced identically.
#
# Because per-site plot counts are small (SJER 8, SOAP 18, TEAK 20 -> validation
# halves 3/8/9), a single split is noisy: SEEDS repeats the split and reports the
# held-out ΔF1 distribution and the fraction of seeds multichm wins, so the
# verdict is not a one-split artifact.
#
# Usage:
#   export CLAUDE_JOB_DIR=/path/to/work
#   Rscript scripts/calval_multichm.R SITES=SJER,SOAP,TEAK SEED=1 FRAC=0.5 \
#           SEEDS=1,2,3,4,5
#
# Reads (cached; no LiDAR recomputation):
#   $CLAUDE_JOB_DIR/neon/<SITE>/sweep_results.csv          (CHM-VWF grid)
#   $CLAUDE_JOB_DIR/neon/<SITE>/multichm_sweep_results.csv (multichm arm, #37)
# Writes:
#   $CLAUDE_JOB_DIR/neon/<SITE>/calval_multichm_metrics.csv (long form)
#   $CLAUDE_JOB_DIR/neon/calval_multichm_addendum.md        (results-doc fragment)

source(.find("calval_lib.R"))        # rung_lab, classes, plot_table, assign_split
source(.find("model_bench_lib.R"))   # pool (canonical sum-counts pooler)

## res-rule: the density-derived CHM resolution both arms use, from a row's OWN
## measured first-return density (identical to analyze_multichm_sweep.R).
res_rule <- function(frdens) ifelse(frdens >= 8, 0.25, 0.5)
fmt <- function(x, k = 2) ifelse(is.na(x), "  --", formatC(x, format = "f", digits = k))

## Ensure a frame carries the columns pool() consumes: a `site` column (pool keys
## n_plots on site::plot::rung) and tp_core (recovered per row if absent).
ensure_cols <- function(df, site) {
  if (is.null(df$site)) df$site <- site
  if (is.null(df$tp_core)) df$tp_core <- round(df$precision * df$n_det)
  df
}

## Pair two single-row-per-(plot,rung) arms on a plot subset: keep only the
## (plot, rung) cells scored by BOTH, so every pooled Δ uses an identical
## population. Returns the restricted mc/cv frames + the dropped-cell count.
pair_cells <- function(mc, cv, plots) {
  mc <- mc[mc$plot %in% plots, ]; cv <- cv[cv$plot %in% plots, ]
  km <- paste(mc$plot, mc$rung, sep = "::"); kc <- paste(cv$plot, cv$rung, sep = "::")
  common <- intersect(km, kc)
  list(mc = mc[km %in% common, ], cv = cv[kc %in% common, ],
       n_drop = length(union(km, kc)) - length(common))
}

## Per-rung pooled head-to-head from already-paired arms: one row per present
## rung with both arms' recall/precision/F1/understory + the deltas.
h2h_by_rung <- function(P) {
  do.call(rbind, lapply(rung_lab, function(rl) {
    a <- P$mc[P$mc$rung == rl, ]; b <- P$cv[P$cv$rung == rl, ]
    if (!nrow(a) || !nrow(b)) return(NULL)
    pa <- pool(a); pb <- pool(b)
    data.frame(rung = rl, n_pair = length(unique(a$plot)), n_ref = pa$n_ref,
               rec_mc = pa$recall, prec_mc = pa$precision, F1_mc = pa$F1,
               und_mc = pa$rec_understory,
               rec_cv = pb$recall, prec_cv = pb$precision, F1_cv = pb$F1,
               und_cv = pb$rec_understory,
               d_recall = pa$recall - pb$recall, d_F1 = pa$F1 - pb$F1)
  }))
}

## Pool both arms over ALL paired cells (collapsing rungs) -> one Δ per site.
h2h_pooled <- function(P) {
  pa <- pool(P$mc); pb <- pool(P$cv)
  data.frame(n_pair = length(unique(paste(P$mc$plot, P$mc$rung))), n_ref = pa$n_ref,
             rec_mc = pa$recall, prec_mc = pa$precision, F1_mc = pa$F1,
             und_mc = pa$rec_understory,
             rec_cv = pb$recall, prec_cv = pb$precision, F1_cv = pb$F1,
             und_cv = pb$rec_understory,
             d_recall = pa$recall - pb$recall, d_F1 = pa$F1 - pb$F1)
}

## Tune CHM-VWF: best (chm_res, vwf_a) per rung by pooled F1 on the calibration
## plots, with the hypothesis-neutral tie-break calval_split.R uses (finest
## chm_res, then lower vwf_a on equal F1). Returns one (rung, chm_res, vwf_a) row.
tune_cv <- function(cv, calib_plots) {
  do.call(rbind, lapply(rung_lab, function(rl) {
    s <- cv[cv$rung == rl & cv$plot %in% calib_plots, ]
    if (!nrow(s)) return(NULL)
    grid <- unique(s[, c("chm_res", "vwf_a")])
    pg <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      ss <- s[s$chm_res == grid$chm_res[i] & s$vwf_a == grid$vwf_a[i], ]
      data.frame(chm_res = grid$chm_res[i], vwf_a = grid$vwf_a[i], F1 = pool(ss)$F1)
    }))
    pg <- pg[order(-pg$F1, pg$chm_res, pg$vwf_a), ][1, ]
    data.frame(rung = rl, chm_res = pg$chm_res, vwf_a = pg$vwf_a)
  }))
}

## Build the single-row-per-(plot,rung) CHM-VWF frame for the tuned selection:
## for each rung take the cv rows at the calibration-selected (chm_res, vwf_a).
apply_tuned <- function(cv, sel) {
  do.call(rbind, lapply(rung_lab, function(rl) {
    p <- sel[sel$rung == rl, ]; if (!nrow(p)) return(NULL)
    cv[cv$rung == rl & cv$chm_res == p$chm_res & cv$vwf_a == p$vwf_a, ]
  }))
}

## long-form rows for the metrics CSV
mk_long <- function(P, site, seed, split, detector) {
  rows <- h2h_by_rung(P)
  if (is.null(rows) || !nrow(rows)) return(NULL)
  mc_l <- data.frame(site = site, seed = seed, split = split, detector = "multichm",
    rung = rows$rung, n_pair = rows$n_pair, n_ref = rows$n_ref,
    recall = rows$rec_mc, precision = rows$prec_mc, F1 = rows$F1_mc,
    rec_understory = rows$und_mc)
  cv_l <- data.frame(site = site, seed = seed, split = split, detector = detector,
    rung = rows$rung, n_pair = rows$n_pair, n_ref = rows$n_ref,
    recall = rows$rec_cv, precision = rows$prec_cv, F1 = rows$F1_cv,
    rec_understory = rows$und_cv)
  rbind(mc_l, cv_l)
}

## ==========================================================================
## Main (wrapped so the helpers above can be sourced for unit tests without
## running the analysis -- mirrors detect_multichm_sweep.R::run_main).
## ==========================================================================
run_main <- function() {
d     <- .job_dir()
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
SEED  <- if (is.null(A$SEED))  1L  else as.integer(A$SEED)
FRAC  <- if (is.null(A$FRAC))  0.5 else as.numeric(A$FRAC)
SEEDS <- if (is.null(A$SEEDS)) 1:5 else as.integer(strsplit(A$SEEDS, ",")[[1]])
all_long <- list()
md <- c("<!-- generated by scripts/calval_multichm.R -->", "")
# accumulate the cross-site held-out summary (one row per site)
hold_summary <- list()

for (SITE in SITES) {
  nd  <- file.path(d, "neon", SITE)
  fmc <- file.path(nd, "multichm_sweep_results.csv")
  fcv <- file.path(nd, "sweep_results.csv")
  if (!file.exists(fmc) || !file.exists(fcv)) {
    cat(sprintf("[skip] %s missing CSV (multichm:%s chmvwf:%s)\n",
                SITE, file.exists(fmc), file.exists(fcv))); next }
  mc <- ensure_cols(read.csv(fmc, stringsAsFactors = FALSE), SITE)
  cv <- ensure_cols(read.csv(fcv, stringsAsFactors = FALSE), SITE)
  # matched-discipline CHM-VWF cell: one row per (plot, rung) at a=0.10 and the
  # row's OWN density-derived res -- the §8 baseline.
  cvm <- cv[cv$vwf_a == 0.10 & abs(cv$chm_res - res_rule(cv$frdens)) < 1e-9, ]
  pt  <- plot_table(cv)

  cat(sprintf("\n############### SITE %s ###############\n", SITE))
  cat(sprintf("plots: %d (%s)\n", nrow(pt),
      paste(names(table(pt$plotType)), table(pt$plotType), sep = "=", collapse = ", ")))

  ## ---- in-sample reference (all plots), per rung + pooled ----------------
  Pall <- pair_cells(mc, cvm, pt$plot)
  ref_rung   <- h2h_by_rung(Pall)
  ref_pooled <- h2h_pooled(Pall)
  cat("\n=== IN-SAMPLE reference (all plots): multichm vs CHM-VWF[matched] ===\n")
  cat(sprintf("  pooled over all cells: F1 mc/vwf = %s / %s  (ΔF1 %+.2f), recall %s / %s\n",
      fmt(ref_pooled$F1_mc), fmt(ref_pooled$F1_cv), ref_pooled$d_F1,
      fmt(ref_pooled$rec_mc), fmt(ref_pooled$rec_cv)))

  ## ---- headline split (SEED) --------------------------------------------
  ptS  <- assign_split(pt, SEED, FRAC)
  cal  <- ptS$plot[ptS$split == "calib"]
  val  <- ptS$plot[ptS$split == "valid"]
  cat(sprintf("\n-- SEED=%d FRAC=%.2f -- calib n=%d, valid n=%d\n",
              SEED, FRAC, length(cal), length(val)))
  if (!length(val)) { cat("  [warn] no validation plots; held-out tables empty\n"); next }

  Pm  <- pair_cells(mc, cvm, val)                       # vs matched CHM-VWF
  sel <- tune_cv(cv, cal)                               # calib-tuned CHM-VWF
  cvt <- ensure_cols(apply_tuned(cv, sel), SITE)
  Pt  <- pair_cells(mc, cvt, val)                       # vs tuned CHM-VWF

  hm  <- h2h_by_rung(Pm)
  ht  <- h2h_by_rung(Pt)
  cat("\n=== HELD-OUT (validation plots): multichm vs CHM-VWF ===\n")
  cat("  [matched]=density-derived res+a0.10 (§8 baseline)  [tuned]=calib-tuned (chm_res,vwf_a)\n")
  cat(sprintf("%-7s %5s | %-22s | %-13s | %-13s | %-13s\n", "rung", "n_pr",
              "multichm rec/prec/F1", "vwf[m] r/F1", "vwf[t] r/F1", "ΔF1 m / t"))
  for (rl in rung_lab) {
    a <- hm[hm$rung == rl, ]; t <- ht[ht$rung == rl, ]
    if (!nrow(a)) next
    dft <- if (nrow(t)) sprintf("%+.2f", t$F1_mc - t$F1_cv) else "  -- "
    f1t <- if (nrow(t)) sprintf("%s/%s", fmt(t$rec_cv), fmt(t$F1_cv)) else "   --     "
    cat(sprintf("%-7s %5d | %s / %s / %s | %s/%s | %-12s | %+.2f / %s\n",
        rl, a$n_pair, fmt(a$rec_mc), fmt(a$prec_mc), fmt(a$F1_mc),
        fmt(a$rec_cv), fmt(a$F1_cv), f1t, a$d_F1, dft))
  }
  # Each ΔF1 must subtract within its OWN pairing: Pm and Pt run complete-case
  # pair_cells independently (vs cvm vs cvt), so their multichm cell sets can
  # differ by a plot. The tuned delta therefore uses ptp$F1_mc (multichm on the
  # tuned pairing), NOT pm$F1_mc -- mixing them would compare unequal populations.
  pm <- h2h_pooled(Pm); ptp <- h2h_pooled(Pt)
  cat(sprintf("\n  POOLED held-out: F1 mc=%s  vwf[matched]=%s (ΔF1 %+.2f)  vwf[tuned]=%s (ΔF1 %+.2f)\n",
      fmt(pm$F1_mc), fmt(pm$F1_cv), pm$d_F1, fmt(ptp$F1_cv), ptp$F1_mc - ptp$F1_cv))
  cat(sprintf("  POOLED held-out understory recall: mc=%s  vwf[matched]=%s\n",
      fmt(pm$und_mc), fmt(pm$und_cv)))

  all_long[[paste0(SITE, "_m")]] <- mk_long(Pm, SITE, SEED, "valid", "chm_vwf_matched")
  all_long[[paste0(SITE, "_t")]] <- mk_long(Pt, SITE, SEED, "valid", "chm_vwf_tuned")

  ## ---- multi-seed robustness --------------------------------------------
  cat(sprintf("\n=== multi-seed robustness (SEEDS=%s, FRAC=%.2f) ===\n",
              paste(SEEDS, collapse = ","), FRAC))
  rb <- do.call(rbind, lapply(SEEDS, function(sd) {
    p   <- assign_split(pt, sd, FRAC)
    va  <- p$plot[p$split == "valid"]; ca <- p$plot[p$split == "calib"]
    if (!length(va)) return(NULL)
    Pms <- pair_cells(mc, cvm, va)
    sl  <- tune_cv(cv, ca)
    Pts <- pair_cells(mc, ensure_cols(apply_tuned(cv, sl), SITE), va)
    hr  <- h2h_by_rung(Pms); tr <- h2h_by_rung(Pts)
    do.call(rbind, lapply(rung_lab, function(rl) {
      a <- hr[hr$rung == rl, ]; tt <- tr[tr$rung == rl, ]
      if (!nrow(a)) return(NULL)
      data.frame(seed = sd, rung = rl, F1_mc = a$F1_mc, d_F1_m = a$d_F1,
                 d_F1_t = if (nrow(tt)) tt$F1_mc - tt$F1_cv else NA_real_) })) }))
  agg <- do.call(rbind, lapply(rung_lab, function(rl) {
    z <- rb[rb$rung == rl, ]; if (!nrow(z)) return(NULL)
    data.frame(rung = rl, n_seeds = nrow(z),
               F1mc_med = median(z$F1_mc), F1mc_min = min(z$F1_mc), F1mc_max = max(z$F1_mc),
               dF1m_med = median(z$d_F1_m), win_m = mean(z$d_F1_m > 0),
               dF1t_med = median(z$d_F1_t, na.rm = TRUE),
               win_t = mean(z$d_F1_t > 0, na.rm = TRUE)) }))
  cat(sprintf("%-7s %14s %18s %18s\n", "rung", "mc F1 med[min..max]",
              "ΔF1[matched] med/win", "ΔF1[tuned] med/win"))
  for (rl in rung_lab) { a <- agg[agg$rung == rl, ]; if (!nrow(a)) next
    cat(sprintf("%-7s %s [%s..%s] %+9.2f / %3.0f%% %+9.2f / %3.0f%%\n",
        rl, fmt(a$F1mc_med), fmt(a$F1mc_min), fmt(a$F1mc_max),
        a$dF1m_med, 100 * a$win_m, a$dF1t_med, 100 * a$win_t)) }

  # per-site held-out verdict (pooled over rungs across seeds: mean ΔF1 + win)
  win_m_all <- mean(rb$d_F1_m > 0); win_t_all <- mean(rb$d_F1_t > 0, na.rm = TRUE)
  cat(sprintf("\n  VERDICT %s: held-out ΔF1[matched] median=%+.2f (multichm wins %.0f%% of seed×rung);\n",
              SITE, median(rb$d_F1_m), 100 * win_m_all))
  cat(sprintf("           held-out ΔF1[tuned]   median=%+.2f (multichm wins %.0f%%); in-sample ΔF1=%+.2f\n",
              median(rb$d_F1_t, na.rm = TRUE), 100 * win_t_all, ref_pooled$d_F1))

  hold_summary[[SITE]] <- data.frame(site = SITE,
    insample_dF1 = ref_pooled$d_F1, heldout_dF1_matched = median(rb$d_F1_m),
    win_matched = win_m_all, heldout_dF1_tuned = median(rb$d_F1_t, na.rm = TRUE),
    win_tuned = win_t_all, heldout_F1_mc = pm$F1_mc, heldout_und_mc = pm$und_mc,
    heldout_und_cv = pm$und_cv)

  ## ---- markdown fragment for the addendum -------------------------------
  md <- c(md, sprintf("### %s held-out (SEED=%d, FRAC=%.2f; valid n=%d)",
                      SITE, SEED, FRAC, length(val)), "",
          "| rung | n_pair | multichm rec/prec/F1 | vwf[matched] rec/F1 | ΔF1[matched] | vwf[tuned] F1 | ΔF1[tuned] |",
          "|---|---|---|---|---|---|---|")
  for (rl in rung_lab) { a <- hm[hm$rung == rl, ]; t <- ht[ht$rung == rl, ]
    if (!nrow(a)) next
    md <- c(md, sprintf("| %s | %d | %s / %s / %s | %s / %s | %+.2f | %s | %s |",
      rl, a$n_pair, fmt(a$rec_mc), fmt(a$prec_mc), fmt(a$F1_mc),
      fmt(a$rec_cv), fmt(a$F1_cv), a$d_F1,
      if (nrow(t)) fmt(t$F1_cv) else "--",
      if (nrow(t)) sprintf("%+.2f", t$F1_mc - t$F1_cv) else "--")) }
  md <- c(md, "",
          sprintf("Pooled held-out: multichm F1 %s vs vwf[matched] %s (ΔF1 %+.2f), vwf[tuned] %s (ΔF1 %+.2f); understory recall mc %s vs vwf %s.",
                  fmt(pm$F1_mc), fmt(pm$F1_cv), pm$d_F1, fmt(ptp$F1_cv),
                  ptp$F1_mc - ptp$F1_cv, fmt(pm$und_mc), fmt(pm$und_cv)), "")

  rows_long <- do.call(rbind, list(all_long[[paste0(SITE, "_m")]],
                                   all_long[[paste0(SITE, "_t")]]))
  write.csv(rows_long, file.path(nd, "calval_multichm_metrics.csv"), row.names = FALSE)
}

if (length(hold_summary)) {
  hs <- do.call(rbind, hold_summary)
  cat("\n==================== CROSS-SITE HELD-OUT SUMMARY ====================\n")
  print(hs, row.names = FALSE, digits = 2)
  md <- c("# multichm cal/val (issue #39) — generated fragment", "",
          "Per-site: in-sample vs held-out ΔF1 (multichm − CHM-VWF), and the",
          "fraction of seed×rung splits multichm wins held-out.", "",
          "| Site | in-sample ΔF1 | held-out ΔF1[matched] | win% | held-out ΔF1[tuned] | win% |",
          "|---|---|---|---|---|---|",
          unlist(lapply(seq_len(nrow(hs)), function(i) { r <- hs[i, ]
            sprintf("| %s | %+.2f | %+.2f | %.0f%% | %+.2f | %.0f%% |", r$site,
                    r$insample_dF1, r$heldout_dF1_matched, 100 * r$win_matched,
                    r$heldout_dF1_tuned, 100 * r$win_tuned) })), "", md)
  writeLines(md, file.path(d, "neon", "calval_multichm_addendum.md"))
  cat(sprintf("\nwrote -> %s\n", file.path(d, "neon", "calval_multichm_addendum.md")))
}
cat("\n==================== DONE ====================\n")
}

if (sys.nframe() == 0L) run_main()
