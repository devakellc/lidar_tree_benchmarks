#!/usr/bin/env Rscript
# Calibration / validation split for the NEON density-ladder sweep (issue #3).
#
# The sweep's "best (chm_res, vwf_a) per density rung" table (analyze_sweep.R,
# section 2) is an IN-SAMPLE optimum: parameters are picked on the SAME plots
# they are scored on, pooled over every plot. This script asks whether that
# choice survives out-of-sample. It splits plots into a calibration subset and a
# held-out validation subset, tunes (chm_res, vwf_a) per density rung on the
# calibration plots ONLY (by pooled F1), then reports recall/precision/F1 on the
# validation plots with those calibration-selected parameters. The full pooled
# in-sample optimum is shown alongside for reference.
#
# The split is deterministic (set.seed(SEED)) and STRATIFIED by
#   (a) plotType  : tower vs distributed (different mapped extent / structure)
#   (b) crown-class mix : overstory-dominated vs understory-present, from
#       n_dominant+n_codominant vs n_intermediate+n_suppressed.
# Both strata are balanced across calib/valid. Per-site plot counts are small
# (SJER 8, SOAP 18, TEAK 20), so a single split is noisy: SEEDS runs several
# seeds and reports the held-out F1 distribution and the modal calibration
# optimum per rung, so the verdict is not a one-split artifact.
#
# Pooling rule is copied from analyze_sweep.R (the pool() function): pooled
# recall = sum(TP)/sum(n_ref); pooled precision = sum(tp_core)/sum(n_det) with
# tp_core = round(precision*n_det); per-class TP = round(rec_class*n_class).
# DO NOT average per-plot rates -- that over-weights small plots.
#
# Usage:
#   export CLAUDE_JOB_DIR=/path/to/work
#   Rscript scripts/calval_split.R SITES=SJER,SOAP,TEAK SEED=1 FRAC=0.5
#   Rscript scripts/calval_split.R SITES=SOAP SEED=1 FRAC=0.5 SEEDS=1,2,3,4,5
#
# Writes (NEW file, never overwrites sweep_results.csv):
#   $CLAUDE_JOB_DIR/neon/<SITE>/calval_metrics.csv
#     long form: site, rung, split, chm_res, vwf_a, n_plots, n_ref, n_det,
#                TP, recall, precision, F1
#
# See analyze_sweep.R section 5 for the run pointer. Self-contained: the pooling
# logic lives here, not duplicated into analyze_sweep.R.

args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER","SOAP","TEAK") else strsplit(A$SITES, ",")[[1]]
SEED  <- if (is.null(A$SEED)) 1L else as.integer(A$SEED)
FRAC  <- if (is.null(A$FRAC)) 0.5 else as.numeric(A$FRAC)   # calibration fraction
SEEDS <- if (is.null(A$SEEDS)) 1:5 else as.integer(strsplit(A$SEEDS, ",")[[1]])

d  <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
rung_lab <- c("native","8","4","2","1")
classes  <- c("dominant","codominant","intermediate","suppressed")

## ---- pooling (copied/adapted from analyze_sweep.R::pool) ------------------
# Adds the per-plot true-positive and per-class TP columns the pooler needs.
add_tp_cols <- function(r) {
  r$tp_core <- round(r$precision * r$n_det)
  for (cl in classes) {
    nc <- r[[paste0("n_", cl)]]; rc <- r[[paste0("rec_", cl)]]
    r[[paste0("tp_", cl)]] <- ifelse(nc > 0, round(rc * nc), 0)
  }
  r
}
pool <- function(df) {
  out <- data.frame(
    n_plots = length(unique(df$plot)),
    n_ref   = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall    = sum(df$TP) / sum(df$n_ref),
    precision = sum(df$tp_core, na.rm = TRUE) / sum(df$n_det))
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  out
}

## ---- per-plot stratification covariates ----------------------------------
# Each plot is one row of a sweep parameter combo; crown-class counts are
# identical across combos for a plot (they describe the ground truth), so we
# read them from the native / res=0.5 / a=0.10 slice for a stable per-plot view.
plot_table <- function(r) {
  # crown-class counts + plotType are constant per plot across all sweep combos
  # (verified empirically), so one row per plot suffices and captures EVERY plot
  # even if a plot is missing a specific (chm_res, vwf_a) slice in the grid.
  ref <- r[!duplicated(r$plot), ]
  os <- ref$n_dominant + ref$n_codominant            # overstory stems
  us <- ref$n_intermediate + ref$n_suppressed         # understory stems
  data.frame(
    plot     = ref$plot,
    plotType = ref$plotType,
    n_over   = os, n_under = us,
    # understory-present iff understory stems strictly outnumber-or-equal a
    # nontrivial share: classify as "understory" when us > 0 AND us >= os/4,
    # else "overstory". (At SJER understory is near-absent -- handled below.)
    crown_mix = ifelse(us > 0 & us >= os / 4, "understory", "overstory"),
    stringsAsFactors = FALSE)
}

## ---- deterministic stratified calib/valid assignment ---------------------
# Strata = plotType x crown_mix. Within each stratum, a deterministic shuffle
# (set.seed) then take ceil(FRAC * n) to calibration. If a stratum has a single
# plot it goes to calibration (so validation never silently borrows tuning data
# from an unrepresented stratum -- documented as a caveat for tiny sites).
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

## ---- pick best (chm_res, vwf_a) per rung by pooled F1 on a plot subset ----
best_per_rung <- function(r, plots) {
  s0 <- r[r$plot %in% plots, ]
  do.call(rbind, lapply(rung_lab, function(rl) {
    s <- s0[s0$rung == rl, ]; if (!nrow(s)) return(NULL)
    grid <- unique(s[, c("chm_res","vwf_a")])
    pg <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
      ss <- s[s$chm_res == grid$chm_res[i] & s$vwf_a == grid$vwf_a[i], ]
      cbind(chm_res = grid$chm_res[i], vwf_a = grid$vwf_a[i], pool(ss)) }))
    # tie-break (hypothesis-NEUTRAL, deterministic): on equal F1 prefer the
    # FINEST chm_res (ascending chm_res), then the lower vwf_a. This does not
    # privilege any particular resolution -- it simply favours the finer CHM and
    # gentler window, so the selection cannot be accused of steering toward the
    # value the doc claims is optimal.
    pg <- pg[order(-pg$F1, pg$chm_res, pg$vwf_a), ][1, ]
    cbind(rung = rl, pg[, c("chm_res","vwf_a","F1")]) }))
}

## ---- pooled metrics for a fixed (chm_res, vwf_a) per rung on a subset -----
apply_params <- function(r, plots, sel) {
  s0 <- r[r$plot %in% plots, ]
  do.call(rbind, lapply(rung_lab, function(rl) {
    p <- sel[sel$rung == rl, ]
    ss <- s0[s0$rung == rl & s0$chm_res == p$chm_res & s0$vwf_a == p$vwf_a, ]
    if (!nrow(ss)) return(NULL)
    cbind(rung = rl, chm_res = p$chm_res, vwf_a = p$vwf_a, pool(ss)) }))
}

## ---- VWF-slope spread at the calib-optimal chm_res, on a held-out subset --
# At each rung, fix chm_res to the CALIBRATION optimum (chm_res comes from the
# calibration plots), then pool F1 across the three vwf_a values on the supplied
# `plots` subset; the spread (max-min) quantifies the VWF main effect. To make
# the VWF claim genuinely OUT-OF-SAMPLE this is called on the VALIDATION plots.
# If a rung's validation subset has no rows (too thin to pool), n_pool is 0 and
# the spread is NA -- reported honestly rather than falling back to calibration.
vwf_spread <- function(r, plots, sel) {
  s0 <- r[r$plot %in% plots, ]
  do.call(rbind, lapply(rung_lab, function(rl) {
    p  <- sel[sel$rung == rl, ]
    s  <- s0[s0$rung == rl & s0$chm_res == p$chm_res, ]
    if (!nrow(s)) return(data.frame(rung = rl, chm_res = p$chm_res,
               f1_min = NA_real_, f1_max = NA_real_, spread = NA_real_,
               n_pool = 0L))
    f1 <- sapply(sort(unique(s$vwf_a)), function(av) pool(s[s$vwf_a == av, ])$F1)
    data.frame(rung = rl, chm_res = p$chm_res,
               f1_min = min(f1), f1_max = max(f1), spread = max(f1) - min(f1),
               n_pool = length(unique(s$plot))) }))
}

fmt <- function(x, k = 3) formatC(x, format = "f", digits = k)

## ==========================================================================
## Main: per-site single-seed report + multi-seed robustness
## ==========================================================================
all_long <- list()       # rows for calval_metrics.csv
robust    <- list()      # per-site multi-seed aggregation for the doc

for (SITE in SITES) {
  nd <- file.path(d, "neon", SITE)
  fp <- file.path(nd, "sweep_results.csv")
  if (!file.exists(fp)) { cat(sprintf("[skip] no sweep_results for %s\n", SITE)); next }
  r  <- add_tp_cols(read.csv(fp, stringsAsFactors = FALSE))
  pt <- plot_table(r)

  cat(sprintf("\n############### SITE %s ###############\n", SITE))
  cat(sprintf("plots: %d (%s)\n", nrow(pt),
      paste(names(table(pt$plotType)), table(pt$plotType), sep="=", collapse=", ")))
  cat("crown-mix x plotType strata (all plots):\n")
  print(table(pt$plotType, pt$crown_mix))

  ## ---- headline single-seed split (SEED arg) ----
  ptS <- assign_split(pt, SEED, FRAC)
  cal <- ptS$plot[ptS$split == "calib"]
  val <- ptS$plot[ptS$split == "valid"]
  cat(sprintf("\n-- SEED=%d FRAC=%.2f --  calib n=%d, valid n=%d\n",
              SEED, FRAC, length(cal), length(val)))
  cat("  calib strata:\n"); print(table(ptS$plotType[ptS$split=="calib"],
                                        ptS$crown_mix[ptS$split=="calib"]))
  cat("  valid strata:\n"); print(table(ptS$plotType[ptS$split=="valid"],
                                        ptS$crown_mix[ptS$split=="valid"]))
  if (!length(val)) cat("  [warn] no validation plots at this FRAC/seed; held-out table empty\n")

  selC  <- best_per_rung(r, cal)                 # calibration-selected params
  selF  <- best_per_rung(r, pt$plot)             # full pooled in-sample reference
  m_cal <- apply_params(r, cal, selC)            # in-sample (calibration) metrics
  m_val <- if (length(val)) apply_params(r, val, selC) else NULL  # held-out
  m_full<- apply_params(r, pt$plot, selF)        # full pooled in-sample optimum

  cat("\n=== calibration-selected (chm_res, vwf_a) per rung; metrics by split ===\n")
  cat("    [CAL]=calibration plots  [VAL]=held-out validation plots\n")
  cat(sprintf("%-7s %7s %6s | %-22s | %-22s\n",
              "rung","chm_res","vwf_a","CAL  rec / prec / F1","VAL  rec / prec / F1"))
  for (rl in rung_lab) {
    c1 <- m_cal[m_cal$rung == rl, ]; v1 <- if (!is.null(m_val)) m_val[m_val$rung==rl,] else NULL
    if (!nrow(c1)) next
    vs <- if (!is.null(v1) && nrow(v1)) {
      sprintf("%s / %s / %s", fmt(v1$recall,2), fmt(v1$precision,2), fmt(v1$F1,2))
    } else "   --       "
    cat(sprintf("%-7s %7.2f %6.2f | %s / %s / %s | %s\n",
        rl, c1$chm_res, c1$vwf_a, fmt(c1$recall,2), fmt(c1$precision,2),
        fmt(c1$F1,2), vs))
  }
  cat("\n=== full-pooled IN-SAMPLE optimum (all plots; reference) ===\n")
  for (rl in rung_lab) { f1 <- m_full[m_full$rung==rl,]; if (!nrow(f1)) next
    cat(sprintf("%-7s chm_res=%.2f vwf_a=%.2f  rec=%s prec=%s F1=%s\n",
        rl, f1$chm_res, f1$vwf_a, fmt(f1$recall,2), fmt(f1$precision,2), fmt(f1$F1,2))) }

  ## ---- verdict on the headline finding (single seed) ----
  # VWF-slope spread is computed OUT-OF-SAMPLE on the held-out validation plots
  # (chm_res still comes from the calibration optimum). If no validation plots
  # exist for a rung the spread is NA and flagged, not silently filled in.
  vs <- if (length(val)) vwf_spread(r, val, selC) else NULL
  cat("\n=== does the finding survive (SEED single split)? ===\n")
  n05 <- sum(selC$chm_res == 0.5); ntot <- nrow(selC)
  cat(sprintf("  calib-optimal chm_res == 0.5 in %d/%d rungs\n", n05, ntot))
  cat("  VWF-slope F1 spread is pooled on HELD-OUT validation plots (out-of-sample)\n")
  for (rl in rung_lab) { z <- selC[selC$rung==rl,]
    w <- if (!is.null(vs)) vs[vs$rung==rl,] else NULL
    if (!nrow(z)) next
    if (is.null(w) || !nrow(w) || is.na(w$spread)) {
      cat(sprintf("    rung %-6s: chm_res=%.2f vwf_a=%.2f | VWF-slope F1 spread=  -- (no/thin valid subset, n_pool=%d)\n",
          rl, z$chm_res, z$vwf_a, if (is.null(w) || !nrow(w)) 0L else w$n_pool))
    } else {
      cat(sprintf("    rung %-6s: chm_res=%.2f vwf_a=%.2f | VWF-slope F1 spread=%s (min %s, max %s, n_pool=%d)\n",
          rl, z$chm_res, z$vwf_a, fmt(w$spread,3), fmt(w$f1_min,3), fmt(w$f1_max,3), w$n_pool)) }
  }
  max_spread <- if (!is.null(vs)) max(vs$spread, na.rm = TRUE) else NA_real_
  cat(sprintf("  max VWF-slope F1 spread across rungs (held-out) = %s\n", fmt(max_spread,3)))

  ## ---- collect long-form rows for the CSV (this SEED) ----
  mk <- function(df, split) if (is.null(df) || !nrow(df)) NULL else
    data.frame(site=SITE, rung=df$rung, split=split, chm_res=df$chm_res,
               vwf_a=df$vwf_a, n_plots=df$n_plots, n_ref=df$n_ref, n_det=df$n_det,
               TP=df$TP, recall=df$recall, precision=df$precision, F1=df$F1)
  long <- rbind(mk(m_cal,"calib"), mk(m_val,"valid"), mk(m_full,"full_insample"))
  long$seed <- SEED
  write.csv(long, file.path(nd, "calval_metrics.csv"), row.names = FALSE)
  cat(sprintf("\nwrote -> %s\n", file.path(nd, "calval_metrics.csv")))
  all_long[[SITE]] <- long

  ## ---- multi-seed robustness ----
  cat(sprintf("\n=== multi-seed robustness (SEEDS=%s, FRAC=%.2f) ===\n",
              paste(SEEDS, collapse=","), FRAC))
  rb <- do.call(rbind, lapply(SEEDS, function(sd) {
    pts <- assign_split(pt, sd, FRAC)
    ca  <- pts$plot[pts$split=="calib"]; va <- pts$plot[pts$split=="valid"]
    sc  <- best_per_rung(r, ca)
    mv  <- if (length(va)) apply_params(r, va, sc) else NULL
    do.call(rbind, lapply(rung_lab, function(rl) {
      z <- sc[sc$rung==rl,]; if (!nrow(z)) return(NULL)
      v <- if (!is.null(mv)) mv[mv$rung==rl,] else NULL
      data.frame(seed=sd, rung=rl, opt_chm=z$chm_res, opt_vwf=z$vwf_a,
                 val_F1 = if (!is.null(v) && nrow(v)) v$F1 else NA_real_,
                 val_rec= if (!is.null(v) && nrow(v)) v$recall else NA_real_,
                 val_prec=if (!is.null(v) && nrow(v)) v$precision else NA_real_) })) }))
  agg <- do.call(rbind, lapply(rung_lab, function(rl) {
    z <- rb[rb$rung==rl,]; if (!nrow(z)) return(NULL)
    tabchm <- sort(table(z$opt_chm), decreasing = TRUE)
    modal_chm <- as.numeric(names(tabchm)[1])
    data.frame(rung=rl,
               modal_chm = modal_chm,
               chm05_frac = mean(z$opt_chm == 0.5),
               valF1_med = median(z$val_F1, na.rm=TRUE),
               valF1_min = min(z$val_F1, na.rm=TRUE),
               valF1_max = max(z$val_F1, na.rm=TRUE),
               n_seeds = sum(!is.na(z$val_F1))) }))
  cat(sprintf("%-7s %9s %10s %10s %s\n",
              "rung","modal_chm","chm0.5frac","valF1_med","[valF1 min..max]"))
  for (rl in rung_lab) { a <- agg[agg$rung==rl,]; if (!nrow(a)) next
    cat(sprintf("%-7s %9.2f %10s %10s [%s .. %s]\n", rl, a$modal_chm,
        fmt(a$chm05_frac,2), fmt(a$valF1_med,3), fmt(a$valF1_min,3),
        fmt(a$valF1_max,3))) }
  robust[[SITE]] <- list(agg=agg, raw=rb)
}

cat("\n==================== DONE ====================\n")
cat("calval_metrics.csv written per site under $CLAUDE_JOB_DIR/neon/<SITE>/\n")
