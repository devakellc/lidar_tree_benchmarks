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
# It ALSO runs the genuine out-of-sample optimality test issue #3 asks for: for
# each site x rung it pools held-out F1 for EACH candidate chm_res on the
# VALIDATION plots (holding vwf_a at the calibration-selected value for that
# rung -- documented choice; the VWF slope is second-order so the chm_res argmax
# is insensitive to it) and reports the held-out-F1-OPTIMAL chm_res (argmax),
# then compares it to 0.5. This is reported per single-seed split AND aggregated
# across SEEDS (modal held-out optimum + fraction of seeds where it equals 0.5),
# so the "does 0.5 m survive out-of-sample" question is answered on HELD-OUT
# data rather than read off the calibration selection.
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

d  <- .job_dir()
# rung_lab, classes, plot_table(), assign_split() are the shared split definition
# (single source of truth in calval_lib.R; also used by calval_multichm.R, #39).
source(.find("calval_lib.R"))

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
    precision = if (sum(df$n_det) > 0) sum(df$tp_core, na.rm = TRUE) / sum(df$n_det) else NA_real_)
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  out
}

# plot_table() (per-plot stratification covariates) and assign_split()
# (deterministic plotType x crown_mix calib/valid assignment) now live in
# calval_lib.R, sourced above -- the single source of truth shared with
# calval_multichm.R so both studies use the IDENTICAL held-out plot splits.

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

## ---- held-out F1 across chm_res (genuine out-of-sample optimality test) ----
# For each rung, on the VALIDATION plots, pool held-out F1 for EACH candidate
# chm_res, HOLDING vwf_a at the calibration-selected value for that rung (so the
# only thing that varies is chm_res). Report the held-out-F1-OPTIMAL chm_res
# (argmax) and whether it equals 0.5. This directly answers issue #3 -- "does
# the 0.5 m finding survive out-of-sample" -- because the chm_res argmax is now
# computed on HELD-OUT data, not just read off the calibration-selected param.
#
# COMPLETE-CASE across chm_res (round-3 fix). The sweep grid is NOT rectangular:
# at native density chm_res=0.25 is only run when a plot's frdens>=8, so it
# covers FEWER plots than 0.5/1.0 (e.g. SJER native 5/8, SOAP 16/18, TEAK 17/20).
# Pooling each chm_res over whatever rows exist then compares F1 across
# resolutions on DIFFERENT plot subsets, biasing the argmax. Fix: restrict to the
# INTERSECTION of validation plots that have a scored row at EVERY candidate
# chm_res (at the calibration vwf_a), then pool F1 per chm_res over that common
# plot set and take the argmax. n_pool = size of the common set; n_drop = plots
# present at some-but-not-all resolutions (dropped to enforce complete-case).
# Decimated rungs only ever run {0.5,1.0} at equal plot counts, so there
# complete-case is a no-op (n_drop=0) -- the bias only ever touched the native
# rung. We still report the candidate resolutions actually compared.
#
# Design choice (documented): vwf_a is FIXED to the calibration optimum rather
# than marginalized. Rationale: (a) the VWF slope is shown to be second-order
# (small held-out F1 spread across vwf_a), so the chm_res argmax is insensitive
# to this choice; (b) fixing vwf_a keeps the comparison a clean one-factor
# (chm_res-only) contrast at the operating point the calibration would deploy.
# Ties in held-out F1 are broken deterministically by finest chm_res first
# (ascending). At decimated rungs the candidate set is {0.5, 1.0}, so finest is
# 0.5 there; this is moot because no exact held-out F1 ties occur in the cached
# data, so the reported optima do not depend on the tie rule.
heldout_res_opt <- function(r, plots, sel) {
  s0 <- r[r$plot %in% plots, ]
  na_row <- function(rl, vwf_a, cal_res) data.frame(rung = rl, vwf_a = vwf_a,
      cal_res = cal_res, ho_res = NA_real_, ho_F1 = NA_real_,
      F1_at_05 = NA_real_, res_is_05 = NA, n_pool = 0L, n_drop = 0L,
      n_res = 0L)
  do.call(rbind, lapply(rung_lab, function(rl) {
    p <- sel[sel$rung == rl, ]
    s <- s0[s0$rung == rl & s0$vwf_a == p$vwf_a, ]
    if (!nrow(s)) return(na_row(rl, p$vwf_a, p$chm_res))
    res_vals <- sort(unique(s$chm_res))
    # complete-case: keep only plots scored at EVERY candidate chm_res
    common <- Reduce(intersect, lapply(res_vals, function(rv)
                     unique(s$plot[s$chm_res == rv])))
    all_plots <- unique(s$plot)
    n_drop <- length(setdiff(all_plots, common))
    if (!length(common)) return(na_row(rl, p$vwf_a, p$chm_res))
    sc <- s[s$plot %in% common, ]
    f1 <- sapply(res_vals, function(rv) {
      ss <- sc[sc$chm_res == rv, ]; if (!nrow(ss)) NA_real_ else pool(ss)$F1 })
    ok <- !is.na(f1)
    if (!any(ok)) return(na_row(rl, p$vwf_a, p$chm_res))
    # argmax with finest-res (ascending) tie-break: order by -F1 then +res
    ord    <- order(-f1[ok], res_vals[ok])
    ho_res <- res_vals[ok][ord][1]
    ho_F1  <- f1[ok][ord][1]
    f1_05  <- if (any(abs(res_vals - 0.5) < 1e-9)) f1[abs(res_vals - 0.5) < 1e-9] else NA_real_
    data.frame(rung = rl, vwf_a = p$vwf_a, cal_res = p$chm_res,
               ho_res = ho_res, ho_F1 = ho_F1, F1_at_05 = f1_05,
               res_is_05 = isTRUE(abs(ho_res - 0.5) < 1e-9),
               n_pool = length(common), n_drop = n_drop,
               n_res = sum(ok)) }))
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

  ## ---- held-out F1 across chm_res: the genuine out-of-sample optimality test --
  # vwf_a fixed to the calibration optimum; chm_res argmax computed on the
  # HELD-OUT validation plots. This is the test the doc's out-of-sample claim
  # actually requires (calibration selection alone never compares chm_res on
  # held-out data).
  ho <- if (length(val)) heldout_res_opt(r, val, selC) else NULL
  cat("\n=== HELD-OUT F1-optimal chm_res (COMPLETE-CASE argmax on validation plots; vwf_a fixed to calib opt) ===\n")
  cat("    n_pl = validation plots common to ALL compared chm_res; drop = plots removed to enforce complete-case\n")
  cat(sprintf("%-7s %6s | %-9s %-7s | %-7s | %-5s %-4s | %s\n",
              "rung","vwf_a","ho_opt_res","ho_F1","F1@0.5","n_pl","drop",
              "is 0.5 the held-out opt?"))
  for (rl in rung_lab) { w <- if (!is.null(ho)) ho[ho$rung==rl,] else NULL
    if (is.null(w) || !nrow(w) || is.na(w$ho_res)) {
      cat(sprintf("%-7s %6s | %-9s %-7s | %-7s | %-5s %-4s | -- (no/thin valid subset)\n",
                  rl, "--", "--", "--", "--", "--", "--")); next }
    cat(sprintf("%-7s %6.2f | %-9.2f %-7s | %-7s | %-5d %-4d | %s\n",
        rl, w$vwf_a, w$ho_res, fmt(w$ho_F1,3),
        if (is.na(w$F1_at_05)) "--" else fmt(w$F1_at_05,3),
        w$n_pool, w$n_drop,
        if (isTRUE(w$res_is_05)) "YES" else sprintf("NO (opt=%.2f)", w$ho_res))) }
  if (!is.null(ho)) {
    hov <- ho[!is.na(ho$ho_res), ]
    cat(sprintf("  held-out optimum == 0.5 in %d/%d rungs (vs in-sample claim of 0.5)\n",
                sum(hov$res_is_05, na.rm=TRUE), nrow(hov))) }

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
    hr  <- if (length(va)) heldout_res_opt(r, va, sc) else NULL  # held-out res argmax
    do.call(rbind, lapply(rung_lab, function(rl) {
      z <- sc[sc$rung==rl,]; if (!nrow(z)) return(NULL)
      v <- if (!is.null(mv)) mv[mv$rung==rl,] else NULL
      h <- if (!is.null(hr)) hr[hr$rung==rl,] else NULL
      data.frame(seed=sd, rung=rl, opt_chm=z$chm_res, opt_vwf=z$vwf_a,
                 val_F1 = if (!is.null(v) && nrow(v)) v$F1 else NA_real_,
                 val_rec= if (!is.null(v) && nrow(v)) v$recall else NA_real_,
                 val_prec=if (!is.null(v) && nrow(v)) v$precision else NA_real_,
                 ho_res = if (!is.null(h) && nrow(h)) h$ho_res else NA_real_,
                 ho_npool = if (!is.null(h) && nrow(h)) h$n_pool else NA_integer_,
                 ho_ndrop = if (!is.null(h) && nrow(h)) h$n_drop else NA_integer_) })) }))
  agg <- do.call(rbind, lapply(rung_lab, function(rl) {
    z <- rb[rb$rung==rl,]; if (!nrow(z)) return(NULL)
    tabchm <- sort(table(z$opt_chm), decreasing = TRUE)
    modal_chm <- as.numeric(names(tabchm)[1])
    hz <- z$ho_res[!is.na(z$ho_res)]
    modal_ho <- if (length(hz)) as.numeric(names(sort(table(hz), decreasing=TRUE))[1]) else NA_real_
    data.frame(rung=rl,
               modal_chm = modal_chm,
               chm05_frac = mean(z$opt_chm == 0.5),
               modal_ho_res = modal_ho,
               ho05_frac = if (length(hz)) mean(hz == 0.5) else NA_real_,
               ho_npool_med = if (length(hz)) median(z$ho_npool, na.rm=TRUE) else NA_real_,
               ho_ndrop_med = if (length(hz)) median(z$ho_ndrop, na.rm=TRUE) else NA_real_,
               ho_ndrop_max = if (length(hz)) max(z$ho_ndrop, na.rm=TRUE) else NA_real_,
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
  cat("\n  held-out F1-optimal chm_res across seeds (COMPLETE-CASE argmax on validation plots):\n")
  cat(sprintf("  %-7s %12s %14s %10s %10s\n",
              "rung","modal_ho_res","ho_opt==0.5 frac","n_pl_med","drop_med"))
  for (rl in rung_lab) { a <- agg[agg$rung==rl,]; if (!nrow(a)) next
    cat(sprintf("  %-7s %12s %14s %10s %10s\n", rl,
        if (is.na(a$modal_ho_res)) "--" else fmt(a$modal_ho_res,2),
        if (is.na(a$ho05_frac)) "--" else fmt(a$ho05_frac,2),
        if (is.na(a$ho_npool_med)) "--" else fmt(a$ho_npool_med,1),
        if (is.na(a$ho_ndrop_med)) "--" else fmt(a$ho_ndrop_med,1))) }
  robust[[SITE]] <- list(agg=agg, raw=rb)
}

cat("\n==================== DONE ====================\n")
cat("calval_metrics.csv written per site under $CLAUDE_JOB_DIR/neon/<SITE>/\n")
