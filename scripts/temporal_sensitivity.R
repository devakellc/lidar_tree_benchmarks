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

# Temporal-sensitivity cut for the NEON density-ladder sweep (issue #5).
#
# Ground truth pairs each field stem with the apparentindividual measurement
# nearest the 2021 LiDAR (within +/-4 yr). This script quantifies how much that
# temporal mismatch inflates or deflates the metrics by comparing, at the modal
# detection parameters (chm_res=0.5, vwf_a=0.10), the +/-4 yr BASELINE
# (sweep_results.csv) against the EXACT-2021 re-score (sweep_results_2021.csv,
# produced by run_sweep.R MEAS_YEAR=2021). For each density rung it pools both
# cuts and reports recall / precision / F1 / height_rmse and the delta
# (exact-2021 minus +/-4 yr), plus n_plots and n_ref per cut (plots can drop
# below the MINTREES=6 floor once restricted to the 2021 subset).
#
# Pooling reuses analyze_sweep.R's rule: pooled recall = sum(TP)/sum(n_ref) (NOT
# a mean of per-plot rates); precision = sum(tp_core)/sum(n_det) with
# tp_core = round(precision * n_det); height_rmse = TP-weighted RMS of per-plot
# RMSEs.
#
# Usage:
#   Rscript scripts/temporal_sensitivity.R [SITES=TEAK,SOAP]
# Writes work/neon/<SITE>/temporal_compare.csv (one row per rung x cut, plus a
# `delta` cut carrying exact-2021 minus baseline).

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- strsplit(if (is.null(A$SITES)) "TEAK,SOAP" else A$SITES, ",")[[1]]
d     <- .job_dir()

RES <- 0.5; AA <- 0.10                     # modal detection parameters
rung_lab <- c("native", "8", "4", "2", "1")

## ---- pooling (mirrors analyze_sweep.R pool(), headline metrics only) -------
pool <- function(df) {
  if (!nrow(df)) return(NULL)
  tp_core <- round(df$precision * df$n_det)   # recover per-plot core TP
  recall    <- sum(df$TP) / sum(df$n_ref)
  precision <- if (sum(df$n_det) > 0) sum(tp_core, na.rm = TRUE) / sum(df$n_det) else NA_real_
  f1 <- if (!is.na(recall) && !is.na(precision) && (recall + precision) > 0)
    2 * recall * precision / (recall + precision) else NA_real_
  hrmse <- if (sum(df$TP) > 0)
    sqrt(weighted.mean(df$height_rmse^2, df$TP, na.rm = TRUE)) else NA_real_
  data.frame(n_plots = length(unique(df$plot)),
             n_ref = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
             recall = recall, precision = precision, F1 = f1,
             height_rmse = hrmse)
}

## ---- per-site comparison --------------------------------------------------
compare_site <- function(site) {
  nd  <- file.path(d, "neon", site)
  f_base <- file.path(nd, "sweep_results.csv")
  f_2021 <- file.path(nd, "sweep_results_2021.csv")
  if (!file.exists(f_base)) { cat(sprintf("[%s] no baseline sweep_results.csv -> skip\n", site)); return(invisible(NULL)) }
  if (!file.exists(f_2021)) { cat(sprintf("[%s] no sweep_results_2021.csv (run MEAS_YEAR=2021 first) -> skip\n", site)); return(invisible(NULL)) }
  rb <- read.csv(f_base, stringsAsFactors = FALSE)
  r2 <- read.csv(f_2021, stringsAsFactors = FALSE)
  mb <- rb[rb$chm_res == RES & rb$vwf_a == AA, ]
  m2 <- r2[r2$chm_res == RES & r2$vwf_a == AA, ]

  rows <- list()
  for (rl in rung_lab) {
    pb <- pool(mb[mb$rung == rl, ])
    p2 <- pool(m2[m2$rung == rl, ])
    if (is.null(pb) && is.null(p2)) next
    add <- function(cut, p) if (!is.null(p))
      cbind(data.frame(site = site, rung = rl, cut = cut), p)
    rows[[length(rows) + 1]] <- add("baseline_pm4yr", pb)
    rows[[length(rows) + 1]] <- add("exact_2021", p2)
    if (!is.null(pb) && !is.null(p2)) {
      del <- data.frame(site = site, rung = rl, cut = "delta",
                        n_plots = p2$n_plots - pb$n_plots,
                        n_ref = p2$n_ref - pb$n_ref,
                        n_det = p2$n_det - pb$n_det, TP = p2$TP - pb$TP,
                        recall = p2$recall - pb$recall,
                        precision = p2$precision - pb$precision,
                        F1 = p2$F1 - pb$F1,
                        height_rmse = p2$height_rmse - pb$height_rmse)
      rows[[length(rows) + 1]] <- del
    }
  }
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out) || !nrow(out)) { cat(sprintf("[%s] no modal-parameter rows -> skip\n", site)); return(invisible(NULL)) }

  cat(sprintf("\n=== %s : +/-4 yr baseline vs exact-2021 (chm_res=%.2f, vwf_a=%.2f) ===\n", site, RES, AA))
  show <- out[, c("rung","cut","n_plots","n_ref","n_det","TP",
                  "recall","precision","F1","height_rmse")]
  print(show, row.names = FALSE, digits = 3)

  op <- file.path(nd, "temporal_compare.csv")
  write.csv(out, op, row.names = FALSE)
  cat(sprintf("-> %s\n", op))
  invisible(out)
}

for (s in SITES) compare_site(s)
