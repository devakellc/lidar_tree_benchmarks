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

# #S1 crown-width + height -> DBH / above-ground biomass allometry (GitHub #74).
#
# The crown issues (#7/#30-#35) score crown DIAMETER vs field crown diameter only;
# none predicts DBH or biomass, yet DBH/AGB are the operational products a
# detection pipeline feeds. This closes the loop: it joins each matched crown's
# d_eq (+ matched apex height) to field `stemDiameter`/`taxonID`, fits
# crown-geometry -> DBH models (general, conifer/broadleaf functional-type), and
# reports how well each segmenter's crowns predict field DBH -- per crown_class
# and per density rung -- then DERIVES above-ground biomass from predicted DBH
# (Jenkins 2003 generic; NEON has no field AGB, so AGB is a labelled derived
# product). A detector whose crowns predict DBH well is more valuable than raw F1
# implies, and we test whether that survives decimation.
#
# Usage:  Rscript scripts/crown_allometry.R SITES=SOAP,SJER,TEAK
# Reads (read-only): the matched crown-metrics CSVs score_crowns_against_field
#   emits ({crown_metrics,crown_metrics_3d,crown_metrics_best,segmentanytree_
#   crown_metrics,forestformer3d_crown_metrics}_results/.csv) + ground_truth_
#   stems.csv (stemDiameter, taxonID, height). Writes work/neon/<SITE>/
#   crown_allometry.csv (one row per matched tree: geometry + field DBH + pred
#   DBH + derived AGB) and prints the skill tables.
suppressMessages({ library(data.table) })
d <- .job_dir()
source(.find("allometry_lib.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else c("SOAP", "SJER", "TEAK")

CM_FILES <- c("crown_metrics_results.csv", "crown_metrics_3d_results.csv",
              "crown_metrics_best_results.csv", "segmentanytree_crown_metrics.csv",
              "forestformer3d_crown_metrics.csv")

## ---- assemble matched crowns + field DBH/taxon/height per site -------------
load_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv"); if (!file.exists(gtf)) return(NULL)
  gt <- read.csv(gtf, stringsAsFactors = FALSE)
  fld <- unique(gt[, c("individualID", "stemDiameter", "taxonID", "height")])
  rows <- list()
  for (f in CM_FILES) {
    p <- file.path(nd, f); if (!file.exists(p)) next
    cm <- read.csv(p, stringsAsFactors = FALSE)
    if (is.null(cm$rung)) cm$rung <- "native"        # the 3-D arm CSV has no rung
    keep <- c("site", "plot", "rung", "algo", "crown_class", "individualID",
              "d_eq", "area")
    cm <- cm[, intersect(keep, names(cm)), drop = FALSE]
    rows[[length(rows) + 1]] <- cm
  }
  if (!length(rows)) return(NULL)
  cm <- rbindlist(rows, fill = TRUE)
  cm <- merge(as.data.frame(cm), fld, by = "individualID", all.x = TRUE)
  cm$site <- site
  cm$field_dbh <- suppressWarnings(as.numeric(cm$stemDiameter))   # cm
  cm$ftype <- functional_type(cm$taxonID)
  cm[is.finite(cm$d_eq) & is.finite(cm$field_dbh) & cm$field_dbh > 0, , drop = FALSE]
}

## ---- R2 of a crown -> DBH lm on a subset (in-sample explained variance) ----
lm_r2 <- function(df, form) {
  df <- df[stats::complete.cases(df[, all.vars(form)]), , drop = FALSE]
  if (nrow(df) < 5) return(data.frame(n = nrow(df), r2 = NA_real_,
                                      rmse = NA_real_, bias = NA_real_))
  fit <- tryCatch(stats::lm(form, data = df), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(n = nrow(df), r2 = NA, rmse = NA, bias = NA))
  fit_stats(stats::predict(fit), df$field_dbh)
}

run_main <- function() {
  cm <- rbindlist(Filter(Negate(is.null), lapply(SITES, load_site)), fill = TRUE)
  if (!nrow(cm)) { cat("No matched crowns with field DBH.\n"); return(invisible()) }
  cm <- as.data.frame(cm)
  nat <- cm[cm$rung == "native", , drop = FALSE]
  cat(sprintf("Matched crowns with field DBH: %d (native %d) across %d sites, %d segmenters\n",
              nrow(cm), nrow(nat), length(unique(cm$site)), length(unique(cm$algo))))

  ## (1) does crown width add over height? pooled native
  cat("\n===== DBH PREDICTABILITY (pooled, native) =====\n")
  cat(sprintf("%-22s %6s %7s %7s %7s\n", "model", "n", "R2", "RMSE", "bias"))
  for (m in list(c("height only", "field_dbh ~ height"),
                 c("crown d_eq only", "field_dbh ~ d_eq"),
                 c("d_eq + height", "field_dbh ~ d_eq + height"))) {
    s <- lm_r2(nat, stats::as.formula(m[2]))
    cat(sprintf("%-22s %6d %7.3f %7.3f %7.3f\n", m[1], s$n, s$r2, s$rmse, s$bias))
  }

  ## (2) per-segmenter crown->DBH skill (native): are this arm's crowns good?
  cat("\n===== CROWN d_eq -> DBH PER SEGMENTER (native, R2 in-sample) =====\n")
  cat(sprintf("%-30s %6s %7s %7s\n", "segmenter", "n", "R2", "RMSE"))
  segs <- sort(unique(nat$algo))
  for (a in segs) {
    s <- lm_r2(nat[nat$algo == a, , drop = FALSE], field_dbh ~ d_eq)
    if (!is.na(s$r2)) cat(sprintf("%-30s %6d %7.3f %7.3f\n", a, s$n, s$r2, s$rmse))
  }

  ## (3) per crown_class (pooled native) + (4) per rung (decimation survival)
  cat("\n===== CROWN d_eq -> DBH by crown_class (pooled native) =====\n")
  for (cl in c("dominant", "codominant", "intermediate", "suppressed")) {
    s <- lm_r2(nat[nat$crown_class == cl, , drop = FALSE], field_dbh ~ d_eq)
    cat(sprintf("  %-13s n=%4d R2=%6.3f RMSE=%6.3f\n", cl, s$n, s$r2, s$rmse))
  }
  cat("\n===== CROWN d_eq -> DBH by density rung (decimation survival) =====\n")
  for (rg in c("native", "8", "4", "2", "1")) {
    s <- lm_r2(cm[cm$rung == rg, , drop = FALSE], field_dbh ~ d_eq)
    if (s$n > 0) cat(sprintf("  rung %-7s n=%4d R2=%6.3f RMSE=%6.3f\n", rg, s$n, s$r2, s$rmse))
  }

  ## (5) DERIVED AGB: fit the general crown->DBH model, predict DBH, derive AGB.
  gen <- stats::lm(field_dbh ~ d_eq, data = nat[stats::complete.cases(nat[, c("field_dbh", "d_eq")]), ])
  cm$pred_dbh <- as.numeric(stats::predict(gen, newdata = cm))
  cm$field_agb <- agb_from_dbh(cm$field_dbh, cm$ftype)        # derived from FIELD dbh
  cm$pred_agb  <- agb_from_dbh(cm$pred_dbh, cm$ftype)         # derived from PRED dbh
  an <- cm[cm$rung == "native" & is.finite(cm$field_agb) & is.finite(cm$pred_agb), ]
  sa <- fit_stats(an$pred_agb, an$field_agb)
  cat(sprintf("\n===== DERIVED AGB (predicted-DBH vs field-DBH AGB, native) =====\n"))
  cat(sprintf("n=%d  AGB R2=%.3f  RMSE=%.1f kg  bias=%+.1f kg  (mean field AGB %.1f kg)\n",
              sa$n, sa$r2, sa$rmse, sa$bias, mean(an$field_agb, na.rm = TRUE)))
  cat("(AGB is a DERIVED product: Jenkins 2003 generic; NEON has no field AGB.)\n")

  for (site in SITES) {
    sub <- cm[cm$site == site, , drop = FALSE]; if (!nrow(sub)) next
    write.csv(sub, file.path(d, "neon", site, "crown_allometry.csv"), row.names = FALSE)
  }
  cat(sprintf("\nDONE: wrote crown_allometry.csv per site (%d rows total)\n", nrow(cm)))
}

if (sys.nframe() == 0L) run_main()
