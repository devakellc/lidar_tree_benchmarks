#!/usr/bin/env Rscript
# Analyse the NEON SOAP density-ladder sweep: pool stems across plots (sum TP /
# sum n_ref -- NOT a mean of per-plot rates, which over-weights small plots),
# produce the density-sensitivity table by crown class, the best-parameter
# table per density rung, and PNG figures. Writes a markdown summary fragment.
args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE <- if (is.null(A$SITE)) "SOAP" else A$SITE
d  <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
nd <- file.path(d, "neon", SITE)
r  <- read.csv(file.path(nd, "sweep_results.csv"), stringsAsFactors = FALSE)

# rung ordering by measured first-return density (native highest)
rung_lab <- c("native","8","4","2","1")
r$rung <- factor(r$rung, levels = rung_lab)
# recover per-plot true-positive counts for correct pooling
r$tp_core <- round(r$precision * r$n_det)
classes <- c("dominant","codominant","intermediate","suppressed")
hbands  <- c("short","mid","tall")          # <8 m, 8-15 m, >=15 m
hb_present <- all(paste0("n_h_", hbands) %in% names(r))
for (cl in classes) {
  nc <- r[[paste0("n_", cl)]]; rc <- r[[paste0("rec_", cl)]]
  r[[paste0("tp_", cl)]] <- ifelse(nc > 0, round(rc * nc), 0)  # 0 when class absent
}
if (hb_present) for (b in hbands) {
  nc <- r[[paste0("n_h_", b)]]; rc <- r[[paste0("rec_h_", b)]]
  r[[paste0("tp_h_", b)]] <- ifelse(nc > 0, round(rc * nc), 0)
} else hbands <- character(0)

pool <- function(df) {
  out <- data.frame(
    n_plots   = length(unique(df$plot)),
    frdens    = round(mean(df$frdens), 2),
    n_ref     = sum(df$n_ref), n_det = sum(df$n_det), TP = sum(df$TP),
    recall    = sum(df$TP) / sum(df$n_ref),
    precision = sum(df$tp_core, na.rm = TRUE) / sum(df$n_det),
    height_rmse = if (sum(df$TP) > 0)
      sqrt(weighted.mean(df$height_rmse^2, df$TP, na.rm = TRUE)) else NA_real_)
  out$F1 <- if (!is.na(out$recall) && !is.na(out$precision) &&
                (out$recall + out$precision) > 0)
    2 * out$recall * out$precision / (out$recall + out$precision) else NA_real_
  for (cl in classes) {
    nref <- sum(df[[paste0("n_", cl)]], na.rm = TRUE)
    tp   <- sum(df[[paste0("tp_", cl)]], na.rm = TRUE)
    out[[paste0("rec_", cl)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_", cl)]]   <- nref
  }
  for (b in hbands) {
    nref <- sum(df[[paste0("n_h_", b)]], na.rm = TRUE)
    tp   <- sum(df[[paste0("tp_h_", b)]], na.rm = TRUE)
    out[[paste0("rec_h_", b)]] <- if (nref) tp / nref else NA_real_
    out[[paste0("n_h_", b)]]   <- nref
  }
  out
}

## ---- 1. headline: density sensitivity at the modal parameters -------------
mod <- r[r$chm_res == 0.5 & r$vwf_a == 0.10, ]
by_rung <- do.call(rbind, lapply(rung_lab, function(rl) {
  s <- mod[mod$rung == rl, ]; if (!nrow(s)) return(NULL)
  cbind(rung = rl, pool(s)) }))
cat("=== Density sensitivity (chm_res=0.5, vwf_a=0.10; pooled over plots) ===\n")
print(by_rung[, c("rung","frdens","n_plots","n_ref","n_det","recall","precision",
                  "F1","rec_dominant","rec_codominant","rec_intermediate",
                  "rec_suppressed","height_rmse")], row.names = FALSE, digits = 2)
if (length(hbands)) {
  cat("\n=== ... same, by HEIGHT band (short<8m, mid 8-15m, tall>=15m) ===\n")
  print(by_rung[, c("rung","frdens","rec_h_tall","n_h_tall","rec_h_mid","n_h_mid",
                    "rec_h_short","n_h_short")], row.names = FALSE, digits = 2)
}

## ---- 2. best parameter set per rung (by pooled F1) ------------------------
cat("\n=== Best (chm_res, vwf_a) per density rung, by pooled F1 ===\n")
best <- do.call(rbind, lapply(rung_lab, function(rl) {
  s <- r[r$rung == rl, ]; if (!nrow(s)) return(NULL)
  grid <- unique(s[, c("chm_res","vwf_a")])
  pg <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    ss <- s[s$chm_res == grid$chm_res[i] & s$vwf_a == grid$vwf_a[i], ]
    cbind(chm_res = grid$chm_res[i], vwf_a = grid$vwf_a[i], pool(ss)) }))
  pg <- pg[order(-pg$F1), ][1, ]
  cbind(rung = rl, pg[, c("frdens","chm_res","vwf_a","recall","precision","F1")]) }))
print(best, row.names = FALSE, digits = 2)

## ---- 3. parameter main effects (pooled over plots & rungs) ----------------
cat("\n=== CHM resolution main effect (pooled, all rungs) ===\n")
for (rs in sort(unique(r$chm_res))) { s <- r[r$chm_res == rs, ]; p <- pool(s)
  cat(sprintf("  res=%.2f  recall=%.2f precision=%.2f F1=%.2f\n", rs, p$recall, p$precision, p$F1)) }
cat("=== VWF slope main effect (pooled, all rungs) ===\n")
for (av in sort(unique(r$vwf_a))) { s <- r[r$vwf_a == av, ]; p <- pool(s)
  cat(sprintf("  a=%.2f   recall=%.2f precision=%.2f F1=%.2f\n", av, p$recall, p$precision, p$F1)) }

## ---- 4. figures -----------------------------------------------------------
fig <- file.path(nd, "figs"); dir.create(fig, showWarnings = FALSE)
xo <- by_rung$frdens; ord <- order(xo)
png(file.path(fig, "density_sensitivity.png"), 1000, 750, res = 130)
par(mar = c(4.2,4.2,2.5,1))
plot(xo[ord], by_rung$recall[ord], type="o", pch=19, ylim=c(0,1), log="x",
     xlab="measured first-return density (pulses/m^2)", ylab="recall",
     main="SOAP density sensitivity by crown class", col="black", lwd=2)
cols <- c(dominant="#1b9e77", codominant="#7570b3",
          intermediate="#d95f02", suppressed="#e7298a")
for (cl in classes) lines(xo[ord], by_rung[[paste0("rec_",cl)]][ord], type="o",
                          pch=1, col=cols[cl], lwd=2)
legend("topleft", c("overall", classes), col=c("black", cols), lwd=2,
       pch=c(19,1,1,1,1), bty="n", cex=0.85)
dev.off()

png(file.path(fig, "precision_recall_f1.png"), 1000, 750, res = 130)
par(mar = c(4.2,4.2,2.5,1))
matplot(xo[ord], by_rung[ord, c("recall","precision","F1")], type="o", pch=19,
        lty=1, log="x", ylim=c(0,1), xlab="first-return density (pulses/m^2)",
        ylab="score", main="SOAP detection vs density (res=0.5, a=0.10)",
        col=c("#1b9e77","#d95f02","#7570b3"), lwd=2)
legend("topleft", c("recall","precision","F1"),
       col=c("#1b9e77","#d95f02","#7570b3"), lwd=2, pch=19, bty="n")
dev.off()

write.csv(by_rung, file.path(nd, "summary_by_rung.csv"), row.names = FALSE)
write.csv(best,    file.path(nd, "summary_best_params.csv"), row.names = FALSE)
cat(sprintf("\nfigures -> %s ; summaries -> %s\n", fig, nd))

## ---- 5. calibration/validation split (issue #3) ---------------------------
# Section 2 above picks "best (chm_res, vwf_a) per rung" IN-SAMPLE: parameters
# are tuned and scored on the SAME pooled plots. To test whether that choice
# (and the headline "chm_res=0.5 m is F1-optimal; VWF slope second-order")
# survives out-of-sample, run the standalone calibration/validation split, which
# tunes on a calibration subset of plots and reports held-out F1 on the rest:
#
#   Rscript scripts/calval_split.R SITES=SJER,SOAP,TEAK SEED=1 FRAC=0.5 \
#           SEEDS=1,2,3,4,5
#
# It reuses this exact pool() rule (sum TP / sum n_ref; tp_core=round(prec*n_det)),
# stratifies the split by plotType x crown-class mix, writes long-form
# calval_metrics.csv per site, and renders results/calibration-validation-results.md.
# The pooling logic is kept self-contained there; it is NOT duplicated here.
