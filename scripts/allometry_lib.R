#!/usr/bin/env Rscript
# #S1 crown-allometry helpers (pure; no I/O, no heavy deps), so the
# taxon classification, the AGB equation, and the error stats are unit-testable
# with no NEON data. crown_allometry.R sources this to turn matched crown
# geometry + field DBH into per-segmenter DBH-prediction skill and a DERIVED AGB.

## ---- functional type from a NEON taxonID ---------------------------------
# conifer vs broadleaf by genus, for the functional-type allometry split. The
# conifer set covers the Sierra mixed-conifer taxa at SOAP/TEAK/SJER (Pinus,
# Abies, Calocedrus, Pseudotsuga, Sequoiadendron, Tsuga, Juniperus, Cupressus);
# every other woody taxon is broadleaf. A blank / NA taxonID returns NA (never a
# fabricated type). Vectorized.
CONIFER_GENERA <- c("PI", "AB", "CA", "PS", "SE", "TS", "JU", "CU")  # 2-letter genus
.genus2 <- function(tx) toupper(substr(tx, 1, 2))
# CADE27 = Calocedrus (conifer) but "CA" also prefixes some hardwoods; the only
# "CA" tree at these sites is CADE (incense cedar), so the 2-letter rule is safe
# here. Explicit hardwood overrides guard the few ambiguous prefixes.
HARDWOOD_OVERRIDE <- c("CANO9", "CASP8")            # hypothetical CA-hardwoods -> broadleaf
functional_type <- function(taxonID) {
  tx <- as.character(taxonID)
  out <- rep(NA_character_, length(tx))
  ok <- !is.na(tx) & nzchar(tx)
  g <- .genus2(tx)
  out[ok] <- ifelse(g[ok] %in% CONIFER_GENERA, "conifer", "broadleaf")
  out[tx %in% HARDWOOD_OVERRIDE] <- "broadleaf"
  out
}

## ---- derived above-ground biomass from DBH (Jenkins et al. 2003) ----------
# National-scale generic AGB estimators: AGB(kg) = exp(b0 + b1*ln(DBH_cm)), with
# the "pine" group for conifers and the "mixed hardwood" group for broadleaves
# (Jenkins, Chojnacky, Heath & Birdsey 2003, For. Sci. 49(1)). DERIVED product:
# NEON provides field DBH, NOT field AGB, so AGB here is a transparent function of
# predicted/observed DBH + functional type, never validated against field AGB.
# NA dbh or NA/unknown type -> NA. Vectorized (dbh and type recycled together).
JENKINS <- list(conifer = c(b0 = -2.5356, b1 = 2.4349),     # pine group
                broadleaf = c(b0 = -2.4800, b1 = 2.4835))   # mixed hardwood group
agb_from_dbh <- function(dbh_cm, type) {
  n <- max(length(dbh_cm), length(type))
  dbh <- rep_len(dbh_cm, n); ty <- rep_len(as.character(type), n)
  out <- rep(NA_real_, n)
  for (g in names(JENKINS)) {
    sel <- !is.na(dbh) & dbh > 0 & !is.na(ty) & ty == g
    if (any(sel)) out[sel] <- exp(JENKINS[[g]]["b0"] + JENKINS[[g]]["b1"] * log(dbh[sel]))
  }
  out
}

## ---- prediction skill: R2 / RMSE / MAE / bias -----------------------------
# Pooled over the finite (pred, obs) pairs: r2 = 1 - SS_res/SS_tot (vs the obs
# mean, so it can go negative for a bad fit), rmse/mae from the residuals, bias =
# mean(pred - obs). Needs >=2 finite pairs, else NA stats with n = #pairs.
fit_stats <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  p <- pred[ok]; o <- obs[ok]; n <- length(p)
  if (n < 2) return(data.frame(n = n, r2 = NA_real_, rmse = NA_real_,
                               mae = NA_real_, bias = NA_real_))
  e <- p - o
  ss_res <- sum((o - p)^2); ss_tot <- sum((o - mean(o))^2)
  data.frame(n = n,
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_,
             rmse = sqrt(mean(e^2)), mae = mean(abs(e)), bias = mean(e))
}
