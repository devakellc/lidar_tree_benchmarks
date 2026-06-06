#!/usr/bin/env Rscript
# Validate detected treetop apex HEIGHT against NEON field-measured tree height.
# Matches stems to treetops on 2D POSITION ONLY (no height gate -> unbiased
# height comparison) within a tight tolerance, at native density, modal params.
# Reports bias, RMSE, MAE, R^2 and a linear fit, pooled and per site / crown
# class, and writes a field-vs-apex scatter. Usage:
#   Rscript scripts/validate_heights.R [SITES=SJER,SOAP,TEAK] [TOL=3] [RES=0.5]
#                                      [A=0.10] [MEAS_YEAR=2021] [OUT=pairs.csv]
# MEAS_YEAR (optional, issue #5): restrict to stems whose nearest field
# measurement is in that exact year, so signed height bias can be compared
# exact-year vs the +/-4 yr baseline. Default unchanged. When MEAS_YEAR is set
# the pairs default to a distinct height_pairs_<YEAR>.csv (never overwriting the
# baseline height_pairs.csv); OUT overrides the path explicitly.
suppressMessages({ library(lidR); library(lasR); library(terra); library(sf) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
source(file.path("scripts", "sweep_lib.R"))

args <- strsplit(commandArgs(TRUE), "="); A <- setNames(lapply(args,`[`,2), sapply(args,`[`,1))
SITES <- strsplit(if (is.null(A$SITES)) "SJER,SOAP,TEAK" else A$SITES, ",")[[1]]
TOL <- as.numeric(if (is.null(A$TOL)) 3 else A$TOL)
RES <- as.numeric(if (is.null(A$RES)) 0.5 else A$RES)
AA  <- as.numeric(if (is.null(A$A)) 0.10 else A$A)
MEAS_YEAR <- if (is.null(A$MEAS_YEAR)) NA_integer_ else as.integer(A$MEAS_YEAR)
OUT <- if (is.null(A$OUT)) {
  if (is.na(MEAS_YEAR)) file.path(d, "neon", "height_pairs.csv")
  else file.path(d, "neon", sprintf("height_pairs_%d.csv", MEAS_YEAR))
} else A$OUT

collect_site <- function(site) {
  nd  <- file.path(d, "neon", site)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"))
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"))
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E) & !is.na(gt$height), ]
  if (!is.na(MEAS_YEAR)) {
    nb <- nrow(gt)
    gt <- gt[!is.na(gt$meas_year) & gt$meas_year == MEAS_YEAR, ]
    cat(sprintf("[%s] MEAS_YEAR=%d : kept %d of %d height stems\n",
                site, MEAS_YEAR, nrow(gt), nb))
  }
  laz <- list.files(file.path(nd, "lidar"), pattern="\\.laz$", recursive=TRUE, full.names=TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  keep <- names(table(gt$plotID))[table(gt$plotID) >= 6]
  keep <- intersect(keep, pc$plotID)
  tmp  <- file.path(tempdir(), "hv"); dir.create(tmp, showWarnings = FALSE)
  pairs <- list()
  for (pid in keep) {
    ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E-cx) <= ph & abs(gt$N-cy) <= ph, ]
    if (!nrow(stems)) next
    prep <- tryCatch(prepare_clip(ctg, cx, cy, NA, tmp, core_half = ph), error=function(e) NULL)
    if (is.null(prep)) next
    det <- tryCatch(detect_lasr(prep$file, RES, AA, prep$frdens), error=function(e) NULL)
    unlink(prep$file)
    if (is.null(det) || !nrow(det)) next
    # POSITION-ONLY match (no height gate) within TOL
    m <- greedy_match(stems$E, stems$N, det$x, det$y, TOL)
    ok <- m > 0
    if (any(ok)) pairs[[pid]] <- data.frame(
      site = site, plot = pid, crown_class = stems$crown_class[ok],
      field_h = stems$height[ok], apex_z = det$z[m[ok]], dbh = stems$stemDiameter[ok])
  }
  do.call(rbind, pairs)
}

stats <- function(df, label) {
  e <- df$apex_z - df$field_h
  fit <- lm(apex_z ~ field_h, df)
  cat(sprintf("%-22s n=%4d  bias=%+5.2f  RMSE=%4.2f  MAE=%4.2f  R2=%.2f  slope=%.2f\n",
      label, nrow(df), mean(e), sqrt(mean(e^2)), mean(abs(e)),
      summary(fit)$r.squared, coef(fit)[2]))
}

all <- do.call(rbind, Filter(Negate(is.null), lapply(SITES, collect_site)))
if (is.null(all) || !nrow(all)) { cat("No height pairs found.\n"); quit(save = "no") }
write.csv(all, OUT, row.names = FALSE)

yr_lab <- if (is.na(MEAS_YEAR)) "all years (+/-4 yr baseline)" else
  sprintf("exact meas_year=%d", MEAS_YEAR)
cat(sprintf("\n=== Detected apex height vs NEON field height [%s] (native density, res=%.2f, a=%.2f, pos-only match <=%gm) ===\n", yr_lab, RES, AA, TOL))
cat("(bias = apex - field; positive = LiDAR apex taller than field-measured)\n\n")
stats(all, "ALL SITES")
for (st in SITES) if (sum(all$site==st)) stats(all[all$site==st, ], paste0("  ", st))
cat("\n-- by crown class (all sites) --\n")
for (cl in c("dominant","codominant","intermediate","suppressed"))
  if (sum(all$crown_class==cl, na.rm=TRUE) >= 5) stats(all[which(all$crown_class==cl), ], paste0("  ", cl))
cat("\n-- tall trees only (field >= 15 m) --\n")
if (sum(all$field_h>=15) >= 5) stats(all[all$field_h>=15, ], "  field>=15m")

## scatter (distinct filename for the exact-year cut so it never overwrites the
## +/-4 yr baseline figure)
fig <- file.path(d, "neon", "figs"); dir.create(fig, showWarnings = FALSE, recursive = TRUE)
fig_png <- if (is.na(MEAS_YEAR)) "height_validation.png" else
  sprintf("height_validation_%d.png", MEAS_YEAR)
png(file.path(fig, fig_png), 1000, 760, res = 130)
par(mar = c(4.2,4.2,2.5,1)); cols <- c(SJER="#1b9e77", SOAP="#d95f02", TEAK="#7570b3")
lim <- c(0, max(all$field_h, all$apex_z, na.rm=TRUE))
plot(all$field_h, all$apex_z, col = cols[all$site], pch = 19, cex = 0.5,
     xlim = lim, ylim = lim, xlab = "NEON field height (m)",
     ylab = "detected CHM apex height (m)", main = "Treetop height vs field truth")
abline(0, 1, lty = 2); abline(lm(apex_z ~ field_h, all), col = "black", lwd = 2)
legend("topleft", c(SITES, "1:1", "fit"), col = c(cols[SITES], "black", "black"),
       pch = c(19,19,19,NA,NA), lty = c(NA,NA,NA,2,1), lwd = 2, bty = "n", cex = 0.8)
dev.off()
cat(sprintf("\npairs -> %s ; scatter -> %s\n", OUT, fig))
