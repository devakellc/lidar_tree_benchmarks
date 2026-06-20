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

# #V4 matcher robustness (GitHub issue #78).
#
# The whole benchmark grades detections with greedy_match: greedy nearest-distance
# 1:1 within a FLAT 4 m radius, gated by a hard height band [0.5*az, az+8]. Two
# weaknesses: (1) the apex<->stem-base offset grows with height/crown size/lean/
# slope, so a flat 4 m is too tight for tall dominants and too loose for dense
# understory; (2) greedy-by-distance is globally suboptimal in dense clusters.
# This arm re-scores the detector on the SAME frozen clips with the hardened
# matchers added in sweep_lib.R and reports whether any pooled rate moves:
#   * baseline       -- greedy, flat tol 4, hard gate (the current benchmark)
#   * scaled         -- greedy, per-stem tol = max(4, k*maxCrownDiameter/2, pos_unc)
#   * optimal        -- Hungarian (clue::solve_LSAP), flat tol 4, hard gate
#   * optimal_scaled -- Hungarian, per-stem tol
#   * soft3d         -- Hungarian, flat tol 4, soft 3-D cost (drops the magic band)
# plus a tol_xy {2,3,4,5} x tol_z_up {5,8,12} sensitivity grid (greedy), and a
# false-positive error-structure split (near a matched stem = over-segmentation
# vs isolated = real understory / field-map gap), per crown_class -- the count
# feeds the #P2 router and #P1 fusion.
#
# Detections are regenerated DETERMINISTICALLY from the cached frozen NORMALIZED
# clips (so every matcher scores identical apexes) via detect_lasr at the
# canonical params (density-derived chm_res, vwf_a = 0.10), matching the benchmark.
# A null result (no pooled movement) confirms greedy-flat-4 is adequate; a delta
# says switch. Pools with the canonical pool() (sum counts, never average rates).
#
# Usage:
#   Rscript scripts/matcher_robustness.R SITE=SOAP
#   Rscript scripts/matcher_robustness.R SITES=SOAP,SJER,TEAK CORES=4 RUNGS=native
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv},
#   the cached field crown widths work/neon/<SITE>/vst/<site>_vst_allyears.rds, and
#   the cached frozen normalized clips work/neon/<SITE>/frozen/<SITE>/<plot>/<rung>/
#   {clip_normalized.laz,manifest.json}.
# Writes: work/neon/<SITE>/matcher_robustness.csv (one row per
#   site x plot x rung x config; pool with pool()).
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))     # pool

## ---- args (KEY=VALUE positional) -----------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES    <- as.integer(if (is.null(A$CORES)) 4 else A$CORES)
RUNGS    <- if (is.null(A$RUNGS)) "native" else strsplit(A$RUNGS, ",")[[1]]
BASE_TOL <- as.numeric(if (is.null(A$BASE_TOL)) 4.0 else A$BASE_TOL)
K        <- as.numeric(if (is.null(A$K)) 1.0 else A$K)
NEAR_TOL <- as.numeric(if (is.null(A$NEAR_TOL)) 4.0 else A$NEAR_TOL)
VWF_A    <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
LAMBDA   <- as.numeric(if (is.null(A$LAMBDA)) 0.5 else A$LAMBDA)
MINTREES <- as.integer(if (is.null(A$MINTREES)) 1 else A$MINTREES)
TX_SET   <- c(2, 3, 4, 5)               # tol_xy sensitivity
TZ_SET   <- c(5, 8, 12)                 # tol_z_up sensitivity

# Headline matcher configs (name, matcher method, scaled per-stem tol?, gate, lambda).
CONFIGS <- list(
  list(name = "baseline",       method = "greedy",  scaled = FALSE, tz = 8, lambda = NULL),
  list(name = "scaled",         method = "greedy",  scaled = TRUE,  tz = 8, lambda = NULL),
  list(name = "optimal",        method = "optimal", scaled = FALSE, tz = 8, lambda = NULL),
  list(name = "optimal_scaled", method = "optimal", scaled = TRUE,  tz = 8, lambda = NULL),
  list(name = "soft3d",         method = "optimal", scaled = FALSE, tz = 8, lambda = LAMBDA))

## ---- field crown diameter per site (cached vst rds) -----------------------
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst",
                   paste0(tolower(site), "_vst_allyears.rds"))
  if (!file.exists(rds)) return(data.frame(individualID = character(),
                                           maxCrownDiameter = numeric()))
  dat <- readRDS(rds)
  ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID), c("individualID", "maxCrownDiameter")]
}

frozen_norm_path <- function(nd, site, pid, rung)
  file.path(nd, "frozen", site, pid, rung, "clip_normalized.laz")
frozen_manifest  <- function(nd, site, pid, rung)
  file.path(nd, "frozen", site, pid, rung, "manifest.json")

## ---- per-class FP error structure for one (plot,rung), baseline matcher ---
# Recomputes the baseline greedy match, finds the core false positives, tags each
# as near a MATCHED stem (over-seg) vs isolated, and attributes it to the nearest
# stem's crown_class. Returns named fp_near_<cls> / fp_iso_<cls> counts.
FP_CLASSES <- c("dominant", "codominant", "intermediate", "suppressed")
fp_by_class <- function(stems, det, cx, cy, ph, near_tol = NEAR_TOL) {
  out <- setNames(rep(0L, 2 * length(FP_CLASSES)),
                  c(paste0("fp_near_", FP_CLASSES), paste0("fp_iso_", FP_CLASSES)))
  in_reg  <- abs(det$x - cx) <= ph + BASE_TOL & abs(det$y - cy) <= ph + BASE_TOL
  detr    <- det[in_reg, , drop = FALSE]
  if (!nrow(detr) || !nrow(stems)) return(out)
  is_core <- abs(detr$x - cx) <= ph & abs(detr$y - cy) <= ph
  m <- greedy_match(stems$E, stems$N, detr$x, detr$y, BASE_TOL,
                    az = stems$height, bz = detr$z)
  matched <- m > 0
  fp_idx <- which(is_core & !(seq_along(detr$x) %in% m[matched]))
  if (!length(fp_idx)) return(out)
  mx <- stems$E[matched]; my <- stems$N[matched]
  for (i in fp_idx) {
    ds <- sqrt((stems$E - detr$x[i])^2 + (stems$N - detr$y[i])^2)
    cls <- stems$crown_class[which.min(ds)]            # nearest stem's class
    if (!cls %in% FP_CLASSES) next
    near <- length(mx) > 0 &&
      any(sqrt((mx - detr$x[i])^2 + (my - detr$y[i])^2) <= near_tol)
    key <- paste0(if (near) "fp_near_" else "fp_iso_", cls)
    out[key] <- out[key] + 1L
  }
  out
}

## ---- per-plot re-score over all matcher configs + the sensitivity grid ----
run_plot <- function(site, pid, pc, gt, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < MINTREES) return(NULL)
  tol_vec <- match_tol(stems$maxCrownDiameter, stems$pos_unc,
                       base_tol = BASE_TOL, k = K)
  rows <- list()
  for (rung in RUNGS) {
    fp <- frozen_norm_path(nd, site, pid, rung)
    if (!file.exists(fp)) next
    frdens <- tryCatch(jsonlite::read_json(frozen_manifest(nd, site, pid, rung),
                                           simplifyVector = TRUE)$frdens,
                       error = function(e) NA_real_)
    if (is.null(frdens) || is.na(frdens)) frdens <- 8        # safe default res rule
    res <- if (frdens >= 8) 0.25 else 0.5
    det <- tryCatch(detect_lasr(fp, res, VWF_A, frdens), error = function(e) NULL)
    if (is.null(det)) next
    fpc <- fp_by_class(stems, det, cx, cy, ph)               # baseline FP split
    score_one <- function(name, tol_xy, method, tz, lambda) {
      sc <- score_plot(stems, det, tol_xy = tol_xy, core_cx = cx, core_cy = cy,
                       core_half = ph, method = method, tol_z_up = tz,
                       lambda = lambda)
      cbind(data.frame(site = site, plot = pid, rung = rung, config = name,
                       matcher = method, tol_z_up = tz,
                       lambda = if (is.null(lambda)) NA_real_ else lambda,
                       frdens = round(frdens, 2), stringsAsFactors = FALSE), sc)
    }
    for (cfg in CONFIGS) {
      r <- score_one(cfg$name, if (cfg$scaled) tol_vec else BASE_TOL,
                     cfg$method, cfg$tz, cfg$lambda)
      r <- cbind(r, as.data.frame(as.list(fpc)))             # carry per-class FP
      rows[[length(rows) + 1]] <- r
    }
    for (tx in TX_SET) for (tz in TZ_SET) {
      r <- score_one(sprintf("sweep_tx%g_tz%g", tx, tz), tx, "greedy", tz, NULL)
      r <- cbind(r, as.data.frame(as.list(fpc)))
      rows[[length(rows) + 1]] <- r
    }
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows, fill = TRUE)
}

## ---- per-site driver ------------------------------------------------------
run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv")
  pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) {
    cat(sprintf("[%s] no ground truth / centroids -- skipped\n", site)); return(NULL) }
  gt <- read.csv(gtf, stringsAsFactors = FALSE)
  pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  gt <- gt[, setdiff(names(gt), "maxCrownDiameter"), drop = FALSE]
  gt <- merge(gt, field_crowns(site), by = "individualID", all.x = TRUE)
  if (is.null(gt$maxCrownDiameter)) gt$maxCrownDiameter <- NA_real_
  if (is.null(gt$pos_unc)) gt$pos_unc <- NA_real_

  plots <- intersect(unique(gt$plotID), pc$plotID)
  cat(sprintf("[%s] re-scoring %d plots over rungs {%s}\n",
              site, length(plots), paste(RUNGS, collapse = ",")))
  if (!length(plots)) return(NULL)
  res_list <- mclapply(plots, function(p)
    tryCatch(run_plot(site, p, pc, gt, nd),
             error = function(e) { message("  ", p, " failed: ",
                                            conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- rbindlist(Filter(Negate(is.null), res_list), fill = TRUE)
  if (!nrow(res)) { cat(sprintf("[%s] no cells scored\n", site)); return(NULL) }
  o <- file.path(nd, "matcher_robustness.csv")
  write.csv(res, o, row.names = FALSE)
  cat(sprintf("[%s] wrote %d cell rows -> %s\n", site, nrow(res), o))
  as.data.frame(res)
}

## ---- pooled report --------------------------------------------------------
print_report <- function(res) {
  cat("\n===== MATCHER ROBUSTNESS (pooled by SUM over plots x rungs) =====\n")
  cat(sprintf("\n%-16s %6s %6s %6s %6s %6s %6s %7s\n",
              "config", "n_ref", "n_det", "recall", "prec", "F1",
              "dF1", "fp_near%"))
  base <- pool(res[res$config == "baseline", , drop = FALSE])
  for (nm in vapply(CONFIGS, `[[`, "", "name")) {
    p <- pool(res[res$config == nm, , drop = FALSE])
    sub <- res[res$config == nm, , drop = FALSE]
    fn <- sum(sub$fp_near, na.rm = TRUE); fi <- sum(sub$fp_isolated, na.rm = TRUE)
    pct <- if (fn + fi > 0) 100 * fn / (fn + fi) else NA_real_
    cat(sprintf("%-16s %6d %6d %6.3f %6.3f %6.3f %+6.3f %6.1f\n",
                nm, p$n_ref, p$n_det, p$recall, p$precision, p$F1,
                p$F1 - base$F1, pct))
  }
  cat("\n--- tol_xy x tol_z_up sensitivity (greedy, pooled F1) ---\n")
  cat(sprintf("%-8s", "tol_xy")); for (tz in TZ_SET) cat(sprintf(" tz=%-4g", tz)); cat("\n")
  for (tx in TX_SET) {
    cat(sprintf("%-8g", tx))
    for (tz in TZ_SET) {
      p <- pool(res[res$config == sprintf("sweep_tx%g_tz%g", tx, tz), , drop = FALSE])
      cat(sprintf(" %6.3f", p$F1))
    }
    cat("\n")
  }
  cat("\n--- baseline FP error structure (core FPs, per crown_class) ---\n")
  cat(sprintf("%-14s %8s %8s\n", "class", "near", "isolated"))
  b <- res[res$config == "baseline", , drop = FALSE]
  for (cl in FP_CLASSES)
    cat(sprintf("%-14s %8d %8d\n", cl,
                sum(b[[paste0("fp_near_", cl)]], na.rm = TRUE),
                sum(b[[paste0("fp_iso_", cl)]],  na.rm = TRUE)))
}

run_main <- function() {
  t0 <- Sys.time(); all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- if (length(all_res)) do.call(rbind, all_res) else NULL
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d cell rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) print_report(res)
}

if (sys.nframe() == 0L) run_main()
