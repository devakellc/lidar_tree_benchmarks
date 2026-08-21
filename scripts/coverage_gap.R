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

# #V5 coverage-gap crediting (GitHub issue #93).
#
# #V4 (matcher_robustness.R) found ~94% of core false positives are ISOLATED --
# real trees the NEON field map never recorded, not over-segmentation -- so every
# precision/F1 number in the benchmark is a biased lower bound. This study makes
# the bias measurable: an isolated core FP of a target arm is reclassified
# `probable-real` when arms from >= MIN_FAM OTHER modality families (chm / pc /
# deep / rgb, see coverage_lib.R::FAMILY_MAP) ALSO leave an isolated FP within
# CRED_R metres of it -- several independent systems all seeing a tree where the
# map has nothing. Credited FPs are removed from the precision DENOMINATOR
# (pool()'s precision_cred/F1_cred); raw precision stays alongside so the bias
# is visible, never hidden. Recall/TP are untouched: an unmapped tree can never
# become a TP, it just stops counting against the arm.
#
# Detections per arm come from the best_treetop_cache written by
# export_best_treetops_geojson.R -- each arm at its best tested configuration
# (one rung per site), i.e. exactly the cells behind the leaderboard -- plus the
# RGB arms' persisted crown boxes (deepforest_boxes/ tile CSVs,
# detectree2_boxes/<plot>.csv; apex z from the native frozen CHM as in
# detect_deepforest_sweep.R). Witness evidence is treated as rung-independent:
# a tree's physical existence does not depend on the decimation rung, so every
# arm testifies at its own best operating point. LADDER=1 additionally
# regenerates the canonical CHM-VWF detector (detect_lasr, density-derived res,
# a=0.10) on every frozen rung so the per-rung bias curve is measured on one arm
# (detector = "chm_vwf_ladder"; witnesses unchanged).
#
# A CRED_R x MIN_FAM sensitivity grid (fp_cred_r*_f* columns) is carried per
# cell so the report can show how the credit moves with the rule's strictness.
#
# Usage:
#   Rscript scripts/coverage_gap.R SITE=SOAP
#   Rscript scripts/coverage_gap.R SITES=SOAP,SJER,TEAK CORES=1 LADDER=1
# (CORES=1 default: the ladder path runs lasR exec, which drops dense cells
#  under fork -- see the repo memory on mclapply + lasR.)
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv},
#   best_treetop_cache/, deepforest_boxes/, detectree2_boxes/, and the frozen
#   normalized clips + manifests.
# Writes: work/neon/<SITE>/coverage_gap.csv (one row per site x plot x cell;
#   pool with pool() -- fp_credited feeds precision_cred/F1_cred).
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))     # pool
source(.find("coverage_lib.R"))        # FAMILY_MAP, readers, credit_isolated

## ---- args (KEY=VALUE positional) -----------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES    <- as.integer(if (is.null(A$CORES)) 1 else A$CORES)
TOL      <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
CRED_R   <- as.numeric(if (is.null(A$CRED_R)) 2.0 else A$CRED_R)
MIN_FAM  <- as.integer(if (is.null(A$MIN_FAM)) 2 else A$MIN_FAM)
VWF_A    <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
CHM_RES  <- as.numeric(if (is.null(A$CHM_RES)) 0.5 else A$CHM_RES)
LADDER   <- is.null(A$LADDER) || A$LADDER != "0"
MINTREES <- as.integer(if (is.null(A$MINTREES)) 1 else A$MINTREES)
SELECTION <- if (is.null(A$SELECTION))
  file.path(d, "neon", "best_treetops_geojson", "best_treetop_selection.csv") else
  A$SELECTION
SENS_R   <- c(1.5, 2.0, 3.0)            # crediting-radius sensitivity
SENS_F   <- c(1L, 2L)                   # min-witness-family sensitivity

LIDAR_ARMS <- names(FAMILY_MAP)[!FAMILY_MAP %in% "rgb"]

frozen_norm_path <- function(nd, site, pid, rung)
  file.path(nd, "frozen", site, pid, rung, "clip_normalized.laz")
frozen_manifest  <- function(nd, site, pid, rung)
  file.path(nd, "frozen", site, pid, rung, "manifest.json")

# native frozen CHM for the optical arms' apex z (detect_deepforest_sweep.R)
plot_chm <- function(nd, site, pid) {
  clip <- frozen_norm_path(nd, site, pid, "native")
  if (!file.exists(clip)) return(NULL)
  las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)
  tryCatch(suppressWarnings(
    lidR::rasterize_canopy(las, res = CHM_RES, algorithm = lidR::p2r())),
    error = function(e) NULL)
}

# site-level DeepForest boxes (tile CSVs; filtered to a plot bbox by caller)
site_deepforest_boxes <- function(nd) {
  fs <- Sys.glob(file.path(nd, "deepforest_boxes", "*.csv"))
  if (!length(fs)) return(NULL)
  b <- rbindlist(lapply(fs, function(f)
    tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)),
    fill = TRUE)
  if (!nrow(b)) NULL else as.data.frame(b)
}

## ---- per-plot: materialize every cell, build the witness pool, credit ------
run_plot <- function(site, pid, pc, gt, nd, df_boxes, sel) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  pstems <- gt[gt$plotID == pid, , drop = FALSE]
  stems <- pstems[abs(pstems$E - cx) <= ph &
                  abs(pstems$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < MINTREES) return(NULL)

  ## 1. cells: every arm's cached best-configuration detections (+ optical).
  ## The selection manifest pins each arm to its selected best rung so stale
  ## extra cache files can never add cells; glob discovery is the fallback.
  cache_dir <- file.path(nd, "best_treetop_cache")
  cells <- list()
  for (arm in LIDAR_ARMS) {
    in_sel <- !is.null(sel) && arm %in% sel$method
    rungs <- if (in_sel) sel$rung[sel$method == arm] else
      cached_rungs(cache_dir, arm, site, pid)$rung
    # exact parameter pinning where the manifest carries it (chm_vwf)
    pr <- NULL
    if (in_sel && all(c("chm_res", "vwf_a") %in% names(sel))) {
      sr <- sel[sel$method == arm, ][1, ]
      if (is.finite(sr$chm_res) && is.finite(sr$vwf_a))
        pr <- c(sprintf("res%s", sr$chm_res), sprintf("a%s", sr$vwf_a))
    }
    for (rung in rungs) {
      det <- read_arm_cache(cache_dir, arm, site, pid, rung, params = pr)
      if (!is.null(det))
        cells[[length(cells) + 1]] <- list(arm = arm, rung = rung, det = det)
    }
  }
  ## Optical cells: a plot with boxes ANYWHERE on the site counts as scored --
  ## an empty in-plot box set is ran-and-found-nothing (0-row detections), so
  ## the plot's stems stay in the recall denominator exactly as in
  ## detect_deepforest_sweep.R. A missing detectree2_boxes/<plot>.csv means the
  ## arm never ran there (cell skipped), an empty one means it found nothing.
  chm <- plot_chm(nd, site, pid)
  if (!is.null(chm)) {
    if (!is.null(df_boxes)) {
      bp <- df_boxes[abs(df_boxes$x - cx) <= ph + TOL &
                     abs(df_boxes$y - cy) <= ph + TOL, , drop = FALSE]
      cells[[length(cells) + 1]] <-
        list(arm = "deepforest", rung = "rgb", det = boxes_to_dets(bp, chm))
    }
    d2f <- file.path(nd, "detectree2_boxes", paste0(pid, ".csv"))
    if (file.exists(d2f)) {
      bp <- tryCatch(read.csv(d2f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(bp))
        cells[[length(cells) + 1]] <-
          list(arm = "detectree2", rung = "rgb", det = boxes_to_dets(bp, chm))
    }
  }
  if (!length(cells)) return(NULL)

  ## 2. witness pool: every cell's core FPs, tagged with its modality family
  ## Eligibility is measured against the plot's FULL mapped stem set (pstems),
  ## not just the scored core: a stem whose position lands outside the core is
  ## still mapped, so a core FP sitting on it is explained by the field map.
  ## Scoring/`isolated` stay on the core stems (score_plot parity, n_ref).
  for (i in seq_along(cells))
    cells[[i]]$fps <- fp_points(stems, cells[[i]]$det, tol_xy = TOL,
                                core_cx = cx, core_cy = cy, core_half = ph,
                                elig_stems = pstems)
  wit <- rbindlist(lapply(cells, function(cl)
    if (nrow(cl$fps)) data.table(cl$fps, fam = arm_family(cl$arm)) else NULL),
    fill = TRUE)
  if (is.null(wit)) wit <- data.table(x = numeric(), y = numeric(),
                                      z = numeric(), isolated = logical(),
                                      fam = character())

  ## 3. score + credit each cell (cached arms + optical)
  score_cell <- function(arm, rung, det, fps) {
    fam <- arm_family(arm)
    if (is.na(fam)) fam <- "chm"                       # ladder rows are CHM-VWF
    sc <- score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                     core_half = ph)
    row <- cbind(
      data.frame(site = site, plot = pid, plotType = ci$plotType,
                 detector = arm, family = fam, rung = rung,
                 n_apex = nrow(det), stringsAsFactors = FALSE),
      sc)
    row$fp_credited <- credit_isolated(fps, wit, fam, r = CRED_R,
                                       min_fam = MIN_FAM)
    row$fp_elig <- sum(fps$credit_eligible)
    row$n_wit_fam <- length(unique(wit$fam[wit$credit_eligible & wit$fam != fam]))
    row$cred_r <- CRED_R; row$min_fam <- MIN_FAM
    for (r in SENS_R) for (f in SENS_F)
      row[[sprintf("fp_cred_r%g_f%d", r, f)]] <-
        credit_isolated(fps, wit, fam, r = r, min_fam = f)
    row
  }
  rows <- lapply(cells, function(cl) score_cell(cl$arm, cl$rung, cl$det, cl$fps))

  ## 4. CHM-VWF density ladder: regenerate the canonical detector per frozen
  ## rung so the bias curve is measured against density on one arm.
  if (LADDER) {
    rungs <- basename(Sys.glob(file.path(nd, "frozen", site, pid, "*")))
    for (rung in rungs) {
      clip <- frozen_norm_path(nd, site, pid, rung)
      if (!file.exists(clip)) next
      frdens <- tryCatch(jsonlite::read_json(frozen_manifest(nd, site, pid, rung),
                                             simplifyVector = TRUE)$frdens,
                         error = function(e) NA_real_)
      if (is.null(frdens) || is.na(frdens)) frdens <- 8
      res <- if (frdens >= 8) 0.25 else 0.5
      det <- tryCatch(detect_lasr(clip, res, VWF_A, frdens), error = function(e) NULL)
      if (is.null(det)) next
      fps <- fp_points(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                       core_half = ph, elig_stems = pstems)
      rows[[length(rows) + 1]] <- score_cell("chm_vwf_ladder", rung, det, fps)
    }
  }
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
  df_boxes <- site_deepforest_boxes(nd)
  sel <- read_selection(SELECTION, site)
  plots <- intersect(unique(gt$plotID), pc$plotID)
  cat(sprintf("[%s] crediting %d plots (rgb boxes: %s; selection: %s; ladder: %s)\n",
              site, length(plots),
              if (is.null(df_boxes)) "none" else nrow(df_boxes),
              if (is.null(sel)) "glob fallback" else sprintf("%d arms", nrow(sel)),
              if (LADDER) "on" else "off"))
  if (!length(plots)) return(NULL)
  res_list <- mclapply(plots, function(p)
    tryCatch(run_plot(site, p, pc, gt, nd, df_boxes, sel),
             error = function(e) { message("  ", p, " failed: ",
                                            conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- rbindlist(Filter(Negate(is.null), res_list), fill = TRUE)
  if (!nrow(res)) { cat(sprintf("[%s] no cells scored\n", site)); return(NULL) }
  o <- file.path(nd, "coverage_gap.csv")
  write.csv(res, o, row.names = FALSE)
  cat(sprintf("[%s] wrote %d cell rows -> %s\n", site, nrow(res), o))
  as.data.frame(res)
}

## ---- pooled report --------------------------------------------------------
print_report <- function(res) {
  lead <- res[res$detector != "chm_vwf_ladder", , drop = FALSE]
  cat(sprintf("\n===== COVERAGE-GAP CREDITING (pooled by SUM; CRED_R=%.1f, MIN_FAM=%d) =====\n",
              CRED_R, MIN_FAM))
  cat(sprintf("%-16s %-5s %5s %6s %6s %6s | %6s %6s | %7s %9s\n",
              "arm", "rung", "n_ref", "recall", "prec", "F1",
              "prec'", "F1'", "dF1", "cred/elig"))
  arms <- unique(lead$detector)
  tab <- list()
  for (a in arms) {
    sub <- lead[lead$detector == a, , drop = FALSE]
    p <- pool(sub)
    elig <- if (!is.null(sub$fp_elig)) sum(sub$fp_elig, na.rm = TRUE) else
      sum(sub$fp_isolated, na.rm = TRUE)
    cred <- sum(sub$fp_credited, na.rm = TRUE)
    tab[[a]] <- p
    cat(sprintf("%-16s %-5s %5d %6.3f %6.3f %6.3f | %6.3f %6.3f | %+6.3f %4d/%-4d\n",
                a, paste(unique(sub$rung), collapse = ","), p$n_ref, p$recall,
                p$precision, p$F1, p$precision_cred, p$F1_cred,
                p$F1_cred - p$F1, cred, elig))
  }
  f1_raw  <- vapply(tab, function(p) p$F1, numeric(1))
  f1_cred <- vapply(tab, function(p) p$F1_cred, numeric(1))
  o_raw  <- names(sort(f1_raw,  decreasing = TRUE))
  o_cred <- names(sort(f1_cred, decreasing = TRUE))
  cat("\nleaderboard raw : ", paste(o_raw,  collapse = " > "), "\n")
  cat("leaderboard cred: ",   paste(o_cred, collapse = " > "), "\n")
  if (!identical(o_raw, o_cred)) cat("** ordering CHANGES under de-biased precision **\n")
  # Arms cover different plot subsets (detectree2, the pc_ladder arms), so the
  # ordering above mixes denominators. Re-derive both orderings on the COMMON
  # plot set (equal_set_guard keyed by site x plot; rungs differ by design).
  ga <- equal_set_guard(lead, arms, key_cols = c("site", "plot"))
  if (nrow(ga)) {
    gtab <- lapply(setNames(arms, arms), function(a)
      pool(ga[ga$detector == a, , drop = FALSE]))
    g_raw  <- names(sort(vapply(gtab, function(p) p$F1, numeric(1)), decreasing = TRUE))
    g_cred <- names(sort(vapply(gtab, function(p) p$F1_cred, numeric(1)), decreasing = TRUE))
    n_common <- length(unique(paste(ga$site, ga$plot)))
    cat(sprintf("\nequal-set (all %d arms on %d common plots):\n", length(arms), n_common))
    cat("  raw : ", paste(g_raw,  collapse = " > "), "\n")
    cat("  cred: ", paste(g_cred, collapse = " > "), "\n")
  }

  lad <- res[res$detector == "chm_vwf_ladder", , drop = FALSE]
  if (nrow(lad)) {
    cat("\n--- CHM-VWF ladder: coverage-gap bias vs density rung ---\n")
    cat(sprintf("%-8s %5s %6s %6s %6s | %6s %6s %7s\n",
                "rung", "n_ref", "recall", "prec", "F1", "prec'", "F1'", "dF1"))
    for (rung in intersect(c("native", "8", "4", "2", "1"), unique(lad$rung))) {
      p <- pool(lad[lad$rung == rung, , drop = FALSE])
      cat(sprintf("%-8s %5d %6.3f %6.3f %6.3f | %6.3f %6.3f %+7.3f\n",
                  rung, p$n_ref, p$recall, p$precision, p$F1,
                  p$precision_cred, p$F1_cred, p$F1_cred - p$F1))
    }
  }

  cat("\n--- crediting-rule sensitivity (pooled dF1 over all leaderboard cells) ---\n")
  cat(sprintf("%-10s", "cred_r")); for (f in SENS_F) cat(sprintf(" min_fam=%d", f)); cat("\n")
  for (r in SENS_R) {
    cat(sprintf("%-10g", r))
    for (f in SENS_F) {
      v <- lead; v$fp_credited <- v[[sprintf("fp_cred_r%g_f%d", r, f)]]
      p <- pool(v)
      cat(sprintf("  %+8.3f", p$F1_cred - p$F1))
    }
    cat("\n")
  }
}

run_main <- function() {
  t0 <- Sys.time(); all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- if (length(all_res)) rbindlist(all_res, fill = TRUE) else NULL
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d cell rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) print_report(as.data.frame(res))
}

if (sys.nframe() == 0L) run_main()
