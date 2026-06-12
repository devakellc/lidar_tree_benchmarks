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

# 3-D instance-segmenter crown-diameter arm for GitHub issue #30.
#
# Issue #7 scores five CHM crown segmenters (dalponte2016, silva2016, watershed,
# lasR region_growing, random_walker) against NEON field crown diameter. The
# model-benchmark DETECTION arms (Li 2012, ptrees, AMS3D) already run on the SAME
# native frozen clips but were never crown-scored. This script closes that gap:
# it runs each 3-D point-instance segmenter at native density, reduces every
# instance to a convex-hull crown diameter (d_eq / d_caliper) via the bridge's
# crown_diameter_table(), matches instance apexes to field stems with the SAME
# greedy_match harness (4 m tol + height gate), and scores pooled RMSE/MAE/bias/R2
# by crown class against the SAME two field diameter definitions as #7.
#
# Segmenters (all on the native NORMALIZED frozen clip, no decimation):
#   - li2012  (lidR Li 2012 point segmentation -> seg@data treeID)
#   - ptrees  (lidRplugins PTrees / Vega 2014  -> seg@data treeID)
#   - ams3d   (crownsegmentr adaptive mean shift -> seg@data crown_id)
# The segmenter invocations are lifted verbatim from the detection arms so the
# crown geometry comes from the identical instance labelling that the detection
# benchmark scored: det_li2012 (detect_li2012_native.R), det_ptrees
# (detect_lidrplugins_sweep.R), det_ams3d (detect_ams3d_sweep.R). Here we keep
# the FULL seg@data (not just the apex set) to derive per-instance diameters.
#
# Two diameter estimates per instance (convex hull of the instance's points,
# crown_diameter_table min_pts=5; same geometric caveat as #7):
#   - d_eq      = 2*sqrt(hull_area/pi)   equivalent-circle -> ninetyCrownDiameter
#   - d_caliper = max pairwise point dist (widest axis)    -> maxCrownDiameter
# These are a POINT-hull estimator, NOT the dissolved-CHM polygon the #7 classical
# arms use -- not identical estimators, though both target the same field column
# (see the results doc caveat). area is the equivalent-circle area implied by
# d_eq (NA below min_pts).
#
# Usage:
#   Rscript scripts/crown_metrics_3d.R SITES=SJER,SOAP,TEAK CORES=8 TOL=4
#
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv}
#   the LiDAR catalog work/neon/<SITE>/lidar/, the cached field-crown widths
#   work/neon/<SITE>/vst/<site>_vst_allyears.rds, and the cached frozen clips
#   work/neon/<SITE>/frozen/<plot>/native/clip_normalized.laz (reused, never
#   regenerated here).
# Writes (NEW file, never overwrites crown_metrics_results.csv):
#   work/neon/<SITE>/crown_metrics_3d_results.csv  (one row per matched tree,
#   canonical cols: site, plot, algo, crown_class, individualID, d_eq,
#   d_caliper, area, field_maxCD, field_ninetyCD; algo in {li2012,ptrees,ams3d}).
suppressMessages({
  library(lidR); library(lidRplugins); library(crownsegmentr)
  library(data.table); library(parallel)
})
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))

## ---- args (KEY=VALUE positional, per repo convention) --------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 8 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL))  4 else A$TOL)
MINTREES <- 6

## ---- per-segmenter instance labelling on a normalized clip ---------------
# Each returns the SAME seg@data table the detection arm consumed, plus the id
# column name, so the caller can derive both crown_diameter_table and the
# per-instance apex from one labelling. NULL = segmenter crashed (skip plot);
# a labelled table with all-NA ids = ran-but-empty (legit no crowns). The
# invocations mirror det_li2012/det_ptrees/det_ams3d verbatim.
seg_li2012 <- function(las, hmin = 2, dt1 = 1.5, dt2 = 2, R = 2) {
  if (sum(las$Z >= hmin) < 1) return(list(data = data.table(), id = "treeID"))
  seg <- tryCatch(lidR::segment_trees(las,
                    lidR::li2012(dt1 = dt1, dt2 = dt2, R = R, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg) || !"treeID" %in% names(seg@data)) return(NULL)
  list(data = seg@data, id = "treeID")
}

# ptrees segfaults (uncatchable) on too-few canopy returns -- skip safely like
# det_ptrees, returning an empty-but-ran labelling so the cell is a legit miss.
seg_ptrees <- function(las, hmin = 2, k = c(30, 15)) {
  if (sum(las$Z >= hmin) < min(k)) return(list(data = data.table(), id = "treeID"))
  seg <- tryCatch(lidR::segment_trees(las,
                    lidRplugins::ptrees(k = k, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg) || !"treeID" %in% names(seg@data)) return(NULL)
  list(data = seg@data, id = "treeID")
}

seg_ams3d <- function(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2) {
  if (is.empty(las) || npoints(las) < 5)
    return(list(data = data.table(), id = "crown_id"))
  seg <- tryCatch(
    crownsegmentr::segment_tree_crowns(
      point_cloud = las,
      crown_diameter_to_tree_height = cd_ratio,
      crown_length_to_tree_height   = cl_ratio,
      segment_crowns_only_above     = min_above,
      ground_height                 = NULL,        # input already normalized
      crown_id_column_name          = "crown_id"),
    error = function(e) NULL)
  if (is.null(seg) || !"crown_id" %in% names(seg@data)) return(NULL)
  list(data = seg@data, id = "crown_id")
}

SEGMENTERS <- list(li2012 = seg_li2012, ptrees = seg_ptrees, ams3d = seg_ams3d)

## ---- field crown diameter per site (join from cached vst rds) ------------
# Identical to crown_metrics_sweep.R::field_crowns: dedup apparentindividual to
# the nearest-to-2021 measurement per individualID; keep max/ninetyCrownDiameter.
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst",
                   paste0(tolower(site), "_vst_allyears.rds"))
  dat <- readRDS(rds)
  ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID),
     c("individualID", "maxCrownDiameter", "ninetyCrownDiameter")]
}

## ---- per-plot 3-D crown benchmark ----------------------------------------
# Reads the cached native frozen normalized clip (NEVER regenerates it), runs
# each segmenter once, and for every segmenter derives crown_diameter_table +
# instance_apex from the SAME labelling, then scores against the plot-core stems
# via the shared score_crowns_against_field glue. Returns canonical-cols rows for
# all segmenters that produced a match, or NULL.
run_plot <- function(site, pid, ctg, pc, gt, fc, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
  if (nrow(stems) < 1) return(NULL)

  prep <- tryCatch(frozen_clip(ctg, site, pid, NA, cx, cy, ph,
                               out_root = file.path(nd, "frozen")),
                   error = function(e) NULL)
  if (is.null(prep)) return(NULL)
  las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)

  rows <- list()
  for (algo in names(SEGMENTERS)) {
    sg <- tryCatch(SEGMENTERS[[algo]](las), error = function(e) NULL)
    if (is.null(sg)) next                       # segmenter crashed -> skip arm
    if (!nrow(sg$data)) next                    # ran-but-empty -> no crowns
    diam <- crown_diameter_table(sg$data, id_col = sg$id, min_pts = 5)
    apex <- instance_apex(sg$data, id_col = sg$id, x = "X", y = "Y", z = "Z")
    r <- score_crowns_against_field(diam, apex, stems, fc, tol = TOL,
                                    site = site, plot = pid, algo = algo)
    if (nrow(r)) rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- per-site driver ------------------------------------------------------
# Mirrors crown_metrics_sweep.R::run_site: live & is_tree stems, authoritative
# rds crown-diameter join (drop any pre-existing CD columns first), >=MINTREES
# plots, mclapply over plots, write the NEW crown_metrics_3d_results.csv.
run_site <- function(site) {
  nd  <- file.path(d, "neon", site)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  gt  <- gt[, setdiff(names(gt),
                      c("maxCrownDiameter", "ninetyCrownDiameter")), drop = FALSE]
  fc  <- field_crowns(site)
  gt  <- merge(gt, fc, by = "individualID", all.x = TRUE)
  gt  <- gt[!is.na(gt$maxCrownDiameter) | !is.na(gt$ninetyCrownDiameter), ]

  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)

  counts <- table(gt$plotID)
  keep <- names(counts)[counts >= MINTREES]
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] plots with >=%d stems w/ field CD: %d (%s)\n",
              site, MINTREES, length(keep), paste(keep, collapse = ",")))
  if (!length(keep)) return(NULL)

  res_list <- mclapply(keep, function(p)
    tryCatch(run_plot(site, p, ctg, pc, gt, fc, nd),
             error = function(e) { message("  plot ", p, " failed: ",
                                            conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res) || !nrow(res)) {
    cat(sprintf("[%s] no crowns matched\n", site)); return(NULL) }
  out <- file.path(nd, "crown_metrics_3d_results.csv")
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("[%s] wrote %d matched-tree rows -> %s\n", site, nrow(res), out))
  res
}

## ---- scoring helpers (mirror crown_metrics_sweep.R) ----------------------
# Pooled error stats over matched trees: RMSE/MAE/bias/R2 (sum of squared errors
# / n, never a mean of per-plot rates). det = detected diameter, fld = field.
err_stats <- function(det, fld) {
  ok <- is.finite(det) & is.finite(fld)
  det <- det[ok]; fld <- fld[ok]; n <- length(det)
  if (n < 2) return(data.frame(n = n, rmse = NA, mae = NA, bias = NA, r2 = NA))
  e <- det - fld
  ss_res <- sum((fld - det)^2)
  ss_tot <- sum((fld - mean(fld))^2)
  data.frame(n = n,
             rmse = sqrt(mean(e^2)), mae = mean(abs(e)),
             bias = mean(e),
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_)
}

print_tables <- function(res) {
  algos <- unique(res$algo)
  cat("\n========= 3-D CROWN-DIAMETER ACCURACY (pooled, all sites) =========\n")
  for (defn in list(c("d_eq", "field_ninetyCD",
                      "equiv-circle d_eq vs ninetyCrownDiameter"),
                    c("d_caliper", "field_maxCD",
                      "max-caliper d_caliper vs maxCrownDiameter"))) {
    cat(sprintf("\n--- %s ---\n", defn[3]))
    cat(sprintf("%-12s %5s %7s %7s %7s %7s\n",
                "algo", "n", "rmse", "mae", "bias", "r2"))
    for (a in algos) {
      s <- err_stats(res[[defn[1]]][res$algo == a], res[[defn[2]]][res$algo == a])
      cat(sprintf("%-12s %5d %7.2f %7.2f %7.2f %7.3f\n",
                  a, s$n, s$rmse, s$mae, s$bias, s$r2))
    }
    cat("  by crown class (rmse / n):\n")
    for (a in algos) {
      cat(sprintf("    %-10s", a))
      for (cc in c("dominant", "codominant", "intermediate", "suppressed")) {
        sel <- res$algo == a & res$crown_class == cc
        s <- err_stats(res[[defn[1]]][sel], res[[defn[2]]][sel])
        cat(sprintf(" %s=%s/%d", substr(cc, 1, 4),
                    ifelse(is.na(s$rmse), "NA", sprintf("%.2f", s$rmse)), s$n))
      }
      cat("\n")
    }
  }
}

## ---- run ------------------------------------------------------------------
run_main <- function() {
  t0 <- Sys.time()
  all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- do.call(rbind, all_res)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d matched-tree rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) {
    cat("\nrows per algorithm:\n"); print(table(res$algo))
    cat("\nrows per crown class:\n")
    print(table(res$crown_class, useNA = "ifany"))
    print_tables(res)
  }
}

if (sys.nframe() == 0L) run_main()
