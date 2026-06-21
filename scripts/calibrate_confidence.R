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

# #P4 per-detection confidence calibration (GitHub issue #66).
#
# Weighted consensus in #P1 needs each detection to carry a TRUSTWORTHY
# probability so a high-confidence SegmentAnyTree mask can outvote a low-
# confidence apex -- but raw scores are not comparable across arms (a CNN's 0.6
# != a watershed prominence of 0.6), and the repo keeps every detection
# unweighted. This arm makes scores cross-arm-comparable: it materializes each
# arm's detections with a raw confidence, labels each TP/FP against field stems,
# builds a per-arm reliability diagram + Expected Calibration Error, fits a
# post-hoc isotonic calibrator (stats::isoreg) mapping raw score -> empirical
# precision, and reports the ensemble payoff (precision at fixed recall when the
# pooled multi-arm detections are ranked by calibrated vs raw score).
#
# Per-arm raw confidence (the issue's prescription):
#   forestformer3d -- the NATIVE per-instance mask score: mean ff3d_score over
#                     the instance's points (ff3d_arm.py writes the ff3d_score
#                     extra dim; the R loaders ignored it -- this exposes it).
#   segmentanytree -- crown point count (run_segmentanytree.py writes no score
#                     field), a size proxy: bigger masks are more often real.
#   chm_vwf/multichm/li2012 -- apex CHM height (AGL): taller apices are more
#                     often real dominant trees (the classical-arm proxy).
# Each arm's raw score is min-max normalized per arm to a predicted probability
# in [0,1] (the uncalibrated p); isotonic recalibration then maps it to precision.
#
# Detections are materialized from the SAME frozen cells as #P1 (native, equal
# set), labelled TP/FP with greedy_match restricted to the plot core (mirroring
# score_plot's precision denominator). A detection is TP if it is the matched
# detection for some core stem, else (in-core, unmatched) FP; buffer detections
# are ignored.
#
# Usage:
#   Rscript scripts/calibrate_confidence.R SITES=SOAP,SJER,TEAK
# Reads: work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv}, the frozen
#   normalized clips + DTMs, and the persisted deep instance clouds.
# Writes: work/neon/<SITE>/confidence_calibration.csv (one row per detection:
#   site,plot,arm,raw,prob,label) and work/neon/<SITE>/confidence_lookup.csv (the
#   per-arm isotonic knots fuse_detectors() consumes to weight votes).
suppressMessages({ library(lidR); library(data.table); library(parallel); library(sf) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))
source(.find("io_bridge.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
CORES <- as.integer(if (is.null(A$CORES)) 1 else A$CORES)
RUNG  <- if (is.null(A$RUNG)) "native" else A$RUNG
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
VWF_A <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
BINS  <- as.integer(if (is.null(A$BINS)) 10 else A$BINS)
MERGE_TOL <- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
SAT_ID_FIELD <- "PredInstance"
ARMS <- c("chm_vwf", "multichm", "li2012", "segmentanytree", "forestformer3d")

fz <- function(nd, site, pid, f) file.path(nd, "frozen", site, pid, RUNG, f)
inst_path <- function(idir, pid) {
  cand <- file.path(idir, paste0(pid, "_", RUNG, c(".laz", ".las")))
  hit <- cand[file.exists(cand)]; if (length(hit)) hit[1] else NA_character_
}

## ---- per-arm detections WITH a raw confidence -> data.frame(x,y,z,score) ---
# All apex z are AGL (normalized clip, or DTM-converted) so greedy_match's height
# gate is apples-to-apples. NULL when the arm has no artefact on this cell.
arm_dets <- function(arm, las, clip, dtm, frdens, res, nd, pid) {
  if (arm == "chm_vwf") {
    det <- tryCatch(detect_lasr(clip, res, VWF_A, frdens), error = function(e) NULL)
    if (is.null(det) || !nrow(det)) return(NULL)
    return(data.frame(x = det$x, y = det$y, z = det$z, score = det$z))  # apex height
  }
  if (arm == "multichm") {
    tt <- tryCatch(lidR::locate_trees(las, lidRplugins::multichm(res = res,
                     ws = ws_factory(VWF_A))), error = function(e) NULL)
    if (is.null(tt) || !nrow(tt)) return(NULL)
    co <- sf::st_coordinates(tt)
    z <- if (ncol(co) >= 3) co[, 3] else tt$Z
    return(data.frame(x = co[, 1], y = co[, 2], z = z, score = z))
  }
  if (arm == "li2012") {
    if (RUNG != "native") return(NULL)
    seg <- tryCatch(lidR::segment_trees(las,
             lidR::li2012(dt1 = 1.5, dt2 = 2, R = 2, hmin = 2)), error = function(e) NULL)
    if (is.null(seg) || !"treeID" %in% names(seg@data)) return(NULL)
    ap <- instance_apex(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
    if (!nrow(ap)) return(NULL)
    return(data.frame(x = ap$x, y = ap$y, z = ap$z, score = ap$z))
  }
  if (arm == "segmentanytree") {
    if (!file.exists(dtm)) return(NULL)
    sp <- inst_path(file.path(nd, "segmentanytree_instances"), pid)
    if (is.na(sp)) return(NULL)
    pts <- tryCatch(read_instance_points_laz(sp, id_field = SAT_ID_FIELD),
                    error = function(e) NULL)
    if (is.null(pts) || !nrow(pts)) return(NULL)
    pts <- pts[!is.na(pts$crown_id), , drop = FALSE]
    if (!nrow(pts)) return(NULL)
    dt <- as.data.table(pts)
    cnt <- dt[, .(npts = .N), by = crown_id]                  # crown point count = score
    ap  <- instance_apex(pts, id_col = "crown_id")
    ap  <- merge(ap, cnt, by.x = "id", by.y = "crown_id")
    agl <- det_to_agl(data.frame(x = ap$x, y = ap$y, z = ap$z), dtm)
    if (!nrow(agl)) return(NULL)
    ap  <- ap[ap$id %in% ap$id, ]                             # keep order
    # det_to_agl drops off-DTM rows; re-join by position
    m <- match(paste(round(agl$x, 3), round(agl$y, 3)),
               paste(round(ap$x, 3), round(ap$y, 3)))
    return(data.frame(x = agl$x, y = agl$y, z = agl$z, score = ap$npts[m]))
  }
  if (arm == "forestformer3d") {
    if (!RUNG %in% c("native", "8") || !file.exists(dtm)) return(NULL)
    fp <- inst_path(file.path(nd, "forestformer3d_instances"), pid)
    if (is.na(fp)) return(NULL)
    las2 <- tryCatch(suppressWarnings(lidR::readLAS(fp)), error = function(e) NULL)
    if (is.null(las2) || lidR::is.empty(las2)) return(NULL)
    dat <- las2@data
    if (!all(c("UserData", "PointSourceID", "ff3d_score") %in% names(dat))) return(NULL)
    pts <- data.frame(block = as.integer(dat$UserData),
                      inst = as.integer(dat$PointSourceID),
                      X = dat$X, Y = dat$Y, Z = dat$Z, score = dat$ff3d_score)
    rel <- dedup_blocks(pts[, c("block", "inst", "X", "Y", "Z")], merge_tol = MERGE_TOL)
    if (!nrow(rel)) return(NULL)
    rel <- as.data.table(rel)
    # mean native ff3d_score per global instance, joined by point identity
    rel[, sc := pts$score[match(paste(round(X, 3), round(Y, 3), round(Z, 3)),
                                paste(round(pts$X, 3), round(pts$Y, 3), round(pts$Z, 3)))]]
    agg <- rel[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z),
                   score = mean(sc, na.rm = TRUE)), by = global_id]
    agl <- det_to_agl(data.frame(x = agg$x, y = agg$y, z = agg$z), dtm)
    if (!nrow(agl)) return(NULL)
    m <- match(paste(round(agl$x, 3), round(agl$y, 3)),
               paste(round(agg$x, 3), round(agg$y, 3)))
    return(data.frame(x = agl$x, y = agl$y, z = agl$z, score = agg$score[m]))
  }
  NULL
}

## ---- label a det set TP/FP against the plot-core stems --------------------
# Mirrors score_plot's precision denominator: greedy_match stems -> dets; a
# CORE detection that is some stem's matched det is TP, else FP. Buffer dets out.
label_dets <- function(det, stems, cx, cy, ph) {
  if (is.null(det) || !nrow(det)) return(NULL)
  is_core <- abs(det$x - cx) <= ph & abs(det$y - cy) <= ph
  if (!any(is_core)) return(NULL)
  m <- greedy_match(stems$E, stems$N, det$x, det$y, TOL,
                    az = stems$height, bz = det$z)
  matched_det <- m[m > 0]
  lab <- integer(nrow(det))
  lab[matched_det] <- 1L
  data.frame(score = det$score[is_core], label = lab[is_core])
}

## ---- per-site driver ------------------------------------------------------
run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv"); pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) return(NULL)
  gt <- read.csv(gtf, stringsAsFactors = FALSE); pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  plots <- intersect(unique(gt$plotID), pc$plotID)
  one_plot <- function(pid) {
    ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
    if (!nrow(stems)) return(NULL)
    clip <- fz(nd, site, pid, "clip_normalized.laz")
    if (!file.exists(clip)) return(NULL)
    las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
    if (is.null(las) || lidR::is.empty(las)) return(NULL)
    frdens <- tryCatch(jsonlite::read_json(fz(nd, site, pid, "manifest.json"),
                       simplifyVector = TRUE)$frdens, error = function(e) NA_real_)
    if (is.null(frdens) || is.na(frdens)) frdens <- 8
    res <- if (frdens >= 8) 0.25 else 0.5
    dtm <- fz(nd, site, pid, "ground_dtm.tif")
    out <- list()
    for (arm in ARMS) {
      det <- tryCatch(arm_dets(arm, las, clip, dtm, frdens, res, nd, pid),
                      error = function(e) NULL)
      lab <- label_dets(det, stems, cx, cy, ph)
      if (!is.null(lab) && nrow(lab))
        out[[arm]] <- data.frame(site = site, plot = pid, arm = arm,
                                 raw = lab$score, label = lab$label)
    }
    if (!length(out)) return(NULL)
    rbindlist(out)
  }
  res <- rbindlist(Filter(Negate(is.null),
            mclapply(plots, function(p) tryCatch(one_plot(p),
              error = function(e) { message("  ", p, ": ", conditionMessage(e)); NULL }),
              mc.cores = CORES, mc.preschedule = FALSE)))
  if (!nrow(res)) { cat(sprintf("[%s] no detections labelled\n", site)); return(NULL) }
  res <- res[is.finite(res$raw), ]
  # per-arm min-max -> predicted prob
  res[, prob := { r <- range(raw); if (diff(r) > 0) (raw - r[1]) / diff(r) else 0.5 }, by = arm]
  write.csv(res, file.path(nd, "confidence_calibration.csv"), row.names = FALSE)
  cat(sprintf("[%s] labelled %d detections across %d arms -> confidence_calibration.csv\n",
              site, nrow(res), length(unique(res$arm))))
  as.data.frame(res)
}

## ---- out-of-sample calibrated probs via k-fold CV ------------------------
# Isotonic is in-sample-perfect (trained ECE -> 0), so report HELD-OUT
# calibration: assign each arm's detections to k folds, fit the calibrator on
# the other folds, predict the held-out fold. Returns res with an added `cal`
# column = out-of-sample calibrated prob (per arm). Deterministic (seeded folds).
oos_calibrated <- function(res, k = 5, seed = 42L) {
  res <- as.data.table(res); res[, cal := NA_real_]
  for (arm in unique(res$arm)) {
    ix <- which(res$arm == arm); n <- length(ix)
    if (n < k) { res$cal[ix] <- mean(res$label[ix]); next }   # too few -> base rate
    set.seed(seed); fold <- sample(rep_len(seq_len(k), n))
    for (f in seq_len(k)) {
      te <- ix[fold == f]; tr <- ix[fold != f]
      cal <- isotonic_calibrate(res$prob[tr], res$label[tr])
      res$cal[te] <- cal(res$prob[te])
    }
  }
  as.data.frame(res)
}

## ---- calibrate + report ---------------------------------------------------
report <- function(res) {
  oos <- oos_calibrated(res)
  cat("\n===== CONFIDENCE CALIBRATION (per arm, native equal-set) =====\n")
  cat("ECE_cal is 5-fold CROSS-VALIDATED (held-out), not in-sample.\n")
  cat(sprintf("%-16s %6s %7s %9s %9s %8s\n",
              "arm", "n", "base_P", "ECE_raw", "ECE_cal_cv", "dECE"))
  lookups <- list()
  for (arm in ARMS) {
    s <- res[res$arm == arm, , drop = FALSE]; if (!nrow(s)) next
    so <- oos[oos$arm == arm, , drop = FALSE]
    baseP <- mean(s$label)
    er <- expected_calibration_error(s$prob, s$label, bins = BINS)
    ec <- expected_calibration_error(so$cal, so$label, bins = BINS)   # held-out
    cat(sprintf("%-16s %6d %7.3f %9.3f %9.3f %+8.3f\n",
                arm, nrow(s), baseP, er, ec, ec - er))
    # the SHIPPED lookup is fit on ALL of this arm's data (max data for #P1)
    cal <- isotonic_calibrate(s$prob, s$label); kx <- sort(unique(s$prob))
    lookups[[arm]] <- data.frame(arm = arm, raw_prob = kx, calibrated = cal(kx))
  }
  cat("\n--- ensemble precision @ fixed recall (pooled multi-arm, held-out) ---\n")
  cat(sprintf("%8s %10s %10s %8s\n", "recall", "P_raw", "P_calib", "gain"))
  for (tr in c(0.5, 0.6, 0.7, 0.8, 0.9)) {
    pr <- precision_at_recall(oos$prob, oos$label, tr)        # raw normalized
    pc <- precision_at_recall(oos$cal,  oos$label, tr)        # held-out calibrated
    cat(sprintf("%8.2f %10.3f %10.3f %+8.3f\n", tr, pr, pc, pc - pr))
  }
  do.call(rbind, lookups)
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
  if (is.null(res) || !nrow(res)) { cat("\nNo detections labelled.\n"); return(invisible()) }
  lk <- report(res)
  # write per-arm calibrated lookup to each site dir (fuse_detectors consumes it)
  for (site in names(all_res))
    write.csv(lk, file.path(d, "neon", site, "confidence_lookup.csv"), row.names = FALSE)
  cat(sprintf("\nDONE: %d labelled detections across %d sites in %.1f min\n",
              nrow(res), length(all_res), dt))
}

if (sys.nframe() == 0L) run_main()
