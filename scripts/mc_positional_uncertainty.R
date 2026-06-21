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

# #V3 Monte-Carlo stem-position uncertainty (GitHub issue #67).
#
# Every leaderboard delta the router (#P2) would switch arms on is a single point
# estimate against stem coordinates that carry real uncertainty -- neon_ground_
# truth.R stores `pos_unc` (NEON coordinate uncertainty + 0.3 m TruPulse) per
# stem, yet score_plot ignores it. This puts confidence bands on recall/precision/
# F1 (and per-crown-class recall) by re-scoring every arm -- AND the #P1 fusion
# union/layered arm -- under K reproducible draws of the field positions, so we can
# see which arm-vs-arm gaps survive ground-truth jitter and whether fusion is MORE
# stable than any single arm.
#
# Efficient design: detections (including the fused sets) are materialized ONCE per
# cell (the expensive step); each of K draws then only perturbs the stems
# (perturb_positions, sigma = pos_unc) and re-runs the cheap greedy_match/
# score_plot, pooling by the canonical summed-count rule. Reproducible: draw k uses
# seed = seed_for(site, plot, "native") + k, so a re-run reproduces the bands.
#
# Usage:
#   Rscript scripts/mc_positional_uncertainty.R SITES=SOAP,SJER,TEAK K=200
#   Rscript scripts/mc_positional_uncertainty.R SITE=SOAP K=300 TOLS=3,4,5
# Reads: ground_truth_stems.csv (with pos_unc), plot_centroids.csv, the frozen
#   normalized clips + DTMs, and the persisted deep instance clouds.
# Writes: work/neon/<SITE>/positional_uncertainty.csv (per arm x tol: median +
#   5th/95th-percentile bands on recall/precision/F1 + per-class recall).
suppressMessages({ library(lidR); library(data.table); library(parallel); library(sf) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))
source(.find("io_bridge.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else c("SOAP", "SJER", "TEAK")
K      <- as.integer(if (is.null(A$K)) 200 else A$K)
RUNG   <- if (is.null(A$RUNG)) "native" else A$RUNG
TOLS   <- as.numeric(if (is.null(A$TOLS)) 4 else strsplit(A$TOLS, ",")[[1]])
VWF_A  <- as.numeric(if (is.null(A$VWF_A)) 0.10 else A$VWF_A)
MERGE_TOL <- as.numeric(if (is.null(A$MERGE_TOL)) 2.0 else A$MERGE_TOL)
Z_TOL  <- as.numeric(if (is.null(A$Z_TOL)) 5.0 else A$Z_TOL)
# CORES parallelizes the per-cell materialization (prep_cell). KEEP IT AT 1:
# prep_cell -> materialize -> detect_lasr uses lasR::exec, which DEADLOCKS under
# mclapply fork (the same hazard matcher_robustness.R / fuse_detectors.R guard).
CORES  <- as.integer(if (is.null(A$CORES)) 1 else A$CORES)
SAT_ID_FIELD <- "PredInstance"
CHM_ARMS <- c("chm_vwf", "multichm")
CLASSES  <- c("dominant", "codominant", "intermediate", "suppressed")

fz <- function(nd, site, pid, f) file.path(nd, "frozen", site, pid, RUNG, f)
inst_path <- function(idir, pid) {
  cand <- file.path(idir, paste0(pid, "_", RUNG, c(".laz", ".las")))
  hit <- cand[file.exists(cand)]; if (length(hit)) hit[1] else NA_character_
}

## ---- materialize the fixed per-cell detections (single arms + fusion) -----
materialize <- function(las, clip, dtm, frdens, res, nd, site, pid) {
  out <- list()
  out$chm_vwf  <- tryCatch(detect_lasr(clip, res, VWF_A, frdens), error = function(e) NULL)
  out$multichm <- tryCatch({
    tt <- lidR::locate_trees(las, lidRplugins::multichm(res = res, ws = ws_factory(VWF_A)))
    co <- sf::st_coordinates(tt)
    data.frame(x = co[, 1], y = co[, 2], z = if (ncol(co) >= 3) co[, 3] else tt$Z)
  }, error = function(e) NULL)
  out$li2012 <- tryCatch({
    seg <- lidR::segment_trees(las, lidR::li2012(dt1 = 1.5, dt2 = 2, R = 2, hmin = 2))
    reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
  }, error = function(e) NULL)
  if (file.exists(dtm)) {
    sp <- inst_path(file.path(nd, "segmentanytree_instances"), pid)
    if (!is.na(sp)) out$segmentanytree <- tryCatch({
      det <- read_instances_laz(sp, id_field = SAT_ID_FIELD)
      if (is.null(det)) NULL else det_to_agl(det, dtm) }, error = function(e) NULL)
    fp <- inst_path(file.path(nd, "forestformer3d_instances"), pid)
    if (!is.na(fp)) out$forestformer3d <- tryCatch({
      det <- ff3d_collapse(fp, merge_tol = MERGE_TOL)
      if (is.null(det)) NULL else agl_guard(det, dtm) }, error = function(e) NULL)
  }
  out <- out[!vapply(out, function(z) is.null(z) || !nrow(z), logical(1))]
  if (length(out) < 2) return(NULL)
  # #P1 fusion union + layered over the single arms
  stack <- rbindlist(lapply(names(out), function(a)
    data.table(arm = a, x = out[[a]]$x, y = out[[a]]$y, z = out[[a]]$z)))
  fp <- fusion_points(stack$arm, stack$x, stack$y, stack$z, n_arms = length(out),
                      merge_tol = MERGE_TOL, z_tol = Z_TOL, chm_arms = CHM_ARMS)
  out[["fusion_union"]]   <- fp$union[, c("x", "y", "z")]
  out[["fusion_layered"]] <- fp$layered[, c("x", "y", "z")]
  out
}

## ---- per-cell fixed detections + stems, for the MC loop -------------------
prep_cell <- function(site, pid, pc, gt, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < 1) return(NULL)
  clip <- fz(nd, site, pid, "clip_normalized.laz"); if (!file.exists(clip)) return(NULL)
  las <- tryCatch(suppressWarnings(lidR::readLAS(clip)), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)
  frdens <- tryCatch(jsonlite::read_json(fz(nd, site, pid, "manifest.json"),
                     simplifyVector = TRUE)$frdens, error = function(e) NA_real_)
  if (is.null(frdens) || is.na(frdens)) frdens <- 8
  res <- if (frdens >= 8) 0.25 else 0.5
  dets <- materialize(las, clip, fz(nd, site, pid, "ground_dtm.tif"), frdens, res, nd, site, pid)
  if (is.null(dets)) return(NULL)
  list(site = site, plot = pid, cx = cx, cy = cy, ph = ph, stems = stems,
       dets = dets, seed = seed_for(site, pid, RUNG))
}

## ---- driver ---------------------------------------------------------------
run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv"); pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) return(NULL)
  gt <- read.csv(gtf, stringsAsFactors = FALSE); pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  if (is.null(gt$pos_unc)) gt$pos_unc <- NA_real_
  plots <- intersect(unique(gt$plotID), pc$plotID)
  cells <- Filter(Negate(is.null), mclapply(plots, function(p)
    tryCatch(prep_cell(site, p, pc, gt, nd), error = function(e) NULL),
    mc.cores = CORES, mc.preschedule = FALSE))
  if (!length(cells)) { cat(sprintf("[%s] no cells\n", site)); return(NULL) }
  arms <- sort(unique(unlist(lapply(cells, function(c) names(c$dets)))))
  cat(sprintf("[%s] %d cells, %d arms, K=%d draws, tol={%s}\n",
              site, length(cells), length(arms), K, paste(TOLS, collapse = ",")))

  # K draws x tol: each draw perturbs stems per cell, scores every arm, pools.
  out <- list()
  for (tol in TOLS) {
    draws <- vector("list", K)
    for (k in seq_len(K)) {
      rows <- list()
      for (cell in cells) {
        st <- cell$stems
        pj <- if (k == 0L) list(E = st$E, N = st$N) else
          perturb_positions(st$E, st$N, st$pos_unc, seed = cell$seed + k)
        sj <- st; sj$E <- pj$E; sj$N <- pj$N
        for (a in names(cell$dets)) {
          sc <- tryCatch(score_plot(sj, cell$dets[[a]], tol_xy = tol,
                  core_cx = cell$cx, core_cy = cell$cy, core_half = cell$ph),
                  error = function(e) NULL)
          if (!is.null(sc)) {
            sc$site <- cell$site; sc$plot <- cell$plot; sc$rung <- RUNG; sc$arm <- a
            rows[[length(rows) + 1]] <- sc
          }
        }
      }
      long <- rbindlist(rows, fill = TRUE)
      draws[[k]] <- rbindlist(lapply(arms, function(a) {
        p <- tryCatch(pool(as.data.frame(long[long$arm == a])), error = function(e) NULL)
        if (is.null(p)) return(NULL)
        data.table(arm = a, draw = k, recall = p$recall, precision = p$precision,
                   F1 = p$F1, rec_dominant = p$rec_dominant, rec_codominant = p$rec_codominant,
                   rec_intermediate = p$rec_intermediate, rec_suppressed = p$rec_suppressed,
                   rec_understory = p$rec_understory)
      }), fill = TRUE)
    }
    dd <- rbindlist(draws, fill = TRUE); dd$tol <- tol
    out[[as.character(tol)]] <- dd
  }
  alld <- rbindlist(out, fill = TRUE)
  # summarize: median + 5/95 bands per (arm, tol, metric)
  metrics <- c("recall", "precision", "F1", "rec_dominant", "rec_codominant",
               "rec_intermediate", "rec_suppressed", "rec_understory")
  summ <- alld[, {
    s <- lapply(metrics, function(m) quantile(get(m), c(.05, .5, .95), na.rm = TRUE))
    setNames(as.list(unlist(s)), paste0(rep(metrics, each = 3), c("_p05", "_p50", "_p95")))
  }, by = .(arm, tol)]
  summ$site <- site
  write.csv(summ, file.path(nd, "positional_uncertainty.csv"), row.names = FALSE)
  cat(sprintf("[%s] wrote %d arm x tol band rows -> positional_uncertainty.csv\n",
              site, nrow(summ)))
  list(summ = summ, draws = alld[, .(site = site, arm, tol, draw, F1, rec_understory)])
}

ARM_ORDER <- c("chm_vwf", "multichm", "li2012", "segmentanytree", "forestformer3d",
               "fusion_union", "fusion_layered")
run_main <- function() {
  t0 <- Sys.time(); res <- list(); drw <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) {
      message("site ", site, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) { res[[site]] <- r$summ; drw[[site]] <- r$draws }
  }
  if (!length(res)) { cat("\nNo bands produced.\n"); return(invisible()) }
  draws <- rbindlist(drw, fill = TRUE)
  cat("\n===== POSITIONAL-UNCERTAINTY BANDS (F1, median [p05,p95], tol=4) =====\n")
  cat(sprintf("%-16s %s\n", "arm", "per site: F1 median [p05, p95]"))
  for (a in ARM_ORDER) {
    line <- sprintf("%-16s", a)
    for (site in names(res)) {
      s <- res[[site]]; row <- s[s$arm == a & s$tol == 4, ]
      if (nrow(row)) line <- paste0(line, sprintf("  %s: %.3f [%.3f,%.3f]",
                                    site, row$F1_p50, row$F1_p05, row$F1_p95))
    }
    cat(line, "\n")
  }
  # stability: which arm has the tightest F1 band (p95-p05), pooled across sites?
  cat("\n--- F1 band WIDTH (p95-p05) at tol=4, smaller = more stable ---\n")
  for (a in ARM_ORDER) {
    w <- sapply(names(res), function(site) { s <- res[[site]]
      row <- s[s$arm == a & s$tol == 4, ]; if (nrow(row)) row$F1_p95 - row$F1_p05 else NA })
    if (any(is.finite(w))) cat(sprintf("%-16s mean width %.3f\n", a, mean(w, na.rm = TRUE)))
  }
  cat(sprintf("\nDONE: %d sites, K=%d in %.1f min\n", length(res), K,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

if (sys.nframe() == 0L) run_main()
