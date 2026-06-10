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

# Point-cloud-detector DENSITY LADDER (issue #38): native + 8 pts/m^2 only.
#
# Issue #6 (detect_pc_sweep.R) compared CHM-VWF vs three point-cloud detectors at
# NATIVE density only. This script adds the missing density-RESPONSE arm for the
# point-cloud detectors, but ONLY at the top two rungs -- native and 8 pts/m^2.
# The sparser rungs (4/2/1) are deliberately OUT OF SCOPE: a point segmenter
# below ~3 first-return/m^2 produces spurious maxima and empty segments
# (treetop-lasr-vs-lidr-comparison.md shows point methods over-detect there), so
# a 1 pts/m^2 Li 2012 number would be noise, not signal. The full CHM-VWF ladder
# (run_sweep.R) already covers 8/4/2/1.
#
# Detectors (identical to #6 -- same windows, same scorer):
#   chm_vwf       : detect_lasr(file, res, a, frdens) -- pit-filled CHM + VWF
#                   local_maximum_raster. res is density-derived per (plot,rung)
#                   (0.25 m if frdens>=8 else 0.5 m), per run_sweep.R's res_set.
#   lidr_lmf_pc   : lidR lmf on the POINT CLOUD.
#   lidr_li2012   : lidR Li 2012 point-cloud segmentation; apex = max-Z/treeID.
#   lasr_lmax_pc  : lasR point-based local_maximum (not local_maximum_raster).
# All three point-cloud extractors live in pc_detect_lib.R, shared with #6.
#
# Clip provider: frozen_clip (model_bench_lib.R). Chosen over the prepare_clip +
# inline-decimate path because its decimation is SEEDED by (site,plot,rung) and
# cached, so (a) the 8-rung clip is reproducible across runs and arms, and (b)
# all four detectors are guaranteed to score the SAME bytes per (plot,rung) --
# they read one cached normalized clip. The shared cache under neon/<SITE>/frozen
# is the same one the model-benchmark arms use (same core_half=plot_half, same
# buffer=25), so the 8-rung clips are reused, not rebuilt.
#
# Pooling: the canonical sum-counts pool() + equal_set_guard() from
# model_bench_lib.R, keyed by (site, plot, rung). recall = sum(TP)/sum(n_ref);
# per-class TP recovered as round(rec_class * n_class); understory =
# intermediate + suppressed. Deltas vs chm_vwf are reported WITHIN each rung over
# an identical plot population (the guard drops any (plot,rung) cell an arm
# failed on).
#
# Usage:
#   Rscript scripts/detect_pc_ladder.R [SITES=SJER,SOAP,TEAK] [CORES=6]
#                                      [TOL=4] [A=0.10] [RUNG=8]
# Output:
#   $CLAUDE_JOB_DIR/neon/<SITE>/pc_detect_ladder_results.csv  (one row per
#       plot x rung x detector)
#   $CLAUDE_JOB_DIR/neon/pc_detect_ladder_pooled.csv          (rung x detector)
suppressMessages({ library(lidR); library(lasR); library(sf)
                   library(data.table); library(parallel) })
options(lidR.progress = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("pc_detect_lib.R"))
source(.find("model_bench_lib.R"))

## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER","SOAP","TEAK") else
           strsplit(A$SITES, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 6 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL))  4.0  else A$TOL)
A_VWF <- as.numeric(if (is.null(A$A))    0.10 else A$A)
RUNG  <- as.numeric(if (is.null(A$RUNG)) 8    else A$RUNG)
MINTREES <- 6
ARMS  <- c("chm_vwf","lidr_lmf_pc","lidr_li2012","lasr_lmax_pc")

## ---- run all 4 arms on ONE frozen (plot,rung) clip, score each -----------
# prep: a frozen_clip() result. stems: field stems in the plot core. Returns up
# to 4 scored rows (one per arm), or NULL if the clip is unusable.
run_cell <- function(prep, stems, cx, cy, ph, SITE, pid, plotType, rung_lbl) {
  las <- tryCatch(readLAS(prep$normalized), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)
  frdens <- prep$frdens; pdens <- prep$pdens
  ws  <- ws_factory(A_VWF)
  # CHM resolution is density-derived, NOT hardcoded (run_sweep.R res_set): finest
  # 0.25 m where first-return density supports it (>=8 pts/m^2), else 0.5 m. The
  # point-cloud arms use the ws window, not a raster res, so this is chm_vwf-only.
  res <- if (frdens >= 8) 0.25 else 0.5

  dets <- list(); tim <- c()
  t0 <- Sys.time()
  dets$chm_vwf <- tryCatch(detect_lasr(prep$normalized, res, A_VWF, frdens),
                           error = function(e) NULL)
  tim["chm_vwf"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lidr_lmf_pc <- tryCatch(det_lidr_lmf_pc(las, ws), error = function(e) NULL)
  tim["lidr_lmf_pc"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lidr_li2012 <- tryCatch(det_lidr_li2012(las), error = function(e) NULL)
  tim["lidr_li2012"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  t0 <- Sys.time()
  dets$lasr_lmax_pc <- tryCatch(det_lasr_lmax_pc(prep$normalized, ws),
                                error = function(e) NULL)
  tim["lasr_lmax_pc"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out <- list()
  for (nm in ARMS) {
    det <- dets[[nm]]
    if (is.null(det)) next
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx,
                              core_cy = cy, core_half = ph),
                   error = function(e) NULL)
    if (is.null(sc)) next
    out[[length(out) + 1]] <- cbind(
      data.frame(site = SITE, plot = pid, plotType = plotType, detector = nm,
                 rung = rung_lbl, frdens = round(frdens, 2),
                 pdens = round(pdens, 2),
                 chm_res = if (nm == "chm_vwf") res else NA_real_,
                 n_apex = nrow(det), secs = round(tim[[nm]], 2)),
      sc)
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

## ---- per-plot worker: native + (8 if not upsampling) ---------------------
run_plot <- function(pid, gt, pc, ctg, SITE, froot) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
  if (nrow(stems) < 1) return(NULL)

  # native first -> gives native all-return density for the no-upsampling guard
  np <- tryCatch(frozen_clip(ctg, SITE, pid, NA, cx, cy, ph, out_root = froot),
                 error = function(e) NULL)
  if (is.null(np)) return(NULL)
  rungs <- c(NA, pc_rungs_for(np$pdens, RUNG))   # native always; RUNG if dense

  rows <- list()
  for (rung in rungs) {
    prep <- if (is.na(rung)) np else
      tryCatch(frozen_clip(ctg, SITE, pid, rung, cx, cy, ph, out_root = froot),
               error = function(e) NULL)
    if (is.null(prep)) next
    lbl <- if (is.na(rung)) "native" else as.character(rung)
    r <- run_cell(prep, stems, cx, cy, ph, SITE, pid, ci$plotType, lbl)
    if (!is.null(r)) rows[[length(rows) + 1]] <- r
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- per-site driver ------------------------------------------------------
run_site <- function(SITE) {
  nd  <- file.path(d, "neon", SITE)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"),
                  stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"),
                  stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  froot <- file.path(nd, "frozen")

  counts <- table(gt$plotID)
  keep   <- intersect(names(counts)[counts >= MINTREES], pc$plotID)
  cat(sprintf("[%s] plots: %d (%s)\n", SITE, length(keep),
              paste(keep, collapse = ",")))
  if (!length(keep)) return(NULL)

  t0 <- Sys.time()
  res_list <- mclapply(keep, function(p)
                tryCatch(run_plot(p, gt, pc, ctg, SITE, froot),
                         error = function(e) {
                           message("plot ", p, " failed: ",
                                   conditionMessage(e)); NULL }),
                mc.cores = CORES, mc.preschedule = FALSE)
  results <- do.call(rbind, Filter(Negate(is.null), res_list))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(results)) { cat(sprintf("[%s] no results\n", SITE)); return(NULL) }

  results$tp_core <- round(results$precision * results$n_det)
  out <- file.path(nd, "pc_detect_ladder_results.csv")
  write.csv(results, out, row.names = FALSE)
  cat(sprintf("[%s] DONE: %d rows (%d plots, rungs %s) in %.1f min -> %s\n",
              SITE, nrow(results), length(unique(results$plot)),
              paste(sort(unique(results$rung)), collapse = "/"), dt, out))
  results
}

## ---- pooled report: rung x detector --------------------------------------
report <- function(all) {
  all <- equal_set_guard(all, arms = ARMS)     # drop cells missing any arm
  dropped <- attr(all, "dropped")
  cat("\n=== Equal-plot-set guard (per site x plot x rung) ===\n")
  cat(sprintf("(site,plot,rung) cells scored by all %d arms: %d\n",
              length(ARMS), length(unique(paste(all$site, all$plot, all$rung)))))
  if (length(dropped)) {
    cat(sprintf("DROPPED %d cell(s) where >=1 arm failed/empty:\n",
                length(dropped)))
    for (k in dropped) cat(sprintf("  %s\n", k))
  } else cat("DROPPED 0 cells: every arm scored every (plot,rung).\n")
  if (!nrow(all)) { cat("no common cells to pool\n"); return(invisible()) }

  pooled_all <- list()
  for (rg in c("native", as.character(RUNG))) {
    sub <- all[all$rung == rg, ]
    if (!nrow(sub)) next
    pr <- do.call(rbind, lapply(ARMS, function(a) {
      s <- sub[sub$detector == a, ]; if (!nrow(s)) return(NULL)
      cbind(rung = rg, detector = a, pool(s),
            secs_med = round(median(s$secs), 2)) }))
    np <- length(unique(paste(sub$site, sub$plot)))
    cat(sprintf("\n=== Rung %s: pooled across all sites (%d plots) ===\n",
                rg, np))
    print(pr[, c("detector","n_ref","n_det","recall","precision","F1",
                 "rec_intermediate","rec_suppressed","rec_understory",
                 "n_understory","secs_med")],
          row.names = FALSE, digits = 3)
    base <- pr[pr$detector == "chm_vwf", ]
    cat(sprintf("  -- delta vs chm_vwf (rung %s) --\n", rg))
    for (i in seq_len(nrow(pr))) {
      p <- pr[i, ]
      cat(sprintf("    %-13s  recall %+0.3f  precision %+0.3f  understory %+0.3f\n",
                  p$detector, p$recall - base$recall,
                  p$precision - base$precision,
                  p$rec_understory - base$rec_understory))
    }
    pooled_all[[rg]] <- pr
  }

  cat("\n=== Understory (intermediate+suppressed) recall by site x rung x detector ===\n")
  for (S in sort(unique(all$site))) for (rg in c("native", as.character(RUNG))) {
    sub <- all[all$site == S & all$rung == rg, ]
    if (!nrow(sub)) next
    ps <- do.call(rbind, lapply(ARMS, function(a) {
      s <- sub[sub$detector == a, ]; if (!nrow(s)) return(NULL)
      data.frame(detector = a, recall = round(pool(s)$recall, 3),
                 rec_understory = round(pool(s)$rec_understory, 3),
                 n_understory = pool(s)$n_understory) }))
    cat(sprintf("-- %s / rung %s --\n", S, rg))
    print(ps, row.names = FALSE)
  }

  pooled <- do.call(rbind, pooled_all)
  outp <- file.path(d, "neon", "pc_detect_ladder_pooled.csv")
  write.csv(pooled, outp, row.names = FALSE)
  cat(sprintf("\npooled -> %s\n", outp))
  invisible(pooled)
}

run_main <- function() {
  all <- do.call(rbind, Filter(Negate(is.null), lapply(SITES, run_site)))
  if (is.null(all) || !nrow(all)) { cat("no results across any site\n"); return(invisible()) }
  report(all)
}

if (sys.nframe() == 0L) run_main()
