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

# #P5 classical Treeiso arm (GitHub issue #73).
#
# TreeisoNet (the LEARNED variant, #M7/#20) is in the repo; the classical
# cut-pursuit Treeiso is not. Fusion needs a strong NON-learned 3-D member to
# test whether the deep arms actually beat a good classical instance segmenter,
# and Treeiso's graph-cut errors are decorrelated from CHM/LM errors, so it adds
# diversity to the #P1 consensus pool. Treeiso was designed for dense TLS/ULS, so
# DEGRADED understory recall on sparse discrete-return ALS is EXPECTED and is the
# point of measuring it here.
#
# Runs the vendored Treeiso (external/treeiso/, MIT) as a CPU subprocess on the
# frozen RAW-with-ground clip per (plot, rung): run_treeiso.py removes ground,
# runs init/intermediate/final cut-pursuit, and writes per-point instance ids in
# a `treeiso` extra dim. This arm persists that labelled cloud to
# work/neon/<SITE>/treeiso_instances/<plot>_<rung>.laz (so #V1 IoU/PQ + #P1
# fusion can consume the per-point labels), reduces it to apexes via the bridge,
# converts the absolute-UTM apex z to AGL with the cached frozen ground_dtm.tif
# (Treeiso keeps absolute Z, like SegmentAnyTree), matches to field stems with
# score_plot, and writes treeiso_results.csv (detector = "treeiso").
#
# Usage:
#   Rscript scripts/detect_treeiso_sweep.R SITE=SOAP RUNGS=native
#   Rscript scripts/detect_treeiso_sweep.R SITES=SOAP,SJER,TEAK RUNGS=native,8,4,2,1 CORES=4
# Env: PYTHON=~/miniconda3/envs/treeiso/bin/python (see external/treeiso/README.md).
# Reads: work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv} + the cached
#   frozen clips work/neon/<SITE>/frozen/<SITE>/<plot>/<rung>/{clip_rawground.laz,
#   ground_dtm.tif}. Writes treeiso_results.csv + persists treeiso_instances/.
suppressMessages({ library(lidR); library(data.table); library(parallel) })
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- .job_dir()
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R")); source(.find("io_bridge.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (!is.null(A$SITES)) strsplit(A$SITES, ",")[[1]] else
  if (!is.null(A$SITE)) A$SITE else "SOAP"
RUNGS <- if (is.null(A$RUNGS)) "native" else strsplit(A$RUNGS, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 4 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL)) 4 else A$TOL)
MINTREES <- as.integer(if (is.null(A$MINTREES)) 6 else A$MINTREES)
PYTHON  <- if (!is.null(A$PYTHON)) A$PYTHON else "~/miniconda3/envs/treeiso/bin/python"
PYTHON  <- path.expand(PYTHON)
RUNNER  <- .find("../external/treeiso/run_treeiso.py")
if (is.null(RUNNER)) RUNNER <- Find(file.exists,
  c("external/treeiso/run_treeiso.py", file.path(getwd(), "external/treeiso/run_treeiso.py")))
ID_FIELD <- "treeiso"

fz <- function(nd, site, pid, rung, f) file.path(nd, "frozen", site, pid, rung, f)

run_plot <- function(site, pid, pc, gt, nd) {
  ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, , drop = FALSE]
  if (nrow(stems) < 1) return(NULL)
  idir <- file.path(nd, "treeiso_instances"); dir.create(idir, showWarnings = FALSE, recursive = TRUE)
  rows <- list()
  for (rung in RUNGS) {
    raw <- fz(nd, site, pid, rung, "clip_rawground.laz")
    dtm <- fz(nd, site, pid, rung, "ground_dtm.tif")
    if (!file.exists(raw) || !file.exists(dtm)) next
    out <- file.path(idir, sprintf("%s_%s.laz", pid, rung))
    if (!file.exists(out)) {                                   # checkpoint: reuse
      ok <- tryCatch(system2(PYTHON, c(shQuote(RUNNER), shQuote(raw), shQuote(out)),
                             stdout = FALSE, stderr = FALSE), error = function(e) 1L)
      if (!identical(as.integer(ok), 0L) || !file.exists(out)) next  # treeiso failed -> skip cell
    }
    det <- tryCatch(read_instances_laz(out, id_field = ID_FIELD), error = function(e) NULL)
    if (is.null(det) || !nrow(det)) next
    n_apex <- nrow(det)
    det <- tryCatch(det_to_agl(det, dtm), error = function(e) NULL)
    if (is.null(det) || !nrow(det)) next
    frdens <- tryCatch(jsonlite::read_json(fz(nd, site, pid, rung, "manifest.json"),
                       simplifyVector = TRUE)$frdens, error = function(e) NA_real_)
    sc <- tryCatch(score_plot(stems, det, tol_xy = TOL, core_cx = cx, core_cy = cy,
                              core_half = ph), error = function(e) NULL)
    if (is.null(sc)) next
    rows[[length(rows) + 1]] <- cbind(
      data.frame(site = site, plot = pid, plotType = ci$plotType, detector = "treeiso",
                 rung = rung, frdens = round(as.numeric(frdens), 2), n_apex = n_apex,
                 stringsAsFactors = FALSE), sc)
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows, fill = TRUE)
}

run_site <- function(site) {
  nd <- file.path(d, "neon", site)
  gtf <- file.path(nd, "ground_truth_stems.csv"); pcf <- file.path(nd, "plot_centroids.csv")
  if (!file.exists(gtf) || !file.exists(pcf)) { cat(sprintf("[%s] no GT\n", site)); return(NULL) }
  gt <- read.csv(gtf, stringsAsFactors = FALSE); pc <- read.csv(pcf, stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), , drop = FALSE]
  counts <- table(gt$plotID); keep <- intersect(names(counts)[counts >= MINTREES], pc$plotID)
  cat(sprintf("[%s] treeiso on %d plots x {%s}\n", site, length(keep), paste(RUNGS, collapse = ",")))
  if (!length(keep)) return(NULL)
  res <- rbindlist(Filter(Negate(is.null), mclapply(keep, function(p)
    tryCatch(run_plot(site, p, pc, gt, nd),
             error = function(e) { message("  ", p, ": ", conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)), fill = TRUE)
  if (!nrow(res)) { cat(sprintf("[%s] no rows\n", site)); return(NULL) }
  o <- file.path(nd, "treeiso_results.csv"); write.csv(res, o, row.names = FALSE)
  cat(sprintf("[%s] wrote %d rows -> %s\n", site, nrow(res), o))
  as.data.frame(res)
}

run_main <- function() {
  if (is.null(RUNNER) || !file.exists(PYTHON))
    cat(sprintf("NOTE: PYTHON=%s exists=%s ; RUNNER=%s\n", PYTHON, file.exists(PYTHON), RUNNER))
  t0 <- Sys.time(); all_res <- list()
  for (site in SITES) {
    r <- tryCatch(run_site(site), error = function(e) { message(site, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) all_res[[site]] <- r
  }
  res <- if (length(all_res)) do.call(rbind, all_res) else NULL
  dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\nDONE: %d rows across %d sites in %.1f min\n",
              if (is.null(res)) 0 else nrow(res), length(all_res), dt))
  if (!is.null(res) && nrow(res)) {
    cat("\npooled recall/precision/F1 by rung (treeiso):\n")
    for (rg in RUNGS) { s <- res[res$rung == rg, , drop = FALSE]; if (!nrow(s)) next
      p <- pool(s); cat(sprintf("  %-7s n_ref=%4d recall=%.3f prec=%.3f F1=%.3f rec_und=%.3f\n",
                                rg, p$n_ref, p$recall, p$precision, p$F1, p$rec_understory)) }
  }
}

if (sys.nframe() == 0L) run_main()
