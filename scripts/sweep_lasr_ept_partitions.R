#!/usr/bin/env Rscript
# Partition-count sweep for lasR EPT acquisition on a LARGE AOI.
#
# Motivation: the original sweep (112.5M pts, see ept-acquisition-sweep-results.md)
# found a fixed target of 32 partitions (~36 chunks) optimal at
# concurrent_files(16), beating 64 (~100 chunks) by ~40%. The optimum was a
# roughly constant *chunk count* across the 31M and 112M AOIs -- not constant
# points-per-chunk. This script tests whether that still holds ~8x larger
# (10k-acre AOI, ~0.8-0.9B pts): does a much bigger AOI finally want more
# chunks, or does the fixed-32 default remain best?
#
# Fixes the worker count, varies LASR_EPT_PARTITIONS. Prefetch / VSI cache are
# left at defaults (the original sweep showed <=1 s sensitivity to both).
#
# Usage:
#   export CLAUDE_JOB_DIR=/tmp/lidar_bench_10k
#   export LASR_EPT_AOI_BBOX="-13549926.53,4944632.44,-13541563.99,4953027.51"
#   Rscript scripts/sweep_lasr_ept_partitions.R

suppressMessages(library(lasR))

url <- Sys.getenv("LASR_EPT_URL",
  unset = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json")

raw_bbox <- Sys.getenv("LASR_EPT_AOI_BBOX", unset = "")
if (!nzchar(raw_bbox)) stop("set LASR_EPT_AOI_BBOX='xmin,ymin,xmax,ymax' (EPSG:3857)")
bb <- as.numeric(strsplit(raw_bbox, ",", fixed = TRUE)[[1]])
stopifnot(length(bb) == 4L, all(is.finite(bb)), bb[1] < bb[3], bb[2] < bb[4])

reps     <- max(1L, as.integer(Sys.getenv("LASR_SWEEP_REPS",     unset = "2")))
seq_reps <- max(1L, as.integer(Sys.getenv("LASR_SWEEP_SEQ_REPS", unset = "1")))
workers  <- max(2L, as.integer(Sys.getenv("LASR_SWEEP_WORKERS",  unset = "16")))
parts    <- as.integer(strsplit(Sys.getenv("LASR_SWEEP_PARTS", unset = "32,64,96,128,192"), ",")[[1]])
do_seq   <- tolower(Sys.getenv("LASR_SWEEP_SEQ",  unset = "true")) %in% c("1","true","yes")
do_par8  <- tolower(Sys.getenv("LASR_SWEEP_PAR8", unset = "true")) %in% c("1","true","yes")

out_dir <- Sys.getenv("CLAUDE_JOB_DIR", unset = tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(out_dir, "sweep_lasr_ept_partitions.csv")

q <- reader_rectangles(bb[1], bb[2], bb[3], bb[4],
                       filter = "-drop_class 7 18 -drop_withheld") + summarise()

read_chunks <- function(logf) {
  if (!file.exists(logf)) return(NA_integer_)
  lines <- readLines(logf, warn = FALSE)
  hit <- grep("Chunks:", lines, value = TRUE)
  if (!length(hit)) return(NA_integer_)
  as.integer(sub(".*Chunks: ([0-9]+).*", "\\1", hit[1]))
}

run_case <- function(label, ncores, wk, partitions, n_reps) {
  prev <- Sys.getenv("LASR_EPT_PARTITIONS", unset = NA_character_)
  on.exit({
    if (is.na(prev)) Sys.unsetenv("LASR_EPT_PARTITIONS")
    else Sys.setenv(LASR_EPT_PARTITIONS = prev)
  }, add = TRUE)
  if (is.na(partitions)) Sys.unsetenv("LASR_EPT_PARTITIONS")
  else Sys.setenv(LASR_EPT_PARTITIONS = as.character(partitions))

  times <- rep(NA_real_, n_reps); npoints <- NA_real_; chunks <- NA_integer_; err <- ""
  for (i in seq_len(n_reps)) {
    gc(full = TRUE)
    logf <- tempfile(fileext = ".log")
    t0 <- Sys.time()
    ans <- tryCatch(
      exec(q, on = url, ncores = ncores, progress = FALSE,
           verbose = (i == 1L), log_file = if (i == 1L) logf else ""),
      error = function(e) { err <<- conditionMessage(e); NULL })
    if (is.null(ans)) { cat(sprintf("    rep %d/%d ERROR: %s\n", i, n_reps, err)); flush.console(); next }
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    npoints  <- ans$npoints
    if (i == 1L) chunks <- read_chunks(logf)
    cat(sprintf("    rep %d/%d: %6.1f s  (npoints=%s, chunks=%s)  [%s]\n",
                i, n_reps, times[i], format(npoints, big.mark = ","),
                chunks, format(Sys.time(), "%H:%M:%S"))); flush.console()
  }
  data.frame(
    label = label, workers = wk,
    partitions = ifelse(is.na(partitions), "default", as.character(partitions)),
    reps = n_reps,
    median_s = suppressWarnings(median(times, na.rm = TRUE)),
    min_s = suppressWarnings(min(times, na.rm = TRUE)),
    max_s = suppressWarnings(max(times, na.rm = TRUE)),
    chunks = chunks, npoints = as.numeric(npoints),
    mpts_s = npoints / suppressWarnings(median(times, na.rm = TRUE)) / 1e6,
    error = err, stringsAsFactors = FALSE)
}

# Build case list: partition sweep at `workers` first (the question), then anchors.
cases <- lapply(parts, function(p) list(
  label = sprintf("par%d_parts%d", workers, p),
  ncores = concurrent_files(workers), wk = workers, parts = p, reps = reps))
if (do_par8) cases <- c(cases, list(list(
  label = "par8_parts32", ncores = concurrent_files(8), wk = 8L, parts = 32L, reps = reps)))
if (do_seq) cases <- c(cases, list(list(
  label = "sequential", ncores = sequential(), wk = 1L, parts = NA_integer_, reps = seq_reps)))

cat(sprintf("lasR %s | available_threads=%d | host_cpus=%d\n",
            as.character(packageVersion("lasR")), lasR:::available_threads(),
            parallel::detectCores()))
cat(sprintf("AOI bbox (3857): %s\n", paste(sprintf("%.1f", bb), collapse = ",")))
cat(sprintf("3857 extent: %.0f x %.0f m\n", bb[3]-bb[1], bb[4]-bb[2]))
cat(sprintf("cases: %d | reps: %d (seq %d) | started %s\n\n",
            length(cases), reps, seq_reps, format(Sys.time(), "%H:%M:%S")))

# Warm-up: prime DNS/TLS + EPT metadata so the first timed run isn't penalized.
cat("[warm-up] priming connection + EPT metadata ...\n"); flush.console()
invisible(tryCatch(exec(reader(depth = 0) + summarise(), on = url, ncores = sequential()),
                   error = function(e) cat("  warm-up note:", conditionMessage(e), "\n")))

t_start <- Sys.time()
all_results <- list()
for (k in seq_along(cases)) {
  cs <- cases[[k]]
  cat(sprintf("[%d/%d] %s (workers=%d, partitions=%s, reps=%d) ...\n",
              k, length(cases), cs$label, cs$wk,
              ifelse(is.na(cs$parts), "default", cs$parts), cs$reps)); flush.console()
  all_results[[k]] <- run_case(cs$label, cs$ncores, cs$wk, cs$parts, cs$reps)
  out <- do.call(rbind, all_results)
  write.csv(out, csv, row.names = FALSE)   # incremental: survive interruption
  cat(sprintf("    -> wrote %s (%d/%d cases)\n\n", csv, k, length(cases))); flush.console()
}

out <- do.call(rbind, all_results)
rownames(out) <- NULL
# npoints consistency check (all configs read the same AOI)
np <- unique(out$npoints[is.finite(out$npoints)])
cat(sprintf("=== DONE in %.1f min | distinct npoints: %s ===\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins")),
            paste(format(np, big.mark = ","), collapse = " | ")))
cat("\n=== results (sorted by median_s) ===\n")
print(out[order(out$median_s), c("label","workers","partitions","chunks",
                                 "median_s","min_s","max_s","mpts_s","error")],
      row.names = FALSE)
cat("\nWrote:", csv, "\n")
