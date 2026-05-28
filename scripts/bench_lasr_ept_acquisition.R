#!/usr/bin/env Rscript
# Benchmark: lasR EPT AOI acquisition sequential vs parallel-files.
# Times reader_rectangles(AOI) + summarise() on a remote ept.json.

suppressMessages(library(lasR))

url <- Sys.getenv(
  "LASR_EPT_URL",
  unset = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json"
)

default_bbox <- c(-13578452, 4988399, -13577792, 4989059)
raw_bbox <- Sys.getenv("LASR_EPT_AOI_BBOX", unset = "")
if (nzchar(raw_bbox)) {
  bb <- suppressWarnings(as.numeric(strsplit(raw_bbox, ",", fixed = TRUE)[[1]]))
  if (length(bb) != 4L || any(is.na(bb)) || bb[1] >= bb[3] || bb[2] >= bb[4]) {
    stop("LASR_EPT_AOI_BBOX must be 'xmin,ymin,xmax,ymax' with xmin<xmax and ymin<ymax")
  }
  bbox <- bb
} else {
  bbox <- default_bbox
}
xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

depth <- as.integer(Sys.getenv("LASR_EPT_DEPTH", unset = "-1"))
reps <- as.integer(Sys.getenv("LASR_BENCH_REPS", unset = "3"))
if (is.na(reps) || reps < 1) reps <- 1L

stopifnot(has_omp_support())

q <- reader_rectangles(
  xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
  depth = depth,
  filter = "-drop_class 7 18 -drop_withheld"
) + summarise()

avail <- lasR:::available_threads()
cat(sprintf("lasR version:          %s\n", as.character(packageVersion("lasR"))))
cat(sprintf("lasR.so built:         %s\n", format(file.info(system.file("libs/lasR.so", package = "lasR"))$mtime)))
cat(sprintf("Endpoint:              %s\n", url))
cat(sprintf("AOI (EPSG:3857):       [%.0f, %.0f, %.0f, %.0f]\n", xmin, ymin, xmax, ymax))
cat(sprintf("Depth:                 %d\n", depth))
cat(sprintf("Host CPUs (detect):    %d\n", parallel::detectCores()))
cat(sprintf("lasR available_threads: %d  (requests above this are capped)\n", avail))
cat(sprintf("Repeats per config:    %d\n", reps))

# Warm-up: prime DNS/TLS and EPT metadata cache.
cat("\n[warm-up]\n")
invisible(exec(reader(depth = 0) + summarise(), on = url, ncores = sequential()))

read_exec_config <- function(log_file) {
  lines <- readLines(log_file, warn = FALSE)
  grep("Concurrent files:|Concurrent points:|Chunks:", lines, value = TRUE)
}

time_run <- function(label, ncores) {
  times <- numeric(reps)
  npoints <- NA_real_
  exec_cfg <- character()
  for (i in seq_len(reps)) {
    gc(full = TRUE)
    cat(sprintf("  %-18s rep=%d/%d ...\n", label, i, reps)); flush.console()
    logf <- tempfile(fileext = ".log")
    t0 <- Sys.time()
    ans <- exec(q, on = url, ncores = ncores, progress = FALSE,
                verbose = (i == 1L), log_file = if (i == 1L) logf else "")
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    npoints <- ans$npoints
    if (i == 1L && file.exists(logf)) exec_cfg <- read_exec_config(logf)
    cat(sprintf("  %-18s rep=%d wall=%.2f s npoints=%d\n", label, i, times[i], ans$npoints))
    flush.console()
  }
  list(label = label, npoints = npoints, times = times, exec_cfg = exec_cfg)
}

cores_files <- as.integer(Sys.getenv("LASR_BENCH_FILE_CORES", unset = as.character(half_cores())))
if (is.na(cores_files) || cores_files < 2) cores_files <- 2L
if (cores_files > avail) {
  cat(sprintf("\nNOTE: LASR_BENCH_FILE_CORES=%d exceeds available_threads=%d; lasR will cap it.\n",
              cores_files, avail))
}

cat("\n[runs]\n")
res_seq <- time_run("sequential", sequential())
res_par <- time_run(sprintf("concurrent_files(%d)", cores_files), concurrent_files(cores_files))

if (!is.na(res_seq$npoints) && !is.na(res_par$npoints)) stopifnot(res_seq$npoints == res_par$npoints)

seq_med <- stats::median(res_seq$times)
par_med <- stats::median(res_par$times)
speedup <- seq_med / par_med

cat("\n[summary]\n")
cat(sprintf("  sequential:            median=%.2f s  (n=%d)\n", seq_med, reps))
cat(sprintf("  concurrent_files(%d):  median=%.2f s  (n=%d)\n", cores_files, par_med, reps))
cat(sprintf("  speedup:               %.2fx\n", speedup))
cat(sprintf("  npoints:               %d\n", as.integer(res_seq$npoints)))
if (length(res_seq$exec_cfg)) {
  cat("\n[exec config: sequential]\n")
  cat(paste0("  ", res_seq$exec_cfg, collapse = "\n"), "\n")
}
if (length(res_par$exec_cfg)) {
  cat("\n[exec config: parallel-files]\n")
  cat(paste0("  ", res_par$exec_cfg, collapse = "\n"), "\n")
  cat("  (EPT auto-partition target defaults to 4 x concurrent_files workers)\n")
}

