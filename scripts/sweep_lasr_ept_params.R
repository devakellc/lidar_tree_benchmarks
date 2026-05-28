#!/usr/bin/env Rscript
# Sweep lasR EPT acquisition tuning knobs on a fixed remote AOI.
# Usage:
#   export CLAUDE_JOB_DIR=/tmp/work
#   Rscript scripts/sweep_lasr_ept_params.R

suppressMessages(library(lasR))

url <- Sys.getenv(
  "LASR_EPT_URL",
  unset = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json"
)
bb <- if (nzchar(Sys.getenv("LASR_EPT_AOI_BBOX", unset = ""))) {
  as.numeric(strsplit(Sys.getenv("LASR_EPT_AOI_BBOX"), ",", fixed = TRUE)[[1]])
} else {
  c(-13579688, 4987253, -13576556, 4990205)  # ~4x prior 5x test AOI
}
reps <- as.integer(Sys.getenv("LASR_SWEEP_REPS", unset = "2"))
if (is.na(reps) || reps < 1) reps <- 1L

q <- reader_rectangles(
  bb[1], bb[2], bb[3], bb[4],
  filter = "-drop_class 7 18 -drop_withheld"
) + summarise()

# Each row: label, workers, partitions (NA = default), prefetch (NA = default),
# vsi_cache_mb (NA = default / unset)
grid <- data.frame(
  label = c(
    "seq",
    "par8_default",
    "par16_default",
    "par16_parts32",
    "par16_parts36",
    "par16_parts64",
    "par8_prefetch8",
    "par8_prefetch16",
    "par8_prefetch32",
    "par16_parts32_prefetch8",
    "par16_parts32_prefetch16",
    "par16_parts32_prefetch32",
    "par8_cache256",
    "par8_cache1024",
    "par16_parts32_cache256",
    "par16_parts32_cache1024"
  ),
  workers = c(1, 8, 16, 16, 16, 16, 8, 8, 8, 16, 16, 16, 8, 8, 16, 16),
  partitions = c(NA, NA, NA, 32, 36, 64, NA, NA, NA, 32, 32, 32, NA, NA, 32, 32),
  prefetch = c(NA, NA, NA, NA, NA, NA, 8, 16, 32, 8, 16, 32, NA, NA, NA, NA),
  cache_mb = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 256, 1024, 256, 1024),
  stringsAsFactors = FALSE
)

restore_env <- function(prev) {
  for (nm in names(prev)) {
    if (is.na(prev[[nm]])) Sys.unsetenv(nm) else Sys.setenv(setNames(prev[[nm]], nm))
  }
}

capture_env <- function(vars) {
  out <- list()
  for (v in vars) {
    raw <- Sys.getenv(v, unset = NA_character_)
    out[[v]] <- if (is.na(raw) || !nzchar(raw)) NA_character_ else raw
  }
  out
}

read_exec_cfg <- function(logf) {
  if (!file.exists(logf)) return(list(chunks = NA_integer_, files = NA_integer_))
  lines <- readLines(logf, warn = FALSE)
  chunks <- sub(".*Chunks: ([0-9]+).*", "\\1",
                grep("Chunks:", lines, value = TRUE)[1])
  files <- sub(".*Concurrent files: ([0-9]+).*", "\\1",
               grep("Concurrent files:", lines, value = TRUE)[1])
  list(
    chunks = as.integer(chunks),
    files = as.integer(files)
  )
}

run_case <- function(label, workers, partitions, prefetch, cache_mb) {
  env_vars <- c("LASR_EPT_PARTITIONS", "LASR_EPT_PREFETCH", "VSI_CACHE_SIZE")
  prev <- capture_env(env_vars)

  if (is.na(partitions)) Sys.unsetenv("LASR_EPT_PARTITIONS") else
    Sys.setenv(LASR_EPT_PARTITIONS = as.character(partitions))
  if (is.na(prefetch)) Sys.unsetenv("LASR_EPT_PREFETCH") else
    Sys.setenv(LASR_EPT_PREFETCH = as.character(prefetch))
  if (is.na(cache_mb)) Sys.unsetenv("VSI_CACHE_SIZE") else
    Sys.setenv(VSI_CACHE_SIZE = as.character(as.integer(cache_mb) * 1024L * 1024L))

  on.exit(restore_env(prev), add = TRUE)

  ncores <- if (workers <= 1L) sequential() else concurrent_files(workers)
  times <- numeric(reps)
  npoints <- NA_real_
  chunks <- NA_integer_
  cfiles <- NA_integer_

  for (i in seq_len(reps)) {
    gc(full = TRUE)
    logf <- tempfile(fileext = ".log")
    t0 <- Sys.time()
    ans <- exec(q, on = url, ncores = ncores, progress = FALSE,
                verbose = (i == 1L), log_file = if (i == 1L) logf else "")
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    npoints <- ans$npoints
    if (i == 1L) {
      cfg <- read_exec_cfg(logf)
      chunks <- cfg$chunks
      cfiles <- cfg$files
    }
  }

  data.frame(
    label = label,
    workers = workers,
    partitions = ifelse(is.na(partitions), "default", partitions),
    prefetch = ifelse(is.na(prefetch), "default", prefetch),
    cache_mb = ifelse(is.na(cache_mb), "default", cache_mb),
    reps = reps,
    median_s = median(times),
    min_s = min(times),
    max_s = max(times),
    chunks = chunks,
    concurrent_files = cfiles,
    npoints = as.integer(npoints),
    mpts_s = npoints / median(times) / 1e6,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("lasR %s | available_threads=%d\n",
            as.character(packageVersion("lasR")), lasR:::available_threads()))
cat(sprintf("AOI bbox: %s\n", paste(bb, collapse = ",")))
cat(sprintf("reps per case: %d\n\n", reps))

results <- lapply(seq_len(nrow(grid)), function(i) {
  row <- grid[i, ]
  cat(sprintf("[%d/%d] %s ...\n", i, nrow(grid), row$label)); flush.console()
  run_case(row$label, row$workers, row$partitions, row$prefetch, row$cache_mb)
})
out <- do.call(rbind, results)
rownames(out) <- NULL

out_dir <- Sys.getenv("CLAUDE_JOB_DIR", unset = tempdir())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv <- file.path(out_dir, "sweep_lasr_ept_params.csv")
write.csv(out, csv, row.names = FALSE)

cat("\n=== results (sorted by median_s) ===\n")
print(out[order(out$median_s), ], row.names = FALSE)
cat("\nWrote:", csv, "\n")
