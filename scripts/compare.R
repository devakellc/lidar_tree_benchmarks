#!/usr/bin/env Rscript
# Compare lasR vs lidR treetop detections.
d  <- Sys.getenv("CLAUDE_JOB_DIR")
ar <- commandArgs(trailingOnly = TRUE)
fa <- if (length(ar) >= 1) ar[1] else file.path(d, "tops_lasr.csv")
fb <- if (length(ar) >= 2) ar[2] else file.path(d, "tops_lidr.csv")
a  <- read.csv(fa)   # lasR
b  <- read.csv(fb)   # lidR

# Greedy nearest-neighbour matching within a tolerance (crown-scale).
match_tops <- function(a, b, tol = 2.0) {
  used <- rep(FALSE, nrow(b)); matched <- 0; dz <- numeric(0)
  for (i in seq_len(nrow(a))) {
    dx <- b$x - a$x[i]; dy <- b$y - a$y[i]
    dd <- sqrt(dx * dx + dy * dy); dd[used] <- Inf
    j  <- which.min(dd)
    if (length(j) && is.finite(dd[j]) && dd[j] <= tol) {
      used[j] <- TRUE; matched <- matched + 1; dz <- c(dz, b$z[j] - a$z[i])
    }
  }
  list(matched = matched, dz = dz)
}

for (tol in c(1.0, 2.0)) {
  m <- match_tops(a, b, tol)
  cat(sprintf(
    "tol=%.1f m | matched %d of (lasR %d, lidR %d) | Jaccard=%.2f | dZ mean=%+.2f sd=%.2f m\n",
    tol, m$matched, nrow(a), nrow(b),
    m$matched / (nrow(a) + nrow(b) - m$matched),
    if (length(m$dz)) mean(m$dz) else NA,
    if (length(m$dz)) sd(m$dz) else NA))
}
cat(sprintf("\nCount: lasR=%d  lidR=%d  diff=%+d (%.1f%%)\n",
            nrow(a), nrow(b), nrow(a) - nrow(b),
            100 * (nrow(a) - nrow(b)) / nrow(b)))
cat(sprintf("Height: lasR mean=%.2f  lidR mean=%.2f m\n", mean(a$z), mean(b$z)))
