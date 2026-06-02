#!/usr/bin/env Rscript
# Controlled test: run BOTH local-maximum implementations on ONE shared CHM
# (the lidR pit-free CHM written to chm_shared.tif). Isolates the
# detection algorithm from CHM-construction differences.
#
# Window function uses the SAME density-derived wfloor as the detect scripts,
# so "same ws" between lasR and lidR is literally true.
d   <- Sys.getenv("CLAUDE_JOB_DIR")
tif <- file.path(d, "chm_shared.tif")
las_f <- system.file("extdata", "MixedConifer.las", package = "lasR")

suppressMessages(library(lasR))
ans0 <- exec(rasterize(1, "count", filter = keep_first()), on = las_f)
v    <- terra::values(ans0); dens <- mean(v[!is.na(v) & v > 0])
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("shared_chm: density=%.2f -> wfloor=%.1f m\n", dens, wfloor))

## lasR LM on the shared CHM
r    <- load_raster(tif)
seed <- local_maximum_raster(r, ws, min_height = 2)
ans  <- exec(r + seed, on = las_f)  # LAS defines the region; CHM is loaded
tops <- ans[[length(ans)]]
cc   <- sf::st_coordinates(tops)
A    <- data.frame(x = cc[, 1], y = cc[, 2], z = cc[, 3])

## lidR LM on the same CHM
suppressMessages(library(lidR))
chm  <- terra::rast(tif)
tt   <- locate_trees(chm, lmf(ws = ws, hmin = 2, shape = "circular"))
xy   <- sf::st_coordinates(tt)
B    <- data.frame(x = xy[, 1], y = xy[, 2], z = tt$Z)

## match
match_n <- function(a, b, tol) {
  used <- rep(FALSE, nrow(b)); m <- 0
  for (i in seq_len(nrow(a))) {
    dd <- sqrt((b$x - a$x[i])^2 + (b$y - a$y[i])^2); dd[used] <- Inf
    j <- which.min(dd)
    if (length(j) && is.finite(dd[j]) && dd[j] <= tol) { used[j] <- TRUE; m <- m + 1 }
  }
  m
}
cat(sprintf("SAME CHM | lasR LM=%d  lidR LM=%d\n", nrow(A), nrow(B)))
for (tol in c(1.0, 2.0)) {
  m <- match_n(A, B, tol)
  cat(sprintf("  tol=%.1f m matched=%d Jaccard=%.2f\n",
              tol, m, m / (nrow(A) + nrow(B) - m)))
}
