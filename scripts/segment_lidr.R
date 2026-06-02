#!/usr/bin/env Rscript
# Step 6 (crown segmentation) + Step 7 (crown metrics) in lidR on the toy tile.
# Mirrors detect_lidr.R for Steps 0/4/5, then runs dalponte2016 and writes
# both per-tree treetops and convex-hull crown polygons (.stdtreemetrics).
suppressMessages({ library(lidR); library(sf) })
f   <- system.file("extdata", "MixedConifer.las", package = "lasR")
d   <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")
las <- readLAS(f, select = "xyzr")

## Step 0
fr   <- filter_poi(las, ReturnNumber == 1L)
g    <- pixel_metrics(fr, ~length(Z), res = 1)
v    <- terra::values(g); dens <- mean(v[!is.na(v) & v > 0])
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lidR seg Step 0: density=%.2f -> res=%.2f m, wfloor=%.1f m\n",
            dens, res, wfloor))

## Steps 4-5 -- pit-free CHM, smooth (if QL2), VWF
t0    <- Sys.time()
chm   <- rasterize_canopy(las, res = res,
           pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5), subcircle = 0.2))
chm_lm <- if (dens < 8) terra::focal(chm, w = matrix(1/9, 3, 3), na.rm = TRUE) else chm
ttops <- locate_trees(chm_lm, lmf(ws = ws, hmin = 2, shape = "circular"))

## Step 6 -- dalponte2016 region growing on the same CHM the LM used
# max_cr in lidR is PIXELS, so 10 m / res for a comparable max diameter to
# the lasR script (which uses max_cr in meters).
max_cr_px <- as.integer(round(10 / res))
seg <- segment_trees(las, dalponte2016(chm_lm, ttops,
                                       th_seed = 0.45, th_cr = 0.55,
                                       max_cr = max_cr_px))
dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## Step 7 -- crown metrics with convex hulls
crowns <- crown_metrics(seg, func = .stdtreemetrics, geom = "convex")

# Save
sf::write_sf(crowns,  file.path(d, "crowns_lidr.gpkg"),  delete_dsn = TRUE)
sf::write_sf(ttops,   file.path(d, "tops_lidr_seg.gpkg"), delete_dsn = TRUE)

cat(sprintf(
  "lidR seg: %d crowns / %d seeds in %.2f s; mean area %.1f m², Z %.1f-%.1f m\n",
  nrow(crowns), nrow(ttops), dt,
  mean(crowns$convhull_area), min(crowns$Z), max(crowns$Z)))
