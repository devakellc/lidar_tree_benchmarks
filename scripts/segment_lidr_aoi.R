#!/usr/bin/env Rscript
# Step 6 (crown segmentation) + Step 7 (crown metrics) in lidR on the real AOI.
# Mirrors detect_lidr_aoi.R for Steps 0-5, then dalponte2016 + crown_metrics.
suppressMessages({ library(lidR); library(sf) })
d <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")
f   <- file.path(d, "aoi.laz")
las <- readLAS(f, filter = "-drop_class 7 18 -drop_withheld")

## Step 0
fr   <- filter_poi(las, ReturnNumber == 1L)
g    <- pixel_metrics(fr, ~length(Z), res = 1)
v    <- terra::values(g); dens <- mean(v[!is.na(v) & v > 0])
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lidR seg AOI Step 0: density=%.2f -> res=%.2f m, wfloor=%.1f m\n",
            dens, res, wfloor))

## Steps 1-5
t0     <- Sys.time()
las    <- normalize_height(las, tin())
chm    <- rasterize_canopy(las, res = res,
            pitfree(thresholds = c(0, 10, 20, 30), max_edge = c(0, 1.5), subcircle = 0.2))
chm_lm <- if (dens < 8) terra::focal(chm, w = matrix(1/9, 3, 3), na.rm = TRUE) else chm
ttops  <- locate_trees(chm_lm, lmf(ws = ws, hmin = 2, shape = "circular"))

## Step 6 -- dalponte2016 region growing
max_cr_px <- as.integer(round(10 / res))   # 10 m diameter, lidR uses pixels
seg <- segment_trees(las, dalponte2016(chm_lm, ttops,
                                       th_seed = 0.45, th_cr = 0.55,
                                       max_cr = max_cr_px))
dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## Step 7
crowns <- crown_metrics(seg, func = .stdtreemetrics, geom = "convex")

sf::write_sf(crowns, file.path(d, "crowns_lidr_aoi.gpkg"), delete_dsn = TRUE)
sf::write_sf(ttops,  file.path(d, "tops_lidr_aoi_seg.gpkg"), delete_dsn = TRUE)

cat(sprintf(
  "lidR seg AOI: %d crowns / %d seeds in %.1f s; mean area %.1f m², Z %.1f-%.1f m\n",
  nrow(crowns), nrow(ttops), dt,
  mean(crowns$convhull_area), min(crowns$Z), max(crowns$Z)))
