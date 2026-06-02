#!/usr/bin/env Rscript
# Full approach on real USGS 3DEP (raw, multi-return) AOI -- lidR.
# drop noise -> Step 0 density -> normalize -> pit-free CHM -> variable-window LM.
suppressMessages(library(lidR))
f   <- file.path(Sys.getenv("CLAUDE_JOB_DIR"), "aoi.laz")
las <- readLAS(f, filter = "-drop_class 7 18 -drop_withheld")

## Step 0 -- first-return density (1 m cells), same rule as lasR
fr   <- filter_poi(las, ReturnNumber == 1L)
g    <- pixel_metrics(fr, ~length(Z), res = 1)
v    <- terra::values(g); dens <- mean(v[!is.na(v) & v > 0])
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lidR Step 0: density=%.2f -> res=%.2f m, window floor=%.1f m\n",
            dens, res, wfloor))

## normalize -> pit-free CHM -> VWF
# Step 5 smoothing branch (approach §2 Step 5): QL2 (dens < 8) gets a 3x3 mean
# pre-LM smooth to suppress noise peaks; QL1/0 skips it.
t0    <- Sys.time()
las   <- normalize_height(las, tin())
chm   <- rasterize_canopy(las, res = res,
           pitfree(thresholds = c(0, 10, 20, 30), max_edge = c(0, 1.5),
                   subcircle = 0.2))
chm_lm <- if (dens < 8) terra::focal(chm, w = matrix(1/9, 3, 3), na.rm = TRUE) else chm
ttops <- locate_trees(chm_lm, lmf(ws = ws, hmin = 2, shape = "circular"))
dt    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

terra::writeRaster(chm, file.path(Sys.getenv("CLAUDE_JOB_DIR"), "chm_aoi.tif"),
                   overwrite = TRUE)
xy  <- sf::st_coordinates(ttops)
out <- data.frame(x = xy[, 1], y = xy[, 2], z = ttops$Z)
write.csv(out, file.path(Sys.getenv("CLAUDE_JOB_DIR"), "tops_lidr_aoi.csv"),
          row.names = FALSE)
cat(sprintf("lidR: %d treetops in %.1f s; height range %.1f-%.1f m\n",
            nrow(out), dt, min(out$z), max(out$z)))