#!/usr/bin/env Rscript
# Step 6 (crown segmentation) + Step 7 (crown metrics) in lasR on the toy tile.
# Mirrors detect_lasr.R for Steps 0/4/5, then adds region_growing (Dalponte
# 2016) and polygonizes the crown raster for per-tree metrics.
suppressMessages(library(lasR))
f <- system.file("extdata", "MixedConifer.las", package = "lasR")
d <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")

## Step 0 -- first-return density
ans0 <- exec(rasterize(1, "count", filter = keep_first()), on = f)
v    <- terra::values(ans0); dens <- mean(v[!is.na(v) & v > 0])
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lasR seg Step 0: density=%.2f -> res=%.2f m, wfloor=%.1f m\n",
            dens, res, wfloor))

## Steps 4-6 -- pit-free CHM, smooth (if QL2), VWF, region_growing
# NOTE: lasR's pit_fill is not the Khosravipour algorithm used by lidR::pitfree.
# max_cr is in DATA UNITS (meters) in lasR vs PIXELS in lidR.
t0   <- Sys.time()
del  <- triangulate(filter = keep_first())
chm  <- rasterize(res, del)
chm2 <- pit_fill(chm)
if (dens < 8) {
  smooth <- focal(chm2, size = 3, fun = "mean")
  seed   <- local_maximum_raster(smooth, ws, min_height = 2)
  crowns <- region_growing(smooth, seed,
                           th_tree = 2, th_seed = 0.45, th_cr = 0.55, max_cr = 10)
  ans    <- exec(del + chm + chm2 + smooth + seed + crowns, on = f)
} else {
  seed   <- local_maximum_raster(chm2, ws, min_height = 2)
  crowns <- region_growing(chm2, seed,
                           th_tree = 2, th_seed = 0.45, th_cr = 0.55, max_cr = 10)
  ans    <- exec(del + chm + chm2 + seed + crowns, on = f)
}
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

## Step 7 -- polygonize + per-tree metrics
# lasR returns SpatRaster proxies without values loaded; reload from file path.
crowns_r <- terra::rast(terra::sources(ans$region_growing))
seeds_sf <- sf::st_as_sf(ans$local_maximum)

# CHM used for the LM (post-smooth if applied) for apex height per crown
chm_src  <- if (dens < 8) ans$focal else ans$pit_fill
chm_lm   <- terra::rast(terra::sources(chm_src))

# Polygons: one per unique treeID label (background value 0 dropped).
polys <- terra::as.polygons(crowns_r, dissolve = TRUE, values = TRUE)
names(polys)[1] <- "treeID"
polys <- polys[polys$treeID > 0, ]
polys$area_m2 <- terra::expanse(polys, unit = "m")

# Apex height per crown via zonal max on the LM CHM (raster aligns 1:1)
zmax <- terra::zonal(chm_lm, crowns_r, fun = "max", na.rm = TRUE)
names(zmax) <- c("treeID", "z_max")
polys$z_max <- zmax$z_max[match(polys$treeID, zmax$treeID)]

# Save crown raster, polygons, and a treetop CSV consistent with detect script
terra::writeRaster(crowns_r, file.path(d, "crowns_lasr.tif"), overwrite = TRUE)
sf::write_sf(sf::st_as_sf(polys), file.path(d, "crowns_lasr.gpkg"), delete_dsn = TRUE)
cc <- sf::st_coordinates(seeds_sf)
write.csv(data.frame(x = cc[, 1], y = cc[, 2], z = cc[, 3]),
          file.path(d, "tops_lasr_seg.csv"), row.names = FALSE)

cat(sprintf("lasR seg: %d crowns / %d seeds in %.2f s; mean area %.1f m², Z %.1f-%.1f m\n",
            nrow(polys), nrow(seeds_sf), dt,
            mean(polys$area_m2), min(polys$z_max, na.rm = TRUE), max(polys$z_max, na.rm = TRUE)))
