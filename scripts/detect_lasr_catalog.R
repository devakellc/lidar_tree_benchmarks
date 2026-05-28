#!/usr/bin/env Rscript
# lasR streaming demo on the retiled AOI tiles from tile_aoi.R.
# lasR automatically buffers between tiles for stages that need it
# (triangulate, normalize, region_growing, etc.) — see approach §3.
suppressMessages(library(lasR))
d <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")
tdir <- file.path(d, "tiles")
if (!dir.exists(tdir) || length(list.files(tdir, "\\.laz$")) == 0L)
  stop("No tiles found; run scripts/tile_aoi.R first")

read <- reader_las(filter = "-drop_class 7 18 -drop_withheld")

## Step 0 -- density via summarise() on the whole tile set
stats <- exec(read + summarise(), on = tdir)
dens <- if (!is.null(stats$pulse_density) && stats$pulse_density > 0) stats$pulse_density else stats$density
if (is.null(dens) || is.na(dens) || dens <= 0)
  stop("Could not compute density")
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")

res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lasR catalog Step 0: density=%.2f -> res=%.2f m, wfloor=%.1f m\n",
            dens, res, wfloor))

## Steps 1-5 -- normalize, pit-free CHM, smooth (if QL2), VWF
# Each stage uses automatic inter-tile buffers as needed (approach §3).
t0   <- Sys.time()
norm <- normalize()
del  <- triangulate(filter = keep_first())
chm  <- rasterize(res, del)
chm2 <- pit_fill(chm)
if (dens < 8) {
  smooth <- focal(chm2, size = 3, fun = "mean")
  seed   <- local_maximum_raster(smooth, ws, min_height = 2)
  ans    <- exec(read + norm + del + chm + chm2 + smooth + seed, on = tdir)
} else {
  seed <- local_maximum_raster(chm2, ws, min_height = 2)
  ans  <- exec(read + norm + del + chm + chm2 + seed, on = tdir)
}
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

tops <- ans$local_maximum
sf::write_sf(sf::st_as_sf(tops), file.path(d, "tops_lasr_catalog.gpkg"), delete_dsn = TRUE)

cc <- sf::st_coordinates(tops)
cat(sprintf("lasR catalog: %d treetops in %.1f s; Z range %.1f-%.1f m\n",
            nrow(cc), dt, min(cc[, 3]), max(cc[, 3])))
