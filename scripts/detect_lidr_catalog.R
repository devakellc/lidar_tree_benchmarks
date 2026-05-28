#!/usr/bin/env Rscript
# LAScatalog streaming demo on the retiled AOI tiles produced by tile_aoi.R.
# Exercises the approach doc's §3 lidR guidance: opt_chunk_buffer >= largest
# crown radius (10-30 m) so trees at seams aren't missed/duplicated.
suppressMessages({ library(lidR); library(sf) })
d <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")

tdir <- file.path(d, "tiles")
if (!dir.exists(tdir) || length(list.files(tdir, "\\.laz$")) == 0L)
  stop("No tiles found; run scripts/tile_aoi.R first")

ctg <- readLAScatalog(tdir)
opt_filter(ctg)       <- "-drop_class 7 18 -drop_withheld"
opt_select(ctg)       <- "xyzrn"
opt_chunk_buffer(ctg) <- 20    # >= largest crown radius (approach §3)

cat(sprintf("Catalog: %d tile files, chunk_buffer=%s, chunk_size=%s\n",
            nrow(ctg@data), opt_chunk_buffer(ctg), opt_chunk_size(ctg)))

## Step 0 -- density on the catalog as a whole
dens <- lidR::density(ctg)
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lidR catalog Step 0: density=%.2f -> res=%.2f m, wfloor=%.1f m\n",
            dens, res, wfloor))

## Steps 3-5 -- normalize, pit-free CHM, smooth (if QL2), VWF
# normalize_height on a catalog needs an output template; write to tiles_norm/.
ndir <- file.path(d, "tiles_norm")
dir.create(ndir, recursive = TRUE, showWarnings = FALSE)
opt_output_files(ctg) <- file.path(ndir, "tile_{XLEFT}_{YBOTTOM}_norm")

t0 <- Sys.time()
ctg_norm <- normalize_height(ctg, tin())
# Subsequent stages return in-memory; clear the output template inherited from
# the parent catalog so locate_trees returns an sf object instead of file paths.
opt_output_files(ctg_norm) <- ""
opt_chunk_buffer(ctg_norm) <- 20
chm  <- rasterize_canopy(ctg_norm, res = res,
          pitfree(thresholds = c(0, 10, 20, 30), max_edge = c(0, 1.5), subcircle = 0.2))
chm_lm <- if (dens < 8) terra::focal(chm, w = matrix(1/9, 3, 3), na.rm = TRUE) else chm
ttops  <- locate_trees(ctg_norm, lmf(ws = ws, hmin = 2, shape = "circular"))
dt     <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

terra::writeRaster(chm, file.path(d, "chm_catalog.tif"), overwrite = TRUE)
sf::write_sf(ttops, file.path(d, "tops_lidr_catalog.gpkg"), delete_dsn = TRUE)

cat(sprintf("lidR catalog: %d treetops in %.1f s; Z range %.1f-%.1f m\n",
            nrow(ttops), dt, min(ttops$Z), max(ttops$Z)))
