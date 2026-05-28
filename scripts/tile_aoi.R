#!/usr/bin/env Rscript
# Retile the 25 ha AOI laz into a small grid of tiles under
# $CLAUDE_JOB_DIR/tiles/. catalog_retile() aligns tiles to its own grid origin,
# so a 505 m AOI with a target ~253 m chunk size produces a 3x3 grid (9 tiles).
# This gives the catalog/streaming scripts inter-tile seams to demonstrate
# buffering on without needing a true wall-to-wall dataset.
suppressMessages(library(lidR))
d  <- Sys.getenv("CLAUDE_JOB_DIR")
if (!nzchar(d)) stop("Set CLAUDE_JOB_DIR to a writable directory")
f  <- file.path(d, "aoi.laz")
stopifnot(file.exists(f))

tdir <- file.path(d, "tiles")
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

ctg <- readLAScatalog(f)
e   <- as.numeric(sf::st_bbox(ctg))   # xmin, ymin, xmax, ymax
side <- max(e[3] - e[1], e[4] - e[2]) / 2 + 1   # 2x2 -> ~half-extent

opt_chunk_size(ctg) <- ceiling(side)
opt_chunk_buffer(ctg) <- 0
opt_output_files(ctg) <- file.path(tdir, "tile_{XLEFT}_{YBOTTOM}")
opt_laz_compression(ctg) <- TRUE

invisible(catalog_retile(ctg))

tiles <- list.files(tdir, pattern = "\\.laz$", full.names = TRUE)
cat(sprintf("tile_aoi: wrote %d tiles to %s\n", length(tiles), tdir))
for (t in tiles) cat("  ", basename(t), "\n")
