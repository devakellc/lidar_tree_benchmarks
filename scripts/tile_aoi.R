#!/usr/bin/env Rscript
.bs_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.bs_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bs <- Find(file.exists, c(
  if (!is.null(.bs_ofile) && length(.bs_ofile) && nzchar(.bs_ofile))
    file.path(dirname(.bs_ofile), "bootstrap.R"),
  if (length(.bs_file)) file.path(dirname(sub("^--file=", "", .bs_file[1])),
                                  "bootstrap.R"),
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs, .bs_ofile, .bs_file)

# Retile the 25 ha AOI laz into a small grid of tiles under
# $CLAUDE_JOB_DIR/tiles/. catalog_retile() aligns tiles to its own grid origin,
# so a 505 m AOI with a target ~253 m chunk size produces a 3x3 grid (9 tiles).
# This gives the catalog/streaming scripts inter-tile seams to demonstrate
# buffering on without needing a true wall-to-wall dataset.
suppressMessages(library(lidR))
d  <- .job_dir()
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
