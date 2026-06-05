#!/usr/bin/env Rscript
# Download NEON discrete-return LiDAR tiles (DP1.30003.001) that overlap the
# live field stems of a site, via byTileAOP. Usage:
#   Rscript scripts/neon_download_lidar.R SITE=TEAK YEAR=2021
# Reads work/neon/<SITE>/ground_truth_stems.csv (built by neon_ground_truth.R).
suppressMessages(library(neonUtilities))
args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
site <- if (is.null(A$SITE)) "SOAP" else A$SITE
year <- if (is.null(A$YEAR)) "2021" else A$YEAR
d  <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
nd <- file.path(d, "neon", site)

gt <- read.csv(file.path(nd, "ground_truth_stems.csv"))
lt <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
cat(sprintf("[%s] live trees: %d  unique 1km tiles: %d\n",
            site, nrow(lt), length(unique(lt$tile))))
savep <- file.path(nd, "lidar"); dir.create(savep, showWarnings = FALSE, recursive = TRUE)
options(timeout = 3600)
# buffer >= the per-plot clip reach (core_half + BUF = 45 m for tower plots in
# sweep_lib.R), so a plot whose clip box crosses a 1 km tile boundary always has
# the neighbouring tile present and the clip is never silently truncated.
byTileAOP(dpID = "DP1.30003.001", site = site, year = year,
          easting = lt$E, northing = lt$N, buffer = 50,
          check.size = FALSE, savepath = savep, include.provisional = FALSE)
laz <- list.files(savep, pattern = "\\.laz$", recursive = TRUE)
cat(sprintf("[%s] downloaded %d laz tiles\n", site, length(laz)))
