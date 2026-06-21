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

# #X1: download NEON high-res RGB camera mosaics (DP3.30010.001) for the tiles
# that overlap a site's live field stems, via byTileAOP -- the RGB the DeepForest
# arm (detect_deepforest_sweep.R) runs predict_tile on. Keyed off the stem E/N
# (the LiDAR downloader's pattern). YEAR defaults to 2021 to match the 2021
# LiDAR/ground-truth epoch the rest of the benchmark uses (RGB exists for SOAP/
# SJER/TEAK at 2021). NEON is reprocessing DP3.30010 for RELEASE-2026, so the
# year is recorded in the savepath.
#   Rscript scripts/neon_download_aop.R SITE=SOAP YEAR=2021
# Reads work/neon/<SITE>/ground_truth_stems.csv; writes work/neon/<SITE>/rgb/.
suppressMessages(library(neonUtilities))
args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
site <- if (is.null(A$SITE)) "SOAP" else A$SITE
year <- if (is.null(A$YEAR)) "2021" else A$YEAR
prov <- isTRUE(as.logical(if (is.null(A$PROVISIONAL)) "FALSE" else A$PROVISIONAL))
d  <- .job_dir(); nd <- file.path(d, "neon", site)
gt <- read.csv(file.path(nd, "ground_truth_stems.csv"))
lt <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
cat(sprintf("[%s] live trees: %d (RGB DP3.30010 %s, provisional=%s)\n",
            site, nrow(lt), year, prov))
savep <- file.path(nd, "rgb"); dir.create(savep, showWarnings = FALSE, recursive = TRUE)
options(timeout = 7200)
byTileAOP(dpID = "DP3.30010.001", site = site, year = year,
          easting = lt$E, northing = lt$N, buffer = 50,
          check.size = FALSE, savepath = savep, include.provisional = prov)
tif <- list.files(savep, pattern = "\\.tif$", recursive = TRUE)
cat(sprintf("[%s] downloaded %d RGB tiles -> %s\n", site, length(tif), savep))
