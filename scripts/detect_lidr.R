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

# Tree-top detection in lidR, following the documented approach:
# Step 0 measure density -> DERIVE params -> pit-free CHM -> variable-window LM.
suppressMessages(library(lidR))
f   <- system.file("extdata", "MixedConifer.las", package = "lasR")
las <- readLAS(f, select = "xyzr")          # height-normalized; keep ReturnNumber

## Step 0 -- measure first-return density the same way as the lasR script:
## mean of first-return counts in 1 m cells (count == pts/m²)
fr   <- filter_poi(las, ReturnNumber == 1L)
g    <- pixel_metrics(fr, ~length(Z), res = 1)
v    <- terra::values(g); dens <- mean(v[!is.na(v) & v > 0])

## derive parameters from density (identical rule to the lasR script)
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens)
wfloor  <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lidR Step 0: density=%.2f -> res=%.2f m, window floor=%.1f m\n",
            dens, res, wfloor))

## Steps 4-5 -- pit-free CHM + variable-window local maximum
# Step 5 smoothing branch (approach §2 Step 5): QL2 (dens < 8) gets a 3x3 mean
# pre-LM smooth to suppress noise peaks; QL1/0 skips it.
t0    <- Sys.time()
chm   <- rasterize_canopy(las, res = res,
           pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5),
                   subcircle = 0.2))
chm_lm <- if (dens < 8) terra::focal(chm, w = matrix(1/9, 3, 3), na.rm = TRUE) else chm
ttops <- locate_trees(chm_lm, lmf(ws = ws, hmin = 2, shape = "circular"))
dt    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

terra::writeRaster(chm, file.path(.job_dir(), "chm_shared.tif"),
                   overwrite = TRUE)   # for the shared-CHM controlled test
xy  <- sf::st_coordinates(ttops)
out <- data.frame(x = xy[, 1], y = xy[, 2], z = ttops$Z)
write.csv(out, file.path(.job_dir(), "tops_lidr.csv"),
          row.names = FALSE)

ref <- tryCatch(length(unique(na.omit(readLAS(f)@data$treeID))), error = function(e) NA)
cat(sprintf("lidR: %d treetops in %.2f s; Z range %.1f-%.1f m; ref treeID count=%s\n",
            nrow(out), dt, min(out$z), max(out$z), ref))
