#!/usr/bin/env Rscript
# Tree-top detection in lasR, following the documented approach:
# Step 0 measure density -> DERIVE params -> pit-free CHM -> variable-window LM.
suppressMessages(library(lasR))
f <- system.file("extdata", "MixedConifer.las", package = "lasR")

## Step 0 -- measure first-return density (1 m cells: count == pts/m²)
ans0 <- exec(rasterize(1, "count", filter = keep_first()), on = f)
v    <- terra::values(ans0); dens <- mean(v[!is.na(v) & v > 0])

## derive parameters from density (rules from the approach doc, §1)
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens)
wfloor  <- max(2, round(2.5 * spacing, 1))   # smallest crown ~2-3x spacing
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lasR Step 0: density=%.2f -> res=%.2f m, window floor=%.1f m\n",
            dens, res, wfloor))

## Steps 4-5 -- pit-free CHM + variable-window local maximum
# NOTE: lasR's `pit_fill` is not the Khosravipour pit-free algorithm used by
# lidR::pitfree(); it is a TIN + post-hoc pit filling. See comparison doc.
# Step 5 smoothing branch (approach §2 Step 5): QL2 (dens < 8) gets a 3x3 mean
# pre-LM smooth to suppress noise peaks; QL1/0 skips it.
t0   <- Sys.time()
del  <- triangulate(filter = keep_first())
chm  <- rasterize(res, del)
chm2 <- pit_fill(chm)
if (dens < 8) {
  smooth <- focal(chm2, size = 3, fun = "mean")
  seed   <- local_maximum_raster(smooth, ws, min_height = 2)
  ans    <- exec(del + chm + chm2 + smooth + seed, on = f)
} else {
  seed <- local_maximum_raster(chm2, ws, min_height = 2)
  ans  <- exec(del + chm + chm2 + seed, on = f)
}
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

tops <- ans[[length(ans)]]
xy   <- sf::st_coordinates(tops)            # 3D points: X, Y, Z(height)
out  <- data.frame(x = xy[, 1], y = xy[, 2], z = xy[, 3])
write.csv(out, file.path(Sys.getenv("CLAUDE_JOB_DIR"), "tops_lasr.csv"),
          row.names = FALSE)
cat(sprintf("lasR: %d treetops in %.2f s; Z range %.1f-%.1f m\n",
            nrow(out), dt, min(out$z), max(out$z)))
