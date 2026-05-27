#!/usr/bin/env Rscript
# Full approach on real USGS 3DEP (raw, multi-return) AOI -- lasR.
# drop noise -> Step 0 density -> normalize -> pit-free CHM -> variable-window LM.
suppressMessages(library(lasR))
f    <- file.path(Sys.getenv("CLAUDE_JOB_DIR"), "aoi.laz")
read <- reader_las(filter = "-drop_class 7 18 -drop_withheld")

## Step 0 -- first-return density (1 m cells), drives res + window floor
ans0 <- exec(read + rasterize(1, "count", filter = keep_first()), on = f)
v    <- terra::values(ans0); dens <- mean(v[!is.na(v) & v > 0])
res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lasR Step 0: density=%.2f -> res=%.2f m, window floor=%.1f m\n",
            dens, res, wfloor))

## Steps 1-5 (ground already classified) -> normalize -> pit-free CHM -> VWF
t0   <- Sys.time()
norm <- normalize()
del  <- triangulate(filter = keep_first())
chm  <- rasterize(res, del)
chm2 <- pit_fill(chm)
seed <- local_maximum_raster(chm2, ws, min_height = 2)
ans  <- exec(read + norm + del + chm + chm2 + seed, on = f)
dt   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

tops <- ans[[length(ans)]]
xy   <- sf::st_coordinates(tops)
out  <- data.frame(x = xy[, 1], y = xy[, 2], z = xy[, 3])
write.csv(out, file.path(Sys.getenv("CLAUDE_JOB_DIR"), "tops_lasr_aoi.csv"),
          row.names = FALSE)
cat(sprintf("lasR: %d treetops in %.1f s; height range %.1f-%.1f m\n",
            nrow(out), dt, min(out$z), max(out$z)))
