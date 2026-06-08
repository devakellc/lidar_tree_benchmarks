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

# Full lasR-native AOI pipeline on remote USGS EPT (no PDAL extraction step).
# Requires lasR pre-devel EPT improvements and variable-window local maximum.
suppressMessages(library(lasR))

# CA_CarrHirzDeltaFires_2_2019 public EPT (EPSG:3857)
url <- Sys.getenv(
  "LASR_EPT_URL",
  unset = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public/CA_CarrHirzDeltaFires_2_2019/ept.json"
)

# Same AOI as scripts/extract.json (Web Mercator coordinates).
# Override via env: "xmin,ymin,xmax,ymax"
default_bbox <- c(-13578452, 4988399, -13577792, 4989059)
raw_bbox <- Sys.getenv("LASR_EPT_AOI_BBOX", unset = "")
if (nzchar(raw_bbox)) {
  bb <- suppressWarnings(as.numeric(strsplit(raw_bbox, ",", fixed = TRUE)[[1]]))
  if (length(bb) != 4L || any(is.na(bb)) || bb[1] >= bb[3] || bb[2] >= bb[4]) {
    stop("LASR_EPT_AOI_BBOX must be 'xmin,ymin,xmax,ymax' with xmin<xmax and ymin<ymax")
  }
  bbox <- bb
} else {
  bbox <- default_bbox
}
xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

# Optional EPT depth cap (default -1 = full depth)
depth <- as.integer(Sys.getenv("LASR_EPT_DEPTH", unset = "-1"))

parallel_mode <- tolower(Sys.getenv("LASR_EPT_PARALLEL", unset = "concurrent-files"))
ncores <- switch(
  parallel_mode,
  "sequential"        = sequential(),
  "concurrent-points" = concurrent_points(half_cores()),
  "concurrent-files"  = concurrent_files(half_cores()),
  stop("LASR_EPT_PARALLEL must be sequential|concurrent-points|concurrent-files")
)

# Query only the AOI; all acquisition is done by lasR's EPT reader.
read <- reader_rectangles(
  xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
  depth = depth,
  filter = "-drop_class 7 18 -drop_withheld"
)

## Step 0 -- use pre-devel native density fields from summarise()
stats <- exec(read + summarise(), on = url, ncores = ncores)
dens <- stats$pulse_density
if (is.null(dens) || is.na(dens) || dens <= 0) {
  dens <- stats$density
}
if (is.null(dens) || is.na(dens) || dens <= 0) {
  stop("Could not compute density from summarise() on EPT AOI")
}
if (dens < 1) stop("density < 1 pts/m²: metrics collapse; do not proceed")

res     <- if (dens >= 8) 0.25 else if (dens >= 4) 0.50 else 1.0
spacing <- 1 / sqrt(dens)
wfloor  <- max(2, round(2.5 * spacing, 1))
ws <- function(h) { y <- 0.1 * h + 3; y[h < 2] <- wfloor; y[h > 20] <- 5; y }
cat(sprintf("lasR EPT Step 0: density=%.2f -> res=%.2f m, window floor=%.1f m\n",
            dens, res, wfloor))

## Steps 1-5 -- normalize -> pit-free CHM -> variable-window local maximum
# NOTE: lasR's `pit_fill` is not the Khosravipour pit-free algorithm used by
# lidR::pitfree(); it is TIN + post-hoc pit filling (see comparison doc).
# Step 5 smoothing branch: QL2 (dens < 8) gets a 3x3 mean pre-LM smooth.
t0   <- Sys.time()
norm <- normalize()
del  <- triangulate(filter = keep_first())
chm  <- rasterize(res, del)
chm2 <- pit_fill(chm)
if (dens < 8) {
  smooth <- focal(chm2, size = 3, fun = "mean")
  seed   <- local_maximum_raster(smooth, ws, min_height = 2)
  ans    <- exec(read + norm + del + chm + chm2 + smooth + seed, on = url, ncores = ncores)
} else {
  seed <- local_maximum_raster(chm2, ws, min_height = 2)
  ans  <- exec(read + norm + del + chm + chm2 + seed, on = url, ncores = ncores)
}
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

tops <- ans[[length(ans)]]
xy   <- sf::st_coordinates(tops)
out  <- data.frame(x = xy[, 1], y = xy[, 2], z = xy[, 3])
write.csv(out, file.path(.job_dir(), "tops_lasr_ept_aoi.csv"),
          row.names = FALSE)
cat(sprintf("lasR EPT: %d treetops in %.1f s; height range %.1f-%.1f m\n",
            nrow(out), dt, min(out$z), max(out$z)))
