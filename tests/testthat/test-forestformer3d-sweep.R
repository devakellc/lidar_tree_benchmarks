suppressMessages({ library(lidR); library(terra) })
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)
source(file.path("..", "..", "scripts", "model_runner.R"), local = TRUE)

# Write a merged labelled LAZ (UTM): UserData = cylinder/block, PointSourceID =
# per-cylinder instance id. Same geometry as synth_block_points().
write_ff3d_laz <- function(path) {
  p <- synth_block_points()
  las <- LAS(data.frame(X = p$X, Y = p$Y, Z = p$Z,
                        UserData = as.integer(p$block),
                        PointSourceID = as.integer(p$inst)))
  lidR::writeLAS(las, path)
}

test_that("ff3d_collapse reads UserData/PointSourceID, dedups across blocks, reduces", {
  f <- tempfile(fileext = ".laz"); write_ff3d_laz(f)
  det <- ff3d_collapse(f, merge_tol = 2.0)
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 4L)                          # A, B, merged-C, D
  c_row <- det[abs(det$x - 40) < 1 & abs(det$y - 40) < 1, ]
  expect_equal(nrow(c_row), 1L); expect_equal(c_row$z, 12)
})

test_that("ff3d_collapse returns NULL on an unreadable file (schema failure -> skip)", {
  expect_null(ff3d_collapse(tempfile(fileext = ".laz")))
})

test_that("agl_guard: empty in -> 0-row; partial off-DTM -> AGL; all off-DTM -> NULL", {
  dtm <- tempfile(fileext = ".tif")
  r <- terra::rast(xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                   resolution = 1, vals = 5)            # flat ground at z=5
  terra::writeRaster(r, dtm, overwrite = TRUE)
  # empty in -> empty out (legit ran-but-empty)
  e <- data.frame(x = numeric(), y = numeric(), z = numeric())
  expect_equal(nrow(agl_guard(e, dtm)), 0L)
  # on-DTM apex -> AGL (z - 5)
  on <- data.frame(x = 50, y = 50, z = 25)
  g <- agl_guard(on, dtm); expect_equal(g$z, 20)
  # all apexes off the raster -> wholesale drop -> NULL (skip the cell)
  off <- data.frame(x = c(1e6, 1e6), y = c(1e6, 1e6), z = c(25, 30))
  expect_null(agl_guard(off, dtm))
})

# Spec §6 gated live smoke: push one real cylinder through the container and prove
# the centering-offset round-trip (the "fragile step", §5) — output coords land in
# the input UTM bbox, not centered near 0. Portable + opt-in: skips unless the
# operator points FF3D_REPO/FF3D_CKPT at a checkout and the ff3d-sm120 image exists.
test_that("LIVE (gated): a cylinder through ff3d_arm.py is restored to input UTM", {
  IMG  <- Sys.getenv("FF3D_IMAGE", "ff3d-sm120")
  REPO <- Sys.getenv("FF3D_REPO"); CKPT <- Sys.getenv("FF3D_CKPT")
  skip_if(REPO == "" || CKPT == "" || !dir.exists(REPO) || !file.exists(CKPT),
          "set FF3D_REPO + FF3D_CKPT to run the live FF3D smoke")
  skip_if(Sys.which("docker") == "", "docker not on PATH")
  skip_if(tryCatch(system2("docker", c("image", "inspect", IMG),
                           stdout = FALSE, stderr = FALSE) != 0,
                   error = function(e) TRUE), paste("image", IMG, "absent"))
  gpu <- file.path("..", "..", "gpu", "forestformer3d-sm120")
  ENTRY <- normalizePath(file.path(gpu, "ff3d_entry.sh"))
  DRIVER <- normalizePath(file.path(gpu, "ff3d_arm.py"))
  PATCH  <- normalizePath(file.path(gpu, "ff3d_repo.patch"))
  # synthetic raw-with-ground cloud at SOAP-ish UTM (easting ~3e5, northing ~4.1e6)
  set.seed(7); cx <- 299680; cy <- 4101160; gz <- 1100
  cone <- function(tx, ty, h, n = 500) { z <- runif(n, 0, h); r <- (1 - z / h) * 3
    data.frame(X = tx + r * cos(runif(n, 0, 2 * pi)),
               Y = ty + r * sin(runif(n, 0, 2 * pi)), Z = gz + z) }
  df <- rbind(
    data.frame(X = runif(4000, cx - 16, cx + 16), Y = runif(4000, cy - 16, cy + 16), Z = gz),
    cone(cx - 6, cy - 4, 18), cone(cx + 5, cy + 6, 15), cone(cx, cy, 20))
  ind <- tempfile("ffin"); dir.create(ind)
  suppressWarnings(lidR::writeLAS(lidR::LAS(df), file.path(ind, "cyl_000.laz")))
  outl <- tempfile(fileext = ".laz")
  det <- run_docker_arm(IMG, ind, outl, cmd = c("bash", ENTRY),
                        extra = c(CKPT, REPO, PATCH, DRIVER),
                        mounts = c(REPO, dirname(CKPT), dirname(ENTRY)),
                        reader = function(p) ff3d_collapse(p, merge_tol = 2),
                        gpus = "all", timeout = 1200)
  expect_false(is.null(det))                       # container ran (NULL = crash/missing)
  out <- suppressWarnings(lidR::readLAS(outl))      # read the merged LAZ directly
  expect_gt(nrow(out@data), 0L)
  expect_true("UserData" %in% names(out@data))      # block id carried
  expect_true(all(out@data$X > cx - 20 & out@data$X < cx + 20))  # UTM restored, not ~0
  expect_true(all(out@data$Y > cy - 20 & out@data$Y < cy + 20))
})
