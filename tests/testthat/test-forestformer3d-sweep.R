suppressMessages({ library(lidR); library(terra) })
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)

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
