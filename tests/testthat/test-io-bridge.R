source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)
suppressMessages({ library(lidR); library(terra) })

test_that("instances_to_det reduces a labeled cloud to per-id apexes, dropping 0", {
  dt <- data.frame(X = c(10,10,40,40, 25), Y = c(10,10,12,12, 25),
                   Z = c(18, 9, 12, 5,  3),
                   pred_itc = c(1L,1L,2L,2L, 0L))   # id 0 = unassigned
  det <- instances_to_det(dt, id_field = "pred_itc")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)                 # two trees; the 0-point dropped
  expect_setequal(det$z, c(18, 12))           # max-Z apex per id
})

test_that("instances_to_det is NA-safe when the id column already has NA", {
  dt <- data.frame(X = c(1,2), Y = c(1,2), Z = c(5,6),
                   pred_itc = c(NA_integer_, 0L))    # NA and 0 both unassigned
  expect_equal(nrow(instances_to_det(dt, id_field = "pred_itc")), 0L)
})

test_that("det_to_agl subtracts the DTM and reports dropped (off-raster) apexes", {
  r <- terra::rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100)
  terra::values(r) <- 50                              # flat ground at 50 m
  f <- tempfile(fileext=".tif"); terra::writeRaster(r, f)
  det <- data.frame(x = c(25, 75, 250), y = c(25, 75, 25), z = c(77, 62, 40))
  agl <- det_to_agl(det, f)                           # 3rd point is off the raster
  expect_equal(agl$z, c(27, 12))                      # 77-50, 62-50
  expect_equal(nrow(agl), 2L)
  expect_equal(attr(agl, "n_dropped"), 1L)
})

test_that("CRS/units round-trip: a known UTM stem reduces to its apex within tol", {
  e0 <- 320000; n0 <- 4100000              # UTM 11N; apex (max-Z) at (e0+0.1, n0+0.1)
  dt <- data.frame(X = c(e0+0.1, e0,   e0-0.1),
                   Y = c(n0+0.1, n0,   n0-0.1),
                   Z = c(27,     14,   6),
                   pred_itc = c(1L, 1L, 1L))
  las <- lidR::LAS(data.frame(X=dt$X, Y=dt$Y, Z=dt$Z)); sf::st_crs(las) <- 32611L
  las <- lidR::add_lasattribute(las, dt$pred_itc, "pred_itc", "instance id")
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  det <- read_instances_laz(f, id_field = "pred_itc")
  expect_equal(nrow(det), 1L)
  expect_lt(sqrt((det$x - (e0+0.1))^2 + (det$y - (n0+0.1))^2), 0.5)  # apex xy in UTM
  expect_equal(det$z, 27)                                            # apex elevation
})
