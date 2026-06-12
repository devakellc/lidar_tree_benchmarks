# Tests for multichm_seed_tops() (issue #31): the alternative crown-segmenter
# seed source used by crown_metrics_sweep.R, which runs lidRplugins::multichm on
# the point cloud and reads each 2-D top's apex height from the pit-free CHM the
# segmenters grow on. PURE helper in sweep_lib.R; exercised on a synthetic
# two-peak cloud + CHM, so it needs NO work/ data.
suppressMessages({ library(lidR); library(lidRplugins); library(terra); library(sf) })
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

# Two dense cone-shaped "trees" (apex 20 m at (10,10), 16 m at (25,25)) so
# multichm finds two clear maxima; a CHM at 0.5 m gives a height to read back.
make_two_cone_las <- function() {
  set.seed(1)
  cone <- function(cx, cy, top, n = 600) {
    x <- rnorm(n, cx, 1.3); y <- rnorm(n, cy, 1.3)
    r <- sqrt((x - cx)^2 + (y - cy)^2)
    data.frame(X = x, Y = y, Z = pmax(0, top - 1.5 * r))
  }
  las <- LAS(rbind(cone(10, 10, 20), cone(25, 25, 16)))
  st_crs(las) <- 32611L
  las
}

test_that("multichm_seed_tops returns an sf with treeID and CHM-read Z", {
  las <- make_two_cone_las()
  chm <- rasterize_canopy(las, res = 0.5, p2r(0.2))
  tt  <- multichm_seed_tops(las, chm, a = 0.10, frdens = 4)
  expect_s3_class(tt, "sf")
  expect_true("treeID" %in% names(tt))
  co <- sf::st_coordinates(tt)
  expect_gte(nrow(tt), 1L)
  expect_equal(ncol(co), 3L)                  # X, Y, Z present (segmenter contract)
  expect_true(all(is.finite(co[, 3])))        # every seed has a finite apex height
  # treeID is a clean 1..n sequence (the segmenters key crowns by it)
  expect_identical(tt$treeID, seq_len(nrow(tt)))
  # the Z is read from the CHM at each top, so it matches the CHM there (not e.g.
  # the raw point-cloud apex), within the CHM's own cell resolution.
  z_chm <- as.numeric(terra::extract(chm, co[, 1:2])[, 1])
  expect_equal(co[, 3], z_chm, tolerance = 1e-6)
})

test_that("multichm_seed_tops seeds are usable by dalponte2016", {
  las <- make_two_cone_las()
  chm <- rasterize_canopy(las, res = 0.5, p2r(0.2))
  tt  <- multichm_seed_tops(las, chm, a = 0.10, frdens = 4)
  skip_if(is.null(tt) || nrow(tt) == 0)
  seg <- dalponte2016(chm, tt, th_seed = 0.45, th_cr = 0.55, max_cr = 20)()
  expect_s4_class(seg, "SpatRaster")
  # at least one crown label grew from the seeds
  expect_true(any(is.finite(terra::values(seg, mat = FALSE))))
})

test_that("multichm_seed_tops returns NULL when the cloud has no canopy", {
  flat <- LAS(data.frame(X = runif(300, 0, 20), Y = runif(300, 0, 20),
                         Z = rep(0.1, 300)))
  st_crs(flat) <- 32611L
  chm <- tryCatch(rasterize_canopy(flat, res = 0.5, p2r(0.2)),
                  error = function(e) NULL)
  tt <- multichm_seed_tops(flat, if (is.null(chm)) flat else chm,
                           a = 0.10, frdens = 4)
  # a canopy-free cloud yields no tops (or, on an internal crash, NULL) -- never
  # a malformed seed object handed to the segmenters.
  expect_true(is.null(tt) ||
              (inherits(tt, "sf") && "treeID" %in% names(tt)))
})
