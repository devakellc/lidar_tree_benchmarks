source(file.path("..", "..", "scripts", "detect_ams3d_sweep.R"), local = TRUE)  # path from tests/testthat/

test_that("det_ams3d returns the x,y,z detection contract on a 2-tree cloud", {
  las <- synth_las_normalized()
  det <- det_ams3d(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
  expect_gte(nrow(det), 1L)          # at least one crown found
  expect_true(max(det$z) > 8)        # apex is a real height, normalized
})

test_that("det_ams3d returns an empty 0-row frame, never NULL, when no crown", {
  empty <- LAS(data.frame(X = 0, Y = 0, Z = 0.1))
  st_crs(empty) <- 32611L
  det <- det_ams3d(empty, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})

test_that("det_ams3d returns NULL when the segmenter errors (crash != empty)", {
  testthat::local_mocked_bindings(
    segment_tree_crowns = function(...) stop("boom"), .package = "crownsegmentr")
  det <- det_ams3d(synth_las_normalized(), cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_null(det)
})

test_that("det_ams3d returns a 0-row frame when the segmenter ran but found no crown", {
  testthat::local_mocked_bindings(
    segment_tree_crowns = function(point_cloud, ...) {
      point_cloud@data$crown_id <- NA_integer_   # ran fine, assigned no crowns
      point_cloud
    }, .package = "crownsegmentr")
  det <- det_ams3d(synth_las_normalized(), cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})
