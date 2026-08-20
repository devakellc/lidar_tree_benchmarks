source(file.path("..", "..", "scripts", "detect_li2012_native.R"), local = TRUE)

test_that("det_li2012 returns the detection contract on a 2-tree clip", {
  las <- synth_las_normalized()
  det <- det_li2012(las, hmin = 2)
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
  expect_gte(nrow(det), 1L)                 # finds at least one apex
})

test_that("det_li2012 returns a 0-row frame (not NULL) on a no-canopy clip", {
  las <- synth_las_normalized()
  las@data$Z <- pmin(las@data$Z, 0.5)        # crush everything below hmin
  det <- det_li2012(las, hmin = 2)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)                # ran-but-empty -> legit recall 0
})

## ---- #V6 instance persistence ----------------------------------------------
test_that("det_li2012 persists the per-point treeID cloud when inst_path is set", {
  las <- synth_las_normalized()
  f <- tempfile(fileext = ".laz")
  on.exit(unlink(f), add = TRUE)
  det <- det_li2012(las, hmin = 2, inst_path = f)
  expect_true(file.exists(f))
  pts <- read_instance_points_laz(f, id_field = "treeID")
  expect_gt(nrow(pts), 0L)
  # reducing the persisted cloud reproduces the returned apex set exactly
  redet <- instances_to_det(pts, id_field = "crown_id")
  expect_equal(redet[order(redet$x), ], det[order(det$x), ],
               ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("det_li2012 writes no artifact when the canopy guard short-circuits", {
  las <- synth_las_normalized()
  las@data$Z <- pmin(las@data$Z, 0.5)
  f <- tempfile(fileext = ".laz")
  det <- det_li2012(las, hmin = 2, inst_path = f)
  expect_equal(nrow(det), 0L)
  expect_false(file.exists(f))
})
