source(file.path("..", "..", "scripts", "detect_lidrplugins_sweep.R"), local = TRUE)

test_that("lmfauto and multichm return the x,y,z contract on a 2-tree cloud", {
  las <- synth_las_normalized()
  for (det in list(det_lmfauto(las, hmin = 2),
                   det_multichm(las, res = 0.5, a = 0.10))) {
    expect_s3_class(det, "data.frame")
    expect_false(inherits(det, "sf"))
    expect_identical(names(det), c("x", "y", "z"))
    expect_true(all(vapply(det, is.numeric, logical(1))))
  }
})

test_that("det_ptrees honors the contract (frame-or-NULL) on the toy cloud", {
  det <- det_ptrees(synth_las_normalized(), hmin = 2)
  expect_true(is.null(det) ||
              (is.data.frame(det) && !inherits(det, "sf") &&
               identical(names(det), c("x", "y", "z"))))
})
