source(file.path("..", "..", "scripts", "detect_multichm_sweep.R"), local = TRUE)

test_that("det_multichm_run returns the x,y,z contract on a 2-tree cloud", {
  det <- det_multichm_run(synth_las_normalized(), res = 0.5, a = 0.10)
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
})

test_that("det_multichm_run honors the contract (frame-or-NULL) on no trees", {
  # A flat near-ground cloud has no canopy; multichm must yield a 0-row frame
  # (or NULL on an internal crash) -- never a malformed or sf table.
  flat <- LAS(data.frame(X = runif(300, 0, 20), Y = runif(300, 0, 20),
                         Z = rep(0.1, 300)))
  st_crs(flat) <- 32611L
  det <- det_multichm_run(flat, res = 0.5, a = 0.10)
  expect_true(is.null(det) ||
              (is.data.frame(det) && !inherits(det, "sf") &&
               identical(names(det), c("x", "y", "z")) && nrow(det) == 0L))
})
