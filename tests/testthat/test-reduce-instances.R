source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

test_that("reduce_instances collapses each id to its max-Z apex as x,y,z", {
  pts <- synth_labelled_points()          # tree 1 apex (10,10,18); tree 2 (40,12,12)
  det <- reduce_instances(pts, id_col = "id")
  expect_s3_class(det, "data.frame")
  expect_false(inherits(det, "sf"))
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)             # NA-id noise point dropped
  det <- det[order(det$x), ]
  expect_equal(det$z, c(18, 12))
  expect_equal(det$x, c(10, 40))
  expect_equal(det$y, c(10, 12))
})

test_that("reduce_instances returns a 0-row x,y,z frame on all-NA ids", {
  pts <- data.table::data.table(X = 1, Y = 1, Z = 1, id = NA_integer_)
  det <- reduce_instances(pts, id_col = "id")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})
