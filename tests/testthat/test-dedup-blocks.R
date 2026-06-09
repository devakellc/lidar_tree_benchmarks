source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

test_that("dedup_blocks merges cross-block dups, keeps same-block neighbors distinct", {
  rel <- data.table::as.data.table(dedup_blocks(synth_block_points(), merge_tol = 2.0))
  gid <- function(bx, by) unique(rel[abs(X - bx) < 0.5 & abs(Y - by) < 0.5]$global_id)
  a <- gid(10, 10); b <- gid(11, 10)
  expect_length(a, 1L); expect_length(b, 1L)
  expect_false(a == b)                       # same block within tol -> distinct
  expect_equal(gid(40, 40), gid(40.1, 40))   # different blocks within tol -> merged
  expect_equal(length(unique(rel$global_id)), 4L)
})

test_that("dedup_blocks + reduce_instances yields the max-Z apex across merged blocks", {
  rel <- dedup_blocks(synth_block_points(), merge_tol = 2.0)
  det <- reduce_instances(rel, id_col = "global_id")
  expect_equal(nrow(det), 4L)
  c_row <- det[abs(det$x - 40) < 1 & abs(det$y - 40) < 1, ]
  expect_equal(nrow(c_row), 1L); expect_equal(c_row$z, 12)   # 12 > 11.8
})

test_that("dedup_blocks drops unassigned (0/NA) and is 0-row safe", {
  empty <- data.frame(block = integer(), inst = integer(),
                      X = numeric(), Y = numeric(), Z = numeric())
  rel <- dedup_blocks(empty)
  expect_equal(nrow(rel), 0L)
  det <- reduce_instances(rel, id_col = "global_id")
  expect_identical(names(det), c("x", "y", "z")); expect_equal(nrow(det), 0L)
})
