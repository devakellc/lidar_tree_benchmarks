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

test_that("dedup_blocks keys on (block,inst): same inst id in different blocks merges", {
  # FF3D numbers instances PER cylinder, so inst id 1 recurs in every block. A tree
  # in the overlap is (block 0, inst 1) AND (block 1, inst 1) at ~the same XY — the
  # composite key + cross-block union must merge them into ONE global_id; a distant
  # instance with the same id must stay its own tree.
  pts <- data.table::data.table(
    block = c(0L, 0L,  1L, 1L,   1L,  1L),
    inst  = c(1L, 1L,  1L, 1L,   2L,  2L),     # inst 1 in BOTH block 0 and block 1
    X     = c(70, 70,  70.1, 70.1,  200, 200),
    Y     = c(70, 70,  70.0, 70.0,  200, 200),
    Z     = c(20, 15,  19,   14,    18,  9))
  rel  <- data.table::as.data.table(dedup_blocks(pts, merge_tol = 2.0))
  near <- unique(rel[abs(X - 70) < 1 & abs(Y - 70) < 1]$global_id)
  far  <- unique(rel[abs(X - 200) < 1]$global_id)
  expect_length(near, 1L)            # (b0,inst1) + (b1,inst1) at (70,70) -> merged
  expect_length(far, 1L)
  expect_false(near == far)          # the distant instance stays distinct
  expect_equal(length(unique(rel$global_id)), 2L)
})

test_that("dedup_blocks drops unassigned (0/NA) and is 0-row safe", {
  empty <- data.frame(block = integer(), inst = integer(),
                      X = numeric(), Y = numeric(), Z = numeric())
  rel <- dedup_blocks(empty)
  expect_equal(nrow(rel), 0L)
  det <- reduce_instances(rel, id_col = "global_id")
  expect_identical(names(det), c("x", "y", "z")); expect_equal(nrow(det), 0L)
})
