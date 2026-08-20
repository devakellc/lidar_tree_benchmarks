# Unit tests for the #P1 cross-arm fusion helpers in model_bench_lib.R:
# fuse_apexes (cross-arm single-linkage clustering within merge_tol + a height
# gate) and fusion_points (union / majority / height-layered operating points).
# All synthetic; no NEON data.
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

# Five apexes from three arms:
#   chm (0,0,20) + sat (0.1,0,20)   -> co-located, same height -> ONE tree (votes 2)
#   li  (0,0,5)                      -> same x,y but 15 m lower -> understory, SEPARATE
#   chm (50,50,18) + sat (50.1,50,18)-> co-located -> ONE tree (votes 2)
fx_arm <- c("chm", "sat", "li", "chm", "sat")
fx_x   <- c(0,     0.1,   0,    50,    50.1)
fx_y   <- c(0,     0,     0,    50,    50)
fx_z   <- c(20,    20,    5,    18,    18)

test_that("fuse_apexes merges co-located cross-arm apexes, never same-arm", {
  f <- fuse_apexes(fx_arm, fx_x, fx_y, fx_z, merge_tol = 2.0, z_tol = 5.0)
  expect_identical(names(f), c("cluster", "x", "y", "z", "votes", "arms"))
  expect_equal(nrow(f), 3L)                            # 3 distinct trees
  top <- f[abs(f$x) < 1 & f$z > 10, ]                  # the (0,0,20) overstory cluster
  expect_equal(top$votes, 2L); expect_equal(top$arms, "chm,sat")
  expect_equal(top$z, 20)                              # representative = max-Z member
})

test_that("fuse_apexes height gate keeps an understory apex separate from the overstory", {
  f <- fuse_apexes(fx_arm, fx_x, fx_y, fx_z, merge_tol = 2.0, z_tol = 5.0)
  under <- f[abs(f$x) < 1 & f$z < 10, ]                # li (0,0,5)
  expect_equal(nrow(under), 1L)
  expect_equal(under$votes, 1L); expect_equal(under$arms, "li")
})

test_that("fuse_apexes does not merge two apexes from the SAME arm", {
  f <- fuse_apexes(c("chm", "chm"), c(0, 0.1), c(0, 0), c(20, 20),
                   merge_tol = 2.0, z_tol = 5.0)
  expect_equal(nrow(f), 2L)                            # same arm -> distinct trees
})

test_that("fuse_apexes is 0-row safe", {
  expect_equal(nrow(fuse_apexes(character(0), numeric(0), numeric(0), numeric(0))), 0L)
})

## ---- fusion_points: union / majority / layered ---------------------------
# n_arms = 3 (chm, sat, li). union = all clusters; majority = votes >= ceil(3/2)=2;
# layered = chm-family apexes in the overstory (z >= 0.5*max) + point/deep apexes
# in the understory (z < 0.5*max), then fused.
test_that("fusion_points emits union, majority, and layered operating points", {
  fp <- fusion_points(fx_arm, fx_x, fx_y, fx_z, n_arms = 3,
                      merge_tol = 2.0, z_tol = 5.0,
                      chm_arms = c("chm", "multichm"), overstory_frac = 0.5)
  expect_equal(nrow(fp$union), 3L)                     # all three trees
  expect_equal(nrow(fp$majority), 2L)                  # only the two votes>=2 clusters
  expect_true(all(fp$majority$votes >= 2))
  # layered: H = 0.5*20 = 10. Keep chm z>=10 (both chm) + non-chm z<10 (li at 5).
  # -> overstory (0,0,20) & (50,50,18) from chm, understory (0,0,5) from li = 3.
  expect_equal(nrow(fp$layered), 3L)
  expect_true(any(fp$layered$z < 10))                  # understory apex retained
  expect_false(any(fp$layered$z > 10 & fp$layered$z != 20 & fp$layered$z != 18))
})
