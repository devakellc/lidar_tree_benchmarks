# Tests for the shared point-cloud-detector library (pc_detect_lib.R) used by
# both detect_pc_sweep.R (#6, native-only) and detect_pc_ladder.R (#38, the
# native + 8 pts/m^2 density ladder). Covers the no-upsampling rung guard
# (pc_rungs_for) and characterization smoke tests for the three point-cloud apex
# extractors on the synthetic normalized LAS fixture (no NEON data required).
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)
source(file.path("..", "..", "scripts", "pc_detect_lib.R"), local = TRUE)

## ---- pc_rungs_for: no-upsampling guard (mirrors run_sweep.R) --------------
# A decimated rung is runnable only if the target density is STRICTLY below the
# plot's native all-return density; decimating to >= native would be upsampling
# (a no-op that aliases the rung onto native), which run_sweep.R skips with
# `rung >= native_pdens -> next`.

test_that("pc_rungs_for keeps a rung strictly below native density", {
  expect_equal(pc_rungs_for(15, 8), 8)
  expect_equal(pc_rungs_for(20, c(8, 4)), c(8, 4))
})

test_that("pc_rungs_for drops a rung at or above native density (no upsampling)", {
  expect_equal(pc_rungs_for(6, 8), numeric(0))   # 8 >= 6 native -> would upsample
  expect_equal(pc_rungs_for(8, 8), numeric(0))   # 8 >= 8 native -> alias of native
})

test_that("pc_rungs_for drops everything when native density is unknown", {
  expect_equal(pc_rungs_for(NA, 8), numeric(0))
})

## ---- extractor characterization smoke tests -------------------------------
# The fixture has two well-separated apexes at (10,10,18) and (40,12,12).

test_that("det_lidr_lmf_pc returns a contract-valid x,y,z frame with apexes", {
  las <- synth_las_normalized()
  det <- det_lidr_lmf_pc(las, ws_factory(0.10))
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
  expect_gte(nrow(det), 1L)
})

test_that("det_lidr_li2012 returns a contract-valid x,y,z frame", {
  las <- synth_las_normalized()
  det <- det_lidr_li2012(las)
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_true(all(vapply(det, is.numeric, logical(1))))
})

test_that("empty / sub-canopy input yields a 0-row x,y,z frame, never NULL", {
  empty <- suppressWarnings(lidR::LAS(data.frame(X = numeric(), Y = numeric(),
                                                 Z = numeric())))
  det <- det_lidr_lmf_pc(empty, ws_factory(0.10))
  expect_s3_class(det, "data.frame")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)
})
