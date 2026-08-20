# Unit tests for the #V4 matcher-hardening helpers in sweep_lib.R:
# match_tol (size/uncertainty-scaled per-stem tolerance), optimal_match
# (Hungarian optimal 1:1 assignment, drop-in for greedy_match), and fp_structure
# (false-positive error-structure split). All synthetic; no NEON data.
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

## ---- match_tol -----------------------------------------------------------
# tol_i = max(base_tol, k * crown_radius_i, pos_unc_i); NA/<=0 terms ignored.
test_that("match_tol scales tolerance by crown radius and positional uncertainty", {
  tol <- match_tol(crown_diam = c(10, 2, NA), pos_unc = c(0.5, 0.3, 6),
                   base_tol = 4, k = 1.0)
  expect_equal(tol, c(5, 4, 6))          # max(4,5,.5)=5; max(4,1,.3)=4; max(4,0,6)=6
})

test_that("match_tol floors at base_tol and honours k", {
  expect_equal(match_tol(crown_diam = c(10, 2), pos_unc = NULL, base_tol = 4, k = 0.5),
               c(4, 4))                   # k*radius = 2.5,0.5 -> both below floor
  expect_equal(match_tol(crown_diam = c(NA, NA), pos_unc = c(NA, NA)), c(4, 4))
  expect_equal(length(match_tol(crown_diam = numeric(0))), 0L)
})

## ---- optimal_match: drop-in agreement on uncontested sets ------------------
test_that("optimal_match matches within tol and leaves out-of-tol stems unmatched", {
  testthat::skip_if_not_installed("clue")
  expect_equal(optimal_match(0, 0, 1, 0, tol = 2), 1L)     # within tol -> det 1
  expect_equal(optimal_match(0, 0, 5, 0, tol = 2), 0L)     # outside tol -> unmatched
  # two cleanly separated pairs: same as greedy
  ax <- c(0, 20); ay <- c(0, 0); bx <- c(0.3, 20.2); by <- c(0, 0)
  expect_equal(optimal_match(ax, ay, bx, by, tol = 3),
               greedy_match(ax, ay, bx, by, tol = 3))
})

## ---- optimal_match beats greedy in a dense cluster -----------------------
# stem b's locally-shortest edge (b->d1, 0.4) steals d1 from stem a, whose only
# other det (d2 at 1.8) is out of tol -> greedy leaves a unmatched (1 match).
# Hungarian assigns a->d1 (0.6) + b->d2 (0.8) for 2 matches.
test_that("optimal_match recovers a match greedy loses to a closer neighbour", {
  testthat::skip_if_not_installed("clue")
  ax <- c(0, 1); ay <- c(0, 0); bx <- c(0.6, 1.8); by <- c(0, 0)
  g <- greedy_match(ax, ay, bx, by, tol = 1.0)
  o <- optimal_match(ax, ay, bx, by, tol = 1.0)
  expect_equal(sum(g > 0), 1L)            # greedy: only one stem matched
  expect_equal(sum(o > 0), 2L)            # optimal: both matched
  expect_equal(o, c(1L, 2L))
})

## ---- optimal_match honours the height gate (same as greedy) ---------------
test_that("optimal_match applies the [0.5*az, az+tol_z_up] height gate", {
  testthat::skip_if_not_installed("clue")
  # d1 is closer (0.5 m) but far too short (z=2 vs stem 20 m); d2 is height-OK.
  o <- optimal_match(0, 0, c(0.5, 1.0), c(0, 0), tol = 3,
                     az = 20, bz = c(2, 18), tol_z_up = 8)
  expect_equal(o, 2L)                     # the short det is gated out
})

## ---- optimal_match: non-square sets ---------------------------------------
test_that("optimal_match handles more detections than stems and vice versa", {
  testthat::skip_if_not_installed("clue")
  # 1 stem, 3 dets -> nearest only
  expect_equal(optimal_match(0, 0, c(2.5, 0.4, 5), c(0, 0, 0), tol = 3), 2L)
  # 3 stems, 1 det -> exactly one stem matched
  m <- optimal_match(c(0, 1, 2), c(0, 0, 0), 0.1, 0, tol = 3)
  expect_equal(sum(m > 0), 1L); expect_equal(m[1], 1L)
})

## ---- optimal_match: per-stem tol vector + soft 3-D cost -------------------
test_that("optimal_match accepts a per-stem tol vector", {
  testthat::skip_if_not_installed("clue")
  # stem1 tol 1 (det at 0.9 ok); stem2 tol 0.5 (det at 0.9 out) -> stem2 unmatched
  m <- optimal_match(c(0, 10), c(0, 0), c(0.9, 10.9), c(0, 0), tol = c(1.0, 0.5))
  expect_equal(m, c(1L, 0L))
})

test_that("optimal_match soft 3-D cost (lambda) prefers the height-consistent det", {
  testthat::skip_if_not_installed("clue")
  # both within horizontal tol; d1 closer in XY but 8 m off in height, d2 aligned.
  o <- optimal_match(0, 0, c(0.5, 0.9), c(0, 0), tol = 3,
                     az = 10, bz = c(2, 10), lambda = 1)
  expect_equal(o, 2L)                     # 3-D cost picks the aligned det
})

test_that("optimal_match soft 3-D cost still rejects height-IMPOSSIBLE matches", {
  testthat::skip_if_not_installed("clue")
  # The soft path must keep a GENEROUS-but-finite height envelope, not drop the
  # gate: a 10 m stem can never own a 50 m apex even at lambda > 0 (#V4 review).
  expect_equal(optimal_match(0, 0, 0.3, 0, tol = 4, az = 10, bz = 50, lambda = 0.5), 0L)
  # but a near-band miss the HARD gate rejects (bz just below 0.5*az) is still
  # recoverable inside the soft envelope:
  expect_equal(optimal_match(0, 0, 0.3, 0, tol = 4, az = 10, bz = 4.5, lambda = 0.5), 1L)
})

test_that("match_tol caps the crown-radius term so a corrupt crown record can't blow up tol", {
  # a 70 m maxCrownDiameter (radius 35) must clamp to tol_cap, mirroring ws_factory.
  expect_equal(match_tol(crown_diam = c(70, 6), pos_unc = NULL,
                         base_tol = 4, k = 1.0, tol_cap = 12),
               c(12, 4))
})

## ---- fp_structure --------------------------------------------------------
test_that("fp_structure splits false positives into near-a-stem vs isolated", {
  s <- fp_structure(fpx = c(1, 50), fpy = c(0, 50),
                    mstemx = c(0, 10), mstemy = c(0, 10), near_tol = 4)
  expect_equal(unname(s["near"]), 1L)     # (1,0) within 4 m of stem (0,0)
  expect_equal(unname(s["isolated"]), 1L) # (50,50) isolated
})

test_that("fp_structure is edge-safe", {
  expect_equal(unname(fp_structure(numeric(0), numeric(0), 0, 0, 4)), c(0L, 0L))
  expect_equal(unname(fp_structure(c(1, 2), c(0, 0), numeric(0), numeric(0), 4)),
               c(0L, 2L))                 # no matched stems -> all isolated
})

## ---- greedy_match accepts a per-stem tol vector (back-compat for scalar) ---
# stem1 (tol 5) reaches the det at 4 m; stem2 (tol 0.3) does not reach it at
# 3.5 m, so only stem1 matches. A scalar-tol recycling bug would instead let
# stem2 (closer) steal the det -> m = c(0,1); the per-stem fix gives m = c(1,0).
test_that("greedy_match applies a per-stem tolerance vector", {
  m <- greedy_match(ax = c(0, 0.5), ay = c(0, 0), bx = 4, by = 0,
                    tol = c(5, 0.3))
  expect_equal(m, c(1L, 0L))
  # scalar path unchanged
  expect_equal(greedy_match(0, 0, 1, 0, tol = 2), 1L)
})

## ---- score_plot: optimal method + FP-structure columns --------------------
test_that("score_plot exposes method='optimal' and FP-structure columns", {
  testthat::skip_if_not_installed("clue")
  stems <- data.frame(E = c(0, 1), N = c(0, 0), height = c(10, 10),
                      crown_class = c("dominant", "dominant"),
                      stringsAsFactors = FALSE)
  det <- data.frame(x = c(0.6, 1.8), y = c(0, 0), z = c(10, 10))
  g <- score_plot(stems, det, tol_xy = 1.0, core_cx = 0, core_cy = 0,
                  core_half = 20, method = "greedy")
  o <- score_plot(stems, det, tol_xy = 1.0, core_cx = 0, core_cy = 0,
                  core_half = 20, method = "optimal")
  expect_equal(g$TP, 1)                   # greedy loses one in the cluster
  expect_equal(o$TP, 2)                   # optimal recovers both
  expect_true(all(c("fp_near", "fp_isolated") %in% names(g)))
})
