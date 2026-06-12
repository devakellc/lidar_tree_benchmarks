# Tests for the random-walker th_cr stop rule (issue #35), the per-crown height
# cutoff that closes the over-grow gap of the pure-argmax random walker (#7).
# apply_thcr_cutoff() lives in sweep_lib.R and is a PURE relabel of the argmax
# label vector -> background below TH_CR * seed apex height; no work/ data needed.
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

test_that("apply_thcr_cutoff drops a short fringe pixel and keeps a tall one", {
  # crown 1 seed apex = 20 m; th_cr 0.55 -> cutoff 11 m.
  #  pixel 1: at the seed apex (20 m)        -> KEEP
  #  pixel 2: tall fringe (15 m, >= 11)      -> KEEP
  #  pixel 3: short fringe (8 m, < 11)        -> DROP to background
  #  pixel 4: background already (label 0)    -> stays background
  out_lab    <- c(1L, 1L, 1L, 0L)
  chm_vals   <- c(20, 15, 8, 4)
  lab_height <- c(20)                       # label 1 -> 20 m apex
  res <- apply_thcr_cutoff(out_lab, chm_vals, lab_height, th_cr = 0.55)
  expect_identical(res, c(1L, 1L, 0L, 0L))
})

test_that("apply_thcr_cutoff applies an independent cutoff per crown", {
  # crown 1 apex 20 (cutoff 11), crown 2 apex 6 (cutoff 3.3). A 5 m pixel is
  # below crown 1's cutoff but above crown 2's, so its fate depends on its label.
  lab_height <- c(20, 6)
  # same 5 m pixel, labelled to crown 1 then to crown 2
  expect_identical(
    apply_thcr_cutoff(c(1L), c(5), lab_height, th_cr = 0.55), 0L)   # < 11 -> drop
  expect_identical(
    apply_thcr_cutoff(c(2L), c(5), lab_height, th_cr = 0.55), 2L)   # > 3.3 -> keep
})

test_that("apply_thcr_cutoff leaves labels with no seed height untouched", {
  # label 2 has no apex height (NA) -> no finite cutoff -> never dropped, even
  # for a very short pixel. A label index beyond lab_height is treated the same.
  out_lab    <- c(1L, 2L, 3L)
  chm_vals   <- c(2, 2, 2)
  lab_height <- c(20, NA_real_)             # label 3 is beyond the vector
  res <- apply_thcr_cutoff(out_lab, chm_vals, lab_height, th_cr = 0.55)
  expect_identical(res, c(0L, 2L, 3L))      # only the seeded crown 1 is cut
})

test_that("apply_thcr_cutoff is a no-op on an all-background vector", {
  expect_identical(
    apply_thcr_cutoff(c(0L, 0L), c(5, 9), c(20), th_cr = 0.55),
    c(0L, 0L))
})

test_that("apply_thcr_cutoff errors on a length mismatch", {
  expect_error(apply_thcr_cutoff(c(1L, 1L), c(5), c(20), th_cr = 0.55))
})
