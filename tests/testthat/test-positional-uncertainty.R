# Unit tests for the #V3 Monte-Carlo positional-uncertainty helper
# perturb_positions in sweep_lib.R: draws each stem's (E,N) from a 2-D Gaussian
# with per-stem sigma = pos_unc, reproducibly. Synthetic; no NEON data.
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

test_that("perturb_positions is deterministic for a given seed", {
  E <- c(0, 10, 20); N <- c(0, 5, 10); s <- c(0.5, 1.0, 0.3)
  a <- perturb_positions(E, N, s, seed = 123L)
  b <- perturb_positions(E, N, s, seed = 123L)
  expect_identical(names(a), c("E", "N"))
  expect_equal(a$E, b$E); expect_equal(a$N, b$N)
})

test_that("a zero-sigma stem never moves", {
  p <- perturb_positions(c(0, 10), c(0, 0), sigma = c(0, 2), seed = 1L)
  expect_equal(p$E[1], 0); expect_equal(p$N[1], 0)     # sigma 0 -> exact
  expect_true(p$E[2] != 10 || p$N[2] != 0)             # sigma 2 -> moved
})

test_that("different seeds give different draws", {
  E <- rep(0, 5); N <- rep(0, 5); s <- rep(1, 5)
  expect_false(isTRUE(all.equal(perturb_positions(E, N, s, 1L)$E,
                                perturb_positions(E, N, s, 2L)$E)))
})

test_that("draws are unbiased and scale with per-stem sigma", {
  set.seed(NULL)
  E <- c(100, 100); N <- c(200, 200); sig <- c(0.5, 3.0)
  K <- 4000
  dE <- matrix(NA_real_, K, 2)
  for (k in seq_len(K)) dE[k, ] <- perturb_positions(E, N, sig, seed = k)$E - E
  expect_lt(abs(mean(dE[, 1])), 0.05)                  # ~unbiased
  expect_lt(abs(mean(dE[, 2])), 0.20)
  expect_lt(abs(sd(dE[, 1]) - 0.5), 0.05)              # empirical sd ~ sigma
  expect_lt(abs(sd(dE[, 2]) - 3.0), 0.20)
})

test_that("perturb_positions recycles a scalar sigma and is length-stable", {
  p <- perturb_positions(c(0, 1, 2), c(0, 1, 2), sigma = 1.0, seed = 7L)
  expect_length(p$E, 3L); expect_length(p$N, 3L)
})

test_that("NA / negative sigma is treated as zero (no move, no NA out)", {
  p <- perturb_positions(c(5, 6), c(5, 6), sigma = c(NA, -1), seed = 3L)
  expect_equal(p$E, c(5, 6)); expect_equal(p$N, c(5, 6))
})
