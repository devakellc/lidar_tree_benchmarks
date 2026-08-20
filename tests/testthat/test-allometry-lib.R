# Unit tests for the #S1 crown-allometry helpers in allometry_lib.R:
# functional_type (taxonID -> conifer/broadleaf), agb_from_dbh (Jenkins 2003
# generic AGB), and fit_stats (R2/RMSE/bias/MAE). Pure; no NEON data.
source(file.path("..", "..", "scripts", "allometry_lib.R"), local = TRUE)

test_that("functional_type classifies Sierra conifers and hardwoods", {
  expect_equal(functional_type(c("PIPO", "CADE27", "ABCO", "PSME", "PILA")),
               rep("conifer", 5))
  expect_equal(functional_type(c("QUKE", "QUCH2", "QUWI2", "ACMA3")),
               rep("broadleaf", 4))
  expect_true(is.na(functional_type("")))           # blank -> NA, not fabricated
})

test_that("agb_from_dbh matches the Jenkins 2003 generic forms", {
  # conifer (pine): exp(-2.5356 + 2.4349*ln(dbh)); dbh=30 cm -> ~313 kg
  expect_equal(agb_from_dbh(30, "conifer"), exp(-2.5356 + 2.4349 * log(30)),
               tolerance = 1e-9)
  expect_equal(round(agb_from_dbh(30, "conifer")), 313)
  # hardwood: exp(-2.48 + 2.4835*ln(dbh))
  expect_equal(agb_from_dbh(30, "broadleaf"), exp(-2.48 + 2.4835 * log(30)),
               tolerance = 1e-9)
  # monotone increasing in DBH
  expect_true(agb_from_dbh(40, "conifer") > agb_from_dbh(20, "conifer"))
})

test_that("agb_from_dbh is vectorized and NA-safe", {
  out <- agb_from_dbh(c(30, NA, 25), c("conifer", "conifer", NA))
  expect_length(out, 3L)
  expect_true(is.na(out[2]))                        # NA dbh -> NA
  expect_true(is.na(out[3]))                        # NA type -> NA
  expect_true(is.finite(out[1]))
})

test_that("fit_stats reports R2/RMSE/bias/MAE and a perfect fit", {
  s <- fit_stats(pred = c(1, 2, 3, 4), obs = c(1, 2, 3, 4))
  expect_identical(names(s), c("n", "r2", "rmse", "mae", "bias"))
  expect_equal(s$n, 4L); expect_equal(s$rmse, 0); expect_equal(s$bias, 0)
  expect_equal(s$r2, 1)
})

test_that("fit_stats bias and RMSE are signed/positive as expected", {
  # pred over-predicts obs by a constant 2 -> bias +2, rmse 2, r2 (vs mean) = ?
  s <- fit_stats(pred = c(3, 4, 5), obs = c(1, 2, 3))
  expect_equal(s$bias, 2); expect_equal(s$rmse, 2); expect_equal(s$mae, 2)
})

test_that("fit_stats needs >=2 finite pairs, else NA stats", {
  s <- fit_stats(pred = c(1, NA), obs = c(NA, 2))
  expect_equal(s$n, 0L); expect_true(is.na(s$rmse))
})
