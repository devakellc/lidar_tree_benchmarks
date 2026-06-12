# Tests for pool_crown_by_rung() (issue #33: crown-diameter density ladder), the
# per-(algo, rung) crown-diameter pooler in sweep_lib.R. It must pool WITHIN each
# rung by the summed-squared-error rule (never a mean of per-plot RMSEs) and must
# default a missing/NA rung to "native" for backward compatibility with pre-#33
# crown_metrics_results.csv files. PURE helper; synthetic rows, no work/ data.
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

# Hand-computed reference RMSE/bias for a set of (det, fld) pairs.
ref <- function(det, fld) {
  e <- det - fld
  list(rmse = sqrt(mean(e^2)), bias = mean(e), n = length(det))
}

test_that("pool_crown_by_rung pools within rung by the summed-error rule", {
  # native rung: two plots' worth of trees pooled into ONE error stat (the rule
  # forbids averaging a per-plot RMSE -- a small plot must not be upweighted).
  # rung 4: a separate cell, larger errors.
  df <- data.frame(
    algo  = "dalponte2016",
    rung  = c("native", "native", "native", "4", "4"),
    d_eq        = c(5, 6, 7,  10, 12),
    field_ninetyCD = c(5, 5, 5,   5,  5),
    d_caliper   = c(6, 7, 8,  11, 13),
    field_maxCD = c(6, 6, 6,   6,  6),
    stringsAsFactors = FALSE)
  pr <- pool_crown_by_rung(df)

  # d_eq vs ninetyCD, native: errors 0,1,2 over the THREE pooled trees.
  nat <- pr[pr$rung == "native" &
            pr$definition == "d_eq vs ninetyCrownDiameter", ]
  expect_equal(nat$n, 3L)
  r <- ref(c(5, 6, 7), c(5, 5, 5))
  expect_equal(nat$rmse, r$rmse, tolerance = 1e-9)
  expect_equal(nat$bias, r$bias, tolerance = 1e-9)

  # rung 4 is pooled in its OWN cell (errors 5,7), not blended into native.
  r4 <- pr[pr$rung == "4" & pr$definition == "d_eq vs ninetyCrownDiameter", ]
  expect_equal(r4$n, 2L)
  rr <- ref(c(10, 12), c(5, 5))
  expect_equal(r4$rmse, rr$rmse, tolerance = 1e-9)
  expect_equal(r4$bias, rr$bias, tolerance = 1e-9)

  # both diameter definitions are present for each (algo, rung)
  expect_setequal(unique(pr$definition),
                  c("d_eq vs ninetyCrownDiameter", "d_caliper vs maxCrownDiameter"))
})

test_that("rows without a rung column default to native", {
  df <- data.frame(
    algo = "lasr_region_growing",
    d_eq = c(4, 5, 6), field_ninetyCD = c(5, 5, 5),
    d_caliper = c(4, 5, 6), field_maxCD = c(5, 5, 5),
    stringsAsFactors = FALSE)             # NO rung column at all
  pr <- pool_crown_by_rung(df)
  expect_true(all(pr$rung == "native"))
  expect_equal(unique(pr$n[pr$definition == "d_eq vs ninetyCrownDiameter"]), 3L)
})

test_that("NA / empty rung values normalize to native and pool together", {
  df <- data.frame(
    algo = "silva2016",
    rung = c(NA, "", "native"),
    d_eq = c(5, 6, 7), field_ninetyCD = c(5, 5, 5),
    d_caliper = c(5, 6, 7), field_maxCD = c(5, 5, 5),
    stringsAsFactors = FALSE)
  pr <- pool_crown_by_rung(df)
  # all three rows collapse into one native cell (n = 3), not three cells
  nat <- pr[pr$definition == "d_eq vs ninetyCrownDiameter", ]
  expect_equal(nrow(nat), 1L)
  expect_equal(nat$rung, "native")
  expect_equal(nat$n, 3L)
})

test_that("rungs are ordered native -> sparse and an empty frame is handled", {
  df <- data.frame(
    algo = "dalponte2016",
    rung = c("1", "native", "4", "8"),
    d_eq = c(3, 4, 5, 6), field_ninetyCD = c(5, 5, 5, 5),
    d_caliper = c(3, 4, 5, 6), field_maxCD = c(5, 5, 5, 5),
    stringsAsFactors = FALSE)
  pr <- pool_crown_by_rung(df)
  # within d_eq, n=1 per rung so RMSE is NA (crown_err_stats needs n>=2); but the
  # ORDER must still be native, 8, 4, 1 (the RUNG_LEVELS order, sparse last).
  oneeq <- pr[pr$definition == "d_eq vs ninetyCrownDiameter", ]
  expect_identical(oneeq$rung, c("native", "8", "4", "1"))

  empty <- pool_crown_by_rung(data.frame())
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("algo", "rung", "definition", "n", "rmse", "bias") %in%
                  names(empty)))
})
