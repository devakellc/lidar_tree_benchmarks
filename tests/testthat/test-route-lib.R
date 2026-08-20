# Unit tests for the #P2 routing helpers in route_lib.R: oracle_pick (per-cell
# argmax-metric arm) and select_policy_rows (materialize a routing policy's chosen
# rows for pooling). Pure data.frame logic; no NEON data.
source(file.path("..", "..", "scripts", "route_lib.R"), local = TRUE)

# Two cells, two arms each. Cell A: arm y wins (0.8 > 0.5). Cell B: arm x wins
# (0.6 > 0.3). Carries the counts pool() would consume.
long_fix <- function() data.frame(
  site = "S", plot = c("A", "A", "B", "B"), rung = "native",
  arm  = c("x", "y", "x", "y"),
  F1   = c(0.5, 0.8, 0.6, 0.3),
  n_ref = c(10, 10, 8, 8), TP = c(4, 7, 5, 2), n_det = c(9, 9, 7, 7),
  stringsAsFactors = FALSE)

test_that("oracle_pick selects the max-metric arm per cell", {
  pk <- oracle_pick(long_fix(), value = "F1")
  expect_identical(names(pk), c("site", "plot", "rung", "arm"))
  expect_equal(nrow(pk), 2L)
  expect_equal(pk$arm[pk$plot == "A"], "y")
  expect_equal(pk$arm[pk$plot == "B"], "x")
})

test_that("oracle_pick breaks ties deterministically (first arm)", {
  d <- data.frame(site = "S", plot = "A", rung = "native",
                  arm = c("x", "y"), F1 = c(0.7, 0.7),
                  n_ref = 1, TP = 1, n_det = 1, stringsAsFactors = FALSE)
  expect_equal(oracle_pick(d, value = "F1")$arm, "x")
})

test_that("select_policy_rows returns exactly the chosen (cell,arm) rows", {
  long <- long_fix()
  picks <- data.frame(site = "S", plot = c("A", "B"), rung = "native",
                      arm = c("x", "y"), stringsAsFactors = FALSE)  # a fixed-ish policy
  rows <- select_policy_rows(long, picks)
  expect_equal(nrow(rows), 2L)
  expect_equal(rows$TP[rows$plot == "A"], 4)   # arm x at A
  expect_equal(rows$TP[rows$plot == "B"], 2)   # arm y at B
})

test_that("select_policy_rows on the oracle picks the per-cell best rows", {
  long <- long_fix()
  rows <- select_policy_rows(long, oracle_pick(long, value = "F1"))
  expect_setequal(rows$F1, c(0.8, 0.6))        # the two cell maxima
})

test_that("a missing pick (cell with no chosen arm) drops out, not errors", {
  long <- long_fix()
  picks <- data.frame(site = "S", plot = "A", rung = "native", arm = "y",
                      stringsAsFactors = FALSE)
  rows <- select_policy_rows(long, picks)
  expect_equal(nrow(rows), 1L); expect_equal(rows$plot, "A")
})

test_that("a pick naming an arm absent from a cell yields no row for that cell", {
  long <- long_fix()
  picks <- data.frame(site = "S", plot = "A", rung = "native", arm = "zzz",
                      stringsAsFactors = FALSE)
  expect_equal(nrow(select_policy_rows(long, picks)), 0L)
})
