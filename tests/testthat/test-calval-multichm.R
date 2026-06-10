# Unit tests for the multichm cal/val helpers (issue #39). Sourcing
# calval_multichm.R defines its helpers without running the analysis (guarded by
# `if (sys.nframe() == 0L) run_main()`); the chained source pulls in calval_lib.R
# (the shared split) and model_bench_lib.R (the canonical pooler). All tests here
# are pure data-frame logic -- no LiDAR.
source(file.path("..", "..", "scripts", "calval_multichm.R"), local = TRUE)

# minimal scored-row builder: one (plot, rung) cell carrying the columns the
# canonical pooler consumes. Crown-class counts default to zero unless supplied.
mk_cell <- function(plot, rung, n_ref, TP, n_det, prec = 0.5,
                    n_int = 0, rec_int = 0, n_sup = 0, rec_sup = 0, site = "X") {
  data.frame(site = site, plot = plot, rung = as.character(rung),
             n_ref = n_ref, n_det = n_det, TP = TP,
             recall = TP / n_ref, precision = prec, tp_core = round(prec * n_det),
             n_dominant = 0, rec_dominant = 0, n_codominant = 0, rec_codominant = 0,
             n_intermediate = n_int, rec_intermediate = rec_int,
             n_suppressed = n_sup, rec_suppressed = rec_sup,
             stringsAsFactors = FALSE)
}

test_that("res_rule picks 0.25 m only at >= 8 first-returns", {
  expect_equal(res_rule(8), 0.25)
  expect_equal(res_rule(7.99), 0.5)
  expect_equal(res_rule(c(3, 12)), c(0.5, 0.25))
})

test_that("assign_split is deterministic and balances strata", {
  pt <- data.frame(plot = paste0("p", 1:8),
                   plotType = rep(c("tower", "distributed"), each = 4),
                   n_over = 10, n_under = 0, crown_mix = "overstory",
                   stringsAsFactors = FALSE)
  a <- assign_split(pt, seed = 1, frac = 0.5)
  b <- assign_split(pt, seed = 1, frac = 0.5)
  expect_identical(a$split, b$split)                 # same seed -> same split
  # two strata of 4 -> ceil(0.5*4)=2 calib each
  expect_equal(sum(a$split == "calib"), 4L)
  expect_equal(sum(a$split == "valid"), 4L)
})

test_that("assign_split sends a singleton stratum to calibration", {
  pt <- data.frame(plot = paste0("p", 1:3),
                   plotType = c("tower", "tower", "distributed"),
                   n_over = 10, n_under = 0,
                   crown_mix = c("overstory", "overstory", "understory"),
                   stringsAsFactors = FALSE)
  a <- assign_split(pt, seed = 3, frac = 0.5)
  expect_equal(a$split[a$plot == "p3"], "calib")     # distributed/understory n=1
  expect_false(any(is.na(a$split)))
})

test_that("plot_table classifies crown_mix by the understory share rule", {
  r <- data.frame(plot = c("a", "b", "c"), plotType = "tower",
                  n_dominant = c(10, 10, 10), n_codominant = 0,
                  n_intermediate = c(0, 3, 2), n_suppressed = 0)
  pt <- plot_table(r)
  expect_equal(pt$crown_mix[pt$plot == "a"], "overstory")    # us = 0
  expect_equal(pt$crown_mix[pt$plot == "b"], "understory")   # us 3 >= os/4 = 2.5
  expect_equal(pt$crown_mix[pt$plot == "c"], "overstory")    # us 2 < 2.5
})

test_that("pair_cells keeps only (plot,rung) cells scored by BOTH arms", {
  mc <- rbind(mk_cell("p1", "2", 10, 6, 12), mk_cell("p2", "2", 10, 6, 12),
              mk_cell("p3", "2", 10, 6, 12))
  cv <- rbind(mk_cell("p1", "2", 10, 4, 8), mk_cell("p2", "2", 10, 4, 8),
              mk_cell("p3", "1", 10, 4, 8))               # p3 only at a diff rung
  P <- pair_cells(mc, cv, c("p1", "p2", "p3"))
  expect_setequal(paste(P$mc$plot, P$mc$rung), c("p1 2", "p2 2"))
  expect_setequal(paste(P$cv$plot, P$cv$rung), c("p1 2", "p2 2"))
  expect_equal(P$n_drop, 2L)                              # mc p3::2 + cv p3::1
})

test_that("h2h_pooled pools by summed counts and signs ΔF1 correctly", {
  # equal precision (0.5), multichm higher recall -> positive ΔF1.
  mc <- rbind(mk_cell("p1", "2", 10, 6, 12), mk_cell("p2", "2", 10, 6, 12))
  cv <- rbind(mk_cell("p1", "2", 10, 4, 8), mk_cell("p2", "2", 10, 4, 8))
  P <- pair_cells(mc, cv, c("p1", "p2"))
  h <- h2h_pooled(P)
  expect_equal(h$rec_mc, 12 / 20)        # sum(TP)/sum(n_ref), not mean of rates
  expect_equal(h$rec_cv, 8 / 20)
  expect_equal(h$prec_mc, 0.5)           # sum(tp_core)/sum(n_det)
  expect_gt(h$d_F1, 0)
})

test_that("h2h_pooled recovers pooled understory recall from int + suppressed", {
  # understory = intermediate + suppressed; pooled by summed counts.
  mc <- rbind(
    mk_cell("p1", "2", 10, 6, 12, n_int = 4, rec_int = 0.5, n_sup = 2, rec_sup = 0.5),
    mk_cell("p2", "2", 10, 6, 12, n_int = 0, rec_int = 0,   n_sup = 0, rec_sup = 0))
  cv <- rbind(
    mk_cell("p1", "2", 10, 4, 8, n_int = 4, rec_int = 0.25, n_sup = 2, rec_sup = 0),
    mk_cell("p2", "2", 10, 4, 8, n_int = 0, rec_int = 0,    n_sup = 0, rec_sup = 0))
  P <- pair_cells(mc, cv, c("p1", "p2"))
  h <- h2h_pooled(P)
  # mc: round(.5*4)+round(.5*2)=2+1=3 over 6 understory stems -> 0.5
  expect_equal(h$und_mc, 0.5)
  # cv: round(.25*4)+round(0*2)=1+0=1 over 6 -> ~0.1667
  expect_equal(h$und_cv, 1 / 6)
})
