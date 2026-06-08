source(file.path("..", "..", "scripts", "analyze_model_benchmark.R"), local = TRUE)

# Two arms with DIFFERENT schemas (one has chm_res+tp_core, one does not),
# mimicking lidrplugins_results.csv vs ams3d_results.csv.
mk_arm <- function(det, rung, plot, n_ref, TP, n_det, with_extra) {
  base <- data.frame(site = "SOAP", plot = plot, plotType = "distributed",
                     detector = det, rung = rung, pdens = 14, frdens = 7,
                     n_apex = n_det, n_ref = n_ref, n_det = n_det, TP = TP,
                     recall = TP / n_ref, precision = TP / n_det, F1 = NA_real_,
                     height_rmse = 1,
                     rec_dominant = TP / n_ref, n_dominant = n_ref,
                     rec_codominant = NA_real_, n_codominant = 0,
                     rec_intermediate = NA_real_, n_intermediate = 0,
                     rec_suppressed = NA_real_, n_suppressed = 0,
                     rec_h_short = NA_real_, n_h_short = 0,
                     rec_h_mid = NA_real_, n_h_mid = 0,
                     rec_h_tall = TP / n_ref, n_h_tall = n_ref)
  if (with_extra) { base$chm_res <- 0.5; base$tp_core <- TP }
  base
}

test_that("harmonize_union fills missing cols and recomputes tp_core for all rows", {
  a <- mk_arm("chm_vwf", "8", "p1", 10, 6, 8, TRUE)
  b <- mk_arm("ams3d",   "8", "p1", 10, 5, 20, FALSE)   # no chm_res/tp_core
  u <- harmonize_union(list(a, b))
  expect_true(all(c("chm_res", "tp_core") %in% names(u)))
  expect_equal(nrow(u), 2L)
  # tp_core recomputed = round(precision * n_det) for BOTH rows
  expect_equal(u$tp_core[u$detector == "ams3d"], round((5 / 20) * 20))   # 5
  expect_true(is.na(u$chm_res[u$detector == "ams3d"]))                   # filled NA
})

test_that("pool_arms guards to common cells and pools per (detector,rung)", {
  u <- harmonize_union(list(
    mk_arm("chm_vwf", "8", "p1", 10, 6, 8, TRUE),
    mk_arm("ams3d",   "8", "p1", 10, 5, 20, FALSE),
    mk_arm("chm_vwf", "8", "p2", 10, 4, 9, TRUE)))   # p2 has no ams3d -> dropped
  pooled <- pool_arms(u, arms = c("chm_vwf", "ams3d"), rungs = "8")
  expect_setequal(unique(pooled$detector), c("chm_vwf", "ams3d"))
  # only p1/8 survives the guard, so chm_vwf recall = 6/10 (NOT pooling p2)
  expect_equal(pooled$recall[pooled$detector == "chm_vwf"], 0.6)
  expect_equal(pooled$n_plots[pooled$detector == "chm_vwf"], 1L)
  expect_true("rec_h_tall" %in% names(pooled))         # height bands flow through
})

test_that("deltas_vs_baseline subtracts the baseline arm per rung", {
  pooled <- data.frame(detector = c("chm_vwf", "ams3d"), rung = c("8", "8"),
                       recall = c(0.6, 0.5), F1 = c(0.5, 0.4),
                       rec_understory = c(0.2, 0.4))
  dl <- deltas_vs_baseline(pooled, baseline = "chm_vwf")
  expect_equal(dl$d_recall[dl$detector == "ams3d"], -0.1)
  expect_equal(dl$d_understory[dl$detector == "ams3d"], 0.2)
  expect_false("chm_vwf" %in% dl$detector)            # baseline itself excluded
})

test_that("model synthesis fails explicitly without the CHM-VWF baseline", {
  u <- harmonize_union(list(mk_arm("ams3d", "8", "p1", 10, 5, 20, FALSE)))
  expect_error(require_chm_vwf_baseline(u),
               "chm_vwf missing - run detect_lidrplugins_sweep.R first",
               fixed = TRUE)
})
