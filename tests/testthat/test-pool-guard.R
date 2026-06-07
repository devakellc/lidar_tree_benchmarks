source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

mk_row <- function(site, plot, rung, det, n_ref, TP, n_det, tp_core) {
  data.frame(site = site, plot = plot, rung = rung, detector = det,
             n_ref = n_ref, TP = TP, n_det = n_det, tp_core = tp_core,
             n_dominant = n_ref, rec_dominant = TP / n_ref,
             n_codominant = 0, rec_codominant = NA_real_,
             n_intermediate = 0, rec_intermediate = NA_real_,
             n_suppressed = 0, rec_suppressed = NA_real_,
             frdens = 10, secs = 1)
}

test_that("pool sums counts: recall = sum(TP)/sum(n_ref)", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p2","8","ams3d",90,45,50,40))
  p <- pool(df)
  expect_equal(p$n_ref, 100); expect_equal(p$TP, 50)
  expect_equal(p$recall, 0.5)            # NOT mean(0.5, 0.5)-by-rate weighting
  expect_equal(p$precision, 44 / 56)
})

test_that("equal_set_guard drops (site,plot,rung) missing any arm", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p1","8","li2012",10,6,7,5),
              mk_row("SOAP","p2","8","ams3d",10,5,6,4))  # li2012 missing on p2/8
  g <- equal_set_guard(df, arms = c("ams3d","li2012"))
  expect_equal(nrow(g), 2L)                       # only p1/8 survives
  expect_true(all(paste(g$site,g$plot,g$rung) == "SOAP p1 8"))
  expect_equal(attr(g, "dropped"), "SOAP::p2::8")
})
