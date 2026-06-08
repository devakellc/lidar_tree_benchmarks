source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

mk_row <- function(site, plot, rung, det, n_ref, TP, n_det, tp_core) {
  data.frame(site = site, plot = plot, rung = rung, detector = det,
             n_ref = n_ref, TP = TP, n_det = n_det,
             precision = tp_core / n_det,            # score_plot emits precision, not tp_core
             recall = TP / n_ref,
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

test_that("pool pools height bands by summed counts when present", {
  hrow <- function(plot, n_ref, TP, tp_core, n_tall, rec_tall) {
    data.frame(site = "SOAP", plot = plot, rung = "native", detector = "x",
               n_ref = n_ref, TP = TP, n_det = TP, precision = tp_core / TP,
               recall = TP / n_ref,
               n_dominant = 0, rec_dominant = NA_real_,
               n_codominant = 0, rec_codominant = NA_real_,
               n_intermediate = 0, rec_intermediate = NA_real_,
               n_suppressed = 0, rec_suppressed = NA_real_,
               n_h_short = 0, rec_h_short = NA_real_,
               n_h_mid = 0, rec_h_mid = NA_real_,
               n_h_tall = n_tall, rec_h_tall = rec_tall)
  }
  df <- rbind(hrow("p1", 10, 5, 5, 4, 0.50),   # tall: 2 of 4
              hrow("p2", 10, 6, 6, 6, 0.50))    # tall: 3 of 6
  p <- pool(df)
  expect_equal(p$n_h_tall, 10)                  # summed n
  expect_equal(p$rec_h_tall, 0.5)               # (2+3)/(4+6), NOT mean(0.5,0.5)
  expect_true(is.na(p$rec_h_short))             # band absent -> NA, n=0
  expect_equal(p$n_h_short, 0)
})

test_that("pool omits height-band columns when n_h_* absent", {
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p2","8","ams3d",90,45,50,40))
  p <- pool(df)
  expect_false("rec_h_tall" %in% names(p))      # additive: no bands in, none out
})

test_that("equal_set_guard requires the named arms, not just a matching count", {
  # cell has TWO detectors but the required li2012 is missing (ams3d + a stray 'foo')
  df <- rbind(mk_row("SOAP","p1","8","ams3d",10,5,6,4),
              mk_row("SOAP","p1","8","foo",   10,5,6,4),
              mk_row("SOAP","p2","8","ams3d",10,5,6,4),
              mk_row("SOAP","p2","8","li2012",10,6,7,5))
  g <- equal_set_guard(df, arms = c("ams3d","li2012"))
  expect_true(all(paste(g$site,g$plot,g$rung) == "SOAP p2 8"))  # only p2 has BOTH required arms
  expect_true("SOAP::p1::8" %in% attr(g, "dropped"))            # p1 dropped despite 2 detectors
})
