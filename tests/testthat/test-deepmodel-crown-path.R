test_that("deep-model crown arm reads DTMs from frozen_clip layout", {
  env <- new.env(parent = globalenv())
  source(file.path("..", "..", "scripts", "crown_metrics_deepmodel.R"),
         local = env)

  nd <- file.path("work", "neon", "SOAP")
  expect_equal(
    env$frozen_dtm_path(nd, "SOAP", "SOAP_001", "native"),
    file.path(nd, "frozen", "SOAP", "SOAP_001", "native", "ground_dtm.tif"))
  expect_equal(
    env$frozen_dtm_path(nd, "SOAP", "SOAP_001", "8"),
    file.path(nd, "frozen", "SOAP", "SOAP_001", "8", "ground_dtm.tif"))
})
