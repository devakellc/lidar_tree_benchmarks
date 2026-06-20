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

test_that("deep-model crown rows preserve their source rung", {
  env <- new.env(parent = globalenv())
  source(file.path("..", "..", "scripts", "crown_metrics_deepmodel.R"),
         local = env)

  env$RUNGS <- c("native", "8")
  env$plot_half <- function(plotType) 20
  env$det_to_agl <- function(apex, dtm) apex
  env$score_crowns_against_field <- function(diam_table, apex, stems, field_cd,
                                             tol, site, plot, algo) {
    data.frame(site = site, plot = plot, algo = algo,
               crown_class = "dominant", individualID = "T1",
               d_eq = 2, d_caliper = 3, area = pi,
               field_maxCD = 4, field_ninetyCD = 2,
               stringsAsFactors = FALSE)
  }

  nd <- tempfile("deepmodel")
  idir <- file.path(nd, "instances")
  dir.create(idir, recursive = TRUE)
  file.create(file.path(idir, "SOAP_001_native.laz"))
  file.create(file.path(idir, "SOAP_001_8.laz"))
  dir.create(file.path(nd, "frozen", "SOAP", "SOAP_001", "native"),
             recursive = TRUE)
  dir.create(file.path(nd, "frozen", "SOAP", "SOAP_001", "8"),
             recursive = TRUE)
  file.create(file.path(nd, "frozen", "SOAP", "SOAP_001", "native",
                        "ground_dtm.tif"))
  file.create(file.path(nd, "frozen", "SOAP", "SOAP_001", "8",
                        "ground_dtm.tif"))

  model <- list(dir = "instances", load = function(path) {
    list(diam = data.frame(id = 1L, n_pts = 6L, d_eq = 2, d_caliper = 3),
         apex = data.frame(id = 1L, x = 0, y = 0, z = 10))
  })
  pc <- data.frame(plotID = "SOAP_001", easting = 0, northing = 0,
                   plotType = "tower")
  gt <- data.frame(plotID = "SOAP_001", E = 0, N = 0, height = 10,
                   crown_class = "dominant", individualID = "T1")
  fc <- data.frame(individualID = "T1", maxCrownDiameter = 4,
                   ninetyCrownDiameter = 2)

  res <- env$run_plot_model("SOAP", "SOAP_001", "segmentanytree", model,
                            pc, gt, fc, nd)
  expect_identical(names(res)[1:3], c("site", "plot", "rung"))
  expect_equal(res$rung, c("native", "8"))
})
