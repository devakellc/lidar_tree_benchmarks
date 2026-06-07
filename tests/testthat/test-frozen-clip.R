source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)
suppressMessages(library(lidR))

test_that("seed_for is deterministic and varies by key", {
  expect_equal(seed_for("SOAP", "SOAP_001", 8), seed_for("SOAP", "SOAP_001", 8))
  expect_false(seed_for("SOAP", "SOAP_001", 8) == seed_for("SOAP", "SOAP_001", 4))
  expect_false(seed_for("SOAP", "SOAP_001", 8) == seed_for("SOAP", "SOAP_002", 8))
  expect_type(seed_for("SOAP", "SOAP_001", 8), "integer")
})

test_that("seeded homogenize decimation is reproducible", {
  set.seed(123)
  big <- LAS(data.frame(X = runif(5000, 0, 50), Y = runif(5000, 0, 50),
                        Z = runif(5000, 0, 30)))
  s <- seed_for("SOAP", "SOAP_001", 8)
  set.seed(s); a <- decimate_points(big, homogenize(density = 8, res = 5))
  set.seed(s); b <- decimate_points(big, homogenize(density = 8, res = 5))
  expect_equal(npoints(a), npoints(b))
  expect_equal(a@data$X, b@data$X)     # identical point subset, not just count
})
