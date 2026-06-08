source(file.path("..", "..", "scripts", "model_runner.R"), local = TRUE)

test_that("run_python_arm returns x,y,z from a valid CSV the arm writes", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "arm.py"); out <- file.path(d, "det.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('x y z\\n1 2 3\\n4 5 6\\n')"), py)
  det <- run_python_arm("python3", py, input = "ignored.laz", out_csv = out)
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L); expect_equal(det$z, c(3, 6))
})

test_that("run_python_arm returns NULL on non-zero exit even if a stale CSV exists", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z\\n9 9 9", out)  # pre-existing
  py <- file.path(d, "boom.py"); writeLines("import sys; sys.exit(3)", py)
  expect_null(run_python_arm("python3", py, input = "x.laz", out_csv = out))
})

test_that("run_python_arm returns NULL on a wrong-schema CSV (contract failure)", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "bad.py"); out <- file.path(d, "bad.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('a b c\\n1 2 3\\n')"), py)  # not x,y,z
  expect_null(run_python_arm("python3", py, input = "x.laz", out_csv = out))
})

test_that("run_python_arm returns a 0-row frame on a valid header with no rows", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "empty.py"); out <- file.path(d, "empty.csv")
  writeLines(c("import sys", "open(sys.argv[2],'w').write('x y z\\n')"), py)
  det <- run_python_arm("python3", py, input = "x.laz", out_csv = out)
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 0L)                 # ran-but-empty -> legit recall 0
})

test_that("run_python_arm logs driver failures when label is set", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z\n9 9 9", out)
  py <- file.path(d, "boom.py"); writeLines("import sys; sys.exit(3)", py)
  expect_message(
    expect_null(run_python_arm("python3", py, input = "x.laz", out_csv = out,
                              label = "SOAP_001/native")),
    "SOAP_001/native: GPU driver failed")
})

test_that("run_python_crown_arm returns filtered labeled crown points", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "crowns.py"); out <- file.path(d, "crowns.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('x y z crown_id\\n1 2 3 7\\nNA 2 3 8\\n')"), py)
  pts <- run_python_crown_arm("python3", py, input = "x.laz", out_csv = out)
  expect_identical(names(pts), c("X", "Y", "Z", "crown_id"))
  expect_equal(nrow(pts), 1L)
  expect_equal(pts$crown_id, 7L)
})

test_that("run_python_crown_arm returns NULL on stale or wrong-schema output", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z crown_id\\n9 9 9 1", out)
  py <- file.path(d, "boom.py"); writeLines("import sys; sys.exit(3)", py)
  expect_null(run_python_crown_arm("python3", py, input = "x.laz", out_csv = out))

  bad <- file.path(d, "bad.py"); bad_out <- file.path(d, "bad.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('x y z\\n1 2 3\\n')"), bad)
  expect_null(run_python_crown_arm("python3", bad, input = "x.laz", out_csv = bad_out))
})
