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

test_that(".read_detection_csv returns NULL on non-finite coordinates", {
  out <- tempfile(fileext = ".csv")
  writeLines("x y z\n1 foo 3\n4 5 Inf", out)
  expect_null(.read_detection_csv(out))
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

## ---- Docker backend (#19) ------------------------------------------------

# A stand-in `docker` CLI: records its full argv (one token per line) to
# $DOCKER_ARGV_FILE when set, then emulates `docker run ...` by stripping the run
# flags + image token and exec-ing the in-container command. Identity mounts
# (-v dir:dir) make in-container paths == host paths, so the command just runs.
.fake_docker <- function(dir) {
  p <- file.path(dir, "fake-docker.sh")
  writeLines(c(
    "#!/bin/sh",
    "[ -n \"$DOCKER_ARGV_FILE\" ] && printf '%s\\n' \"$@\" > \"$DOCKER_ARGV_FILE\"",
    "shift",                                  # drop 'run'
    "while [ $# -gt 0 ]; do",
    "  case \"$1\" in",
    "    --rm) shift ;;",
    "    --gpus) shift 2 ;;",
    "    -v) shift 2 ;;",
    "    *) break ;;",
    "  esac",
    "done",
    "shift",                                  # drop the image token
    "exec \"$@\""), p)
  Sys.chmod(p, "0755")
  p
}

test_that("run_docker_arm requires cmd so image CMD is not overwritten", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z\\n9 9 9", out)
  expect_error(run_docker_arm("img", input = file.path(d, "in.laz"),
                              out_csv = out, docker = .fake_docker(d),
                              cmd = character()),
               "cmd must be supplied")
  expect_true(file.exists(out))                 # API misuse must not unlink output
})

test_that("run_docker_arm runs a container and returns x,y,z; argv is well-formed", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "arm.py"); out <- file.path(d, "det.csv")
  writeLines(c("import sys",
               "open(sys.argv[2],'w').write('x y z\\n1 2 3\\n4 5 6\\n')"), py)
  argv <- file.path(d, "argv.txt"); Sys.setenv(DOCKER_ARGV_FILE = argv)
  on.exit(Sys.unsetenv("DOCKER_ARGV_FILE"), add = TRUE)
  det <- run_docker_arm("img:tag", input = file.path(d, "in.laz"), out_csv = out,
                        cmd = c("python3", py), docker = .fake_docker(d))
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L); expect_equal(det$z, c(3, 6))
  in_abs <- normalizePath(file.path(d, "in.laz"), mustWork = FALSE)
  out_abs <- normalizePath(out, mustWork = FALSE)
  mdir <- normalizePath(dirname(in_abs), mustWork = FALSE)   # == dirname(out_abs)
  expect_identical(readLines(argv),
                   c("run", "--rm", "--gpus", "all", "-v", paste0(mdir, ":", mdir),
                     "img:tag", "python3", py, in_abs, out_abs))
})

test_that("run_docker_arm: gpus=NULL omits --gpus and mounts= adds -v entries", {
  d <- tempfile(); dir.create(d); m1 <- tempfile(); dir.create(m1)
  py <- file.path(d, "arm.py"); out <- file.path(d, "det.csv")
  writeLines(c("import sys", "open(sys.argv[2],'w').write('x y z\\n1 2 3\\n')"), py)
  argv <- file.path(d, "argv.txt"); Sys.setenv(DOCKER_ARGV_FILE = argv)
  on.exit(Sys.unsetenv("DOCKER_ARGV_FILE"), add = TRUE)
  det <- run_docker_arm("img", input = file.path(d, "in.laz"), out_csv = out,
                        cmd = c("python3", py), mounts = m1, gpus = NULL,
                        docker = .fake_docker(d))
  expect_equal(nrow(det), 1L)
  rec <- readLines(argv)
  expect_false("--gpus" %in% rec)
  vidx <- which(rec == "-v")
  expect_length(vidx, 2L)                                    # data dir + m1
  expect_true(any(grepl(normalizePath(m1), rec[vidx + 1], fixed = TRUE)))
})

test_that("run_docker_arm handles input/output paths with spaces (shQuote)", {
  d <- file.path(tempdir(), "dir with spaces"); dir.create(d, showWarnings = FALSE)
  py <- file.path(d, "arm.py"); out <- file.path(d, "det out.csv")
  writeLines(c("import sys", "open(sys.argv[2],'w').write('x y z\\n7 8 9\\n')"), py)
  Sys.unsetenv("DOCKER_ARGV_FILE")
  det <- run_docker_arm("img", input = file.path(d, "in put.laz"), out_csv = out,
                        cmd = c("python3", py), docker = .fake_docker(d))
  expect_equal(det$z, 9)
})

test_that("run_docker_arm returns NULL on non-zero container exit (stale CSV present)", {
  d <- tempfile(); dir.create(d)
  out <- file.path(d, "stale.csv"); writeLines("x y z\n9 9 9", out)
  py <- file.path(d, "boom.py"); writeLines("import sys; sys.exit(3)", py)
  Sys.unsetenv("DOCKER_ARGV_FILE")
  expect_null(run_docker_arm("img", input = file.path(d, "in.laz"), out_csv = out,
                             cmd = c("python3", py), docker = .fake_docker(d)))
})

test_that("run_docker_arm returns NULL when the container writes no output", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "noop.py"); writeLines("pass", py)
  Sys.unsetenv("DOCKER_ARGV_FILE")
  expect_null(run_docker_arm("img", input = file.path(d, "in.laz"),
                             out_csv = file.path(d, "none.csv"),
                             cmd = c("python3", py), docker = .fake_docker(d)))
})

test_that("run_docker_arm returns a 0-row frame on a valid header with no rows", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "empty.py"); out <- file.path(d, "empty.csv")
  writeLines(c("import sys", "open(sys.argv[2],'w').write('x y z\\n')"), py)
  Sys.unsetenv("DOCKER_ARGV_FILE")
  det <- run_docker_arm("img", input = file.path(d, "in.laz"), out_csv = out,
                        cmd = c("python3", py), docker = .fake_docker(d))
  expect_identical(names(det), c("x", "y", "z")); expect_equal(nrow(det), 0L)
})

test_that("run_docker_arm honors a custom reader= (NULL skip, throw skip, valid det)", {
  d <- tempfile(); dir.create(d)
  py <- file.path(d, "cloud.py"); out <- file.path(d, "labeled.bin")
  writeLines(c("import sys", "open(sys.argv[2],'w').write('LABELED')"), py)  # non-CSV
  Sys.unsetenv("DOCKER_ARGV_FILE")
  base <- list(image = "img", input = file.path(d, "in.laz"), out_csv = out,
               cmd = c("python3", py), docker = .fake_docker(d))
  expect_null(do.call(run_docker_arm, c(base, list(reader = function(p) NULL))))
  expect_null(do.call(run_docker_arm, c(base, list(reader = function(p) stop("x")))))
  expect_null(do.call(run_docker_arm,
                      c(base, list(reader = function(p) data.frame(a = 1)))))
  expect_null(do.call(run_docker_arm,
                      c(base, list(reader = function(p) data.frame(x = 1, y = NaN, z = 3)))))
  det <- do.call(run_docker_arm,
                 c(base, list(reader = function(p) data.frame(x = 1, y = 2, z = 3))))
  expect_equal(det$z, 3)
})
