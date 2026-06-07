#!/usr/bin/env Rscript
# Runs the benchmark library unit tests. Usage: Rscript tests/run_tests.R
suppressMessages(library(testthat))
# testthat 3.3.2 raises a hard error (not a warning) when no test files exist
# yet. Suppress that specific case so the harness exits 0 on an empty suite;
# all other errors are re-thrown.
tryCatch(
  test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE),
  error = function(e) {
    if (grepl("No test files found", conditionMessage(e))) {
      message("No test files found — harness is ready.")
    } else {
      stop(e)
    }
  }
)
