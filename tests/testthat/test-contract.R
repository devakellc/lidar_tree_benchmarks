source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

test_that("assert_detection_contract passes a valid x,y,z frame", {
  ok <- data.frame(x = 1, y = 2, z = 3)
  expect_true(assert_detection_contract(ok))
  expect_true(assert_detection_contract(data.frame(x = numeric(), y = numeric(),
                                                    z = numeric())))
})

test_that("assert_detection_contract rejects contract violations", {
  expect_error(assert_detection_contract(NULL), "NULL")
  expect_error(assert_detection_contract(data.frame(X = 1, Y = 2, Z = 3)),
               "x, y, z")                                   # capital cols
  expect_error(assert_detection_contract(list(x = 1, y = 2, z = 3)),
               "data.frame")                                # not a data.frame
  expect_error(assert_detection_contract(
    data.frame(x = "a", y = 2, z = 3)), "numeric")          # non-numeric
})
