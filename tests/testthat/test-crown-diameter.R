source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

test_that("crown_diameter_table reports d_eq, d_caliper, n_pts per instance", {
  # one wide instance (id=1, ~4 m caliper) and one sparse instance (id=2, 2 pts)
  pts <- data.table::data.table(
    X = c(0, 4, 0, 4,  100, 100.5),
    Y = c(0, 0, 4, 4,  100, 100),
    Z = c(5, 5, 5, 5,  3,   3),
    crown_id = c(1, 1, 1, 1,  2, 2))
  tab <- crown_diameter_table(pts, id_col = "crown_id", min_pts = 4)
  expect_identical(sort(names(tab)),
                   sort(c("id", "n_pts", "d_eq", "d_caliper")))
  r1 <- tab[tab$id == 1, ]
  expect_equal(r1$n_pts, 4L)
  expect_equal(round(r1$d_caliper, 2), round(sqrt(32), 2))  # diagonal of 4x4
  # sparse instance below the floor: diameters NA, but row + n_pts still present
  r2 <- tab[tab$id == 2, ]
  expect_equal(r2$n_pts, 2L)
  expect_true(is.na(r2$d_eq) && is.na(r2$d_caliper))
})
