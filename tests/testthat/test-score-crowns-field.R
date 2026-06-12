source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)        # greedy_match
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)  # under test

# Two well-separated instances on a normalized clip:
#   id=1: a 4x4 m square of canopy points at Z~5, apex (0,0,5) -> 4 points,
#         hull caliper = diagonal sqrt(32), so >= min_pts (5? no -> 4 < 5).
# We want id=1 ABOVE the floor so its diameters carry through, so give it 6 pts.
#   id=2: a sparse 2-point instance far away (below min_pts -> diameters NA).
synth_instances <- function() {
  data.table::data.table(
    X = c(0, 4, 0, 4, 2, 2,    100, 100.5),
    Y = c(0, 0, 4, 4, 0, 4,    100, 100),
    Z = c(5, 5, 5, 5, 5, 5,    3,   3),
    crown_id = c(1, 1, 1, 1, 1, 1,  2, 2))
}

# Field stems: one sits on instance 1's apex (0,0) at field height 6 m (the apex
# z=5 is height-consistent: within [0.5*6, 6+8]); one on instance 2 (100,100) at
# field height 3.2 m. Carries crown_class + field crown diameters.
synth_field <- function() {
  list(
    stems = data.frame(
      individualID = c("S1", "S2"),
      E = c(0.1, 100.1), N = c(0.0, 100.0),
      height = c(6.0, 3.2),
      crown_class = c("dominant", "suppressed"),
      stringsAsFactors = FALSE),
    field_cd = data.frame(
      individualID = c("S1", "S2"),
      maxCrownDiameter = c(7.0, 1.5),
      ninetyCrownDiameter = c(5.0, 1.2),
      stringsAsFactors = FALSE))
}

test_that("instance_apex keeps the id and reports the max-Z apex per instance", {
  ap <- instance_apex(synth_instances(), id_col = "crown_id")
  expect_identical(names(ap), c("id", "x", "y", "z"))
  expect_setequal(ap$id, c(1, 2))
  expect_equal(ap$z[ap$id == 2], 3)            # both id=2 pts at z=3
})

test_that("score_crowns_against_field matches 1:1 and carries d_eq/d_caliper", {
  fx   <- synth_field()
  inst <- synth_instances()
  diam <- crown_diameter_table(inst, id_col = "crown_id", min_pts = 5)
  apex <- instance_apex(inst, id_col = "crown_id")
  out  <- score_crowns_against_field(diam, apex, fx$stems, fx$field_cd, tol = 4,
                                     site = "TEST", plot = "P1", algo = "li2012")
  expect_identical(names(out), CROWN_COLS)
  expect_equal(nrow(out), 2L)                  # both stems matched 1:1

  r1 <- out[out$individualID == "S1", ]
  expect_equal(r1$algo, "li2012")
  expect_equal(r1$crown_class, "dominant")
  expect_equal(r1$field_maxCD, 7.0)
  expect_equal(r1$field_ninetyCD, 5.0)
  # instance 1 (6 pts >= min_pts): caliper = diagonal of the 4x4 square
  expect_equal(round(r1$d_caliper, 2), round(sqrt(32), 2))
  expect_true(is.finite(r1$d_eq) && r1$d_eq > 0)
  # area is the equivalent-circle area implied by d_eq
  expect_equal(round(r1$area, 4), round(pi * (r1$d_eq / 2)^2, 4))

  # instance 2 (2 pts < min_pts): diameters NA, but the row + field CD persist
  r2 <- out[out$individualID == "S2", ]
  expect_true(is.na(r2$d_eq) && is.na(r2$d_caliper) && is.na(r2$area))
  expect_equal(r2$field_maxCD, 1.5)
})

test_that("score_crowns_against_field returns a 0-row canonical frame on no match", {
  fx   <- synth_field()
  # apex far from every stem -> no pair within tol
  apex <- data.frame(id = 1, x = 9000, y = 9000, z = 5)
  diam <- data.frame(id = 1, n_pts = 6L, d_eq = 4, d_caliper = 5)
  out  <- score_crowns_against_field(diam, apex, fx$stems, fx$field_cd, tol = 4,
                                     site = "TEST", plot = "P1", algo = "ptrees")
  expect_identical(names(out), CROWN_COLS)
  expect_equal(nrow(out), 0L)
})
