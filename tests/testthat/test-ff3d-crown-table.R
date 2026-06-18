source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)        # greedy_match
source(file.path("..", "..", "scripts", "model_bench_lib.R"), local = TRUE)

# FF3D dedup -> crown-diameter glue (#34). ff3d_crown_table composes the #M8
# dedup_blocks() cross-block merge with crown_diameter_table + instance_apex,
# both keyed by global_id, so the driver derives per-TREE diameters (one row per
# global_id) -- never per-cylinder. Uses the shared synth_block_points() fixture:
#   A (b0,inst1 @10,10), B (b0,inst2 @11,10): same block, distinct trees;
#   C (b0,inst3 @40,40) + C' (b1,inst1 @40.1,40): cross-block dup -> ONE tree;
#   D (b1,inst2 @60,60): isolated. -> 4 global_ids, each 3 points except the
#   merged C which has 6 (>= crown_diameter_table min_pts=5).
test_that("ff3d_crown_table: per-global_id diameters after cross-block dedup", {
  ct <- ff3d_crown_table(synth_block_points(), merge_tol = 2.0, min_pts = 5)
  diam <- ct$diam; apex <- ct$apex

  # one diameter row + one apex row per merged tree (4 global_ids, not 5 cylinders)
  expect_identical(sort(names(diam)),
                   sort(c("id", "n_pts", "d_eq", "d_caliper")))
  expect_identical(names(apex), c("id", "x", "y", "z"))
  expect_equal(nrow(diam), 4L)
  expect_equal(nrow(apex), 4L)
  expect_setequal(diam$id, apex$id)

  # the merged-C tree (the only one with >= min_pts points) carries a finite
  # diameter; its apex is the cross-block max-Z (12 > 11.8).
  ap_c <- apex[abs(apex$x - 40) < 1 & abs(apex$y - 40) < 1, ]
  expect_equal(nrow(ap_c), 1L)
  expect_equal(ap_c$z, 12)
  d_c <- diam[diam$id == ap_c$id, ]
  expect_equal(d_c$n_pts, 6L)                 # 3 (b0,inst3) + 3 (b1,inst1)
  expect_true(is.finite(d_c$d_eq) && d_c$d_eq > 0)
  expect_true(is.finite(d_c$d_caliper) && d_c$d_caliper > 0)

  # the isolated 3-point trees stay below min_pts: row + n_pts present, NA diam.
  d_sparse <- diam[diam$n_pts < 5, ]
  expect_equal(nrow(d_sparse), 3L)
  expect_true(all(is.na(d_sparse$d_eq) & is.na(d_sparse$d_caliper)))
})

test_that("ff3d_crown_table feeds score_crowns_against_field end to end", {
  ct <- ff3d_crown_table(synth_block_points(), merge_tol = 2.0, min_pts = 5)
  # a field stem sitting on the merged-C apex (40,40) at a height-consistent 14 m
  # (apex z=12 is within [0.5*14, 14+8]); carries crown_class + field diameters.
  stems <- data.frame(individualID = "C1", E = 40.1, N = 40.0, height = 14,
                      crown_class = "dominant", stringsAsFactors = FALSE)
  field_cd <- data.frame(individualID = "C1", maxCrownDiameter = 6,
                        ninetyCrownDiameter = 4, stringsAsFactors = FALSE)
  out <- score_crowns_against_field(ct$diam, ct$apex, stems, field_cd, tol = 4,
                                    site = "T", plot = "P", algo = "forestformer3d")
  expect_identical(names(out), CROWN_COLS)
  expect_equal(nrow(out), 1L)
  expect_equal(out$algo, "forestformer3d")
  expect_equal(out$individualID, "C1")
  expect_true(is.finite(out$d_eq) && is.finite(out$d_caliper))
  expect_equal(out$field_maxCD, 6)
})

test_that("ff3d_crown_table is 0-row safe on an empty block table", {
  empty <- data.frame(block = integer(), inst = integer(),
                      X = numeric(), Y = numeric(), Z = numeric())
  ct <- ff3d_crown_table(empty)
  expect_equal(nrow(ct$diam), 0L)
  expect_equal(nrow(ct$apex), 0L)
  expect_identical(names(ct$apex), c("id", "x", "y", "z"))
})
