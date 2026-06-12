# Tests for the marker-controlled (seeded) watershed helpers (issue #32):
# seeds_to_marker_raster() rasterizes treetops into a seed-label image and
# priority_flood_watershed() floods one basin per marker over a CHM. Both live in
# sweep_lib.R and are PURE (no work/ data); a tiny synthetic CHM exercises them.
suppressMessages(library(terra))
source(file.path("..", "..", "scripts", "sweep_lib.R"), local = TRUE)

# A 1 m-resolution CHM with two clear peaks separated by a low valley, so a
# seeded watershed must split the canopy into exactly two basins (one per seed).
#  - left half peaks at col 1-2, right half peaks at col 5-6
#  - col 4 is a valley (still canopy, >= hmin) where the two basins meet
make_two_peak_chm <- function() {
  m <- rbind(
    c(20, 18, 12,  9, 14, 19),
    c(19, 17, 11,  8, 13, 18),
    c(18, 16, 10,  8, 12, 17))
  r <- terra::rast(nrows = nrow(m), ncols = ncol(m),
                   xmin = 0, xmax = ncol(m), ymin = 0, ymax = nrow(m))
  terra::values(r) <- m                       # row-major == terra cell order
  r
}

test_that("seeds_to_marker_raster labels one cell per seed by treeID", {
  chm <- make_two_peak_chm()
  # two seeds: one over the left peak, one over the right peak (cell centres)
  seed_xy <- cbind(c(0.5, 5.5), c(2.5, 2.5))  # x in cols 1 & 6; y in top row
  lab <- seeds_to_marker_raster(chm, seed_xy, treeID = c(7L, 9L))
  expect_length(lab, terra::ncell(chm))
  expect_equal(sum(lab != 0L), 2L)            # exactly two marker cells
  expect_setequal(lab[lab != 0L], c(7L, 9L))  # carry the treeIDs, not 1..L
  # the marker cells are the ones the seed xy fell into
  expect_identical(lab[terra::cellFromXY(chm, seed_xy)], c(7L, 9L))
})

test_that("seeds_to_marker_raster first seed wins a shared cell, drops off-extent", {
  chm <- make_two_peak_chm()
  # two seeds in the SAME cell (col 1) + one seed off the CHM extent
  seed_xy <- cbind(c(0.5, 0.6, 99), c(2.5, 2.5, 2.5))
  lab <- seeds_to_marker_raster(chm, seed_xy, treeID = c(3L, 4L, 5L))
  expect_equal(sum(lab != 0L), 1L)            # one occupied cell, off-extent dropped
  expect_equal(lab[lab != 0L], 3L)            # first seed (treeID 3) wins
})

test_that("priority_flood_watershed grows exactly one basin per marker", {
  chm <- make_two_peak_chm()
  seed_xy <- cbind(c(0.5, 5.5), c(2.5, 2.5))
  markers <- seeds_to_marker_raster(chm, seed_xy, treeID = c(7L, 9L))
  lab <- priority_flood_watershed(chm, markers, hmin = 2)
  # every canopy pixel (all >= 8 m here, so all >= hmin=2) is claimed
  expect_equal(sum(lab != 0L), terra::ncell(chm))
  # exactly the two seed labels, no spurious basins
  expect_setequal(unique(lab[lab != 0L]), c(7L, 9L))
  # the seed cells keep their own label
  expect_identical(lab[terra::cellFromXY(chm, seed_xy)], c(7L, 9L))
  # left columns -> left seed (7), right columns -> right seed (9); the valley
  # (col 4) splits, so both basins are non-trivial in size
  m2 <- matrix(lab, nrow = nrow(chm), byrow = TRUE)
  expect_true(all(m2[, 1:2] == 7L))           # left peak basin
  expect_true(all(m2[, 5:6] == 9L))           # right peak basin
  expect_gt(sum(lab == 7L), 1L)
  expect_gt(sum(lab == 9L), 1L)
})

test_that("priority_flood_watershed leaves seed-less canopy islands as background", {
  # two disconnected canopy blocks; only the left block has a marker. The right
  # block is unreachable over the canopy mask -> stays background (not forced).
  m <- rbind(
    c(20, 18,  0, 15, 14),
    c(19, 17,  0, 13, 12))                     # col 3 is non-canopy (0 < hmin)
  r <- terra::rast(nrows = nrow(m), ncols = ncol(m),
                   xmin = 0, xmax = ncol(m), ymin = 0, ymax = nrow(m))
  terra::values(r) <- m
  markers <- seeds_to_marker_raster(r, cbind(0.5, 1.5), treeID = 4L)  # left block
  lab <- priority_flood_watershed(r, markers, hmin = 2)
  m2  <- matrix(lab, nrow = nrow(r), byrow = TRUE)
  expect_true(all(m2[, 1:2] == 4L))            # connected block flooded
  expect_true(all(m2[, 3] == 0L))              # non-canopy gap stays background
  expect_true(all(m2[, 4:5] == 0L))            # seed-less island stays background
})

test_that("priority_flood_watershed returns all background when no marker on canopy", {
  chm <- make_two_peak_chm()
  markers <- integer(terra::ncell(chm))        # no markers at all
  lab <- priority_flood_watershed(chm, markers, hmin = 2)
  expect_true(all(lab == 0L))
})
