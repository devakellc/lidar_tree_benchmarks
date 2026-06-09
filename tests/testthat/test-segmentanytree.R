# SegmentAnyTree (#M6) host-side contract: the merged LAS carries PredInstance
# (0 = non-tree, 1..N = tree id) at ABSOLUTE UTM Z. The arm reads it with
# read_instances_laz(id_field="PredInstance") -> reduce_instances (max-Z per id,
# 0/NA dropped) -> det_to_agl(ground_dtm.tif). This exercises that whole host-
# side chain GPU-free; the container inference itself is a Docker smoke step.
source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)
suppressMessages({ library(lidR); library(terra) })

test_that("a PredInstance-labeled UTM LAS reduces to per-tree apexes above ground", {
  e0 <- 320000; n0 <- 4100000                  # UTM 11N, flat ground at 100 m
  # two trees (ids 1,2) + non-tree points (id 0); absolute UTM Z.
  dt <- data.frame(
    X = c(e0 + 0.1, e0,       e0 + 40,  e0 + 40.2, e0 + 20),
    Y = c(n0 + 0.1, n0,       n0 + 12,  n0 + 12,   n0 + 20),
    Z = c(127,      109,      118,      105,       101),
    PredInstance = c(1L, 1L, 2L, 2L, 0L))       # id 0 = non-tree, dropped
  las <- lidR::LAS(data.frame(X = dt$X, Y = dt$Y, Z = dt$Z))
  sf::st_crs(las) <- 32611L
  las <- lidR::add_lasattribute(las, dt$PredInstance, "PredInstance", "instance id")
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)

  det_abs <- read_instances_laz(f, id_field = "PredInstance")
  expect_identical(names(det_abs), c("x", "y", "z"))
  expect_equal(nrow(det_abs), 2L)               # two trees; non-tree id 0 dropped
  expect_setequal(round(det_abs$z), c(127, 118)) # max-Z apex per id, absolute

  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = e0 - 50, xmax = e0 + 50,
                   ymin = n0 - 50, ymax = n0 + 50)
  terra::values(r) <- 100                        # flat ground at 100 m
  dtm <- tempfile(fileext = ".tif"); terra::writeRaster(r, dtm)
  det <- det_to_agl(det_abs, dtm)                # absolute Z -> height above ground
  expect_setequal(round(det$z), c(27, 18))       # 127-100, 118-100
  expect_equal(attr(det, "n_dropped"), 0L)
})

test_that("a merged LAS without PredInstance is a schema failure (NULL, skip cell)", {
  las <- lidR::LAS(data.frame(X = c(1, 2), Y = c(1, 2), Z = c(5, 6)))
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  expect_null(read_instances_laz(f, id_field = "PredInstance"))
})
