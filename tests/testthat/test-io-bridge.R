source(file.path("..", "..", "scripts", "io_bridge.R"), local = TRUE)
suppressMessages({ library(lidR); library(terra) })

test_that("instances_to_det reduces a labeled cloud to per-id apexes, dropping 0", {
  dt <- data.frame(X = c(10,10,40,40, 25), Y = c(10,10,12,12, 25),
                   Z = c(18, 9, 12, 5,  3),
                   pred_itc = c(1L,1L,2L,2L, 0L))   # id 0 = unassigned
  det <- instances_to_det(dt, id_field = "pred_itc")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)                 # two trees; the 0-point dropped
  expect_setequal(det$z, c(18, 12))           # max-Z apex per id
})

test_that("instances_to_det is NA-safe when the id column already has NA", {
  dt <- data.frame(X = c(1,2), Y = c(1,2), Z = c(5,6),
                   pred_itc = c(NA_integer_, 0L))    # NA and 0 both unassigned
  expect_equal(nrow(instances_to_det(dt, id_field = "pred_itc")), 0L)
})

test_that("det_to_agl subtracts the DTM and reports dropped (off-raster) apexes", {
  r <- terra::rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100)
  terra::values(r) <- 50                              # flat ground at 50 m
  f <- tempfile(fileext=".tif"); terra::writeRaster(r, f)
  det <- data.frame(x = c(25, 75, 250), y = c(25, 75, 25), z = c(77, 62, 40))
  agl <- det_to_agl(det, f)                           # 3rd point is off the raster
  expect_equal(agl$z, c(27, 12))                      # 77-50, 62-50
  expect_equal(nrow(agl), 2L)
  expect_equal(attr(agl, "n_dropped"), 1L)
})

test_that("CRS/units round-trip: a known UTM stem reduces to its apex within tol", {
  e0 <- 320000; n0 <- 4100000              # UTM 11N; apex (max-Z) at (e0+0.1, n0+0.1)
  dt <- data.frame(X = c(e0+0.1, e0,   e0-0.1),
                   Y = c(n0+0.1, n0,   n0-0.1),
                   Z = c(27,     14,   6),
                   pred_itc = c(1L, 1L, 1L))
  las <- lidR::LAS(data.frame(X=dt$X, Y=dt$Y, Z=dt$Z)); sf::st_crs(las) <- 32611L
  las <- lidR::add_lasattribute(las, dt$pred_itc, "pred_itc", "instance id")
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  det <- read_instances_laz(f, id_field = "pred_itc")
  expect_equal(nrow(det), 1L)
  expect_lt(sqrt((det$x - (e0+0.1))^2 + (det$y - (n0+0.1))^2), 0.5)  # apex xy in UTM
  expect_equal(det$z, 27)                                            # apex elevation
})

test_that("read_instances_laz is strict: a missing id field returns NULL, not 0-row", {
  las <- lidR::LAS(data.frame(X = c(1, 2), Y = c(1, 2), Z = c(5, 6)))
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  expect_null(read_instances_laz(f, id_field = "treeID"))   # schema fail -> skip
})

test_that("read_instances_laz checks schema before accepting an empty LAS", {
  las <- lidR::LAS(data.frame(X = numeric(), Y = numeric(), Z = numeric()))
  f <- tempfile(fileext = ".laz")
  ok <- tryCatch({ lidR::writeLAS(las, f); TRUE }, error = function(e) FALSE)
  skip_if(!ok, "lidR cannot write empty LAS fixture")
  expect_null(read_instances_laz(f, id_field = "treeID"))
})

## ---- full labelled-table reader (#34 deep-model crown arm) ---------------
test_that("read_instance_points_laz keeps the full table, maps id 0 -> NA", {
  # two instances + one unassigned (PredInstance 0) point; absolute UTM coords
  e0 <- 320000; n0 <- 4100000
  dt <- data.frame(X = c(e0, e0 + 0.1, e0 + 4, n0 * 0 + e0 + 50),
                   Y = c(n0, n0 + 0.1, n0,     n0 + 50),
                   Z = c(27, 14, 26, 3),
                   PredInstance = c(1L, 1L, 2L, 0L))
  las <- lidR::LAS(data.frame(X = dt$X, Y = dt$Y, Z = dt$Z))
  sf::st_crs(las) <- 32611L
  las <- lidR::add_lasattribute(las, dt$PredInstance, "PredInstance", "instance")
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  pts <- read_instance_points_laz(f, id_field = "PredInstance")
  expect_identical(names(pts), c("X", "Y", "Z", "crown_id"))
  expect_equal(nrow(pts), 4L)                       # full table, nothing dropped
  expect_true(is.na(pts$crown_id[pts$Z == 3]))      # the 0 label became NA
  # downstream reduce_instances drops the NA, yielding the two real apexes
  det <- reduce_instances(pts, id_col = "crown_id")
  expect_equal(nrow(det), 2L)
  expect_setequal(det$z, c(27, 26))
})

test_that("read_instance_points_laz is strict: missing id field -> NULL", {
  las <- lidR::LAS(data.frame(X = c(1, 2), Y = c(1, 2), Z = c(5, 6)))
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f)
  expect_null(read_instance_points_laz(f, id_field = "PredInstance"))
})

## ---- LAZ <-> PLY bridge (#19) --------------------------------------------

# Build a LAZ on disk with optional extra per-point attributes, return its path.
.laz_fixture <- function(dt, attrs = list()) {
  las <- lidR::LAS(dt)
  for (nm in names(attrs))
    las <- lidR::add_lasattribute(las, attrs[[nm]], nm, nm)
  f <- tempfile(fileext = ".laz"); lidR::writeLAS(las, f); f
}

test_that("laz_to_ply (double) round-trips XYZ and carries LAS fields", {
  set.seed(1); n <- 40
  dt <- data.frame(X = runif(n, 320000, 320020), Y = runif(n, 4100000, 4100020),
                   Z = runif(n, 0, 30),
                   Intensity = as.integer(runif(n, 0, 600)),
                   ReturnNumber = 1L, NumberOfReturns = 1L)
  f <- .laz_fixture(dt)
  ply <- tempfile(fileext = ".ply")
  info <- laz_to_ply(f, ply, fields = c("Intensity", "ReturnNumber"))
  expect_equal(info$offset, c(0, 0, 0))        # double -> absolute, no offset
  expect_equal(info$n, n)
  p <- read_ply(ply)
  r <- lidR::readLAS(f)
  expect_true(all(c("x", "y", "z", "Intensity", "ReturnNumber") %in% names(p)))
  expect_equal(nrow(p), n)
  expect_lt(max(abs(p$x - r$X)), 1e-9)         # double is bit-exact
  expect_lt(max(abs(p$z - r$Z)), 1e-9)
  expect_equal(p$Intensity, r$Intensity)       # ushort carried exactly
})

test_that("PLY double coords preserve absolute UTM losslessly (float32 would not)", {
  # float32 cannot hold a UTM northing ~4.1e6: its ulp there is well over 0.1 m.
  con <- rawConnection(raw(0), "wb"); writeBin(4100987.65, con, size = 4)
  back <- readBin(rawConnectionValue(con), "double", size = 4); close(con)
  expect_gt(abs(back - 4100987.65), 0.1)       # the reason the offset exists
  e0 <- 320123.456; n0 <- 4100987.654
  f <- .laz_fixture(data.frame(X = c(e0, e0 + 1), Y = c(n0, n0 + 1), Z = c(10, 20)))
  ply <- tempfile(fileext = ".ply"); laz_to_ply(f, ply)   # default double
  p <- read_ply(ply); r <- lidR::readLAS(f)
  expect_lt(max(abs(p$x - r$X)), 1e-6)         # sub-mm, lossless
  expect_lt(max(abs(p$y - r$Y)), 1e-6)
})

test_that("coord_type='float' subtracts an offset that read_ply restores", {
  e0 <- 320123.45; n0 <- 4100987.65
  f <- .laz_fixture(data.frame(X = c(e0, e0 + 5), Y = c(n0, n0 + 5), Z = c(10, 20)))
  ply <- tempfile(fileext = ".ply")
  info <- laz_to_ply(f, ply, coord_type = "float")
  expect_false(all(info$offset == 0))          # a local offset was chosen
  raw_local <- read_ply(ply)                   # offset defaults to c(0,0,0)
  expect_lt(max(raw_local$x), 100)             # stored coords are local/small
  p <- read_ply(ply, offset = info$offset)     # restore the absolute frame
  r <- lidR::readLAS(f)
  expect_lt(max(abs(p$x - r$X)), 0.05)         # float32 near 0..5 is fine
  expect_lt(max(abs(p$y - r$Y)), 0.05)
})

test_that("props emits a model-shaped PLY (lowercase names, dtype, const cols)", {
  f <- .laz_fixture(data.frame(X = c(1, 2), Y = c(3, 4), Z = c(5, 6),
                               Intensity = c(10L, 20L)))
  ply <- tempfile(fileext = ".ply")
  laz_to_ply(f, ply, props = list(
    intensity    = list(from = "Intensity", type = "float"),
    semantic_seg = list(const = 1L, type = "int"),
    treeID       = list(const = 0L, type = "int")))
  p <- read_ply(ply)
  pt <- attr(p, "properties")
  expect_identical(pt[["x"]], "double")
  expect_identical(pt[["intensity"]], "float")
  expect_identical(pt[["semantic_seg"]], "int")
  expect_identical(pt[["treeID"]], "int")
  expect_equal(p$intensity, c(10, 20))
  expect_equal(p$semantic_seg, c(1L, 1L))
  expect_equal(p$treeID, c(0L, 0L))
})

test_that("read_ply accepts CRLF headers and rejects truncated bodies", {
  crlf <- tempfile(fileext = ".ply")
  con <- file(crlf, "wb")
  writeBin(charToRaw(paste0("ply\r\nformat binary_little_endian 1.0\r\n",
                            "element vertex 0\r\nproperty double x\r\n",
                            "property double y\r\nproperty double z\r\n",
                            "end_header\r\n")), con)
  close(con)
  p <- read_ply(crlf)
  expect_identical(names(p), c("x", "y", "z"))
  expect_equal(nrow(p), 0L)

  truncated <- tempfile(fileext = ".ply")
  con <- file(truncated, "wb")
  writeBin(charToRaw(paste0("ply\nformat binary_little_endian 1.0\n",
                            "element vertex 1\nproperty double x\n",
                            "property double y\nproperty double z\n",
                            "property int treeID\nend_header\n")), con)
  writeBin(as.raw(c(0, 0)), con)
  close(con)
  expect_error(read_ply(truncated), "truncated body")
})

test_that("uint PLY properties round-trip as unsigned numeric values", {
  f <- .laz_fixture(data.frame(X = 1, Y = 2, Z = 3))
  ply <- tempfile(fileext = ".ply")
  laz_to_ply(f, ply, props = list(bigID = list(const = 4294967295, type = "uint")))
  p <- read_ply(ply)
  expect_identical(attr(p, "properties")[["bigID"]], "uint")
  expect_equal(p$bigID, 4294967295)
})

test_that("read_instances_ply reduces a labeled PLY to apexes; strict on id_field", {
  dt <- data.frame(X = c(10, 10, 40, 40, 25), Y = c(10, 10, 12, 12, 25),
                   Z = c(18, 9, 12, 5, 3))
  f <- .laz_fixture(dt, attrs = list(treeID = c(1L, 1L, 2L, 2L, 0L)))
  ply <- tempfile(fileext = ".ply")
  laz_to_ply(f, ply, props = list(treeID = list(from = "treeID", type = "int")))
  det <- read_instances_ply(ply, id_field = "treeID")
  expect_identical(names(det), c("x", "y", "z"))
  expect_equal(nrow(det), 2L)                  # two trees; the id-0 point dropped
  expect_setequal(det$z, c(18, 12))            # max-Z apex per id
  expect_null(read_instances_ply(ply, id_field = "instance_pred"))  # strict skip
})

test_that("ply_to_laz inverts laz_to_ply (XYZ round-trip through LAS)", {
  dt <- data.frame(X = c(320001.5, 320002.5), Y = c(4100001.5, 4100002.5),
                   Z = c(10, 20))
  f <- .laz_fixture(dt)
  ply <- tempfile(fileext = ".ply"); laz_to_ply(f, ply)
  laz2 <- tempfile(fileext = ".laz"); ply_to_laz(ply, laz2)
  r1 <- lidR::readLAS(f); r2 <- lidR::readLAS(laz2)
  expect_equal(nrow(r2@data), 2L)
  expect_lt(max(abs(sort(r2$X) - sort(r1$X))), 0.01)
  expect_lt(max(abs(sort(r2$Z) - sort(r1$Z))), 0.01)
})

test_that("ply_to_laz preserves carried fields and extra properties", {
  dt <- data.frame(X = c(320001.5, 320002.5), Y = c(4100001.5, 4100002.5),
                   Z = c(10, 20), Intensity = c(10L, 20L),
                   ReturnNumber = c(1L, 2L), NumberOfReturns = c(2L, 2L),
                   Classification = c(2L, 5L))
  f <- .laz_fixture(dt, attrs = list(treeID = c(7L, 8L)))
  ply <- tempfile(fileext = ".ply")
  laz_to_ply(f, ply, fields = c("Intensity", "ReturnNumber", "NumberOfReturns",
                                "Classification"),
             props = list(treeID = list(from = "treeID", type = "int")))
  laz2 <- tempfile(fileext = ".laz"); ply_to_laz(ply, laz2)
  r <- lidR::readLAS(laz2)
  expect_equal(r$Intensity, dt$Intensity)
  expect_equal(r$ReturnNumber, dt$ReturnNumber)
  expect_equal(r$NumberOfReturns, dt$NumberOfReturns)
  expect_equal(r$Classification, dt$Classification)
  expect_true("treeID" %in% names(r@data))
  expect_equal(r$treeID, c(7L, 8L))
})

test_that("emitted binary PLY is readable by Python plyfile (non-circular check)", {
  py <- Sys.which("python3")
  skip_if(!nzchar(py), "python3 not available")
  has_plyfile <- suppressWarnings(system2(py, c("-c", "import plyfile"),
                                          stdout = FALSE, stderr = FALSE)) == 0L
  skip_if(!isTRUE(has_plyfile), "python plyfile not installed")
  f <- .laz_fixture(data.frame(X = c(5, 6, 7), Y = c(8, 9, 10), Z = c(1, 2, 3),
                               Intensity = c(11L, 22L, 33L)))
  ply <- tempfile(fileext = ".ply")
  laz_to_ply(f, ply, fields = "Intensity")
  p <- read_ply(ply)
  code <- paste0("from plyfile import PlyData; v=PlyData.read(r'", ply,
                 "')['vertex']; ",
                 "print(len(v['x']), float(v['x'][0]), int(v['Intensity'][0]))")
  out <- suppressWarnings(system2(py, c("-c", code), stdout = TRUE, stderr = TRUE))
  vals <- as.numeric(strsplit(trimws(out[length(out)]), " ")[[1]])
  expect_equal(vals[1], nrow(p))               # plyfile sees the same vertices
  expect_lt(abs(vals[2] - p$x[1]), 1e-6)       # ... at the same coords
  expect_equal(vals[3], p$Intensity[1])        # ... with the field intact
})

## ---- write_instances_laz (#V6) ---------------------------------------------
test_that("write_instances_laz persists an id column readable by read_instance_points_laz", {
  las <- synth_las_normalized()
  ids <- rep(c(1L, 2L), each = 60)
  ids[c(1, 61)] <- NA_integer_           # unassigned -> 0 on disk -> NA on read
  las@data$treeID <- ids
  f <- tempfile(fileext = ".laz")
  on.exit(unlink(f), add = TRUE)
  out <- write_instances_laz(las, f, id_col = "treeID")
  expect_equal(out, f)
  pts <- read_instance_points_laz(f, id_field = "treeID")
  expect_equal(nrow(pts), lidR::npoints(las))
  expect_equal(sum(is.na(pts$crown_id)), 2L)
  expect_setequal(unique(stats::na.omit(pts$crown_id)), c(1L, 2L))
})

test_that("write_instances_laz preserves a non-treeID id column name", {
  las <- synth_las_normalized()
  las@data$crown_id <- rep(c(3L, 7L), each = 60)
  f <- tempfile(fileext = ".laz")
  on.exit(unlink(f), add = TRUE)
  write_instances_laz(las, f, id_col = "crown_id")
  pts <- read_instance_points_laz(f, id_field = "crown_id")
  expect_setequal(unique(pts$crown_id), c(3L, 7L))
})

test_that("write_instances_laz declines an empty LAS (no artifact = cell skipped)", {
  las <- synth_las_normalized()
  las@data$treeID <- 1L
  empty <- lidR::filter_poi(las, Z > 1000)
  f <- tempfile(fileext = ".laz")
  expect_true(is.na(write_instances_laz(empty, f, id_col = "treeID")))
  expect_false(file.exists(f))
})
