#!/usr/bin/env Rscript
# #I4 I/O bridge. instances_to_det: labeled point table/LAS (instance-id field;
# 0 OR NA = unassigned) -> apex (x,y,z) via the bridge's reduce_instances
# (max-Z per id). det_to_agl: apex absolute elevation -> height above ground via
# the frozen clip's ground_dtm.tif, so z matches score_plot's height gate; it
# drops apexes that fall off the DTM and records the count in attr "n_dropped".
# LAZ<->PLY bridge (laz_to_ply / read_ply / ply_to_laz / read_instances_ply)
# feeds the PLY-ingesting SAT (#M6) / FF3D (#M8) Docker arms; see #19.
.find <- function(rel) Find(file.exists, c(file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("model_bench_lib.R"))
suppressMessages(library(terra))

instances_to_det <- function(pts, id_field = "pred_itc",
                             x = "X", y = "Y", z = "Z") {
  pts <- as.data.frame(pts)
  if (!id_field %in% names(pts))
    return(data.frame(x = numeric(), y = numeric(), z = numeric()))
  ids <- pts[[id_field]]
  pts[[id_field]][!is.na(ids) & ids == 0] <- NA      # 0 or NA = unassigned -> dropped
  det <- reduce_instances(pts, id_col = id_field, x = x, y = y, z = z)
  assert_detection_contract(det)
  det
}

read_instances_laz <- function(path, id_field = "pred_itc") {
  las <- lidR::readLAS(path)
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  if (is.null(las)) return(NULL)
  if (!id_field %in% names(las@data)) return(NULL)   # schema failure -> skip cell
  if (lidR::is.empty(las)) return(empty)              # legit empty with schema
  instances_to_det(las@data, id_field = id_field)
}

# Like read_instances_laz, but returns the FULL labelled point table (not just
# the apex) so the crown-diameter arm can derive both crown_diameter_table and
# instance_apex from the SAME labelling (#34). The instance id is normalized to
# an integer `crown_id` column with the unassigned label (0) mapped to NA, so
# the bridge's crown_diameter_table/instance_apex drop it. Strict on schema: a
# missing id_field returns NULL (-> skip the cell), never a fabricated 0-row; a
# valid-but-empty cloud returns a 0-row data.frame carrying the X,Y,Z,crown_id
# columns. Returns data.frame(X, Y, Z, crown_id).
read_instance_points_laz <- function(path, id_field = "pred_itc") {
  empty <- data.frame(X = numeric(), Y = numeric(), Z = numeric(),
                      crown_id = integer())
  las <- tryCatch(lidR::readLAS(path), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  if (!id_field %in% names(las@data)) return(NULL)   # schema failure -> skip cell
  if (lidR::is.empty(las)) return(empty)              # legit empty with schema
  d <- las@data
  ids <- as.integer(d[[id_field]])
  ids[!is.na(ids) & ids == 0L] <- NA_integer_         # 0 = unassigned -> dropped
  data.frame(X = d$X, Y = d$Y, Z = d$Z, crown_id = ids)
}

# Subtract ground elevation at each apex's (x,y). Apexes off the DTM are dropped;
# attr(.,"n_dropped") carries the count so callers can refuse a wholesale drop.
det_to_agl <- function(det, dtm_path) {
  if (!nrow(det)) { attr(det, "n_dropped") <- 0L; return(det) }
  g <- terra::extract(terra::rast(dtm_path), cbind(det$x, det$y))[, 1]
  ok <- !is.na(g)
  out <- det[ok, , drop = FALSE]; out$z <- out$z - g[ok]
  rownames(out) <- NULL
  attr(out, "n_dropped") <- as.integer(sum(!ok))
  out
}

## ---- ForestFormer3D arm helpers (#M8) ------------------------------------
# Reader for the merged per-plot LAZ ff3d_arm.py writes: UserData = cylinder
# (block), PointSourceID = per-cylinder instance id, coords UTM. Cross-block
# dedup -> canonical reduce -> apex det(x,y,z). NULL on an unreadable file
# (schema failure -> run_docker_arm skips the cell); a valid empty cloud -> a
# legit 0-row frame.
ff3d_collapse <- function(out_laz, merge_tol = 2.0) {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric())
  las <- tryCatch(lidR::readLAS(out_laz), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  if (lidR::is.empty(las)) return(empty)
  d <- las@data
  if (!all(c("UserData", "PointSourceID") %in% names(d))) return(NULL)
  pts <- data.frame(block = as.integer(d$UserData),
                    inst  = as.integer(d$PointSourceID),
                    X = d$X, Y = d$Y, Z = d$Z)
  reduce_instances(dedup_blocks(pts, merge_tol = merge_tol), id_col = "global_id")
}

# z -> AGL with a wholesale-off-DTM backstop. The FF3D centering-offset restore
# is the fragile step; a frame bug puts every apex off the DTM. Empty in -> empty
# out (legit ran-but-empty, recall 0). Non-empty in but ALL apexes dropped
# off-DTM -> NULL so the driver SKIPS the cell (a frame bug must not masquerade as
# a valid 0-recall row). Partial drops (edge apexes) are legitimate.
agl_guard <- function(det_abs, dtm_path) {
  if (!nrow(det_abs)) return(det_abs)
  det <- det_to_agl(det_abs, dtm_path)
  if (!nrow(det) && isTRUE(attr(det, "n_dropped") > 0)) return(NULL)
  det
}

## ---- binary little-endian PLY bridge (#19) -------------------------------
# Hand-rolled vertex-scalar PLY I/O (no ascii, no faces, no list properties) so
# the LAZ<->PLY conversions the SAT (#M6) / FF3D (#M8) Docker arms need stay
# unit-testable in pure R. Coordinates default to absolute UTM as `double`
# (lossless); `coord_type="float"` forces a local `offset` (returned on write,
# re-added on read) because float32 cannot hold a UTM northing to < ~0.5 m.

.PLY_SIZES <- c(char = 1L, uchar = 1L, short = 2L, ushort = 2L,
                int = 4L, uint = 4L, float = 4L, double = 8L)

.ply_default_type <- function(col) switch(col,
  Intensity = "ushort",
  ReturnNumber = , NumberOfReturns = , Classification = "uchar",
  gpstime = "double", "float")

.ply_uint_raw <- function(v) {
  dv <- as.double(v)
  bad <- !is.finite(dv) | dv < 0 | dv > 4294967295 | dv != floor(dv)
  if (any(bad))
    stop("laz_to_ply: uint values must be whole numbers in [0, 2^32 - 1]")
  cbind(as.raw(dv %% 256),
        as.raw(floor(dv / 256) %% 256),
        as.raw(floor(dv / 65536) %% 256),
        as.raw(floor(dv / 16777216) %% 256))
}

# values -> an N x size raw matrix (row i = the little-endian bytes of value i).
.ply_col_raw <- function(v, type) switch(type,
  double = matrix(writeBin(as.double(v), raw(), size = 8, endian = "little"),
                  ncol = 8, byrow = TRUE),
  float  = matrix(writeBin(as.double(v), raw(), size = 4, endian = "little"),
                  ncol = 4, byrow = TRUE),
  int = matrix(writeBin(as.integer(v), raw(), size = 4, endian = "little"),
               ncol = 4, byrow = TRUE),
  uint = .ply_uint_raw(v),
  short = , ushort = {
    iv <- as.integer(v)
    cbind(as.raw(bitwAnd(iv, 255L)), as.raw(bitwAnd(iv %/% 256L, 255L)))
  },
  char = , uchar = matrix(as.raw(bitwAnd(as.integer(v), 255L)), ncol = 1),
  stop("laz_to_ply: unsupported PLY type ", type))

# an N x size raw matrix -> the property numeric/integer vector.
.ply_read_col <- function(mat, type) {
  n <- nrow(mat); bytes <- as.vector(t(mat))
  switch(type,
    double = readBin(bytes, "double", n, size = 8, endian = "little"),
    float  = readBin(bytes, "double", n, size = 4, endian = "little"),
    int    = readBin(bytes, "integer", n, size = 4, endian = "little", signed = TRUE),
    uint   = as.numeric(mat[, 1]) + 256 * as.numeric(mat[, 2]) +
      65536 * as.numeric(mat[, 3]) + 16777216 * as.numeric(mat[, 4]),
    short  = readBin(bytes, "integer", n, size = 2, endian = "little", signed = TRUE),
    ushort = as.integer(mat[, 1]) + 256L * as.integer(mat[, 2]),
    char   = readBin(bytes, "integer", n, size = 1, endian = "little", signed = TRUE),
    uchar  = as.integer(mat[, 1]),
    stop("read_ply: unsupported PLY type ", type))
}

.find_bytes <- function(all, marker, required = TRUE) {
  k <- length(marker); lim <- min(length(all), 65536L)
  hit <- which(all[seq_len(lim)] == marker[1])
  for (i in hit)
    if (i + k - 1 <= length(all) && all(all[i:(i + k - 1)] == marker)) return(i)
  if (required) stop("read_ply: end_header marker not found (not a binary PLY?)")
  NA_integer_
}

# LAZ -> binary PLY. `fields`: LAS columns carried at LAS-native widths. `props`:
# named list(<out> = list(from = <lascol> | const = <value>, type = <ply type>))
# for model-shaped output (renames, dtype, constant columns). Returns
# list(n, offset) invisibly.
laz_to_ply <- function(laz_path, ply_path, fields = character(), props = NULL,
                       coord_type = c("double", "float"), offset = NULL) {
  coord_type <- match.arg(coord_type)
  las <- lidR::readLAS(laz_path)
  d <- las@data; n <- nrow(d)
  xyz <- cbind(as.double(d$X), as.double(d$Y), as.double(d$Z))
  if (is.null(offset))
    offset <- if (n > 0 && coord_type == "float")
      floor(apply(xyz, 2, min)) else c(0, 0, 0)
  xyz <- sweep(xyz, 2, offset)
  cols <- list(x = list(v = xyz[, 1], t = coord_type),
               y = list(v = xyz[, 2], t = coord_type),
               z = list(v = xyz[, 3], t = coord_type))
  for (f in fields) {
    if (!f %in% names(d)) stop("laz_to_ply: field '", f, "' not in the LAS")
    cols[[f]] <- list(v = d[[f]], t = .ply_default_type(f))
  }
  for (nm in names(props)) {
    spec <- props[[nm]]
    ty <- if (is.null(spec$type)) "float" else spec$type
    if (!is.null(spec$const)) {
      v <- rep(spec$const, length.out = n)
    } else {
      from <- if (is.null(spec$from)) nm else spec$from
      if (!from %in% names(d))
        stop("laz_to_ply: prop source '", from, "' not in the LAS")
      v <- d[[from]]
    }
    cols[[nm]] <- list(v = v, t = ty)
  }
  pnames <- names(cols); ptypes <- vapply(cols, `[[`, "", "t")
  body <- if (n > 0)
    as.vector(t(do.call(cbind, lapply(cols, function(c) .ply_col_raw(c$v, c$t)))))
  else raw(0)
  header <- paste0("ply\n", "format binary_little_endian 1.0\n",
                   "comment generated by io_bridge.R laz_to_ply\n",
                   "element vertex ", n, "\n",
                   paste0("property ", ptypes, " ", pnames, collapse = "\n"), "\n",
                   "end_header\n")
  con <- file(ply_path, "wb"); on.exit(close(con))
  writeBin(charToRaw(header), con); writeBin(body, con)
  invisible(list(n = n, offset = offset))
}

# Binary vertex-scalar PLY -> data.frame (one column per property, named by the
# header), with `offset` added back to x/y/z. attr(.,"properties") = name->type.
read_ply <- function(ply_path, offset = c(0, 0, 0)) {
  all <- readBin(ply_path, "raw", n = file.size(ply_path))
  marker <- charToRaw("end_header\n")
  pos <- .find_bytes(all, marker, required = FALSE)
  if (is.na(pos)) {
    marker <- charToRaw("end_header\r\n")
    pos <- .find_bytes(all, marker, required = FALSE)
  }
  if (is.na(pos)) stop("read_ply: end_header marker not found (not a binary PLY?)")
  lines <- strsplit(rawToChar(all[seq_len(pos - 1)]), "\n", fixed = TRUE)[[1]]
  lines <- sub("\r$", "", lines)
  fmt <- grep("^format ", lines, value = TRUE)
  if (!length(fmt) || !grepl("binary_little_endian", fmt))
    stop("read_ply: only binary_little_endian is supported")
  vline <- grep("^element vertex ", lines, value = TRUE)
  n <- as.integer(sub("^element vertex ", "", vline[1]))
  if (!length(vline) || is.na(n) || n < 0)
    stop("read_ply: invalid element vertex line")
  pp <- strsplit(grep("^property ", lines, value = TRUE), " ", fixed = TRUE)
  ptypes <- vapply(pp, `[`, "", 2L); pnames <- vapply(pp, `[`, "", 3L)
  sizes <- as.integer(.PLY_SIZES[ptypes])
  if (anyNA(sizes)) stop("read_ply: unknown property type(s)")
  row_size <- sum(sizes); expected <- n * row_size
  body_start <- pos + length(marker)
  body <- if (body_start <= length(all)) all[body_start:length(all)] else raw(0)
  if (length(body) < expected)
    stop("read_ply: truncated body; expected ", expected, " bytes, got ", length(body))
  if (n > 0) {
    rows <- matrix(body[seq_len(expected)], ncol = row_size, byrow = TRUE)
    off <- 1L; out <- vector("list", length(pnames))
    for (i in seq_along(pnames)) {
      out[[i]] <- .ply_read_col(rows[, off:(off + sizes[i] - 1), drop = FALSE],
                                ptypes[i])
      off <- off + sizes[i]
    }
    df <- as.data.frame(setNames(out, pnames), check.names = FALSE,
                        stringsAsFactors = FALSE)
  } else {
    df <- as.data.frame(setNames(rep(list(numeric(0)), length(pnames)), pnames),
                        check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (all(c("x", "y", "z") %in% names(df))) {
    df$x <- df$x + offset[1]; df$y <- df$y + offset[2]; df$z <- df$z + offset[3]
  }
  attr(df, "properties") <- setNames(ptypes, pnames)
  df
}

# Inverse converter: binary PLY -> LAZ (XYZ + recognized fields). Minimal.
ply_to_laz <- function(ply_path, laz_path, offset = c(0, 0, 0)) {
  p <- read_ply(ply_path, offset = offset)
  out <- data.frame(X = p$x, Y = p$y, Z = p$z)
  if (!is.null(p$Intensity)) out$Intensity <- as.integer(round(p$Intensity))
  else if (!is.null(p$intensity)) out$Intensity <- as.integer(round(p$intensity))
  for (nm in c("ReturnNumber", "NumberOfReturns", "Classification"))
    if (!is.null(p[[nm]])) out[[nm]] <- as.integer(round(p[[nm]]))
  if (!is.null(p$gpstime)) out$gpstime <- as.numeric(p$gpstime)
  las <- lidR::LAS(out)
  recognized <- c("x", "y", "z", "Intensity", "intensity", "ReturnNumber",
                  "NumberOfReturns", "Classification", "gpstime")
  for (nm in setdiff(names(p), recognized))
    las <- lidR::add_lasattribute(las, p[[nm]], nm, nm)
  lidR::writeLAS(las, laz_path)
  invisible(laz_path)
}

# Instance-labeled PLY output -> apex det(x,y,z). Strict: a missing id_field is a
# schema failure (-> NULL, the runner skips the cell), never a fake 0-row.
read_instances_ply <- function(ply_path, id_field = "treeID",
                               offset = c(0, 0, 0),
                               x = "x", y = "y", z = "z") {
  p <- read_ply(ply_path, offset = offset)
  if (!id_field %in% names(p)) return(NULL)
  instances_to_det(p, id_field = id_field, x = x, y = y, z = z)
}
