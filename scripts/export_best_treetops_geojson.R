#!/usr/bin/env Rscript
.bs_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.bs_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
bs <- Find(file.exists, c(
  if (!is.null(.bs_ofile) && length(.bs_ofile) && nzchar(.bs_ofile))
    file.path(dirname(.bs_ofile), "bootstrap.R"),
  if (length(.bs_file)) file.path(dirname(sub("^--file=", "", .bs_file[1])),
                                  "bootstrap.R"),
  file.path("scripts", "bootstrap.R"),
  file.path("..", "..", "scripts", "bootstrap.R"),
  file.path(getwd(), "scripts", "bootstrap.R")))
if (!length(bs)) stop("bootstrap.R not found", call. = FALSE)
source(bs[1]); rm(bs, .bs_ofile, .bs_file)

# Export one WGS84 GeoJSON point layer per tested treetop detector, using the
# best-scoring tested configuration for each site. "Best" means pooled F1 from
# the cached benchmark metric CSVs. For CHM-VWF we use the richer sweep grid
# (rung x chm_res x vwf_a); for the other arms we use each detector's tested
# rung(s). Detections are regenerated from the cached frozen clips and cached as
# per-cell CSVs so GPU-backed arms can resume.
#
# Usage:
#   Rscript scripts/export_best_treetops_geojson.R \
#     [SITES=SJER,SOAP,TEAK] [METHODS=ALL|chm_vwf,multichm,...] \
#     [OUT=work/neon/best_treetops_geojson] [EXTENT=core|recall|all] \
#     [FORCE=0] [SKIP_GPU=0]
#
# Writes:
#   OUT/<detector>.geojson
#   OUT/best_treetop_selection.csv
suppressMessages({
  library(data.table)
  library(sf)
  library(lidR)
  library(lasR)
  library(lidRplugins)
  library(crownsegmentr)
  library(terra)
})
options(lidR.progress = FALSE)

d <- .job_dir()
source(.find("sweep_lib.R"))
source(.find("model_bench_lib.R"))
source(.find("pc_detect_lib.R"))
source(.find("model_runner.R"))
source(.find("io_bridge.R"))

## ---- args ----------------------------------------------------------------
args <- strsplit(commandArgs(TRUE), "=")
A <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
OUT <- if (is.null(A$OUT)) file.path(d, "neon", "best_treetops_geojson") else A$OUT
EXTENT <- if (is.null(A$EXTENT)) "core" else A$EXTENT
if (!EXTENT %in% c("core", "recall", "all"))
  stop("EXTENT must be one of core, recall, all", call. = FALSE)
FORCE <- !is.null(A$FORCE) && A$FORCE != "0"
SKIP_GPU <- !is.null(A$SKIP_GPU) && A$SKIP_GPU != "0"
TOL <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
A_VWF <- as.numeric(if (is.null(A$A)) 0.10 else A$A)
TREEISONET_CONF <- if (is.null(A$TREEISONET_CONF)) "0.22" else A$TREEISONET_CONF
TREEISONET_VOXEL <- if (is.null(A$TREEISONET_VOXEL)) "0.8,0.8,2.0" else A$TREEISONET_VOXEL
SAT_IMAGE <- if (is.null(A$SAT_IMAGE)) "sat-sm120-test" else A$SAT_IMAGE
SAT_BATCH <- is.null(A$SAT_BATCH) || A$SAT_BATCH != "0"
SAT_CORES <- max(1L, as.integer(if (is.null(A$SAT_CORES)) 1L else A$SAT_CORES))
FF3D_IMAGE <- if (is.null(A$FF3D_IMAGE)) "ff3d-sm120" else A$FF3D_IMAGE
TIMEOUT <- as.numeric(if (is.null(A$TIMEOUT)) 1800 else A$TIMEOUT)
FF3D_TIMEOUT <- as.numeric(if (is.null(A$FF3D_TIMEOUT)) 3600 else A$FF3D_TIMEOUT)
FF3D_SPACING <- as.numeric(if (is.null(A$FF3D_SPACING)) 24 else A$FF3D_SPACING)
FF3D_MERGE_TOL <- as.numeric(if (is.null(A$FF3D_MERGE_TOL)) 2.0 else A$FF3D_MERGE_TOL)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

METHOD_REGISTRY <- data.frame(
  method = c("chm_vwf", "lmfauto", "multichm", "ptrees", "ams3d",
             "li2012", "lidr_lmf_pc", "lidr_li2012", "lasr_lmax_pc",
             "treeisonet", "segmentanytree", "forestformer3d"),
  family = c("chm_sweep", "lidrplugins", "lidrplugins", "lidrplugins", "ams3d",
             "li2012", "pc_ladder", "pc_ladder", "pc_ladder",
             "treeisonet", "segmentanytree", "forestformer3d"),
  stringsAsFactors = FALSE
)
METHODS <- if (is.null(A$METHODS) || A$METHODS == "ALL")
  METHOD_REGISTRY$method else strsplit(A$METHODS, ",")[[1]]
bad_methods <- setdiff(METHODS, METHOD_REGISTRY$method)
if (length(bad_methods))
  stop("unknown METHODS: ", paste(bad_methods, collapse = ", "), call. = FALSE)

## ---- helpers -------------------------------------------------------------
zone_epsg <- function(z) {
  num <- as.integer(sub("[NnSs]$", "", z))
  hemi <- toupper(sub("^[0-9]+", "", z))
  (if (hemi == "S") 32700L else 32600L) + num
}

sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
empty_det <- function() data.frame(x = numeric(), y = numeric(), z = numeric())

tops_to_det <- function(tt) {
  if (is.null(tt) || nrow(tt) == 0) return(empty_det())
  co <- sf::st_coordinates(tt)
  z <- if (ncol(co) >= 3) co[, 3] else if (!is.null(tt$Z)) tt$Z else rep(0, nrow(co))
  det <- data.frame(x = co[, 1], y = co[, 2], z = z)
  assert_detection_contract(det)
  det
}

det_lmfauto_export <- function(las, hmin = 2) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::lmfauto(plot = TRUE, hmin = hmin)),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  tops_to_det(tt)
}

det_multichm_export <- function(las, res = 0.5, a = 0.10) {
  tt <- tryCatch(lidR::locate_trees(las, lidRplugins::multichm(res = res,
                                                               ws = ws_factory(a))),
                 error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  tops_to_det(tt)
}

det_ptrees_export <- function(las, hmin = 2, k = c(30, 15)) {
  if (sum(las$Z >= hmin) < min(k)) return(empty_det())
  seg <- tryCatch(lidR::segment_trees(las, lidRplugins::ptrees(k = k, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg) || !"treeID" %in% names(seg@data)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

det_ams3d_export <- function(las, cd_ratio = 0.4, cl_ratio = 0.8, min_above = 2) {
  if (lidR::is.empty(las) || lidR::npoints(las) < 5) return(empty_det())
  seg <- tryCatch(
    crownsegmentr::segment_tree_crowns(
      point_cloud = las,
      crown_diameter_to_tree_height = cd_ratio,
      crown_length_to_tree_height = cl_ratio,
      segment_crowns_only_above = min_above,
      ground_height = NULL,
      crown_id_column_name = "crown_id"),
    error = function(e) NULL)
  if (is.null(seg)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "crown_id", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

det_li2012_export <- function(las, hmin = 2, dt1 = 1.5, dt2 = 2, R = 2) {
  if (sum(las$Z >= hmin) < 1) return(empty_det())
  seg <- tryCatch(lidR::segment_trees(las,
                    lidR::li2012(dt1 = dt1, dt2 = dt2, R = R, hmin = hmin)),
                  error = function(e) NULL)
  if (is.null(seg) || !"treeID" %in% names(seg@data)) return(NULL)
  det <- reduce_instances(seg@data, id_col = "treeID", x = "X", y = "Y", z = "Z")
  assert_detection_contract(det)
  det
}

metric_file_for <- function(site, method) {
  nd <- file.path(d, "neon", site)
  family <- METHOD_REGISTRY$family[match(method, METHOD_REGISTRY$method)]
  switch(family,
    chm_sweep = file.path(nd, "sweep_results.csv"),
    lidrplugins = file.path(nd, "lidrplugins_results.csv"),
    ams3d = file.path(nd, "ams3d_results.csv"),
    li2012 = file.path(nd, "li2012_results.csv"),
    pc_ladder = file.path(nd, "pc_detect_ladder_results.csv"),
    treeisonet = file.path(nd, "treeisonet_results.csv"),
    segmentanytree = file.path(nd, "segmentanytree_results.csv"),
    forestformer3d = file.path(nd, "forestformer3d_results.csv"),
    stop("unhandled method family: ", family)
  )
}

metric_rows <- function(site, method) {
  f <- metric_file_for(site, method)
  if (!file.exists(f)) return(NULL)
  r <- read.csv(f, stringsAsFactors = FALSE)
  if (!nrow(r)) return(NULL)
  family <- METHOD_REGISTRY$family[match(method, METHOD_REGISTRY$method)]
  if (family == "chm_sweep") {
    r$site <- site
    r$detector <- method
    r$source_family <- family
  } else {
    if (!"detector" %in% names(r)) return(NULL)
    r <- r[r$detector == method, , drop = FALSE]
    if (!nrow(r)) return(NULL)
    if (!"site" %in% names(r)) r$site <- site
    r$source_family <- family
  }
  r$source_file <- f
  if (!"tp_core" %in% names(r)) r$tp_core <- round(r$precision * r$n_det)
  r$rung <- as.character(r$rung)
  r
}

group_key <- function(df, cols) {
  do.call(paste, c(lapply(cols, function(nm) df[[nm]]), sep = "\r"))
}

best_selection <- function(site, method) {
  r <- metric_rows(site, method)
  if (is.null(r) || !nrow(r)) return(NULL)
  group_cols <- if (method == "chm_vwf") c("rung", "chm_res", "vwf_a") else "rung"
  k <- group_key(r, group_cols)
  pooled <- do.call(rbind, lapply(split(r, k, drop = TRUE), function(s) {
    p <- pool(s)
    vals <- s[1, group_cols, drop = FALSE]
    cbind(data.frame(site = site, method = method, detector = method,
                     source_family = s$source_family[1],
                     source_file = s$source_file[1],
                     stringsAsFactors = FALSE),
          vals, p)
  }))
  if (is.null(pooled) || !nrow(pooled)) return(NULL)
  for (nm in c("chm_res", "vwf_a")) if (!nm %in% names(pooled)) pooled[[nm]] <- NA_real_
  pooled <- pooled[order(-pooled$F1, -pooled$recall, -pooled$precision,
                         pooled$rung), , drop = FALSE]
  pooled[1, , drop = FALSE]
}

rows_for_selection <- function(sel) {
  r <- metric_rows(sel$site, sel$method)
  if (is.null(r) || !nrow(r)) return(NULL)
  keep <- r$rung == as.character(sel$rung)
  if (sel$method == "chm_vwf") {
    keep <- keep & r$chm_res == sel$chm_res & r$vwf_a == sel$vwf_a
  }
  r[keep, , drop = FALSE]
}

.site_cache <- new.env(parent = emptyenv())
site_context <- function(site) {
  if (exists(site, .site_cache, inherits = FALSE)) return(get(site, .site_cache))
  nd <- file.path(d, "neon", site)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(laz)) stop("no LAZ files under ", file.path(nd, "lidar"), call. = FALSE)
  ctg <- lidR::readLAScatalog(laz, progress = FALSE)
  ctx <- list(site = site, nd = nd, pc = pc, ctg = ctg,
              epsg = zone_epsg(pc$utmZone[1]))
  assign(site, ctx, .site_cache)
  ctx
}

rung_value <- function(rung) {
  if (length(rung) != 1 || is.na(rung) || rung == "native") NA_real_ else as.numeric(rung)
}

cell_cache_path <- function(ctx, method, row) {
  rung <- as.character(row$rung)
  parts <- c(method, ctx$site, row$plot, rung)
  if (method == "chm_vwf")
    parts <- c(parts, sprintf("res%s", row$chm_res), sprintf("a%s", row$vwf_a))
  if (method == "treeisonet")
    parts <- c(parts, sprintf("conf%s", TREEISONET_CONF), sprintf("voxel%s", TREEISONET_VOXEL))
  if (method == "segmentanytree")
    parts <- c(parts, sprintf("image%s", SAT_IMAGE))
  if (method == "forestformer3d")
    parts <- c(parts, sprintf("image%s", FF3D_IMAGE), sprintf("spacing%s", FF3D_SPACING),
               sprintf("merge%s", FF3D_MERGE_TOL))
  dir.create(file.path(ctx$nd, "best_treetop_cache"), recursive = TRUE, showWarnings = FALSE)
  file.path(ctx$nd, "best_treetop_cache", paste0(sanitize(paste(parts, collapse = "__")), ".csv"))
}

read_det_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  det <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(det) || !identical(names(det), c("x", "y", "z"))) return(NULL)
  det <- data.frame(x = as.numeric(det$x), y = as.numeric(det$y), z = as.numeric(det$z))
  assert_detection_contract(det)
  det
}

write_det_cache <- function(det, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(det[, c("x", "y", "z"), drop = FALSE], path, row.names = FALSE)
}

ensure_gpu_allowed <- function(method) {
  if (!SKIP_GPU) return(TRUE)
  if (method %in% c("treeisonet", "segmentanytree", "forestformer3d")) return(FALSE)
  TRUE
}

treeisonet_det <- function(prep, pid, rung) {
  venv <- file.path(.ROOT, "gpu/.venv/bin/python")
  drv <- file.path(.ROOT, "gpu/run_treeisonet.py")
  loc <- file.path(.ROOT, "gpu/store/treeaibox/als_treeloc.pth")
  cfg <- Sys.glob(file.path(.ROOT, "gpu/store/treeaibox/*reclamation*treeloc*.json"))[1]
  if (!file.exists(venv) || !file.exists(drv) || !file.exists(loc) || is.na(cfg))
    stop("TreeisoNet assets are missing", call. = FALSE)
  ocsv <- file.path(tempdir(), sprintf("export_ti_%s_%s.csv", pid, rung))
  run_python_arm(venv, drv, prep$normalized, ocsv,
                 extra = c(loc, cfg, TREEISONET_VOXEL, TREEISONET_CONF),
                 timeout = TIMEOUT, label = sprintf("%s/%s", pid, rung))
}

persisted_segmentanytree_det <- function(ctx, prep, pid, rung) {
  stem <- sprintf("%s_%s", pid, rung)
  candidates <- c(
    file.path(ctx$nd, "segmentanytree_instances", paste0(stem, ".laz")),
    file.path(ctx$nd, "segmentanytree_instances", paste0(stem, ".las")),
    file.path(ctx$nd, "sat_batch_salvage", pid, paste0("sat_", stem, "_out.laz")),
    file.path(ctx$nd, "sat_batch_salvage", pid, paste0("sat_", stem, "_out.las"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) return(NULL)
  det_abs <- read_instances_laz(hit[1], id_field = "PredInstance")
  if (is.null(det_abs)) return(NULL)
  det_to_agl(det_abs, prep$dtm)
}

.find_las_export <- function(dir, stem) {
  direct <- file.path(dir, paste0(stem, c(".las", ".laz")))
  hit <- direct[file.exists(direct)]
  if (length(hit)) return(hit[1])
  esc <- gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", stem)
  hit <- list.files(dir, pattern = paste0("^", esc, "(_out)?\\.la[sz]$"),
                    full.names = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

run_sat_batch_export <- function(image, input_dir, output_dir, timeout, label = NULL) {
  batch_driver <- file.path(.ROOT, "gpu", "run_segmentanytree_batch.py")
  if (!file.exists(batch_driver)) stop("SegmentAnyTree batch driver missing", call. = FALSE)
  if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  in_abs <- normalizePath(input_dir, mustWork = TRUE)
  out_abs <- normalizePath(output_dir, mustWork = FALSE)
  mdirs <- unique(c(in_abs, out_abs, normalizePath(dirname(batch_driver), mustWork = TRUE)))
  vol <- as.vector(rbind("-v", paste0(mdirs, ":", mdirs)))
  args <- c("run", "--rm", "--gpus", "all", "--shm-size=8g", "--ipc=host",
            vol, image, "python3", batch_driver, in_abs, out_abs)
  out <- tryCatch(suppressWarnings(system2("docker", shQuote(args), stdout = TRUE,
                                           stderr = TRUE, timeout = timeout)),
                  error = function(e) NULL)
  st <- if (is.null(out)) 1L else attr(out, "status")
  if (!is.null(st) && st != 0) {
    if (!is.null(label)) message(sprintf("  %s: SAT batch failed (status %s)",
                                         label, st))
    if (length(out)) message(paste(tail(out, 8), collapse = "\n"))
    return(FALSE)
  }
  TRUE
}

segmentanytree_det <- function(ctx, prep, pid, rung) {
  persisted <- persisted_segmentanytree_det(ctx, prep, pid, rung)
  if (!is.null(persisted)) return(persisted)
  driver <- file.path(.ROOT, "gpu", "run_segmentanytree.py")
  if (!file.exists(driver)) stop("SegmentAnyTree driver missing", call. = FALSE)
  ply <- file.path(tempdir(), sprintf("export_sat_%s_%s.ply", pid, rung))
  olas <- file.path(tempdir(), sprintf("export_sat_%s_%s.las", pid, rung))
  laz_to_ply(prep$rawground, ply)
  det_abs <- run_docker_arm(SAT_IMAGE, ply, olas,
    cmd = c("python3", driver),
    mounts = dirname(driver),
    extra_docker = c("--shm-size=8g", "--ipc=host"),
    reader = function(p) read_instances_laz(p, id_field = "PredInstance"),
    gpus = "all", timeout = TIMEOUT,
    label = sprintf("%s/%s", pid, rung))
  if (is.null(det_abs)) return(NULL)
  det_to_agl(det_abs, prep$dtm)
}

ff3d_cyl_centers <- function(cx, cy, ph, spacing) {
  k <- max(1L, ceiling((2 * ph) / spacing) + 1L)
  off <- seq(-ph, ph, length.out = k)
  g <- expand.grid(dx = off, dy = off)
  data.frame(cx = cx + g$dx, cy = cy + g$dy)
}

forestformer3d_det <- function(prep, pid, rung, cx, cy, ph) {
  repo <- file.path(.ROOT, "gpu/store/forestformer3d/ForestFormer3D")
  ckpt <- file.path(repo, "work_dirs/clean_forestformer/epoch_3000_fix.pth")
  entry <- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_entry.sh")
  patch <- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_repo.patch")
  driver <- file.path(.ROOT, "gpu/forestformer3d-sm120/ff3d_arm.py")
  if (!file.exists(repo) || !file.exists(ckpt) || !file.exists(entry) || !file.exists(driver))
    stop("ForestFormer3D assets are missing", call. = FALSE)
  raw <- lidR::readLAS(prep$rawground)
  if (is.null(raw) || lidR::is.empty(raw)) return(empty_det())
  in_dir <- file.path(tempdir(), sprintf("export_ff3d_%s_%s", pid, rung))
  unlink(in_dir, recursive = TRUE, force = TRUE)
  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
  cc <- ff3d_cyl_centers(cx, cy, ph, FF3D_SPACING)
  n_cyl <- 0L
  for (i in seq_len(nrow(cc))) {
    cyl <- lidR::clip_circle(raw, cc$cx[i], cc$cy[i], 16)
    if (lidR::is.empty(cyl) || lidR::npoints(cyl) < 50) next
    lidR::writeLAS(cyl, file.path(in_dir, sprintf("cyl_%03d.laz", n_cyl)))
    n_cyl <- n_cyl + 1L
  }
  if (n_cyl == 0L) return(empty_det())
  out_laz <- file.path(tempdir(), sprintf("export_ff3d_%s_%s.laz", pid, rung))
  det_abs <- run_docker_arm(FF3D_IMAGE, in_dir, out_laz,
    cmd = c("bash", entry),
    extra = c(ckpt, repo, patch, driver),
    mounts = c(repo, dirname(ckpt), dirname(entry)),
    reader = function(p) ff3d_collapse(p, merge_tol = FF3D_MERGE_TOL),
    gpus = "all", timeout = FF3D_TIMEOUT,
    label = sprintf("%s/%s", pid, rung))
  if (is.null(det_abs)) return(NULL)
  agl_guard(det_abs, prep$dtm)
}

generate_cell_det <- function(ctx, row) {
  method <- as.character(row$detector)
  if (!ensure_gpu_allowed(method)) {
    message(sprintf("  %s/%s/%s skipped (SKIP_GPU=1)", method, row$plot, row$rung))
    return(NULL)
  }
  cache <- cell_cache_path(ctx, method, row)
  if (!FORCE) {
    cached <- read_det_cache(cache)
    if (!is.null(cached)) return(cached)
  }

  ci <- ctx$pc[ctx$pc$plotID == row$plot, ][1, ]
  cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  rung_lbl <- as.character(row$rung)
  prep <- frozen_clip(ctx$ctg, ctx$site, row$plot, rung_value(rung_lbl), cx, cy, ph,
                      out_root = file.path(ctx$nd, "frozen"))
  if (is.null(prep)) return(NULL)
  det <- switch(method,
    chm_vwf = detect_lasr(prep$normalized, as.numeric(row$chm_res),
                          as.numeric(row$vwf_a), prep$frdens),
    lmfauto = {
      las <- lidR::readLAS(prep$normalized)
      det_lmfauto_export(las, hmin = 2)
    },
    multichm = {
      las <- lidR::readLAS(prep$normalized)
      res <- if ("chm_res" %in% names(row) && is.finite(as.numeric(row$chm_res)))
        as.numeric(row$chm_res) else if (prep$frdens >= 8) 0.25 else 0.5
      det_multichm_export(las, res = res, a = A_VWF)
    },
    ptrees = {
      las <- lidR::readLAS(prep$normalized)
      det_ptrees_export(las, hmin = 2)
    },
    ams3d = {
      las <- lidR::readLAS(prep$normalized)
      det_ams3d_export(las)
    },
    li2012 = {
      las <- lidR::readLAS(prep$normalized)
      det_li2012_export(las)
    },
    lidr_lmf_pc = {
      las <- lidR::readLAS(prep$normalized)
      det_lidr_lmf_pc(las, ws_factory(A_VWF))
    },
    lidr_li2012 = {
      las <- lidR::readLAS(prep$normalized)
      det_lidr_li2012(las)
    },
    lasr_lmax_pc = det_lasr_lmax_pc(prep$normalized, ws_factory(A_VWF)),
    treeisonet = treeisonet_det(prep, row$plot, rung_lbl),
    segmentanytree = segmentanytree_det(ctx, prep, row$plot, rung_lbl),
    forestformer3d = forestformer3d_det(prep, row$plot, rung_lbl, cx, cy, ph),
    stop("unhandled method: ", method)
  )
  if (is.null(det)) return(NULL)
  assert_detection_contract(det)
  det <- det[is.finite(det$x) & is.finite(det$y) & is.finite(det$z), , drop = FALSE]
  rownames(det) <- NULL
  write_det_cache(det, cache)
  det
}

filter_extent <- function(det, cx, cy, ph) {
  in_core <- abs(det$x - cx) <= ph & abs(det$y - cy) <= ph
  in_recall <- abs(det$x - cx) <= ph + TOL & abs(det$y - cy) <= ph + TOL
  det$in_core <- in_core
  det$in_recall_region <- in_recall
  keep <- switch(EXTENT, core = in_core, recall = in_recall, all = rep(TRUE, nrow(det)))
  det[keep, , drop = FALSE]
}

row_value <- function(row, name, default = NA_real_) {
  if (name %in% names(row)) row[[name]][1] else default
}

det_to_layer <- function(ctx, row, sel, det) {
  if (is.null(det)) return(NULL)
  ci <- ctx$pc[ctx$pc$plotID == row$plot, ][1, ]
  cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
  det <- filter_extent(det, cx, cy, ph)
  if (!nrow(det)) return(NULL)
  props <- data.frame(
    prediction_id = sprintf("%s_%s_%s_%s_%04d", row$detector, ctx$site,
                            row$plot, row$rung, seq_len(nrow(det))),
    detector = as.character(row$detector),
    site = ctx$site,
    plot = as.character(row$plot),
    plotType = as.character(row$plotType),
    selected_rung = as.character(row$rung),
    selected_site_F1 = round(sel$F1, 6),
    selected_site_recall = round(sel$recall, 6),
    selected_site_precision = round(sel$precision, 6),
    selected_site_n_plots = as.integer(sel$n_plots),
    selected_site_n_ref = as.integer(sel$n_ref),
    source_family = as.character(sel$source_family),
    source_file = basename(as.character(sel$source_file)),
    chm_res = as.numeric(row_value(row, "chm_res")),
    vwf_a = as.numeric(row_value(row, "vwf_a")),
    metric_pdens = as.numeric(row_value(row, "pdens")),
    metric_frdens = as.numeric(row_value(row, "frdens")),
    x_utm = det$x,
    y_utm = det$y,
    z_agl_m = det$z,
    in_core = det$in_core,
    in_recall_region = det$in_recall_region,
    stringsAsFactors = FALSE
  )
  sf::st_transform(
    sf::st_as_sf(props, coords = c("x_utm", "y_utm", "z_agl_m"),
                 crs = ctx$epsg, dim = "XYZ", remove = FALSE),
    4326)
}

cell_features <- function(ctx, row, sel) {
  det <- generate_cell_det(ctx, row)
  det_to_layer(ctx, row, sel, det)
}

segmentanytree_batch_features <- function(ctx, rows, sel) {
  layers <- list()
  cells <- list()
  batch_in <- file.path(tempdir(), sprintf("export_sat_batch_in_%s", ctx$site))
  batch_out <- file.path(tempdir(), sprintf("export_sat_batch_out_%s", ctx$site))
  unlink(c(batch_in, batch_out), recursive = TRUE, force = TRUE)
  dir.create(batch_in, recursive = TRUE, showWarnings = FALSE)

  for (j in seq_len(nrow(rows))) {
    row <- rows[j, , drop = FALSE]
    cache <- cell_cache_path(ctx, "segmentanytree", row)
    det <- if (!FORCE) read_det_cache(cache) else NULL
    ci <- ctx$pc[ctx$pc$plotID == row$plot, ][1, ]
    cx <- ci$easting; cy <- ci$northing; ph <- plot_half(ci$plotType)
    rung_lbl <- as.character(row$rung)
    prep <- frozen_clip(ctx$ctg, ctx$site, row$plot, rung_value(rung_lbl), cx, cy, ph,
                        out_root = file.path(ctx$nd, "frozen"))
    if (is.null(prep)) next

    if (is.null(det)) {
      det <- persisted_segmentanytree_det(ctx, prep, row$plot, rung_lbl)
      if (!is.null(det)) {
        det <- det[is.finite(det$x) & is.finite(det$y) & is.finite(det$z), , drop = FALSE]
        rownames(det) <- NULL
        write_det_cache(det, cache)
      }
    }
    if (!is.null(det)) {
      lyr <- det_to_layer(ctx, row, sel, det)
      if (!is.null(lyr) && nrow(lyr)) layers[[length(layers) + 1L]] <- lyr
      next
    }

    stem <- sprintf("sat_%s_%s", row$plot, rung_lbl)
    ply <- file.path(batch_in, paste0(stem, ".ply"))
    laz_to_ply(prep$rawground, ply)
    cells[[length(cells) + 1L]] <- list(row = row, cache = cache, prep = prep, stem = stem)
  }

  if (length(cells)) {
    ok <- FALSE
    if (SAT_BATCH) {
      ok <- run_sat_batch_export(SAT_IMAGE, batch_in, batch_out,
                                 timeout = TIMEOUT * length(cells),
                                 label = sprintf("%s segmentanytree", ctx$site))
    }
    if (ok) {
      for (cell in cells) {
        las <- .find_las_export(batch_out, cell$stem)
        if (is.na(las)) {
          message(sprintf("    %s/%s: SAT batch output missing",
                          cell$row$plot, cell$row$rung))
          next
        }
        det_abs <- read_instances_laz(las, id_field = "PredInstance")
        if (is.null(det_abs)) next
        det <- det_to_agl(det_abs, cell$prep$dtm)
        det <- det[is.finite(det$x) & is.finite(det$y) & is.finite(det$z), , drop = FALSE]
        rownames(det) <- NULL
        write_det_cache(det, cell$cache)
        lyr <- det_to_layer(ctx, cell$row, sel, det)
        if (!is.null(lyr) && nrow(lyr)) layers[[length(layers) + 1L]] <- lyr
      }
    } else {
      message(sprintf("    %s: running single-cell SAT for %d cell(s)",
                      ctx$site, length(cells)))
      run_single <- function(cell) {
        det <- tryCatch(segmentanytree_det(ctx, cell$prep, cell$row$plot,
                                           as.character(cell$row$rung)),
                        error = function(e) {
                          message(sprintf("    %s/%s single-cell SAT failed: %s",
                                          cell$row$plot, cell$row$rung,
                                          conditionMessage(e)))
                          NULL
                        })
        if (is.null(det)) return(NULL)
        det <- det[is.finite(det$x) & is.finite(det$y) & is.finite(det$z), , drop = FALSE]
        rownames(det) <- NULL
        write_det_cache(det, cell$cache)
        det_to_layer(ctx, cell$row, sel, det)
      }
      one_layers <- if (SAT_CORES > 1L && length(cells) > 1L)
        parallel::mclapply(cells, run_single,
                           mc.cores = min(SAT_CORES, length(cells)),
                           mc.preschedule = FALSE) else lapply(cells, run_single)
      for (lyr in one_layers)
        if (!is.null(lyr) && nrow(lyr)) layers[[length(layers) + 1L]] <- lyr
    }
  }
  if (length(layers)) do.call(rbind, layers) else NULL
}

write_geojson <- function(x, path) {
  if (is.null(x) || !nrow(x)) {
    writeLines('{"type":"FeatureCollection","features":[]}', path)
    return(invisible(0L))
  }
  sf::st_write(x, path, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE,
               layer_options = c("RFC7946=YES", "COORDINATE_PRECISION=7"))
  invisible(nrow(x))
}

## ---- main ----------------------------------------------------------------
run_main <- function() {
  selections <- do.call(rbind, Filter(Negate(is.null), lapply(METHODS, function(m)
    do.call(rbind, Filter(Negate(is.null), lapply(SITES, best_selection, method = m))))))
  if (is.null(selections) || !nrow(selections))
    stop("no method/site selections found; check cached benchmark CSVs", call. = FALSE)
  write.csv(selections, file.path(OUT, "best_treetop_selection.csv"), row.names = FALSE)
  cat(sprintf("Selected %d method/site configurations -> %s\n",
              nrow(selections), file.path(OUT, "best_treetop_selection.csv")))

  manifest <- selections
  manifest$features_written <- 0L
  manifest$status <- "selected"
  for (method in METHODS) {
    sel_m <- selections[selections$method == method, , drop = FALSE]
    if (!nrow(sel_m)) {
      out_path <- file.path(OUT, paste0(method, ".geojson"))
      write_geojson(NULL, out_path)
      cat(sprintf("  %-16s %6d features -> %s (no selection)\n", method, 0L, out_path))
      next
    }
    layers <- list()
    total <- 0L
    cat(sprintf("[%s]\n", method))
    for (i in seq_len(nrow(sel_m))) {
      sel <- sel_m[i, , drop = FALSE]
      ctx <- site_context(sel$site)
      rows <- rows_for_selection(sel)
      if (is.null(rows) || !nrow(rows)) next
      rows$detector <- method
      cat(sprintf("  %s best rung=%s F1=%.3f plots=%d\n",
                  sel$site, sel$rung, sel$F1, length(unique(rows$plot))))
      if (method == "segmentanytree") {
        site_sf <- tryCatch(segmentanytree_batch_features(ctx, rows, sel),
                            error = function(e) {
                              message(sprintf("    %s segmentanytree batch failed: %s",
                                              sel$site, conditionMessage(e)))
                              NULL
                            })
      } else {
        site_layers <- list()
        for (j in seq_len(nrow(rows))) {
          lyr <- tryCatch(cell_features(ctx, rows[j, , drop = FALSE], sel),
                          error = function(e) {
                            message(sprintf("    %s/%s failed: %s",
                                            rows$plot[j], rows$rung[j],
                                            conditionMessage(e)))
                            NULL
                          })
          if (!is.null(lyr) && nrow(lyr)) site_layers[[length(site_layers) + 1L]] <- lyr
        }
        site_sf <- if (length(site_layers)) do.call(rbind, site_layers) else NULL
      }
      if (!is.null(site_sf) && nrow(site_sf)) layers[[length(layers) + 1L]] <- site_sf
      nsite <- if (!is.null(site_sf)) nrow(site_sf) else 0L
      total <- total + nsite
      idx <- manifest$method == method & manifest$site == sel$site
      manifest$features_written[idx] <- nsite
      manifest$status[idx] <- if (nsite > 0L) "written" else "empty_or_failed"
      cat(sprintf("    %s features: %d\n", sel$site, nsite))
    }
    out_path <- file.path(OUT, paste0(method, ".geojson"))
    n <- write_geojson(if (length(layers)) do.call(rbind, layers) else NULL, out_path)
    cat(sprintf("  %-16s %6d features -> %s\n", method, n, out_path))
  }
  write.csv(manifest, file.path(OUT, "best_treetop_export_manifest.csv"), row.names = FALSE)
  cat(sprintf("Manifest -> %s\n", file.path(OUT, "best_treetop_export_manifest.csv")))
}

if (sys.nframe() == 0L) run_main()
