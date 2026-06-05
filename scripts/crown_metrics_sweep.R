#!/usr/bin/env Rscript
# Crown-segmentation benchmark for GitHub issue #7.
#
# The density-ladder sweep scores DETECTION (treetops) only. This script closes
# the loop on CROWN DELINEATION: it seeds five crown segmenters from the SAME
# detected treetops, computes per-tree crown area -> diameter, and scores the
# diameter against NEON field crown diameter (maxCrownDiameter /
# ninetyCrownDiameter), pooled by crown class.
#
# Segmenters (all on the per-plot pit-free CHM, native density):
#   - lidR dalponte2016   (seeded region growing)
#   - lidR silva2016      (seeded Voronoi-like)
#   - lidR watershed      (marker-FREE; crowns matched to stems by containment)
#   - lasR region_growing (seeded from local_maximum_raster, mirrors segment_lasr.R)
#   - random_walker       (Grady 2006, seeded; sparse Dirichlet solve via Matrix)
#
# Two diameter estimates per crown, with a geometric caveat:
#   - d_eq      = 2*sqrt(area/pi)          equivalent-circle  -> ninetyCrownDiameter
#   - d_caliper = max pairwise vertex dist polygon max-axis   -> maxCrownDiameter
#
# Usage:
#   Rscript scripts/crown_metrics_sweep.R SITES=SJER,SOAP,TEAK CORES=4 \
#           TOL=4 RES=0.5 A=0.10
#
# Reads (read-only): work/neon/<SITE>/{ground_truth_stems.csv,plot_centroids.csv}
#   the LiDAR catalog work/neon/<SITE>/lidar/, and the cached field-crown widths
#   from work/neon/<SITE>/vst/<site>_vst_allyears.rds.
# Writes (NEW file, never overwrites a shared input):
#   work/neon/<SITE>/crown_metrics_results.csv  (one row per matched tree)
suppressMessages({
  library(lidR); library(lasR); library(terra); library(sf)
  library(data.table); library(Matrix); library(parallel)
})
options(lidR.progress = FALSE, lidR.verbose = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
source(file.path("scripts", "sweep_lib.R"))

## ---- args (KEY=VALUE positional, per repo convention) --------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
CORES <- as.integer(if (is.null(A$CORES)) 4 else A$CORES)
TOL   <- as.numeric(if (is.null(A$TOL))  4    else A$TOL)
RES   <- as.numeric(if (is.null(A$RES))  0.5  else A$RES)
VWF_A <- as.numeric(if (is.null(A$A))    0.10 else A$A)
MINTREES <- 6
RW_TIMEOUT <- 120          # seconds; random-walker solve is timeboxed per plot

## ---- random walker (Grady 2006) on a CHM raster --------------------------
# Builds the 4-neighbour pixel graph over canopy pixels (Z>=hmin), Gaussian edge
# weights w_ij = exp(-beta*(chm_i-chm_j)^2) on the normalized CHM, assigns each
# seed a marker label, and solves the combinatorial Dirichlet problem
# Lu x = -B^T m for marker probabilities; each pixel is argmax-labelled.
# eps regularizes the Laplacian so seed-less canopy islands stay solvable.
# Returns a SpatRaster of crown labels (1..L) aligned to `chm`, with the seed
# order preserved (label k = k-th seed cell), or NULL if it fails.
random_walker <- function(chm, seeds_xy, beta = 1.0, hmin = 2, eps = 1e-6) {
  nr <- nrow(chm); nc <- ncol(chm)
  vals <- terra::values(chm, mat = FALSE)
  valid <- is.finite(vals) & vals >= hmin
  npix <- length(vals)
  if (!any(valid)) return(NULL)
  rng <- range(vals[valid]); span <- diff(rng); if (span <= 0) span <- 1
  vn <- (vals - rng[1]) / span                  # normalize for stable beta
  cells <- terra::cellFromXY(chm, seeds_xy)
  ok <- !is.na(cells) & valid[cells]
  cells <- cells[ok]; keep_seed <- which(ok)
  if (!length(cells)) return(NULL)
  # one label per distinct seed cell (first seed wins a shared cell)
  seed_lab <- integer(npix); k <- 0L; cell_lab <- integer(length(cells))
  for (m in seq_along(cells)) {
    if (seed_lab[cells[m]] == 0L) { k <- k + 1L; seed_lab[cells[m]] <- k }
    cell_lab[m] <- seed_lab[cells[m]]
  }
  L <- k
  vidx <- which(valid); nv <- length(vidx)
  pos <- integer(npix); pos[vidx] <- seq_len(nv)
  rc <- terra::rowColFromCell(chm, vidx); r <- rc[, 1]; c <- rc[, 2]
  buildEdges <- function(dr, dc) {                # right + down avoids dup edges
    nr2 <- r + dr; nc2 <- c + dc
    inb <- nr2 >= 1 & nr2 <= nr & nc2 >= 1 & nc2 <= nc
    cell_a <- vidx[inb]
    cell_b <- terra::cellFromRowCol(chm, nr2[inb], nc2[inb])
    kp <- valid[cell_b]
    list(a = cell_a[kp], b = cell_b[kp])
  }
  e1 <- buildEdges(0, 1); e2 <- buildEdges(1, 0)
  ea <- c(e1$a, e2$a); eb <- c(e1$b, e2$b)
  ia <- pos[ea]; ib <- pos[eb]
  w  <- exp(-beta * (vn[ea] - vn[eb])^2) + 1e-6
  W  <- sparseMatrix(i = c(ia, ib), j = c(ib, ia), x = c(w, w), dims = c(nv, nv))
  dg <- Matrix::rowSums(W)
  Lap <- Diagonal(x = dg + eps) - W
  seed_local <- pos[which(seed_lab > 0)]
  seed_label_local <- seed_lab[which(seed_lab > 0)]
  is_seed <- logical(nv); is_seed[seed_local] <- TRUE
  uidx <- which(!is_seed); sidx <- which(is_seed)
  if (!length(uidx)) return(NULL)
  Lu <- Lap[uidx, uidx]; B <- Lap[uidx, sidx]
  M  <- sparseMatrix(i = seq_along(sidx), j = seed_label_local,
                     x = 1, dims = c(length(sidx), L))
  rhs <- as(-(B %*% M), "matrix")
  X <- tryCatch(as.matrix(Matrix::solve(Lu, rhs)), error = function(e) NULL)
  if (is.null(X)) return(NULL)
  lab_u <- max.col(X, ties.method = "first")
  out_lab <- integer(npix)
  out_lab[vidx[sidx]] <- seed_label_local
  out_lab[vidx[uidx]] <- lab_u
  r_out <- chm; terra::values(r_out) <- ifelse(out_lab > 0, out_lab, NA_integer_)
  list(raster = r_out, seed_cell = cells, seed_label = cell_lab)
}

## ---- crown geometry from a label raster ----------------------------------
# Polygonize a treeID label raster; per crown return area, d_eq, d_caliper.
# d_caliper = max pairwise distance over the polygon's exterior vertices.
crown_geom <- function(label_r) {
  polys <- tryCatch(terra::as.polygons(label_r, dissolve = TRUE, values = TRUE),
                    error = function(e) NULL)
  if (is.null(polys) || nrow(polys) == 0) return(NULL)
  names(polys)[1] <- "treeID"
  polys <- polys[!is.na(polys$treeID) & polys$treeID > 0, ]
  if (nrow(polys) == 0) return(NULL)
  sfp <- sf::st_as_sf(polys)
  area <- as.numeric(sf::st_area(sfp))
  d_eq <- 2 * sqrt(area / pi)
  d_cal <- vapply(seq_len(nrow(sfp)), function(i) {
    cc <- sf::st_coordinates(sfp[i, ])[, c("X", "Y"), drop = FALSE]
    if (nrow(cc) < 2) return(NA_real_)
    max(dist(cc))                                # max pairwise vertex distance
  }, numeric(1))
  cen <- sf::st_coordinates(sf::st_centroid(sfp))
  data.frame(treeID = sfp$treeID, area = area, d_eq = d_eq, d_caliper = d_cal,
             cx = cen[, 1], cy = cen[, 2], stringsAsFactors = FALSE)
}

## ---- field crown diameter per site (join from cached vst rds) ------------
# Dedup apparentindividual to the nearest-to-2021 measurement per individualID;
# keep maxCrownDiameter / ninetyCrownDiameter. Merged onto ground_truth_stems by
# individualID (the rds is the canonical source; ground_truth_stems.csv predates
# the crown-diameter columns added to neon_ground_truth.R for reproducibility).
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst",
                   paste0(tolower(site), "_vst_allyears.rds"))
  dat <- readRDS(rds)
  ai <- as.data.frame(dat$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]
  ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai1 <- ai[!duplicated(ai$individualID),
            c("individualID", "maxCrownDiameter", "ninetyCrownDiameter")]
  ai1
}

## ---- per-plot crown benchmark --------------------------------------------
run_plot <- function(site, pid, ctg, pc, gt, tmpdir) {
  ci <- pc[pc$plotID == pid, ][1, ]
  cx <- ci$easting; cy <- ci$northing
  ph <- plot_half(ci$plotType)
  stems <- gt[gt$plotID == pid &
              abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
  if (nrow(stems) < 1) return(NULL)

  prep <- tryCatch(prepare_clip(ctg, cx, cy, rung = NA, tmpdir, core_half = ph),
                   error = function(e) NULL)
  if (is.null(prep)) return(NULL)
  on.exit(unlink(prep$file), add = TRUE)
  las <- tryCatch(readLAS(prep$file, select = "xyzr"), error = function(e) NULL)
  if (is.null(las) || lidR::is.empty(las)) return(NULL)

  ## pit-free CHM (lidR, mirrors segment_lidr.R) at RES
  chm <- tryCatch(
    rasterize_canopy(las, res = RES,
      pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5), subcircle = 0.2)),
    error = function(e) NULL)
  if (is.null(chm)) return(NULL)
  ws <- ws_factory(VWF_A)

  ## shared seeds: detect treetops ONCE on the CHM
  ttops <- tryCatch(
    locate_trees(chm, lmf(ws = ws, hmin = 2, shape = "circular")),
    error = function(e) NULL)
  if (is.null(ttops) || nrow(ttops) < 1) return(NULL)
  tc <- sf::st_coordinates(ttops)
  seed_x <- tc[, 1]; seed_y <- tc[, 2]; seed_z <- tc[, 3]
  seed_id <- ttops$treeID
  max_cr_px <- as.integer(round(10 / RES))

  ## ----- seeded segmenters: crown label rasters keyed by ttops treeID -----
  crowns_by_algo <- list()

  d_r <- tryCatch(dalponte2016(chm, ttops, th_seed = 0.45, th_cr = 0.55,
                               max_cr = max_cr_px)(), error = function(e) NULL)
  if (!is.null(d_r)) crowns_by_algo[["dalponte2016"]] <-
      list(geom = crown_geom(d_r), by = "treeID")

  s_r <- tryCatch(silva2016(chm, ttops, max_cr_factor = 0.6, exclusion = 0.3)(),
                  error = function(e) NULL)
  if (!is.null(s_r)) crowns_by_algo[["silva2016"]] <-
      list(geom = crown_geom(s_r), by = "treeID")

  ## lasR region_growing seeded from local_maximum_raster on the lasR CHM
  lasr_g <- tryCatch({
    del <- lasR::triangulate(filter = lasR::keep_first())
    chm2 <- lasR::rasterize(RES, del)
    pf  <- lasR::pit_fill(chm2)
    seed <- lasR::local_maximum_raster(pf, ws, min_height = 2)
    cr  <- lasR::region_growing(pf, seed, th_tree = 2, th_seed = 0.45,
                                th_cr = 0.55, max_cr = 10)
    ans <- lasR::exec(del + chm2 + pf + seed + cr, on = prep$file,
                      with = list(progress = FALSE, ncores = 1L))
    cr_r <- terra::rast(terra::sources(ans$region_growing))
    seeds_sf <- sf::st_as_sf(ans$local_maximum)
    list(geom = crown_geom(cr_r), seeds = sf::st_coordinates(seeds_sf))
  }, error = function(e) NULL)
  if (!is.null(lasr_g) && !is.null(lasr_g$geom))
    crowns_by_algo[["lasr_region_growing"]] <-
      list(geom = lasr_g$geom, by = "own_seed", seeds = lasr_g$seeds)

  ## marker-FREE watershed: crowns matched to seeds by polygon containment
  w_r <- tryCatch(lidR::watershed(chm, th_tree = 2)(), error = function(e) NULL)
  if (!is.null(w_r)) crowns_by_algo[["watershed_markerfree"]] <-
      list(geom = crown_geom(w_r), by = "contain")

  ## random walker, timeboxed
  rw_g <- tryCatch({
    setTimeLimit(elapsed = RW_TIMEOUT, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
    rw <- random_walker(chm, cbind(seed_x, seed_y), beta = 1.0, hmin = 2)
    if (is.null(rw)) NULL else {
      g <- crown_geom(rw$raster)
      # rw label = index of the seed cell (in CHM order); map label -> ttop row
      list(geom = g, seed_cell = rw$seed_cell, seed_label = rw$seed_label)
    }
  }, error = function(e) NULL)
  if (!is.null(rw_g) && !is.null(rw_g$geom))
    crowns_by_algo[["random_walker"]] <-
      list(geom = rw_g$geom, by = "rw", seed_cell = rw_g$seed_cell,
           seed_label = rw_g$seed_label, chm = chm)

  ## ----- match SEED treetops to field stems (position + height gate) ------
  m <- greedy_match(stems$E, stems$N, seed_x, seed_y, TOL,
                    az = stems$height, bz = seed_z)
  matched <- which(m > 0)
  if (!length(matched)) return(NULL)

  rows <- list()
  for (algo in names(crowns_by_algo)) {
    cg <- crowns_by_algo[[algo]]
    geom <- cg$geom
    if (is.null(geom) || nrow(geom) == 0) next
    for (si in matched) {
      det <- m[si]                         # index into seed_* of matched detection
      sx <- seed_x[det]; sy <- seed_y[det]
      crow <- NULL
      if (identical(cg$by, "treeID")) {
        j <- which(geom$treeID == seed_id[det])
        if (length(j)) crow <- geom[j[1], ]
      } else if (identical(cg$by, "rw")) {
        # rw label of this detection's cell -> the crown carrying that label
        cell <- cg$seed_cell                # all seed cells in CHM order
        lab  <- cg$seed_label
        # find the rw seed cell nearest this detection's CHM cell
        dcell <- terra::cellFromXY(cg$chm, cbind(sx, sy))
        idx <- which(cell == dcell)
        if (length(idx)) {
          j <- which(geom$treeID == lab[idx[1]])
          if (length(j)) crow <- geom[j[1], ]
        }
      } else {
        # own_seed (lasR) or contain (watershed): crown containing the seed point
        dd <- (geom$cx - sx)^2 + (geom$cy - sy)^2
        # prefer containment via centroid proximity; fall back to nearest crown
        j <- which.min(dd)
        # containment check: nearest crown must be within ~ its own radius
        if (length(j) && is.finite(dd[j]) &&
            sqrt(dd[j]) <= max(geom$d_caliper[j], 6)) crow <- geom[j, ]
      }
      if (is.null(crow)) next
      rows[[length(rows) + 1]] <- data.frame(
        site = site, plot = pid, algo = algo,
        individualID = stems$individualID[si],
        crown_class = stems$crown_class[si],
        d_eq = crow$d_eq, d_caliper = crow$d_caliper, area = crow$area,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

## ---- per-site driver ------------------------------------------------------
run_site <- function(site) {
  nd  <- file.path(d, "neon", site)
  gt  <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc  <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt  <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  ## join field crown diameter from the rds; keep stems with non-NA field CD
  fc  <- field_crowns(site)
  gt  <- merge(gt, fc, by = "individualID", all.x = TRUE)
  gt  <- gt[!is.na(gt$maxCrownDiameter) | !is.na(gt$ninetyCrownDiameter), ]

  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)

  counts <- table(gt$plotID)
  keep <- names(counts)[counts >= MINTREES]
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] plots with >=%d stems w/ field CD: %d (%s)\n",
              site, MINTREES, length(keep), paste(keep, collapse = ",")))
  if (!length(keep)) return(NULL)

  tmpdir <- file.path(tempdir(), paste0("crown_", site))
  dir.create(tmpdir, showWarnings = FALSE)
  res_list <- mclapply(keep, function(p)
    tryCatch(run_plot(site, p, ctg, pc, gt, tmpdir),
             error = function(e) { message("  plot ", p, " failed: ",
                                            conditionMessage(e)); NULL }),
    mc.cores = CORES, mc.preschedule = FALSE)
  res <- do.call(rbind, Filter(Negate(is.null), res_list))
  if (is.null(res) || !nrow(res)) { cat(sprintf("[%s] no crowns matched\n", site)); return(NULL) }

  ## attach field crown diameter to each matched tree
  fld <- gt[, c("individualID", "maxCrownDiameter", "ninetyCrownDiameter")]
  res <- merge(res, fld, by = "individualID", all.x = TRUE)
  names(res)[names(res) == "maxCrownDiameter"]    <- "field_maxCD"
  names(res)[names(res) == "ninetyCrownDiameter"] <- "field_ninetyCD"
  res <- res[, c("site", "plot", "algo", "crown_class", "individualID",
                 "d_eq", "d_caliper", "area", "field_maxCD", "field_ninetyCD")]
  out <- file.path(nd, "crown_metrics_results.csv")
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("[%s] wrote %d matched-tree rows -> %s\n", site, nrow(res), out))
  res
}

## ---- scoring helpers ------------------------------------------------------
# Pooled error stats over matched trees: RMSE/MAE/bias/R2 (sum of squared errors
# / n, never a mean of per-plot rates). det = detected diameter, fld = field.
err_stats <- function(det, fld) {
  ok <- is.finite(det) & is.finite(fld)
  det <- det[ok]; fld <- fld[ok]; n <- length(det)
  if (n < 2) return(data.frame(n = n, rmse = NA, mae = NA, bias = NA, r2 = NA))
  e <- det - fld
  ss_res <- sum((fld - det)^2)
  ss_tot <- sum((fld - mean(fld))^2)
  data.frame(n = n,
             rmse = sqrt(mean(e^2)), mae = mean(abs(e)),
             bias = mean(e),
             r2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_)
}

print_tables <- function(res) {
  algos <- unique(res$algo)
  cat("\n================= CROWN-DIAMETER ACCURACY (pooled, all sites) ============\n")
  for (defn in list(c("d_eq", "field_ninetyCD", "equiv-circle d_eq vs ninetyCrownDiameter"),
                    c("d_caliper", "field_maxCD", "max-caliper d_caliper vs maxCrownDiameter"))) {
    cat(sprintf("\n--- %s ---\n", defn[3]))
    cat(sprintf("%-22s %5s %7s %7s %7s %7s\n",
                "algo", "n", "rmse", "mae", "bias", "r2"))
    for (a in algos) {
      s <- err_stats(res[[defn[1]]][res$algo == a], res[[defn[2]]][res$algo == a])
      cat(sprintf("%-22s %5d %7.2f %7.2f %7.2f %7.3f\n",
                  a, s$n, s$rmse, s$mae, s$bias, s$r2))
    }
    cat("  by crown class (rmse / n):\n")
    for (a in algos) {
      cat(sprintf("    %-20s", a))
      for (cc in c("dominant", "codominant", "intermediate", "suppressed")) {
        sel <- res$algo == a & res$crown_class == cc
        s <- err_stats(res[[defn[1]]][sel], res[[defn[2]]][sel])
        cat(sprintf(" %s=%s/%d", substr(cc, 1, 4),
                    ifelse(is.na(s$rmse), "NA", sprintf("%.2f", s$rmse)), s$n))
      }
      cat("\n")
    }
  }
}

## ---- run ------------------------------------------------------------------
t0 <- Sys.time()
all_res <- list()
for (site in SITES) {
  r <- tryCatch(run_site(site), error = function(e) {
    message("site ", site, " failed: ", conditionMessage(e)); NULL })
  if (!is.null(r)) all_res[[site]] <- r
}
res <- do.call(rbind, all_res)
dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nDONE: %d matched-tree rows across %d sites in %.1f min\n",
            if (is.null(res)) 0 else nrow(res), length(all_res), dt))
if (!is.null(res) && nrow(res)) {
  cat("\nrows per algorithm:\n"); print(table(res$algo))
  cat("\nrows per crown class:\n"); print(table(res$crown_class, useNA = "ifany"))
  print_tables(res)
}
