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

# Export three WGS84 (EPSG:4326, RFC 7946) GeoJSON layers describing the NEON
# density-ladder benchmark geography:
#   sites.geojson  - one Polygon per site (convex hull of its plot centroids),
#                    site-level attributes (plot counts, mean native density,
#                    stem/live-tree totals, mean live-tree height, structure rank)
#   plots.geojson  - one Polygon per plot = the scoring core box (+/-20 m tower,
#                    +/-10 m distributed), with plotType, the `swept` flag, native
#                    first-return/all-return density (null when not swept), and
#                    live-tree / stem counts
#   stems.geojson  - one Point per ground-truth stem (is_tree == TRUE), carrying
#                    species, height, stem diameter, crown class, status, year
#
# Geometry is built in each site's NEON UTM zone (from plot_centroids$utmZone)
# then transformed to WGS84, the GeoJSON standard CRS. Density is read from the
# cached sweep_results.csv native rung, so only swept plots carry a measured
# density; unswept plots get null (density was never measured there - see the
# >=6-live-tree sweep gate in run_sweep.R, NOT a data-quality difference).
#
# Usage:
#   Rscript scripts/export_geojson.R [SITES=SJER,SOAP,TEAK] [OUT=<dir>]
# Reads work/neon/<SITE>/{plot_centroids,ground_truth_stems,sweep_results}.csv.
# Writes <OUT>/{sites,plots,stems}.geojson. OUT defaults to the repo-tracked
# data/ dir, so a plain run regenerates the committed artifacts in place.
suppressMessages({ library(sf) })

d <- .job_dir()
source(.find("sweep_lib.R"))                 # plot_half()

## ---- args ----------------------------------------------------------------
args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
OUT   <- if (is.null(A$OUT)) file.path(.ROOT, "data") else A$OUT
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Cross-site structure gradient (docs/CLAUDE.md: SJER -> SOAP -> TEAK).
STRUCT_RANK <- c(SJER = 1L, SOAP = 2L, TEAK = 3L)
# NEON UTM zone string ("11N") -> EPSG. D17 sites are all UTM 11N (EPSG:32611).
zone_epsg <- function(z) {
  num <- as.integer(sub("[NnSs]$", "", z))
  hemi <- toupper(sub("^[0-9]+", "", z))
  (if (hemi == "S") 32700L else 32600L) + num
}

## ---- per-site assembly ---------------------------------------------------
plot_rows <- list(); stem_rows <- list(); site_rows <- list()

for (site in SITES) {
  nd <- file.path(d, "neon", site)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  swf <- file.path(nd, "sweep_results.csv")
  sw <- if (file.exists(swf)) read.csv(swf, stringsAsFactors = FALSE) else
          data.frame(plot = character(), rung = character(),
                     frdens = numeric(), pdens = numeric())

  epsg <- zone_epsg(pc$utmZone[1])
  domain <- {                                  # parse D17 from an individualID
    m <- regmatches(gt$individualID, regexpr("D[0-9]+", gt$individualID))
    if (length(m)) m[1] else NA_character_
  }

  # live trees = the population the sweep gate counts (run_sweep.R:64)
  gt_live <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  swept <- unique(sw$plot)
  # native-rung density per plot (constant across the chm_res x vwf_a grid)
  natv <- sw[sw$rung == "native", c("plot", "frdens", "pdens")]
  natv <- natv[!duplicated(natv$plot), ]
  dens_of <- function(p, col) { v <- natv[[col]][match(p, natv$plot)]; v }

  ## --- plots: scoring-box polygons ---
  pgon <- function(cx, cy, h)                  # closed CCW square ring
    st_polygon(list(rbind(c(cx-h, cy-h), c(cx+h, cy-h),
                          c(cx+h, cy+h), c(cx-h, cy+h), c(cx-h, cy-h))))
  geoms <- vector("list", nrow(pc))
  pl <- data.frame(plotID = pc$plotID, site = site, plotType = pc$plotType,
                   stringsAsFactors = FALSE)
  pl$swept    <- pc$plotID %in% swept
  pl$frdens   <- dens_of(pc$plotID, "frdens")
  pl$pdens    <- dens_of(pc$plotID, "pdens")
  pl$n_live_core <- NA_integer_; pl$n_live_plot <- NA_integer_
  pl$n_stems_total <- NA_integer_
  for (i in seq_len(nrow(pc))) {
    cx <- pc$easting[i]; cy <- pc$northing[i]; ph <- plot_half(pc$plotType[i])
    geoms[[i]] <- pgon(cx, cy, ph)
    lp <- gt_live[gt_live$plotID == pc$plotID[i], ]
    pl$n_live_plot[i]   <- nrow(lp)
    pl$n_live_core[i]   <- sum(abs(lp$E - cx) <= ph & abs(lp$N - cy) <= ph)
    pl$n_stems_total[i] <- sum(gt$plotID == pc$plotID[i])
  }
  pl_sf <- st_transform(st_sf(pl, geometry = st_sfc(geoms, crs = epsg)), 4326)
  plot_rows[[site]] <- pl_sf

  ## --- stems: points (is_tree == TRUE) ---
  st <- gt[gt$is_tree %in% TRUE & !is.na(gt$E) & !is.na(gt$N), ]
  st_df <- data.frame(individualID = st$individualID, plotID = st$plotID,
                      site = site, taxonID = st$taxonID,
                      scientificName = st$scientificName, height = st$height,
                      stemDiameter = st$stemDiameter, crown_class = st$crown_class,
                      canopyPosition = st$canopyPosition, plantStatus = st$plantStatus,
                      live = st$live, meas_year = st$meas_year,
                      stringsAsFactors = FALSE)
  st_sf <- st_transform(
    st_sf(st_df, geometry = st_sfc(
      lapply(seq_len(nrow(st)), function(j) st_point(c(st$E[j], st$N[j]))),
      crs = epsg)), 4326)
  stem_rows[[site]] <- st_sf

  ## --- site: convex hull of plot centroids ---
  cent <- st_sfc(lapply(seq_len(nrow(pc)),
                        function(i) st_point(c(pc$easting[i], pc$northing[i]))),
                 crs = epsg)
  hull <- st_convex_hull(st_union(cent))
  sr <- data.frame(
    site = site, domain = domain, utm_zone = pc$utmZone[1], epsg = epsg,
    n_plots = nrow(pc), n_plots_swept = sum(pl$swept),
    n_stems_total = nrow(gt), n_live_total = nrow(gt_live),
    frdens_native_mean = round(mean(natv$frdens), 2),
    pdens_native_mean  = round(mean(natv$pdens), 2),
    mean_live_height_m = round(mean(gt_live$height, na.rm = TRUE), 2),
    structure_rank = unname(STRUCT_RANK[site]),
    stringsAsFactors = FALSE)
  site_rows[[site]] <- st_transform(st_sf(sr, geometry = hull), 4326)
}

## ---- write ---------------------------------------------------------------
sites_sf <- do.call(rbind, site_rows)
plots_sf <- do.call(rbind, plot_rows)
stems_sf <- do.call(rbind, stem_rows)

opts <- c("RFC7946=YES", "COORDINATE_PRECISION=7")
write_gj <- function(x, name) {
  f <- file.path(OUT, name)
  st_write(x, f, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE,
           layer_options = opts)
  cat(sprintf("  %-14s %4d features -> %s\n", name, nrow(x), f))
}
cat(sprintf("Writing GeoJSON (WGS84) for %s:\n", paste(SITES, collapse = ", ")))
write_gj(sites_sf, "sites.geojson")
write_gj(plots_sf, "plots.geojson")
write_gj(stems_sf, "stems.geojson")
cat("DONE\n")
