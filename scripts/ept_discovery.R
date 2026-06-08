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

# EPT discovery helper for the native-QL2 cross-check (issue #4).
# Fetches the entwine/USGS 3DEP public boundaries (resources.geojson), reprojects
# each NEON site centroid (UTM 11N / EPSG:32611) into the index CRS (EPSG:4326),
# and finds every public EPT project footprint that covers the site. Candidates
# are written to work/neon/<SITE>/ql2/ept_candidates.csv, ranked with CA QL2
# projects preferred. Run standalone or source()'d for `discover_ept()`.
#
# Usage:
#   Rscript scripts/ept_discovery.R SITES=SOAP,SJER,TEAK [INDEX=/path/resources.geojson]
suppressMessages({ library(sf); library(jsonlite) })

# Public 3DEP/entwine boundary index. usgs.entwine.io/data/boundaries.json is
# gone (404 as of 2026-06); the hobuinc mirror is the live source. We also accept
# a local cached copy via INDEX= to avoid re-downloading the ~8.7 MB file.
INDEX_URLS <- c(
  "https://raw.githubusercontent.com/hobuinc/usgs-lidar/master/boundaries/resources.geojson",
  "https://usgs.entwine.io/boundaries/resources.geojson"
)

fetch_index <- function(local = NULL) {
  if (!is.null(local) && nzchar(local) && file.exists(local)) {
    message("EPT index: using local copy ", local)
    return(local)
  }
  dest <- file.path(tempdir(), "ept_resources.geojson")
  if (file.exists(dest) && file.info(dest)$size > 1e6) return(dest)
  for (u in INDEX_URLS) {
    ok <- tryCatch({
      utils::download.file(u, dest, quiet = TRUE, mode = "wb")
      file.exists(dest) && file.info(dest)$size > 1e6
    }, error = function(e) FALSE)
    if (ok) { message("EPT index fetched from ", u); return(dest) }
    message("EPT index source failed: ", u)
  }
  stop("Could not fetch any EPT boundary index")
}

# Score a candidate project name for "recent CA QL2" preference. Higher = better.
score_candidate <- function(name) {
  s <- 0
  if (grepl("^CA[_]", name)) s <- s + 100               # California project
  yr <- suppressWarnings(as.integer(
    sub(".*?([12][0-9]{3}).*", "\\1", name)))
  if (!is.na(yr) && yr >= 1990 && yr <= 2100) s <- s + (yr - 1990)  # newer better
  if (grepl("QL2|_2_", name)) s <- s + 10               # QL2 hint in name
  s
}

discover_ept <- function(sites, job_dir, index_local = NULL,
                          site_epsg = 32611) {
  idx_path <- fetch_index(index_local)
  bnd <- st_read(idx_path, quiet = TRUE)
  if (is.na(st_crs(bnd))) st_crs(bnd) <- 4326           # resources.geojson = WGS84
  bnd <- st_make_valid(bnd)
  out <- list()
  for (site in sites) {
    pc_path <- file.path(job_dir, "neon", site, "plot_centroids.csv")
    if (!file.exists(pc_path)) { message("no centroids for ", site); next }
    pc <- read.csv(pc_path)
    # site footprint = all plot centroids (UTM), buffered 1 km then unioned, so a
    # project covering the plot margins (not just the exact centroid) still
    # registers as a candidate. Buffer in the metric site CRS (km is meaningless
    # in degrees), then reproject the buffered hull to the WGS84 index CRS.
    pts <- st_as_sf(pc, coords = c("easting", "northing"), crs = site_epsg)
    site_union <- st_transform(st_union(st_buffer(pts, 1000)), 4326)
    pts4326 <- st_transform(pts, 4326)            # raw points for per-project cov
    # candidate projects: those whose footprint intersects the site point cloud
    hit <- bnd[st_intersects(bnd, site_union, sparse = FALSE)[, 1], ]
    if (nrow(hit) == 0) {
      message(site, ": NO covering EPT project found in index")
      cand <- data.frame(site = site, name = NA_character_, url = NA_character_,
                         point_count = NA_real_, n_plots_covered = NA_integer_,
                         score = NA_real_, preferred = NA)
    } else {
      # per-project coverage = how many plot centroids fall inside it
      cov <- st_intersects(hit, pts4326)
      ncov <- vapply(cov, length, integer(1))
      sc <- vapply(as.character(hit$name), score_candidate, numeric(1))
      cand <- data.frame(
        site = site, name = as.character(hit$name),
        url = as.character(hit$url),
        point_count = suppressWarnings(as.numeric(hit$count)),
        n_plots_covered = ncov, score = sc,
        stringsAsFactors = FALSE)
      cand <- cand[order(-cand$score, -cand$n_plots_covered), ]
      cand$preferred <- FALSE
      cand$preferred[1] <- TRUE
    }
    dir.create(file.path(job_dir, "neon", site, "ql2"),
               showWarnings = FALSE, recursive = TRUE)
    write.csv(cand,
              file.path(job_dir, "neon", site, "ql2", "ept_candidates.csv"),
              row.names = FALSE)
    cat(sprintf("\n=== %s: %d candidate project(s) ===\n", site,
                sum(!is.na(cand$name))))
    show <- cand[!is.na(cand$name), c("name", "point_count",
                                      "n_plots_covered", "score", "preferred")]
    if (nrow(show)) print(show, row.names = FALSE)
    if (any(!is.na(cand$url)))
      cat("preferred URL:", cand$url[cand$preferred %in% TRUE][1], "\n")
    out[[site]] <- cand
  }
  invisible(out)
}

if (sys.nframe() == 0L) {
  args <- strsplit(commandArgs(TRUE), "=")
  A <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
  job_dir <- .job_dir()
  sites <- if (is.null(A$SITES)) c("SOAP", "SJER", "TEAK") else
    strsplit(A$SITES, ",")[[1]]
  index_local <- if (is.null(A$INDEX)) NULL else A$INDEX
  discover_ept(sites, job_dir, index_local = index_local)
}
