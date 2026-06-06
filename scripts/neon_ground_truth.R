#!/usr/bin/env Rscript
# Build field-stem ground truth for NEON SOAP woody-vegetation structure
# (DP1.10098.001), geolocated to UTM 11N, paired with the 2021 high-density
# LiDAR (DP1.30003.001, ~20 pts/m^2).
#
# Reimplements the relevant part of geoNEON::getLocTOS() using the public NEON
# locations API (no untrusted code): each mapped stem records a (pointID,
# stemDistance, stemAzimuth) offset from a named grid point inside its base
# plot; we fetch each unique point's UTM coordinate from the API and apply the
# polar offset. Output: ground_truth_stems.csv + tiles_needed.csv.
#
# Ground truth pairs each stem with the apparentindividual measurement nearest
# the 2021 LiDAR acquisition (within +/-4 yr). Two columns support the
# exact-year temporal-sensitivity filter consumed by run_sweep.R /
# validate_heights.R (issue #5): `meas_year` = calendar year of that chosen
# measurement, and `dist21` = |meas_year - 2021|. Selecting meas_year==2021
# re-scores against only stems measured in the LiDAR year. These columns already
# exist in the cached ground_truth_stems.csv -- do NOT re-run this script for the
# temporal cut; it would overwrite the shared ground truth.
#
# Env: CLAUDE_JOB_DIR (working dir; default ./work)
suppressMessages({ library(neonUtilities); library(jsonlite) })

## args: SITE=SOAP (default), used for both veg + tile-assignment
args <- strsplit(commandArgs(TRUE), "=")
A    <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
site <- if (is.null(A$SITE)) "SOAP" else A$SITE

d   <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
nd  <- file.path(d, "neon", site); dir.create(nd, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load woody-veg structure (all years, RELEASE-2026) ----------------
rds_all <- file.path(nd, "vst", paste0(tolower(site), "_vst_allyears.rds"))
if (!file.exists(rds_all)) {
  options(timeout = 1200)
  dat <- loadByProduct(dpID = "DP1.10098.001", site = site,
                       package = "basic", release = "RELEASE-2026",
                       check.size = FALSE, progress = FALSE)
  dir.create(dirname(rds_all), showWarnings = FALSE, recursive = TRUE)
  saveRDS(dat, rds_all)
} else dat <- readRDS(rds_all)

mt <- dat$vst_mappingandtagging
ai <- dat$vst_apparentindividual
pp <- dat$vst_perplotperyear
cat(sprintf("[%s] mappingandtagging: %d rows; apparentindividual: %d rows (years %s)\n",
            site, nrow(mt), nrow(ai),
            paste(range(substr(ai$date, 1, 4), na.rm = TRUE), collapse = "-")))

## plot centroids (most-recent row per plot)
pp <- pp[!is.na(pp$easting), ]
pp <- pp[order(pp$plotID, -as.integer(substr(pp$date, 1, 4))), ]
pp1 <- pp[!duplicated(pp$plotID),
          c("plotID", "plotType", "easting", "northing", "utmZone")]
write.csv(pp1, file.path(nd, "plot_centroids.csv"), row.names = FALSE)

## ---- 2. Mappable stems: geolocate via NEON locations API ------------------
map <- mt[!is.na(mt$stemDistance) & !is.na(mt$stemAzimuth) & !is.na(mt$pointID), ]
# A few individuals carry >1 mapping/tagging record (re-mapped in a later year);
# keep one row per stem so the merge below cannot double-count the recall
# denominator (sibling apparentindividual is already deduped).
map <- map[!duplicated(map$individualID), ]
map$ptloc <- paste0(map$namedLocation, ".", map$pointID)
uniq <- unique(map$ptloc)
cat(sprintf("mappable stems: %d across %d plots; %d unique point locations\n",
            nrow(map), length(unique(map$plotID)), length(uniq)))

cache_f <- file.path(nd, "point_locations.csv")
cache <- if (file.exists(cache_f)) read.csv(cache_f, stringsAsFactors = FALSE) else
         data.frame(ptloc = character(), easting = numeric(),
                    northing = numeric(), zone = character(), unc = numeric())
need <- setdiff(uniq, cache$ptloc)
if (length(need)) {
  cat(sprintf("fetching %d point locations from NEON API...\n", length(need)))
  fetch_one <- function(nm) {
    url <- paste0("https://data.neonscience.org/api/v0/locations/", nm)
    j <- tryCatch(fromJSON(url), error = function(e) NULL)
    if (is.null(j)) return(NULL)
    dd <- j$data
    pr <- dd$locationProperties
    unc <- NA_real_
    if (!is.null(pr) && "locationPropertyValue" %in% names(pr)) {
      k <- grep("Coordinate uncertainty", pr$locationPropertyName)
      if (length(k)) unc <- as.numeric(pr$locationPropertyValue[k[1]])
    }
    data.frame(ptloc = nm, easting = dd$locationUtmEasting,
               northing = dd$locationUtmNorthing,
               zone = as.character(dd$locationUtmZone), unc = unc)
  }
  add <- do.call(rbind, lapply(need, fetch_one))
  cache <- rbind(cache, add)
  write.csv(cache, cache_f, row.names = FALSE)
}
map <- merge(map, cache, by = "ptloc", all.x = TRUE)
ok <- !is.na(map$easting)
cat(sprintf("geolocated %d / %d stems\n", sum(ok), nrow(map)))
map <- map[ok, ]

# Polar offset: azimuth measured clockwise from grid north.
az <- map$stemAzimuth * pi / 180
map$E <- map$easting  + map$stemDistance * sin(az)
map$N <- map$northing + map$stemDistance * cos(az)
map$pos_unc <- ifelse(is.na(map$unc), 0.3, map$unc) + 0.3  # point + rangefinder

## ---- 3. Join nearest-to-2021 apparentindividual measurement ---------------
ai$year <- as.integer(substr(ai$date, 1, 4))
ai <- ai[!is.na(ai$year), ]
ai$dist21 <- abs(ai$year - 2021)
ai <- ai[order(ai$individualID, ai$dist21), ]
ai1 <- ai[!duplicated(ai$individualID),
          c("individualID", "year", "dist21", "height", "stemDiameter",
            "maxCrownDiameter", "ninetyCrownDiameter",
            "plantStatus", "canopyPosition", "growthForm")]
g <- merge(map, ai1, by = "individualID", all.x = TRUE)

# Live trees only for the recall denominator (nearest meas. within 4 yr).
g$live <- grepl("Live", g$plantStatus) & (is.na(g$dist21) | g$dist21 <= 4)
# Trees only (drop shrubs/lianas where growthForm says so).
g$is_tree <- is.na(g$growthForm) | grepl("tree", g$growthForm, ignore.case = TRUE)

# Crown class from canopyPosition; fall back to within-plot height quantiles.
cc_map <- c("Open grown" = "dominant", "Full sun" = "dominant",
            "Partially shaded" = "codominant",
            "Mostly shaded" = "intermediate", "Full shade" = "suppressed")
g$crown_class <- cc_map[g$canopyPosition]
need_cc <- is.na(g$crown_class) & !is.na(g$height)
if (any(need_cc)) {
  for (p in unique(g$plotID[need_cc])) {
    idx <- which(g$plotID == p & need_cc)
    if (length(idx) >= 4) {
      q <- quantile(g$height[idx], c(.33, .66, .9), na.rm = TRUE)
      g$crown_class[idx] <- as.character(  # as.character: a bare factor would
        cut(g$height[idx], c(-Inf, q[1], q[2], q[3], Inf),  # coerce to int codes
            labels = c("suppressed", "intermediate", "codominant", "dominant")))
    }
  }
}

## ---- 4. Tile assignment (1 km NEON tiles, SW origin) ----------------------
g$tile_e <- floor(g$E / 1000) * 1000
g$tile_n <- floor(g$N / 1000) * 1000
g$tile   <- sprintf("%d_%d", g$tile_e, g$tile_n)

out <- g[, c("individualID", "plotID", "tile", "tile_e", "tile_n", "E", "N",
             "pos_unc", "taxonID", "scientificName", "height", "stemDiameter",
             "maxCrownDiameter", "ninetyCrownDiameter",
             "plantStatus", "canopyPosition", "crown_class", "live", "is_tree",
             "year", "dist21")]
names(out)[names(out) == "year"] <- "meas_year"
write.csv(out, file.path(nd, "ground_truth_stems.csv"), row.names = FALSE)

live_tree <- out[out$live & out$is_tree, ]
tiles <- aggregate(individualID ~ tile + tile_e + tile_n, data = live_tree,
                   FUN = length)
names(tiles)[4] <- "n_live_stems"
tiles <- tiles[order(-tiles$n_live_stems), ]
write.csv(tiles, file.path(nd, "tiles_needed.csv"), row.names = FALSE)

cat(sprintf("\n=== ground truth: %d geolocated stems; %d live trees ===\n",
            nrow(out), nrow(live_tree)))
cat("crown-class distribution (live trees):\n")
print(table(live_tree$crown_class, useNA = "ifany"))
cat(sprintf("\ntiles with live stems: %d (top by count)\n", nrow(tiles)))
print(head(tiles, 12), row.names = FALSE)
