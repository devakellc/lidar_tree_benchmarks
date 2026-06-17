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

# Export a WGS84 stem reference layer for the treetop benchmark. Unlike the
# original data/stems.geojson, this keeps the stems that were filtered out of the
# detection ground-truth denominator and records why.
#
# Ground-truth inclusion mirrors the detector scripts:
#   live & is_tree & coordinates present
#   plot has >= 6 live tree stems before core clipping
#   stem lies inside the plot-type scoring core
#
# Usage:
#   Rscript scripts/export_stems_ground_truth_geojson.R \
#     [SITES=SJER,SOAP,TEAK] [OUT=work/neon/best_treetops_geojson]
# Writes:
#   OUT/stems.geojson
#   OUT/stems_ground_truth_used.geojson
#   OUT/stems_filtered_out.geojson
suppressMessages(library(sf))

d <- .job_dir()
source(.find("sweep_lib.R")) # plot_half()

args <- strsplit(commandArgs(TRUE), "=")
A <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITES <- if (is.null(A$SITES)) c("SJER", "SOAP", "TEAK") else strsplit(A$SITES, ",")[[1]]
OUT <- if (is.null(A$OUT)) file.path(d, "neon", "best_treetops_geojson") else A$OUT
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

zone_epsg <- function(z) {
  num <- as.integer(sub("[NnSs]$", "", z))
  hemi <- toupper(sub("^[0-9]+", "", z))
  (if (hemi == "S") 32700L else 32600L) + num
}

add_reason <- function(reasons, idx, label) {
  if (!any(idx, na.rm = TRUE)) return(reasons)
  reasons[idx] <- Map(function(x) c(x, label), reasons[idx])
  reasons
}

site_stems <- function(site) {
  nd <- file.path(d, "neon", site)
  gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"), stringsAsFactors = FALSE)
  epsg <- zone_epsg(pc$utmZone[1])

  coord_ok <- !is.na(gt$E) & !is.na(gt$N)
  live_tree_coord <- (gt$live %in% TRUE) & (gt$is_tree %in% TRUE) & coord_ok
  live_plot_count <- table(gt$plotID[live_tree_coord])
  n_live_plot <- as.integer(live_plot_count[match(gt$plotID, names(live_plot_count))])
  n_live_plot[is.na(n_live_plot)] <- 0L

  mi <- match(gt$plotID, pc$plotID)
  has_centroid <- !is.na(mi)
  plot_type <- ifelse(has_centroid, pc$plotType[mi], NA_character_)
  cx <- ifelse(has_centroid, pc$easting[mi], NA_real_)
  cy <- ifelse(has_centroid, pc$northing[mi], NA_real_)
  core_half <- ifelse(has_centroid, plot_half(plot_type), NA_real_)
  in_core <- has_centroid & coord_ok &
    abs(gt$E - cx) <= core_half & abs(gt$N - cy) <= core_half
  pass_plot_gate <- n_live_plot >= 6L
  ground_truth_used <- live_tree_coord & has_centroid & pass_plot_gate & in_core

  reasons <- rep(list(character()), nrow(gt))
  reasons <- add_reason(reasons, !(gt$live %in% TRUE), "not_live_or_stale_measurement")
  reasons <- add_reason(reasons, !(gt$is_tree %in% TRUE), "not_tree_growth_form")
  reasons <- add_reason(reasons, !coord_ok, "missing_coordinates")
  reasons <- add_reason(reasons, !has_centroid, "missing_plot_centroid")
  reasons <- add_reason(reasons, live_tree_coord & has_centroid & !pass_plot_gate,
                        "plot_below_6_live_tree_gate")
  reasons <- add_reason(reasons, live_tree_coord & has_centroid & pass_plot_gate & !in_core,
                        "outside_scoring_core")
  filter_reason <- vapply(reasons, paste, character(1), collapse = ";")
  filter_reason[ground_truth_used] <- "used_for_ground_truth"

  keep <- coord_ok
  df <- data.frame(
    individualID = gt$individualID[keep],
    site = site,
    plotID = gt$plotID[keep],
    plotType = plot_type[keep],
    taxonID = gt$taxonID[keep],
    scientificName = gt$scientificName[keep],
    height = gt$height[keep],
    stemDiameter = gt$stemDiameter[keep],
    plantStatus = gt$plantStatus[keep],
    canopyPosition = gt$canopyPosition[keep],
    crown_class = gt$crown_class[keep],
    live = gt$live[keep],
    is_tree = gt$is_tree[keep],
    meas_year = gt$meas_year[keep],
    dist21 = gt$dist21[keep],
    pos_unc = gt$pos_unc[keep],
    n_live_tree_plot = n_live_plot[keep],
    scoring_core_half_m = core_half[keep],
    in_scoring_core = in_core[keep],
    ground_truth_used = ground_truth_used[keep],
    filtered_out = !ground_truth_used[keep],
    filter_reason = filter_reason[keep],
    E = gt$E[keep],
    N = gt$N[keep],
    stringsAsFactors = FALSE)

  sf::st_transform(sf::st_as_sf(df, coords = c("E", "N"), crs = epsg,
                                remove = FALSE), 4326)
}

layers <- lapply(SITES, site_stems)
stems <- do.call(rbind, layers)
used <- stems[stems$ground_truth_used %in% TRUE, ]
filtered <- stems[stems$filtered_out %in% TRUE, ]

write_gj <- function(x, file) {
  sf::st_write(x, file.path(OUT, file), driver = "GeoJSON", delete_dsn = TRUE,
               quiet = TRUE,
               layer_options = c("RFC7946=YES", "COORDINATE_PRECISION=7"))
  cat(sprintf("  %-26s %5d features -> %s\n", file, nrow(x), file.path(OUT, file)))
}

cat(sprintf("Writing stem GeoJSON for %s:\n", paste(SITES, collapse = ", ")))
write_gj(stems, "stems.geojson")
write_gj(used, "stems_ground_truth_used.geojson")
write_gj(filtered, "stems_filtered_out.geojson")
cat(sprintf("ground_truth_used=%d filtered_out=%d total=%d\n",
            sum(stems$ground_truth_used), sum(stems$filtered_out), nrow(stems)))
