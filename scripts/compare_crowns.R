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

# Compare lasR vs lidR crown polygons:
#   - counts and area distribution
#   - centroid-based greedy matching within a tolerance
#   - IoU between matched crown pairs
# Usage:
#   Rscript scripts/compare_crowns.R [lasR.gpkg] [lidR.gpkg]
# Defaults to AOI outputs if no args provided.
suppressMessages({ library(sf); library(terra) })

d  <- .job_dir()
ar <- commandArgs(trailingOnly = TRUE)
fa <- if (length(ar) >= 1) ar[1] else file.path(d, "crowns_lasr_aoi.gpkg")
fb <- if (length(ar) >= 2) ar[2] else file.path(d, "crowns_lidr_aoi.gpkg")
stopifnot(file.exists(fa), file.exists(fb))

A <- sf::st_read(fa, quiet = TRUE)
B <- sf::st_read(fb, quiet = TRUE)

# Normalize CRS (assume both in same CRS already; warn if not)
if (!is.na(sf::st_crs(A)) && !is.na(sf::st_crs(B)) && sf::st_crs(A) != sf::st_crs(B)) {
  warning("CRS mismatch; reprojecting B to A")
  B <- sf::st_transform(B, sf::st_crs(A))
}

# Centroids and per-polygon area
ca <- sf::st_coordinates(sf::st_centroid(A))
cb <- sf::st_coordinates(sf::st_centroid(B))
area_a <- as.numeric(sf::st_area(A))
area_b <- as.numeric(sf::st_area(B))

cat(sprintf("Crowns: lasR=%d  lidR=%d  diff=%+d (%.1f%%)\n",
            nrow(A), nrow(B), nrow(A) - nrow(B),
            100 * (nrow(A) - nrow(B)) / nrow(B)))
cat(sprintf("Area m²: lasR mean=%.1f median=%.1f  |  lidR mean=%.1f median=%.1f\n",
            mean(area_a), median(area_a), mean(area_b), median(area_b)))

# Greedy centroid matching, IoU per matched pair (skip empty intersection).
greedy_match <- function(ca, cb, tol = 3) {
  used <- rep(FALSE, nrow(cb)); idx <- integer(0); j_idx <- integer(0)
  for (i in seq_len(nrow(ca))) {
    dd <- sqrt((cb[, 1] - ca[i, 1])^2 + (cb[, 2] - ca[i, 2])^2)
    dd[used] <- Inf
    j <- which.min(dd)
    if (length(j) && is.finite(dd[j]) && dd[j] <= tol) {
      used[j] <- TRUE; idx <- c(idx, i); j_idx <- c(j_idx, j)
    }
  }
  data.frame(a = idx, b = j_idx)
}

for (tol in c(2.0, 3.0)) {
  mm <- greedy_match(ca, cb, tol = tol)
  if (nrow(mm) == 0) {
    cat(sprintf("tol=%.1f m | matched=0\n", tol)); next
  }
  iou <- vapply(seq_len(nrow(mm)), function(k) {
    gi <- suppressWarnings(sf::st_intersection(A[mm$a[k], ], B[mm$b[k], ]))
    if (nrow(gi) == 0) return(0)
    gu <- suppressWarnings(sf::st_union(A[mm$a[k], ], B[mm$b[k], ]))
    ai <- sum(as.numeric(sf::st_area(gi)))
    au <- sum(as.numeric(sf::st_area(gu)))
    if (au == 0) 0 else ai / au
  }, numeric(1))
  cat(sprintf("tol=%.1f m | matched=%d (%.0f%% of lidR) | IoU mean=%.2f median=%.2f\n",
              tol, nrow(mm), 100 * nrow(mm) / nrow(B),
              mean(iou), median(iou)))
}
