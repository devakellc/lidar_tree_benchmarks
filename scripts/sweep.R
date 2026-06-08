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

# How close can we get to the 205-tree reference, and at what cost?
# 1) Extract reference tree positions (apex per treeID) as pseudo-truth.
# 2) Sweep the sensitivity knobs (window scale, hmin, CHM res) in lidR.
# 3) Report count, recall, precision vs the reference -- not just count.
suppressMessages(library(lidR))
d <- .job_dir()
f <- system.file("extdata", "MixedConifer.las", package = "lasR")

las <- readLAS(f)                      # has treeID from the bundled segmentation
dt  <- las@data
ref <- do.call(rbind, lapply(split(dt[!is.na(dt$treeID)], dt$treeID[!is.na(dt$treeID)]),
         function(g) { i <- which.max(g$Z); data.frame(x = g$X[i], y = g$Y[i], z = g$Z[i]) }))
cat(sprintf("reference trees (apex per treeID): %d\n\n", nrow(ref)))

# Greedy match: returns matched count (ref<->det within tol).
match_n <- function(a, b, tol = 2) {
  used <- rep(FALSE, nrow(b)); m <- 0
  for (i in seq_len(nrow(a))) {
    dd <- sqrt((b$x - a$x[i])^2 + (b$y - a$y[i])^2); dd[used] <- Inf
    j <- which.min(dd)
    if (length(j) && is.finite(dd[j]) && dd[j] <= tol) { used[j] <- TRUE; m <- m + 1 }
  }
  m
}

# Window function with a scale k (k<1 => smaller windows => more, closer tops).
mk_ws <- function(k) function(h) {
  y <- k * (0.1 * h + 3); y[h < 2] <- 3 * k; y[h > 20] <- 5 * k; pmax(y, 0.5)
}

run <- function(res, k, hmin) {
  chm <- rasterize_canopy(las, res = res,
           pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5), subcircle = 0.2))
  tt  <- locate_trees(chm, lmf(ws = mk_ws(k), hmin = hmin, shape = "circular"))
  xy  <- sf::st_coordinates(tt)
  det <- data.frame(x = xy[, 1], y = xy[, 2], z = tt$Z)
  m   <- match_n(ref, det, 2)
  data.frame(res = res, k = k, hmin = hmin, n = nrow(det),
             recall = m / nrow(ref), precision = m / nrow(det))
}

grid <- expand.grid(res = c(0.5, 0.25), k = c(1.0, 0.7, 0.5), hmin = c(2, 1))
res  <- do.call(rbind, Map(function(r, k, h) run(r, k, h), grid$res, grid$k, grid$hmin))
res$F1 <- with(res, 2 * recall * precision / (recall + precision))
res <- res[order(-res$F1), ]
print(within(res, { recall <- round(recall, 2); precision <- round(precision, 2)
                    F1 <- round(F1, 2) }), row.names = FALSE)
