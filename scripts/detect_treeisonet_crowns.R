#!/usr/bin/env Rscript
# TreeisoNet treeOff crown arm (#20). Runs the treeOff crown driver
# (gpu/run_treeisonet_crowns.py) on the native normalized frozen clip per SOAP
# plot (SERIAL, single GPU), reduces the per-point instances to apexes + crown
# diameters (the bridge's reduce_instances + crown_diameter_table on the labeled
# CANOPY points), matches apexes to field stems (greedy_match with a height
# gate), and joins NEON field crown diameter. Output is schema-compatible with
# crown_metrics_results.csv (algo = "treeisonet") so analyze_crown_metrics.R can
# union it with the five CHM segmenters from #7.
#
# Crown geometry caveat: the CHM arms take diameters from a dissolved label-
# raster polygon; here it is the convex hull of the crown's canopy points
# (crown_diameter_table). d_eq -> ninetyCrownDiameter, d_caliper -> maxCrownDiameter.
#
# Usage:
#   Rscript scripts/detect_treeisonet_crowns.R [SITE=SOAP] [PLOTS=ALL]
#       [CONF=0.22] [TOL=4] [HMIN=2]
# Output: $CLAUDE_JOB_DIR/neon/<SITE>/treeisonet_crown_metrics.csv
suppressMessages({ library(lidR); library(data.table); library(grDevices) })
options(lidR.progress = FALSE)
d <- Sys.getenv("CLAUDE_JOB_DIR", file.path(getwd(), "work"))
.script_path <- function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile) && length(ofile) && nzchar(ofile)) return(ofile)
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(hit)) return(sub("^--file=", "", hit[1]))
  NA_character_
}
.sp <- .script_path()
.ROOT <- if (!is.na(.sp)) normalizePath(file.path(dirname(.sp), ".."),
                                        mustWork = FALSE) else getwd()
.find <- function(rel) Find(file.exists, c(file.path(.ROOT, "scripts", rel),
                                           file.path("scripts", rel),
                                           file.path("..", "..", "scripts", rel),
                                           file.path(getwd(), "scripts", rel)))
source(.find("sweep_lib.R")); source(.find("model_bench_lib.R"))
source(.find("model_runner.R"))

args  <- strsplit(commandArgs(TRUE), "=")
A     <- setNames(lapply(args, `[`, 2), sapply(args, `[`, 1))
SITE  <- if (is.null(A$SITE))  "SOAP" else A$SITE
PLOTS <- if (is.null(A$PLOTS) || A$PLOTS == "ALL") NULL else strsplit(A$PLOTS, ",")[[1]]
CONF  <- if (is.null(A$CONF))  "0.22" else A$CONF
TOL   <- as.numeric(if (is.null(A$TOL)) 4.0 else A$TOL)
HMIN  <- if (is.null(A$HMIN)) "2" else A$HMIN
MINTREES <- 6
VENV  <- file.path(.ROOT, "gpu/.venv/bin/python")
DRV   <- file.path(.ROOT, "gpu/run_treeisonet_crowns.py")
LOC   <- file.path(.ROOT, "gpu/store/treeaibox/als_treeloc.pth")
LCFG  <- Sys.glob(file.path(.ROOT, "gpu/store/treeaibox/*reclamation*treeloc*.json"))[1]
OFF   <- file.path(.ROOT, "gpu/store/treeaibox/als_treeoff.pth")
OCFG  <- Sys.glob(file.path(.ROOT, "gpu/store/treeaibox/*reclamation*treeoff*.json"))[1]

# Field crown diameter per individualID, nearest-to-2021 measurement (mirrors
# crown_metrics_sweep.R::field_crowns -- the vst rds is the canonical source).
field_crowns <- function(site) {
  rds <- file.path(d, "neon", site, "vst", paste0(tolower(site), "_vst_allyears.rds"))
  ai  <- as.data.frame(readRDS(rds)$vst_apparentindividual)
  ai$year <- as.integer(substr(ai$date, 1, 4))
  ai <- ai[!is.na(ai$year), ]; ai$dist21 <- abs(ai$year - 2021)
  ai <- ai[order(ai$individualID, ai$dist21), ]
  ai[!duplicated(ai$individualID),
     c("individualID", "maxCrownDiameter", "ninetyCrownDiameter")]
}

run_main <- function() {
  stopifnot(file.exists(VENV), file.exists(DRV), file.exists(LOC), !is.na(LCFG),
            file.exists(OFF), !is.na(OCFG))
  nd <- file.path(d, "neon", SITE)
  gt <- read.csv(file.path(nd, "ground_truth_stems.csv"), stringsAsFactors = FALSE)
  pc <- read.csv(file.path(nd, "plot_centroids.csv"),     stringsAsFactors = FALSE)
  gt <- gt[gt$live & gt$is_tree & !is.na(gt$E), ]
  gt <- gt[, setdiff(names(gt), c("maxCrownDiameter", "ninetyCrownDiameter")), drop = FALSE]
  gt <- merge(gt, field_crowns(SITE), by = "individualID", all.x = TRUE)
  gt <- gt[!is.na(gt$maxCrownDiameter) | !is.na(gt$ninetyCrownDiameter), ]
  laz <- list.files(file.path(nd, "lidar"), pattern = "\\.laz$",
                    recursive = TRUE, full.names = TRUE)
  ctg <- readLAScatalog(laz, progress = FALSE)
  counts <- table(gt$plotID); keep <- names(counts)[counts >= MINTREES]
  if (!is.null(PLOTS)) keep <- intersect(keep, PLOTS)
  keep <- intersect(keep, pc$plotID)
  cat(sprintf("[%s] treeisonet crowns: %d plots w/ field CD (conf=%s)\n",
              SITE, length(keep), CONF))

  rows <- list()
  for (pid in keep) {                       # SERIAL -- single GPU
    ci <- pc[pc$plotID == pid, ][1, ]; cx <- ci$easting; cy <- ci$northing
    ph <- plot_half(ci$plotType)
    stems <- gt[gt$plotID == pid & abs(gt$E - cx) <= ph & abs(gt$N - cy) <= ph, ]
    if (nrow(stems) < 1) next
    prep <- tryCatch(frozen_clip(ctg, SITE, pid, NA, cx, cy, ph,
                                 out_root = file.path(nd, "frozen")),
                     error = function(e) NULL)
    if (is.null(prep)) next
    ocsv <- file.path(tempdir(), sprintf("ticr_%s.csv", pid))
    pts <- run_python_crown_arm(VENV, DRV, prep$normalized, ocsv,
      extra = c(LOC, LCFG, OFF, OCFG, "0", CONF, HMIN), timeout = 900)
    if (is.null(pts) || !nrow(pts)) next
    dt  <- as.data.table(pts)
    ap  <- dt[, .(x = X[which.max(Z)], y = Y[which.max(Z)], z = max(Z)), by = crown_id]
    cd  <- crown_diameter_table(pts, id_col = "crown_id", min_pts = 5)
    ap  <- merge(ap, cd, by.x = "crown_id", by.y = "id", all.x = TRUE)
    m <- greedy_match(stems$E, stems$N, ap$x, ap$y, TOL,
                      az = stems$height, bz = ap$z)
    for (si in which(m > 0)) {
      cr <- ap[m[si], ]
      rows[[length(rows) + 1]] <- data.frame(site = SITE, plot = pid,
        algo = "treeisonet", crown_class = stems$crown_class[si],
        individualID = stems$individualID[si], d_eq = cr$d_eq,
        d_caliper = cr$d_caliper, area = if (is.na(cr$d_eq)) NA_real_
          else pi * (cr$d_eq / 2)^2,
        field_maxCD = stems$maxCrownDiameter[si],
        field_ninetyCD = stems$ninetyCrownDiameter[si], stringsAsFactors = FALSE)
    }
    cat(sprintf("  %s: %d crowns, %d matched\n", pid, nrow(ap), sum(m > 0)))
  }
  res <- do.call(rbind, rows)
  if (is.null(res) || !nrow(res)) { cat("no treeisonet crowns matched\n"); return(invisible()) }
  out <- file.path(nd, "treeisonet_crown_metrics.csv")
  write.csv(res, out, row.names = FALSE)
  cat(sprintf("[%s] treeisonet crowns DONE: %d matched-tree rows -> %s\n",
              SITE, nrow(res), out))
}

if (sys.nframe() == 0L) run_main()
