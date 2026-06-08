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

suppressMessages({library(lidR); library(parallel)})
options(lidR.progress = FALSE)
J <- .job_dir()
tiledir <- file.path(J, "tiles"); dir.create(tiledir, showWarnings = FALSE)

## retile the block into 150 m tiles on disk
ctg <- readLAScatalog(file.path(J, "block_big.laz"))
opt_chunk_size(ctg)    <- 150
opt_chunk_buffer(ctg)  <- 0
opt_output_files(ctg)  <- file.path(tiledir, "tile_{XLEFT}_{YBOTTOM}")
opt_progress(ctg)      <- FALSE
rt <- catalog_retile(ctg)
tiles <- list.files(tiledir, pattern = "\\.la[sz]$", full.names = TRUE)
area_m2 <- as.numeric(sf::st_area(sf::st_as_sfc(sf::st_bbox(ctg))))
cat(sprintf("block: %.1f ha, %d tiles (150 m)\n", area_m2/1e4, length(tiles)))

## per-tile: drop noise -> normalize -> Li 2012 segmentation -> tree count
one <- function(f) {
  las <- readLAS(f, filter = "-drop_class 7 18 -drop_withheld")
  if (is.empty(las) || npoints(las) < 100) return(c(n = 0, pts = npoints(las)))
  las <- normalize_height(las, tin())
  seg <- segment_trees(las, li2012(dt1=1.5, dt2=2, R=2, Zu=15, hmin=2, speed_up=10))
  c(n = length(unique(na.omit(seg$treeID))), pts = npoints(las))
}

## run on 16 cores, time the whole parallel segmentation
t <- system.time(res <- mclapply(tiles, one, mc.cores = 16))[3]
ok  <- !sapply(res, function(x) inherits(x, "try-error"))
M   <- do.call(rbind, res[ok])
ntrees <- sum(M[, "n"]); npts <- sum(M[, "pts"])
ha <- area_m2/1e4
cat(sprintf("Li2012 on 16 cores: %.1f s for %.1f ha (%.2fM pts) -> %d trees\n",
            t, ha, npts/1e6, ntrees))
cat(sprintf("throughput: %.2f s/ha (16 cores)\n", t/ha))
ac <- 10000 * 0.404686                      # 10,000 acres in ha
cat(sprintf("=> 10,000 acres (%.0f ha) extrapolated: %.1f min on 16 cores\n",
            ac, (t/ha)*ac/60))
