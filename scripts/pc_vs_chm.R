suppressMessages(library(lidR))
las <- readLAS("aoi.laz",
  filter="-drop_class 7 18 -drop_withheld -keep_xy 586418 4521371 586568 4521521")
las <- normalize_height(las, tin())
ws  <- function(h){y<-0.1*h+3; y[h<2]<-2; y[h>20]<-5; y}

## 1) CHM-based lmf (the method used in the AOI run)
chm <- rasterize_canopy(las, 0.25,
         pitfree(c(0,10,20,30), c(0,1.5), subcircle=0.2))
n_chm <- nrow(locate_trees(chm, lmf(ws=ws, hmin=2, shape="circular")))

## 2) point-cloud lmf (lmf directly on points)
tt_pc <- locate_trees(las, lmf(ws=ws, hmin=2, shape="circular"))
n_pc  <- nrow(tt_pc)

## 3) point-cloud SEGMENTATION (Li 2012) -- true 3D, can split understory
t0 <- Sys.time()
seg <- segment_trees(las, li2012(dt1=1.5, dt2=2, R=2, Zu=15, hmin=2, speed_up=10))
li_dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
ids <- unique(seg$treeID); ids <- ids[!is.na(ids)]
# apex height per Li segment
ap <- tapply(seg$Z[!is.na(seg$treeID)], seg$treeID[!is.na(seg$treeID)], max)
n_li <- length(ids)

cat(sprintf("CHM-lmf trees      : %d\n", n_chm))
cat(sprintf("point-cloud lmf    : %d\n", n_pc))
cat(sprintf("Li2012 segments    : %d   (%.0f s)\n", n_li, li_dt))
cat(sprintf("Li2012 apex height : median %.1f m, q25 %.1f, q75 %.1f\n",
            median(ap), quantile(ap,.25), quantile(ap,.75)))
cat(sprintf("Li2012 trees < 5 m : %d (%.0f%%)  [sub-dominant/regen]\n",
            sum(ap<5), 100*mean(ap<5)))
