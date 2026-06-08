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

suppressMessages(library(lidR))
options(lidR.progress = FALSE)
J <- .job_dir()
las0 <- readLAS(file.path(J,"aoi.laz"),
  filter="-drop_class 7 18 -drop_withheld -keep_xy 586418 4521371 586568 4521521")
area <- 150*150
ws <- function(h){y<-0.1*h+3; y[h<2]<-2; y[h>20]<-5; y}

run_one <- function(las){
  fr  <- sum(las$ReturnNumber==1L)/area              # first-return density
  res <- if(fr>=8) 0.25 else if(fr>=4) 0.50 else 1.0 # density-derived CHM res
  lasn<- normalize_height(las, tin())
  tc  <- system.time({chm <- rasterize_canopy(lasn,res,
            pitfree(c(0,10,20,30),c(0,1.5),subcircle=0.2))
          n_chm <- nrow(locate_trees(chm, lmf(ws=ws,hmin=2,shape="circular")))})[3]
  tp  <- system.time(n_pc <- nrow(locate_trees(lasn, lmf(ws=ws,hmin=2,shape="circular"))))[3]
  tl  <- system.time(seg <- segment_trees(lasn,
            li2012(dt1=1.5,dt2=2,R=2,Zu=15,hmin=2,speed_up=10)))[3]
  n_li<- length(unique(na.omit(seg$treeID)))
  data.frame(fr_density=round(fr,1), res=res, kpts=round(npoints(las)/1e3),
             chm_n=n_chm, chm_s=round(tc,1),
             pc_n=n_pc,   pc_s=round(tp,1),
             li_n=n_li,   li_s=round(tl,1))
}

rows <- list()
for(d in c(1,2,4,8,16)){
  las <- decimate_points(las0, random(d))
  rows[[as.character(d)]] <- run_one(las)
}
rows[["native"]] <- run_one(las0)
out <- do.call(rbind, rows)
print(out, row.names=FALSE)
