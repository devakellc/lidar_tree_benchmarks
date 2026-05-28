suppressMessages(library(lasR))
d   <- Sys.getenv("CLAUDE_JOB_DIR")
tif <- file.path(d, "chm_aoi.tif")
laz <- file.path(d, "aoi.laz")

# Same density-derived wfloor as detect_*_aoi.R so "same ws" is literal.
ans0 <- exec(reader_las(filter = "-drop_class 7 18 -drop_withheld") +
               rasterize(1, "count", filter = keep_first()), on = laz)
v    <- terra::values(ans0); dens <- mean(v[!is.na(v) & v > 0])
spacing <- 1 / sqrt(dens); wfloor <- max(2, round(2.5 * spacing, 1))
ws  <- function(h){y<-0.1*h+3;y[h<2]<-wfloor;y[h>20]<-5;y}
cat(sprintf("shared_chm_aoi: density=%.2f -> wfloor=%.1f m\n", dens, wfloor))

r   <- load_raster(tif)
ans <- exec(r + local_maximum_raster(r, ws, min_height=2), on=laz)
tops<- ans[[length(ans)]]; cc<-sf::st_coordinates(tops)
A   <- data.frame(x=cc[,1], y=cc[,2])
suppressMessages(library(lidR))
chm <- terra::rast(tif)
tt  <- locate_trees(chm, lmf(ws=ws, hmin=2, shape="circular"))
xy  <- sf::st_coordinates(tt); B<-data.frame(x=xy[,1], y=xy[,2])
mn  <- function(a,b,tol){u<-rep(FALSE,nrow(b));m<-0
  for(i in seq_len(nrow(a))){dd<-sqrt((b$x-a$x[i])^2+(b$y-a$y[i])^2);dd[u]<-Inf
    j<-which.min(dd); if(length(j)&&is.finite(dd[j])&&dd[j]<=tol){u[j]<-TRUE;m<-m+1}}; m}
cat(sprintf("SAME CHM (AOI) | lasR LM=%d  lidR LM=%d\n", nrow(A), nrow(B)))
for(tol in c(1,2)){m<-mn(A,B,tol)
  cat(sprintf("  tol=%.1f m matched=%d Jaccard=%.2f\n",tol,m,m/(nrow(A)+nrow(B)-m)))}
