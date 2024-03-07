setwd('/homes/male7736/Desktop/Research/MRF-covariance-modeling/')
library("spam64")
library("geosphere")
ncobj<-ncdf4::nc_open('/homes/male7736/Desktop/Research/Data/waccmx_tonga_O_1min_above400km.nc') #creating nc object
oxy<-ncdf4::ncvar_get( ncobj,varid='O', start = c(1,1,12,50), count=c(-1,-1,1,200) )
lev<-ncdf4::ncvar_get( ncobj,varid='lev')
lat<-ncdf4::ncvar_get(ncobj,varid='lat')
lon<-ncdf4::ncvar_get(ncobj,varid='lon')
lonNdcs = seq(1,720,by=2)
latNdcs = seq(1,361,by=2)
gridPts<-expand.grid(oce::angleRemap(lon[lonNdcs]),lat[latNdcs])
oxy<-oxy[lonNdcs,latNdcs,]
oxyVector <- matrix(oxy, ncol=200)
centers <- R.matlab::readMat('needletCenters.mat') #Longitude,latitude
r = 6378100

for(i in 1:3){
  centerPts <- matrix(unlist(centers$tps[[i]]),nrow=2)
  dists <- distm(t(centerPts),gridPts,fun=distCosine)/r
  distsCenters <- distm(t(centerPts),fun=distCosine)/r
  minNdcs <- apply(dists, 1, FUN = which.min)
  dataPts = oxyVector[minNdcs,]
  corrMat = cor(t(dataPts))
}