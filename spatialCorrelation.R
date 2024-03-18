
library("spam64")
library("geosphere")
a = Sys.info()["sysname"]
if(grepl("ndows",as.character(a["sysname"]))){
  setwd('C:/Users/Test/Documents/Research/MRF-covariance-modeling')
  #library("spam64")
  ncobj<-ncdf4::nc_open('C:/Users/Test/Documents/Research/202201010000.nc') #creating nc object
  oxy<-ncdf4::ncvar_get( ncobj,varid='ETA')#, start = c(1,1,12,85), count=c(-1,-1,1,31) )
  lev<-ncdf4::ncvar_get( ncobj,varid='ZZZ')
}else{
  setwd('/homes/male7736/Desktop/Research/MRF-covariance-modeling/')
  ncobj<-ncdf4::nc_open('/homes/male7736/Desktop/Research/Data/waccmx_tonga_O_1min_above400km.nc') #creating nc object
  oxy<-ncdf4::ncvar_get( ncobj,varid='O', start = c(1,1,12,50), count=c(-1,-1,1,200) )
  lev<-ncdf4::ncvar_get( ncobj,varid='lev')
}
source("NeedletFunctions.R")
lat<-ncdf4::ncvar_get(ncobj,varid='lat')
lon<-ncdf4::ncvar_get(ncobj,varid='lon')
lonNdcs = seq(1,length(lon),by=2)
latNdcs = seq(1,length(lat),by=2)
gridPts<-expand.grid(oce::angleRemap(lon[lonNdcs]),lat[latNdcs])
oxy<-oxy[lonNdcs,latNdcs,]
oxyVector <- matrix(oxy, ncol=dim(oxy)[3])
centers <- R.matlab::readMat('needletCenters.mat') #Longitude,latitude
r = 6378100
lastNdx = 0
#Precision and SAR matrices
nLvls = 4
alphaVals = 1.1^seq(1,nLvls)
lkinfo<-LatticeKrig::LKrigSetup(gridPts, startingLevel=1, nlevel=nLvls,alpha=alphaVals,LKGeometry='LKSphere')
Q<-Needlet.Precision(lkinfo)
B<-fields::spam2full(Q[[2]])
Q<-fields::spam2full(Q[[1]])
Qi <- solve(Q)
for(i in 1:1){
  centerPts <- matrix(unlist(centers$tps[[i]]),nrow=2)
  dists <- distm(t(centerPts),gridPts,fun=distCosine)/r
  distsCenters <- distm(t(centerPts),fun=distCosine)/r
  minNdcs <- apply(dists, 1, FUN = which.min)
  dataPts = oxyVector[minNdcs,]
  corrMat = cor(t(dataPts))
  locs = Q[lastNdx+1,(lastNdx+1):(lastNdx+dim(centerPts)[2])]!=0
  fields::quilt.plot( t(centerPts),Q[lastNdx+1,(lastNdx+1):(lastNdx+dim(centerPts)[2])],main=paste('Precision matrix, j=',i-1))
  points(t(centerPts[,locs]), pch=19)
  fields::quilt.plot( t(centerPts), diag(Qi[seq(lastNdx+1,lastNdx+dim(centerPts)[2]),seq(lastNdx+1,lastNdx+dim(centerPts)[2])]),main=paste('Marginal variances, j=',i-1))
 # fields::quilt.plot( t(centerPts), diag(corrMat),main='empirical variances' )
  lastNdx = lastNdx+dim(centerPts)[2]
  }
# plot3D::scatter2D(centerPts[1,]+180,centerPts[2,],colvar = corrMat[700,],pch=19)
                                                                   


