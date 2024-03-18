## Lattice Krig testing with GLOW and WACCMX output data
## Code for using the LatticeKrig package on WACCMX output from the Tonga eruption

a = Sys.info()["sysname"]
if(grepl("ndows",as.character(a["sysname"]))){
  setwd('C:/Users/Test/Documents/Research/MRF-covariance-modeling')
  #library("spam64")
  ncobj<-ncdf4::nc_open('C:/Users/Test/Documents/Research/202201010000.nc') #creating nc object
  oxy<-ncdf4::ncvar_get( ncobj,varid='ETA')#, start = c(1,1,12,85), count=c(-1,-1,1,31) )
  lev<-ncdf4::ncvar_get( ncobj,varid='ZZZ')
  PHI <- R.matlab::readMat('C:/Users/Test/Documents/Research/needlet covariance/A 2 deg res.mat')
  PHI <- as.matrix(PHI$A)
  
}else{
  setwd('/homes/male7736/Desktop/Research/MRF-covariance-modeling/')
  source('NeedletFunctions.R')
  library("spam64")
  ncobj<-ncdf4::nc_open('/homes/male7736/Desktop/Research/Data/waccmx_tonga_O_1min_above400km.nc') #creating nc object
  oxy<-ncdf4::ncvar_get( ncobj,varid='O', start = c(1,1,12,85), count=c(-1,-1,1,31) )
  lev<-ncdf4::ncvar_get( ncobj,varid='lev')
  #PHI <- R.matlab::readMat('needletA_deg.mat')
  
}
source("NeedletFunctions.R")
lat<-ncdf4::ncvar_get(ncobj,varid='lat')
lon<-ncdf4::ncvar_get(ncobj,varid='lon')

gridPts<-expand.grid(oce::angleRemap(lon),lat)
# y<-matrix(oxy[,,12,],ncol=1)
# results<-vector(mode='list',length=8)
lvls = seq(1,7)
# for(x in lvls)
# {

nLvls = 7
alphaVals = 1.1^seq(1,nLvls)
lkinfo<-LatticeKrig::LKrigSetup(gridPts, startingLevel=1, nlevel=nLvls,alpha=alphaVals,LKGeometry='LKSphere')
centers<-LatticeKrig::LKrigLatticeCenters(lkinfo, Level=nLvls,physicalCoordinates=TRUE)
nFns = c(0,12,42,162,642,2562,10242,40962)
# 
# distMat = fields::rdist.earth(gridPts,centers[1:nFns[min(5,nLvls)+1],] ,R=1)
# ndcs = apply(distMat, 2, which.min)
# for(x in seq(1,min(5,nLvls))){
#   thisLvl = ndcs[1:nFns[x+1]]
#   thisData = array(oxy,dim=c(361*720,31))
#   thisData = thisData[thisLvl,]
# }
# <-LatticeKrig::LKrig.precision(lkinfo)
#   tmp<-LatticeKrig::LKrig(gridPts, y,LKinfo=lkinfo, lambda=0.0001)
#   results[[x]]<-tmp
# }
# 
# coeffs <- Needlet.Simulate(lkinfo)
# field <- PHI%*%coeffs
# field<-matrix(field, ncol=180)
# fields::image.plot(lon[seq(1,360,2)],lat[seq(1,181,2)],field)
for(x in lvls)
{
  png(filename=paste('Level_',x,'with RBF centers.png'))
  fields::image.plot(lon,lat,oxy[,,12],main=paste("Up to level",x))
  points(centers[(nFns[x]+1):nFns[x+1],1]+180, centers[(nFns[x]+1):nFns[x+1],2], pch=19)
  dev.off()
}
