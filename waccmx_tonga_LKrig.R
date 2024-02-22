## Code for using the LatticeKrig package on WACCMX output from the Tonga eruption
setwd('/glade/work/mleduc/Ionosphere_inversion/MRF-covariance-model')
library("spam64")
ncobj<-ncdf4::nc_open('/glade/work/mleduc/Ionosphere_inversion/waccmx_tonga_O_1min_above400km.nc') #creating nc object
oxy<-ncdf4::ncvar_get( ncobj,varid='O', start = c(1,1,1,100), count=c(-1,-1,-1,1) )
lev<-ncdf4::ncvar_get( ncobj,varid='lev')
lat<-ncdf4::ncvar_get(ncobj,varid='lat')
lon<-ncdf4::ncvar_get(ncobj,varid='lon')

gridPts<-expand.grid(oce::angleRemap(lon),lat)


y<-matrix(oxy[,,12],ncol=1)
results<-vector(mode='list',length=8)
lvls = seq(1,7)
for(x in lvls)
{
  alphaVals = 1.1^seq(0,x-1)
  lkinfo<-LatticeKrig::LKrigSetup(gridPts, startingLevel=1, nlevel=x,alpha=alphaVals,LKGeometry='LKSphere')
  #Q     <-LatticeKrig::LKrig.precision(lkinfo)
  tmp<-LatticeKrig::LKrig(gridPts, y,LKinfo=lkinfo, lambda=0.0001)
  results[[x]]<-tmp
  png(filename=paste("lvl 1 to", x, ".png"))
  fields::image.plot(lon,lat,matrix(results[[x]]$residuals,ncol=361),main=paste("Residuals, LKSphere, Levels 1 to", x))
  dev.off()
}
# png(filename='ground truth oxygen concentration data.png')
# fields::image.plot(lon,lat,oxy[,,12],main='Ground truth')
# dev.off()