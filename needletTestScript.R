## Lattice Krig testing with GLOW and WACCMX output data
#Step 1: Read in the test data
#This data set is only over a patch on the sphere but w/e
# testData <- R.matlab::readMat('C:/Users/Test/Documents/Research/needlet covariance/testData.mat')
# #Step 2: Read in the needlets 
needletA <- R.matlab::readMat('C:/Users/Test/Documents/Research/needlet covariance/needletA.mat')
# gridPts <- R.matlab::readMat('C:/Users/Test/Documents/Research/needlet covariance/grid points.mat')
# needletA <- matrix(unlist(needletA),ncol=2901)
# gridPts <- matrix(unlist(gridPts),ncol=2)
# print('Doing something')
# 
# alt140Data <- array(unlist(testData),dim=c(41,71,5))
# y<-matrix(alt140Data[1:41, 1:71, 3],ncol=1)
# lkinfo <- LatticeKrig::LKrigSetup(gridPts, nlevel=3, alpha=c(1,1,1),NC=71)
# results<-LatticeKrig::LKrig( y , LKInfo = lkinfo)

ncObj <- ncdf4::nc_open( 'C:/Users/Test/Documents/Research/waccmx_tonga_2Dionosphere_fields.nc' )

eef <- ncdf4::ncvar_get(ncObj, varid='ED1',start=c(1,1,250),count=c(-1,-1,1))

mlat<-ncdf4::ncvar_get(ncObj, varid='mlat')
mlon<-ncdf4::ncvar_get(ncObj,varid='mlon')

LON<-rep(mlon, 385)+180
LAT<-rep(mlat, 320)
LAT<-matrix(LAT,ncol=320)
LAT<-matrix(t(LAT),ncol=1)
gridPts<-list(LON,LAT)
gridPts<-matrix(unlist(gridPts),ncol=2)
eef<-matrix(eef,ncol=1)
#lkInfo<-LatticeKrig::LKrigSetup(gridPts, startingLevel=1,nlevel=5,alpha=c(1,1,1,1,1),LKGeometry='LKSphere')



#mleResult<-LatticeKrig::LKrig.MLE(gridPts, eef, LKinfo=lkInfo)

#resLKrig<-LatticeKrig::LKrig(gridPts,eef,LKinfo=lkInfo,lambda=mleResult$lambda.MLE)

#resLatticeKrig<-LatticeKrig::LatticeKrig(gridPts, eef , LKGeometry='LKSphere')


