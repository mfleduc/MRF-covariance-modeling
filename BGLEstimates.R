## For estimates with the Basis Graphical Lasso
library("BasisGraphicalLasso")
library("LatticeKrig")

set.seed(684161)
a = Sys.info()["sysname"]
# if(grepl("ndows",as.character(a["sysname"]))){
#   setwd('C:/Users/Test/Documents/Research/MRF-covariance-modeling')
#   print('working')
#   ncobj <- ncdf4::nc_open('C:/Users/Test/Documents/Research/202201010000.nc') #creating nc object
#   oxy<-ncdf4::ncvar_get( ncobj,varid='ETA', start = c(1,1,1), count=c(-1,-1,-1) )
#   lev<-ncdf4::ncvar_get( ncobj,varid='ZZZ')
#   PHI <- R.matlab::readMat('C:/Users/Test/Documents/Research/needlet covariance/A 2 deg res.mat')
#   PHI <- as.matrix(PHI$A)
#   Phi<-PHI
# }else{
#   setwd('/homes/male7736/Desktop/Research/MRF-covariance-modeling/')
#   library("spam64")
#   ncobj<-ncdf4::nc_open('/homes/male7736/Desktop/Research/Data/waccmx_tonga_O_1min_above400km.nc') #creating nc object
#   oxy<-ncdf4::ncvar_get( ncobj,varid='O', start = c(1,1,12,85), count=c(-1,-1,1,21) )
#   lev<-ncdf4::ncvar_get( ncobj,varid='lev')
#   # PHI <- R.matlab::readMat('needletA_2deg.mat')
#   needletData <- R.matlab::readMat('needletA_fullswath.mat')
#   needletCenters <- R.matlab::readMat('needletcenters_fullswath.mat')
#   PHI <- as.matrix(needletData$A)
#   numNeedlets <- needletData$numj
#   Phi<-PHI
#   PHI <- as.spam(PHI, eps = 10^-1)
#   rm(needletData)
# }
source("NeedletFunctions.R")
source("NeedletBGL.R")
lat<-ncdf4::ncvar_get(ncobj,varid='lat')
lon<-ncdf4::ncvar_get(ncobj,varid='lon')

latsSwath = (seq(-90,90,1))
lonsSwath = 0.5*round(2*rbind(5/9*(latsSwath+90), 150+5/9*(latsSwath+90)),digits=0)

gridPts<-expand.grid(lon,lat)
mask = gridPts[,1]>1e6;
for(ii in 1:length(latsSwath)){
  m2 = abs(gridPts[,2]-latsSwath[ii])<10^-1&(gridPts[,1]>=lonsSwath[1,ii])&(gridPts[,1]<=lonsSwath[2,ii])
  mask = mask | m2
}

truth <- matrix(oxy[,,1],ncol=1)
y <- truth[mask]
Z <- cbind(gridPts[mask,2],gridPts[mask,1])
leaveOut <- sample(1:length(latsSwath), 50)
latsToMask = sort(latsSwath[leaveOut])
mask2 = Z[,1]>1e6
for( kk in 1:length(leaveOut)){
  mask2 = mask2 | abs(Z[,1]-latsToMask[kk])<0.01
}
basisFns <- PHI[!mask2,]
data = as.matrix(y[!mask2])
locs = Z[!mask2,]

for( ii in 1:5){
  leaveOut <- sample(1:length(latsSwath), 50)
  latsToMask = sort(latsSwath[leaveOut])
  mask2 = Z[,1]>1e6
  for( kk in 1:25){
    mask2 = mask2 | abs(Z[,1]-latsToMask[ii])<0.01
  }
  truth <- matrix(oxy[,,ii+1],ncol=1) #Use different data for each round of training
  y <- truth[mask]
  data = rbind(data, as.matrix(y[!mask2]));
  locs = rbind(locs, Z[!mask2,]);
  basisFns <- rbind(basisFns,PHI[!mask2,])
}
# 
maxIter = 50
# 
locs = cbind(locs[,2],locs[,1])
centers180 <- WrapTo180(needletCenters$centers)
precision.fitl1<-Needlet.BGL(data, locs, basisFns[,1:5], basisCenters = centers180[,1:5], distance.penalty=TRUE,lambda=7)
precision.fitl2<-Needlet.BGL(data, locs, basisFns[,6:25], basisCenters = centers180[,6:25], distance.penalty=TRUE,lambda=7)
precision.fitl3<-Needlet.BGL(data, locs, basisFns[,26:107], basisCenters = centers180[,26:107], distance.penalty=TRUE,lambda=7)
precision.fitl4<-Needlet.BGL(data, locs, basisFns[,108:431], basisCenters = centers180[,108:431], distance.penalty=TRUE,lambda=7)
precision.fitl5<-Needlet.BGL(data, locs, basisFns[,432:1728], basisCenters = centers180[,432:1728], distance.penalty=TRUE,lambda=7)
# PHI <- as.spam(PHI, eps = 10^-2)
# precision.fitl1 <- BGL(y=as.matrix(data), locs=as.matrix(locs), lambda=7, basis="LatticeKrig",
#                        distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, Phi = as.matrix(basisFns[,1:5]))
# precision.fitl2 <- BGL(y=data, locs=locs, lambda=7, basis="LatticeKrig",
#                        distance.penalty=FALSE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, Phi = as.matrix(basisFns[,6:25]))
# precision.fitl2 <- BGL(y=data, locs=locs, lambda=7, basis="LatticeKrig",
#                        distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, NC=30, startingLevel=2,nlevel=1,LKGeometry = 'LKSphere')
# 
# precision.fitl3 <- BGL(y=data, locs=locs, lambda=7, basis="LatticeKrig",
#                        distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, NC=30, startingLevel=3,nlevel=1,LKGeometry = 'LKSphere')
# 
# precision.fitl4 <- BGL(y=data, locs=locs, lambda=7, basis="LatticeKrig",
#                        distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, NC=30, startingLevel=4,nlevel=1,LKGeometry = 'LKSphere')
# precision.fitl5 <- BGL(y=data, locs=locs, lambda=7, basis="LatticeKrig",
#                        distance.penalty=TRUE,outer_tol=1e-2, MAX_ITER=maxIter,
#                        MAX_RUNTIME_SECONDS=86400, NC=30, startingLevel=5,nlevel=1,LKGeometry = 'LKSphere')
lkInfo <- LKrigSetup(locs, a.wght=4.05, nu=0.5,NC=30, startingLevel=1,nlevel=5,LKGeometry = 'LKSphere')
# 
Q <- bdiag.spam(precision.fitl1$Q,precision.fitl2$Q,precision.fitl3$Q,precision.fitl4$Q,precision.fitl5$Q)
# # Should probably do a different fit for each level desired so we preserve the block structure
maskPred <- abs(Z[,1])<=5
Zmasked <- Z[!maskPred,]
Zpred <- sin(pi/180*Zmasked)
truth <- matrix(oxy[,,12],ncol=1) #Use different data for each round of training and then estimate on a different set
y <- truth[mask]
ypred <- y[!maskPred]
onesVec <- vector(mode='numeric',length=dim(Zpred)[1])+1
Zpred <- cbind(onesVec, Zpred)
dhat <- Needlet.FixedEffects(lkInfo, ypred,Zpred, PHI[!maskPred,],1,Q=Q)

r <- ypred-Zpred%*%dhat
#
chat <- Needlet.CoeffEstimate(lkInfo, r, PHI[!maskPred,], 1,Q=Q)
# # 
Zbg = sin(pi/180*Z)
onesBG <- vector(mode='numeric',length=dim(Zbg)[1])+1
Zbg <- cbind(onesBG, Zbg)
fieldEst <- Zbg%*%dhat + PHI[,]%*%chat
fields::quilt.plot(Z,fieldEst,main='Estimate',zlim=c(0.93,0.99))
# fields::quilt.plot(Z,y,main='truth',zlim=c(0.93,0.99))
# fields::quilt.plot(Z,y-fieldEst,main='truth-estimate',zlim=c(0.93,0.99))