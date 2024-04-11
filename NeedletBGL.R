## Library for implementation of the Basis Graphical Lasso algorithm as described in 
## https://github.com/mlkrock/BasisGraphicalLasso/tree/master for use with other basis functions
## The code as in that github repo can work with a user supplied basis, however 
## there is currently no way to enforce the distance penalty if you wish to use 
## a non-LatticeKrig basis. This code will lean on the BGL package for the 
## computations it can handle, and any other work will be done here.
## This code will need the BasisGraphicalLasso library (and its dependencies) 
## to run properly. If you want to do computations with the LatticeKrig basis, 
## I strongly suggest you use Mitch's code
## One key difference: This code uses the spam format for the basis functions, which 
## are assumed to be locally supported.
## -Matt
library("BasisGraphicalLasso")
library("LatticeKrig")
library("spam")
library("geosphere")
Needlet.BGL <- function(y,locs,Phi,basisCenters=NULL,lambda=1,tau_sq=NULL,
            zero.diagonal.penalty=TRUE,crossValidation=FALSE,
            guess=NULL,outer_tol=NULL,MAX_ITER=NULL,MAX_RUNTIME_SECONDS=NULL,
            verbose=TRUE,distance.penalty=FALSE,distanceFn='cosine'){
  ## This function adds a wrinkle: The distances between basis centers can be 
  ## calculated by using either cosine distance if the centers are lat/lon or 
  ## euclidean distance if desired. The default is cosine distance, and the 
  ## corresponding flags are 'cosine' and 'l2' respectively
  basis.setup <- Needlet.BGLBasisSetup(y,locs,Phi,basisCenters=basisCenters,
            crossValidation=crossValidation,distance.penalty=distance.penalty,
            distanceFn=distanceFn)
                        
  Phi_Phi <- basis.setup$Phi_Phi
  Phi_S_Phi <- basis.setup$Phi_S_Phi
  trS <- basis.setup$trS
  basisDistMat <- basis.setup$basisdistancematrix
  if(is.null(tau_sq))
  {
    cat('Estimating nugget effect variance: ')#Done using BGL package
    tau_sq <- nugget_estimate(Phi_Phi,Phi_S_Phi,trS,n=dim(locs)[1])
    cat(tau_sq, '\n', sep="")
  }
  l <- dim(Phi_Phi)[1]
  if(is.null(guess))
  {
    guess <- diag(l)
  }
  # Set up the penalty matrix using BGL
  penaltymatrix <- penaltymatrixsetup(lambda=lambda,zero.diagonal.penalty=TRUE,
                                l=l,basisdistancematrix=basisDistMat)
  # Execute DC algorithm
  #Again: This uses the original BGL package. 
  cat('Running basis graphical lasso algorithm:\n')
  finalguess <- BGL_DC(penaltymatrix,Phi_Dinv_Phi=Phi_Phi/tau_sq,Phi_Dinv_S_Dinv_Phi=Phi_S_Phi/(tau_sq^2),guess=guess,outer_tol,MAX_ITER,MAX_RUNTIME_SECONDS)
  
  # Return estimated precision matrix and nugget effect variance
  mylist <- list(Q=finalguess,nugget_variance=tau_sq,Phi=basis.setup$Phi)
  
  return(mylist)
}
Needlet.BGLBasisSetup <- function(y,locs,Phi,basisCenters = NULL,
                          crossValidation=FALSE,distance.penalty=FALSE,
                          distanceFn = 'cosine'){
  Phi_Phi <- crossprod.spam(Phi)
  Phi_y <- crossprod.spam(Phi,y)
  trS <- norm(y, "F")^2/ifelse(ncol(y)==1,1,ncol(y)-1)
  if(crossValidation){
    mylist <- list(Phi_Phi=Phi_Phi, Phi_y=Phi_y, trS=trS)
  } else {
    Phi_S_Phi <- tcrossprod.spam(Phi_y)/ ifelse(ncol(Phi_y)==1,1,ncol(Phi_y)-1)
    mylist <- list(Phi_Phi=Phi_Phi, Phi_S_Phi=Phi_S_Phi, trS=trS)
  }
  if(distance.penalty){
    if(distanceFn=='l2'){ #Euclidean distance
      mylist$basisDistanceMatrix <- rdist(basisCenters)
    }else if(distanceFn=='cosine'){ #Cosine distance for data on the sphere
      r = 6378100
      mylist$basisDistanceMatrix <- distm(t(basisCenters),fun=distCosine)/r
    }else{
      stop("So far only cosine ('cosine') and euclidean ('l2') distances are implemented")
    }
  }
  return(mylist)
}
WrapTo180 <- function(lons){
  mask = lons>180;
  lons[mask] <- lons[mask]-360
  return(lons)
}