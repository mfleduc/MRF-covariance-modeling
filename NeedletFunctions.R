Needlet.Precision <- function(lkinfo){
  ## Currently can only handle up to j=4
  #Have to reindex the needlets to have the smallest j=1 inside lkinfo.
  #For now ignoring effects across levels
  jmin = lkinfo$setupArgs$startingLevel-1;
  jmax = jmin+lkinfo$nlevel-1;
  jList <- seq(jmin,jmax)
  adjMat <- R.matlab::readMat('needletAdjacencyMatrix.mat',sparseMatrixClass='SparseM' ) 
  adjList = list()
  for(j in 1:length(jList)){
    tmp<-unlist(adjMat$distances[[j]] )
    mask <- (tmp!=0)
    tmp[mask] <- 1/tmp[mask]
    adj <-spam::spam(tmp ,nrow = 12*4^jList[j])
    invDists <-  spam::diag.spam(-1/spam::rowSums(adj))
    #Adjusting the SAR matrix for the number of nearest neighbors and their distances
    adj2 <- invDists%*%adj
    adj3<-adj2+(unlist(lkinfo$alpha[j]))*(spam::diag.spam(12*4^jList[j]))
    adjList[[j]] = adj3
    }
  B <- adjList[[1]]
  for(j in 2:length(jList)){
    B = spam::bdiag.spam(B,adjList[[j]])
  }
  return(list(Q=t(B)%*%B,B=B))
}
Needlet.LnLike <- function(y,PHI,lkinfo,tau,r.decay=TRUE,look=FALSE,Q=NULL){
  #PHI is the matrix defining the needlet basis fns. Ideally sparse in some manner.
  #Define a cutoff radius? Unclear how to handle
  nSamples <- length(y)
  if(is.null(Q)){
    Q <- Needlet.Precision(lkinfo)[[1]]
  }
  if(is.na(lkinfo$lambda)){
    lambda <- 1
  }else{
    lambda <- lkinfo$lambda
  }
  G = lambda*Q+1/tau^2*t(PHI)%*%PHI #Assuming W = Id : Errors uncorrelated with the same variance. 
  Gchol <- spam::chol.spam( G ) #Bottleneck, but prevents slow inversion of G
  PtDinvy <- t(PHI) %*% y/tau^2 
  rightPc <- backsolve(Gchol ,forwardsolve(Gchol , transpose=TRUE, PtDinvy, upper.tri=TRUE)) #Follows from Sec 3.1
  # Nychka, Douglas, et al. “A Multiresolution Gaussian Process Model for the Analysis of Large Spatial Datasets.” Journal of Computational and Graphical Statistics, vol. 24, no. 2, 2015, pp. 579–99. JSTOR, http://www.jstor.org/stable/24737282. 
  quadform <- sum(y^2/tau^2)-sum(PtDinvy)*rightPc
  Qchol <- spam::chol.spam(Q)
  logdet <- nSamples*2*log(tau) + 2*sum(log(spam::diag(Gchol))) - 2*sum(log(spam::diag(Qchol)))
    
  lnlike <- -1*(logdet + quadform)
  return(lnlike)
  }
Needlet.FixedEffects<-function(lkinfo,y,Z,PHI,tau,Q=NULL){
  #y: data
  #Z: Matrix of size (n datapoints x n predictors), with each column corresponding to a predictor and each row a lattice point
  #phi: Matrix of basis functions
  #tau: Noise stdev (nugget effect)
  if(is.na(lkinfo$lambda)){
    lambda <- 0.01
  }else{
    lambda <- lkinfo$lambda
  }
  #First: calculate M^_{-1}y
  if(is.null(Q)){
    Q <- Needlet.Precision(lkinfo)[[1]]
  }
  G = lambda*Q+1/tau^2*t(PHI)%*%PHI
  Gchol <-spam::chol.spam(G)
  Pty <- t(PHI)%*%y/tau^2
  v <- spam::backsolve.spam(Gchol,spam::forwardsolve.spam(Gchol, Pty, upper.tri=TRUE))
  Miy <- 1/lambda*(y/tau^2-1/tau^2*PHI%*%v)
  ZTMiy <- t(Z)%*%Miy
  # Now: Calculate (Z^TM^{-1}Z)^{-1} in the same manner
  PtZ <- t(PHI)%*%Z/tau^2
  GiPtZ <- spam::backsolve.spam(Gchol,spam::forwardsolve.spam(Gchol, PtZ, upper.tri=TRUE))
  MiZ <- 1/lambda*(Z-PHI%*%GiPtZ)/tau^2
  ZTMiZ <- t(Z)%*%MiZ
  dhat <- solve(ZTMiZ, ZTMiy)
  return(dhat)
}
Needlet.Simulate <-function(lkinfo,Q=NULL){
  #Return the coefficients to save some small amount of trouble
  jmin = lkinfo$setupArgs$startingLevel-1;
  jmax = jmin+lkinfo$nlevel-1;
  if(is.null(Q)){
    Q <- Needlet.Precision(lkinfo)[[1]]
  }
  Qchol <- chol(Q)
  m <- 12*(4^((jmin):(jmax))) # total number of basis functions
  E <- matrix(rnorm(m), nrow = sum(m), ncol = 1) #Random, normally distributed coefficients
  A <- as.matrix(spam::backsolve(Qchol, E))
  return(A)
}
Needlet.CoeffEstimate <- function(lkinfo, y, PHI, tau,Q=NULL){
  #Following Nychka(2015) Section 3.1
  if(is.null(Q)){
    Q <- Needlet.Precision(lkinfo)[[1]]
  }
  if(is.na(lkinfo$lambda)){
    lambda <- 0.01
  }else{
    lambda <- lkinfo$lambda
  }
  G = lambda*Q+1/tau^2*t(PHI)%*%PHI #Assuming W = Id : Errors uncorrelated with the same variance. 
  Gchol <- spam::chol.spam( G ) #Bottleneck, but prevents slow inversion of G
  v <- t(PHI)%*%y/tau^2
  chat <- spam::backsolve.spam(Gchol,spam::forwardsolve.spam(Gchol, v, upper.tri=TRUE))
  return(chat)
}





