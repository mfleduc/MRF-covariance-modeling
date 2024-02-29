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
    tmp<-unlist(adjMat$adjacencies[[j]] )
    adj <-spam::spam(tmp ,nrow = 12*4^jList[j])
    rowAdjust <-  spam::diag.spam(-1/rowSums(adj))
    #Adjusting the SAR matrix for the number of nearest neighbors
    adj2 <- rowAdjust%*%adj
    adj3<-adj2+(unlist(lkinfo$alpha[j])+0.001)*(spam::diag.spam(12*4^jList[j]))
    adjList[[j]] = adj3
    }
  B <- adjList[[1]]
  for(j in 2:length(jList)){
    B = spam::bdiag.spam(B,adjList[[j]])
  }
  return(t(B)%*%B)
}
Needlet.LnLike <- function(p,y,PHI,lkinfo,tau,r.decay=TRUE,look=FALSE){
  #PHI is the matrix defining the needlet basis fns. Ideally sparse in some manner.
  #Define a cutoff radius? Unclear how to handle.
  Q <- Needlet.Precision(lkinfo)
  if(is.na(lkinfo$lambda)){
    lambda <- 1
  }else{
    lambda <- lkinfo$lambda
  }
  G = lambda*Q+1/tau^2*t(PHI)%*%PHI #Assuming W = Id : Errors uncorrelated with the same variance. 
  Gchol <- chol.spam( G ) #Bottleneck, but prevents slow inversion of G
}
Needlet.FixedEffects<-function(lkinfo,y,Z,PHI){
  #y: data
  #Z: Matrix of size (n datapoints x n predictors), with each column corresponding to a predictor and each row a lattice point
  #phi: Matrix of basis functions
}
Needlet.Simulate <-function(lkinfo){
  #Return the coefficients to save some small amount of trouble
  jmin = lkinfo$setupArgs$startingLevel-1;
  jmax = jmin+lkinfo$nlevel-1;
  Q <- Needlet.Precision(lkinfo)
  Qchol <- chol.spam(Q)
  m <- 12*(4^((jmin):(jmax))) # total number of basis functions
  E <- matrix(rnorm(m), nrow = m, ncol = 1) #Random, normally distributed coefficients
  
}
Q<-Needlet.Precision(lkinfo)