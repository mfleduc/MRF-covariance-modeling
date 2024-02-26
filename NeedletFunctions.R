Needlet.Precision.Marginal<-function(lkinfo){
  ## Currently can only handle up to j=4
  #Have to reindex the needlets to have the smallest j=1 inside lkinfo.
  jmin = lkinfo$startingLevel-1;
  jmax = jmin+lkinfo$nlevel-1;
  
  adjMat <- R.matlab::readMat('needletAdjacencyMatrix.mat',sparseMatrixClass='SparseM' ) 
  adjList = list()
  for(j in jmin:jmax){
    tmp=unlist(adjMat$adjacencies[[j+1]] )
    adj <-spam::spam(tmp ,nrow = 12*4^j)
    rowAdjust =  spam::diag.spam(-1/rowSums(adj))
    #Adjusting the SAR matrix for the number of nearest neighbors
    adj2 <- rowAdjust%*%adj
    adj3<-adj2+unlist(lkinfo$alpha[j+1])*diag(12*4^j)
    }
  
}