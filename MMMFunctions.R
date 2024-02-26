################################################################################################
## Multivariate multiresolution model (MMM) functions
## WARNING: Two fields must have same basis structure
## WARNING: Two fields must have same a.wght
################################################################################################

MMM.precision.marginal <- function(LKinfo, r){
    L <- LKinfo$nlevel
    if (length(r) != L) {
      stop("Mismatch between resolution and cross-dependence vector")
    }
    offset <- LKinfo$latticeInfo$offset
    LKinfoCheck(LKinfo)
    ind <- NULL
    ra <- NULL
    da <- rep(0, 2)
    for (j in 1:L) {
        tempB <- LKrigSAR(LKinfo, Level = j)
        alpha.level <- (LKinfo$alpha)[[j]]
        tempra <- 1/sqrt(alpha.level[1] * (1-r[j]^2)) * tempB$ra #
        ra <- c(ra, tempra)
        ind <- rbind(ind, tempB$ind + offset[j])
        da[1] <- da[1] + tempB$da[1]
        da[2] <- da[2] + tempB$da[2]
    }
    if ((da[1] != offset[L + 1]) | (da[2] != offset[L + 1])) {
        stop("Mismatch of dimension with size in LKinfo")
    }
    tempB <- list(ind = ind, ra = ra, da = da)
    tempB <- LKrig.spind2spam(tempB)
    return(t(tempB) %*% (tempB))
}

MMM.precision.cross <- function(LKinfo1, LKinfo2, r){
    if (LKinfo1$nlevel != LKinfo2$nlevel) {
      stop("Mismatch in resolution between processes")
    }
    L <- LKinfo1$nlevel
    if (length(r) != L) {
      stop("Mismatch between resolution and cross-dependence vector")
    }
    # Variable 1
    offset <- LKinfo1$latticeInfo$offset
    ind <- NULL
    ra <- NULL
    da <- rep(0, 2)
    for (j in 1:L) {
        tempB <- LKrigSAR(LKinfo1, Level = j)

        alpha1.level <- (LKinfo1$alpha)[[j]]
        alpha2.level <- (LKinfo2$alpha)[[j]]
        
        tempra <- tempB$ra / (alpha1.level[1] * alpha2.level[1])^0.25
        tempra <- tempra * sqrt(r[j]) / sqrt(1 - r[j]^2)
        
        ra <- c(ra, tempra)
        ind <- rbind(ind, tempB$ind + offset[j])
        da[1] <- da[1] + tempB$da[1]
        da[2] <- da[2] + tempB$da[2]
    }
    if ((da[1] != offset[L + 1]) | (da[2] != offset[L + 1])) {
        stop("Mismatch of dimension with size in LKinfo")
    }
    tempB <- list(ind = ind, ra = ra, da = da)
    tempB1 <- LKrig.spind2spam(tempB)
    # Variable 2
    offset <- LKinfo1$latticeInfo$offset
    ind <- NULL
    ra <- NULL
    da <- rep(0, 2)
    for (j in 1:L) {
        tempB <- LKrigSAR(LKinfo2, Level = j)

        alpha1.level <- (LKinfo1$alpha)[[j]]
        alpha2.level <- (LKinfo2$alpha)[[j]]
        
        tempra <- tempB$ra / (alpha1.level[1] * alpha2.level[1])^0.25
        tempra <- tempra * sqrt(r[j]) / sqrt(1 - r[j]^2)
        
        ra <- c(ra, tempra)
        ind <- rbind(ind, tempB$ind + offset[j])
        da[1] <- da[1] + tempB$da[1]
        da[2] <- da[2] + tempB$da[2]
    }
    if ((da[1] != offset[L + 1]) | (da[2] != offset[L + 1])) {
        stop("Mismatch of dimension with size in LKinfo")
    }
    tempB <- list(ind = ind, ra = ra, da = da)
    tempB2 <- LKrig.spind2spam(tempB)
    return(t(tempB1) %*% (tempB2))
}

MMM.precision <- function(LKinfo1, LKinfo2, r){
  Q11 <- MMM.precision.marginal(LKinfo=LKinfo1,r=r)
  Q22 <- MMM.precision.marginal(LKinfo=LKinfo2,r=r)
  Q12 <- MMM.precision.cross(LKinfo1, LKinfo2, r)
  Q <- rbind(cbind(Q11,t(Q12)),cbind(Q12,Q22))
  return(Q)
}

MMM.sim <- function(x1, LKinfo1, LKinfo2, r, M=1){
  Q <- MMM.precision(LKinfo1,LKinfo2,r=r)
  Qc <- chol(Q)
  m1 <- LKinfo1$latticeInfo$m
  m2 <- LKinfo2$latticeInfo$m
  E <- matrix(rnorm(M * (m1+m2)), nrow = m1+m2, ncol = M)
  A <- as.matrix(backsolve(Qc, E))
  #A <- cbind(A[1:m1],A[(m1+1):(m1+m2)])
  PHI1 <- LKrig.basis(x1, LKinfo1) # bottleneck step
  PHI2 <- LKrig.basis(x1, LKinfo2) # bottleneck step
  process1 <- PHI1 %*% A[1:m1,]
  process2 <- PHI2 %*% A[(m1+1):(m1+m2),]
  #return(cbind(PHI %*% A[,1],PHI %*% A[,2]))
  return(list(v1=process1,v2=process2))
}

# problem here is y.c needs to be a matrix -- have yet to fix this...
MMM.lnlik <- function(p,y.c,PHI1,PHI2,LKinfo1,LKinfo2,tau1,tau2,r.decay=TRUE,look=FALSE){
  if(r.decay){
    r=p[1]*exp(-p[2] * c(0:(LKinfo1$nlevel-1)))
  }else{
    r=p
  }
  # p = r_coefficients
  if(is.vector(y.c)){
    stop("y.c must be a matrix of dimension ns x 2")
  }
  nsamp <- length(y.c[,1])
  ## Quadratic form
  Q <- MMM.precision(LKinfo1,LKinfo2,r=r)
  BIGPHItau <- rbind(cbind(PHI1/tau1,spam(0,nrow=dim(PHI1)[1],ncol=dim(PHI2)[2])),
    cbind(spam(0,nrow=dim(PHI2)[1],ncol=dim(PHI1)[2]),PHI2/tau2))

  G <- Q + t(BIGPHItau) %*% BIGPHItau
  PtDiy <- c(t(PHI1) %*% y.c[,1]/tau1^2, t(PHI2) %*% y.c[,2]/tau2^2) # Phi^T D^-1 y
  GCholesky <- chol(G) # bottleneck
  right.piece <- backsolve(GCholesky,forwardsolve(GCholesky, transpose=TRUE, PtDiy, upper.tri=TRUE))
  quad.form <- sum(c(y.c)^2/c(rep(c(tau1^2,tau2^2),each=nsamp))) - sum(PtDiy * right.piece)
  # checked numerically

  ## Log determinant
  lnDetCov <- nsamp*log(tau1^2) + nsamp*log(tau2^2) + 2*sum(log(diag(GCholesky))) -
    2*sum(log(diag(chol(Q)))) # sanity checked numerically

  ## log-likelihood
  lnlik <- -nsamp * log(2*pi) - 0.5 * lnDetCov - 0.5 * quad.form
  if(look){
    print(p)
    print(lnlik)
  }
  return(lnlik)
}

# See Temperature/CoKriging.R
MMM.cokrig <- function(p1.grd,p2.grd,o1.grd,o2.grd,y.c,LKinfo1,LKinfo2,tau1,tau2,r,M=1){
  PHI1 <- LKrig.basis(o1.grd,LKinfo1)
  PHI2 <- LKrig.basis(o2.grd,LKinfo2)
  PHI01 <- LKrig.basis(p1.grd,LKinfo1)
  PHI02 <- LKrig.basis(p2.grd,LKinfo2)
  nsamp1 <- dim(o1.grd)[1]
  nsamp2 <- dim(o2.grd)[1]
  npred1 <- dim(p1.grd)[1]
  npred2 <- dim(p2.grd)[1]

  ## Setting up multivariate PHI matrices
  BIGPHItau <- rbind(cbind(PHI1/tau1,spam(0,dim(PHI1)[1],dim(PHI2)[2])),
    cbind(spam(0,dim(PHI2)[1],dim(PHI1)[2]),PHI2/tau2)) # D^-1/2 PHI
  BIGPHI0 <- rbind(cbind(PHI01,spam(0,dim(PHI01)[1],dim(PHI02)[2])),
    cbind(spam(0,dim(PHI02)[1],dim(PHI01)[2]),PHI02)) # D^-1/2 PHI
  ## Setting up PHI^t D^-1 y
  PtDiy <- rbind(t(PHI1) %*% y.c$v1/tau1^2, t(PHI2) %*% y.c$v2/tau2^2)
  ## Quadratic form
  Q <- MMM.precision(LKinfo1,LKinfo2,r=r)
  QCholesky <- chol(Q)
  G <- Q + t(BIGPHItau) %*% BIGPHItau
  GCholesky <- chol(G) # bottleneck
  # (GiPtDiy = solve(G) %*% Phi^t D^-1 y)
  GiPtDiy <- backsolve(GCholesky,forwardsolve(GCholesky, transpose=TRUE, PtDiy, upper.tri=TRUE))
  right.piece <- t(BIGPHItau) %*% BIGPHItau %*% GiPtDiy

  ## co-kriged values
  cokrig <- BIGPHI0 %*% backsolve(QCholesky,
    forwardsolve(QCholesky, transpose=TRUE, PtDiy - right.piece, upper.tri=TRUE))
  return(list(v1=as.matrix(cokrig[1:npred1,]),v2=as.matrix(cokrig[(npred1+1):(npred1+npred2),])))
}

MMM.conditional.sim <- function(p1.grd,p2.grd,o1.grd,o2.grd,y,LKinfo1,LKinfo2,tau1,tau2,
  r,Msamp,seed=388){
    nsamp1 <- dim(o1.grd)[1]
    nsamp2 <- dim(o2.grd)[1]
    npred1 <- dim(p1.grd)[1]
    npred2 <- dim(p2.grd)[1]

    set.seed(seed)
    sim <- MMM.sim(x1=rbind(o1.grd,o2.grd,p1.grd,p2.grd),LKinfo1=LKinfo1,LKinfo2=LKinfo2,r=r,
      M=Msamp)

    sim.obs <- list()
    sim.obs$v1 <- sim$v1[1:nsamp1,] + rnorm(Msamp*nsamp1,0,sd=tau1)
    sim.obs$v2 <- sim$v2[(nsamp1+1):(nsamp1+nsamp2),] + rnorm(Msamp*nsamp2,0,sd=tau2)
    sim$v1 <- sim$v1[(nsamp1+nsamp2+1):(nsamp1+nsamp2+npred1),]
    sim$v2 <- sim$v2[(nsamp1+nsamp2+npred1+1):(nsamp1+nsamp2+npred1+npred2),]
  
    sim.ck <- MMM.cokrig(p1.grd=p1.grd,p2.grd=p2.grd,o1.grd=o1.grd,o2.grd=o2.grd,y.c=sim.obs,
      LKinfo1=LKinfo1,LKinfo2=LKinfo2,tau1=tau1,tau2=tau2,r=r,M=Msamp)
    obs.ck <- MMM.cokrig(p1.grd=p1.grd,p2.grd=p2.grd,o1.grd=o1.grd,o2.grd=o2.grd,y.c=y,
      LKinfo1=LKinfo1,LKinfo2=LKinfo2,tau1=tau1,tau2=tau2,r=r,M=1)
    return(list(v1=c(obs.ck$v1) + (sim$v1 - sim.ck$v1),v2=c(obs.ck$v2) + (sim$v2 - sim.ck$v2)))
}
