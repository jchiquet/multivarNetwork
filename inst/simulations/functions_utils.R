library(pbmcapply)
library(mvtnorm)

rgraph <- function(p,nedges) {

  A <- matrix(0,p,p)
  A[sample(which(upper.tri(A)), nedges)] <- 1

  return(A + t(A))
}

rmultivar <- function(n, K, G, agreement = c("full", "indep", "anti"), hardness = 1, eps=.Machine$double.eps, scale=TRUE) {

  S <- switch(match.arg(agreement),
         'full'  = matrix(1,K,K),
         'indep' = diag(K),
         'anti'  = matrix(1,K,K) - diag(K)
  )
  p <- ncol(G) # number of nodes per graph

  ## split the adjacency matrix to multivariate space with the kronecker product
  A <- G %x% S + diag(p*K)

  ## Remove axes from with null eigen values to make it invertible
  eigenA <- eigen(A)
  D <- eigenA$values
  U <- eigenA$vectors
  D[D < eps] <- eps
  B <- U %*% diag(D) %*% t(U)

  ## The concentration matrix is built from B
  ## with hardness manage through diagonal smoothing

  ## Then, generate the associate multivariate Gaussian sample
  X <- rmvnorm(n,sigma=solve(B + hardness * diag(diag(B))), method="chol")

  Xs <- sapply(1:K,function(k) {X[,seq(k,(K*p),by=K)]},simplify=FALSE)

  return(list(X = lapply(Xs, scale,scale,scale), A = A, S=S))
}

perf.roc <- function(theta.hat, theta.star) {

  roc <- function(theta) {
    theta <- as.matrix(theta)
    nzero <- which(theta != 0)
    zero  <- which(theta == 0)

    true.nzero <- which(theta.star != 0)
    true.zero  <- which(theta.star == 0)

    TP <- sum(nzero %in% true.nzero)
    TN <- sum(zero %in%  true.zero)
    FP <- sum(nzero %in% true.zero)
    FN <- sum(zero %in%  true.nzero)

    recall    <- TP/(TP+FN) ## also recall and sensitivity
    fallout   <- FP/(FP+TN) ## also 1 - specificit

    res <-  round(c(fallout,recall),3)
    res[is.nan(res)] <- 0
    names(res) <- c("fallout","recall")
    return(res)
  }

  if (is.list(theta.hat)) {
    return(as.data.frame(do.call(rbind, lapply(theta.hat, roc))))
  } else {
    return(roc(theta.hat))
  }
}

perf.auc <- function(roc) {
  fallout <- c(0,roc$fallout,1)
  recall  <- c(0,roc$recall, 1)
  dx <- diff(fallout)
  return(sum(c(recall[-1]*dx, recall[-length(recall)]*dx))/2)
}
