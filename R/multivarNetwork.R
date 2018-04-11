#' @title Fit a collection of networks from multivariate experiments
#'
#' @description The inference procedure is multivariate neighborhood selection with group Lasso penalty.
#' If an univariate experiment is provided, the usual neighborhood selection with lasso penalty applies.
#'
#' @param X a list of matrices with the same dimension, corresponding to multiple experiments related
#' to the same individuals and variables, e.g., proteomics and transcriptomics data.
#' @param sym.rule a character, either "AND" or "OR" to define the  post-symmetrization rule applied to the coefficients in neighborhood selection
#' @param select a character for defining the cross-validation rule used to choose the penalty level (i.e., the number of edges in the the network).
#' If "none", the whole regulariazation paths is sent back. If "min" or "1se", penalties are cross-validation for each variable
#' with the corresponding rule.
#' @param nlambda integer defining the number of penalty levels used on the grid
#' @param min.ratio the grid of penalties starts by lambda.max, that is,  minimal penalty level for selecting no edge in the network, then decrease on a
#' log scale until min.ratio* lambda.max.
#' @param mc.cores for distributing the computation over the variables. Default is 1 (no parallelization).
#'
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply
#' @import Matrix
#' @import glmnet
#' @import gglasso
#' @import stabsel
#' @export
multivarNetwork <- function(X, sym.rule="AND", select=c("none", "1se", "min" ,git"stabsel"),
                            nlambda=50, min.ratio=1e-3, mc.cores=1, nfold=5, cutoff=0.75) {

  select <- match.arg(select)
  LAPPLY <- ifelse(mc.cores > 1, pbmclapply, mclapply)

  stopifnot(is.list(X) | is.matrix(X))
  if (is.matrix(X)) X <- list(X)

  ## scaling to make everything comparable
  X <- lapply(X, scale)

  ## problem dimension
  K <- length(X)
  p <- ncol(X[[1]])
  n <- nrow(X[[1]])

  ## define group indexes
  group <- rep(rep(1:(p-1),K), K)

  ## create the list of data sets corresponding to each variable
  training_sets <- lapply(1:p, function(j) {
    y <- Reduce("cbind", lapply(X, function(x) x[, j]))
    x <- Reduce("cbind", lapply(X, function(x) x[,-j]))
    list(
      Y = as.vector(y),
      X = bdiag(rep(list(x), K)),
      j = j
    )
  })

  ## Get the largest lambda value across all variable that produce a null vector as estimate
  lambda.max <- max(sapply(training_sets, function(data) {
    return(max(sqrt(rowsum(as.numeric(crossprod(data$X, data$Y)^2), group)/(2*n))))
  }))

  ## Compute an common grid of penalty levels
  lambda <- 10^seq(from=log10(lambda.max),log10(min.ratio*lambda.max),len=nlambda)

  ## extrac the appropriate inference function taking
  f_infer <- neighborhood.selection(K, select, group, lambda)

  ## compute all coefficients
  coefficients <- do.call(rbind,LAPPLY(training_sets, function(d) {
    coef <- f_infer(d$X, d$Y)
    beta  <- matrix(0, p, ncol(coef))
    beta[-d$j, ] <- as.matrix(coef)
    beta
  }, mc.cores=mc.cores))

  networks <- apply(coefficients, 2, function(beta) {
    net <- matrix(beta != 0,p,p)
    net <- Matrix(switch(sym.rule,
                         "AND" = pmin(net,t(net)),
                         "OR"  = pmax(net,t(net))))
    net
  })

  return(list(coefficients = coefficients, networks = networks, penalties = lambda))
}

neighborhood.selection <- function(K, select, group, lambda) {
  ## return a function that take x and y as argument with fixed lambda and group
  if (K == 1) {
    ## single attribute: lasso based learning
    f <- switch(select,
      "none"     = lasso_path(lambda),
      "stabsel"  = lasso_stabsel(lambda), ## add cutoff and PFER
      # default: crossval with min or 1se rule
      lasso_crossval(lambda, select)
    )
  } else {
  ## multiattribute: group-lasso based learning
    f <- switch(select,
   "none"     = group_lasso_path(lambda, group),
   "stabsel"  =
      function(x, y) {
        NULL # TODO
      },
   # default: crossval with min or 1se rule
   group_lasso_crossval(lambda, group, select)
    )
  }
  f
}

lasso_path <- function(lambda) {
  function(x, y) {
    out <- glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- coef.glmnet(out)[-1,]
    abs(beta)
  }
}

lasso_crossval <- function(lambda, select) {
  function(x, y) {
    out <- cv.glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- matrix(coef.cv.glmnet(out, s=paste("lambda",select, sep="."))[-1], ncol=1)
    abs(beta)
  }
}

lasso_stabsel <- function(lambda) {
  function(x, y) {
    out <- stabsel.matrix(x, y, fitfun = glmnet.lasso, args.fitfun = list(lambda = lambda, intercept=FALSE, standardize=FALSE), cutoff = 0.75, PFER = 1)
    beta <- rep(0,ncol(x))
    beta[stab.glmnet$selected] <- 1
    beta
  }
}

group_lasso_path <- function(lambda, group) {
  function(x, y) {
    o <- order(group)
    x <- as.matrix(x)[, o]
    out <- gglasso(x, y, group[o], lambda=lambda, intercept=FALSE)
    beta <- coef(out)[-1,]
    beta <- apply(beta[order(o), , drop=FALSE], 2, function(b) sqrt(rowsum(b^2,group)))
    beta
  }
}

group_lasso_crossval <- function(lambda, group, select) {
  function(x, y) {
    o <- order(group)
    x <- as.matrix(x)[, o]
    out <- cv.gglasso(x, y, group[o], lambda=lambda, intercept=FALSE)
    beta <- matrix(coef(out, s=paste("lambda",select, sep="."))[-1], ncol=1)
    beta <- apply(beta[order(o), , drop=FALSE], 2, function(b) sqrt(rowsum(b^2,group)))
    beta
  }

}

lasso <- function(x, y, group, lambda, cv.choice) {
  if (cv.choice == "none") {
    out <- glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- coef.glmnet(out)[-1,]
  } else {
    out <- cv.glmnet(x, y, lambda=lambda, intercept=FALSE, standardize=FALSE)
    beta <- matrix(coef.cv.glmnet(out, s=paste("lambda",cv.choice, sep="."))[-1], ncol=1)
  }
  abs(beta)
}

grplasso <- function(x, y, group, lambda, cv.choice) {
  o <- order(group)
  x <- as.matrix(x)[, o]
  if (cv.choice == "none") {
    out <- gglasso(x, y, group[o], lambda=lambda, intercept=FALSE)
    beta <- coef(out)[-1,]
  } else {
    out <- cv.gglasso(x, y, group[o], lambda=lambda, intercept=FALSE)
    beta <- matrix(coef(out, s=paste("lambda",cv.choice, sep="."))[-1], ncol=1)
  }
  beta <- apply(beta[order(o), , drop=FALSE], 2, function(b) sqrt(rowsum(b^2,group)))
}

