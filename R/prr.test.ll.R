#' New Permutation of Regression Residuals (PRR) test
#'
#' This is part of the new implementation.
#'
#' @param formula Formula that specifies the model
#' @param data Data frame with the data
#' @param family Specify the distribution to fit
#' @param subset Take a given subset of the observations
#' @param weights Weights to use for the observations
#' @param na.action what to do with NA values
#' @param offset offset to use for each observation
#' @param test.var variable in the formula to permute
#' @param test specify whether to test "both"/"count"/"zero"
#' @param n.rep number of permutation repetitions
#' @param control fitdist.control object to adjust the optimization
#' @param xyz return implicit model matrices T/F
#' @return list of likelihood and permutation based p-values
#' @export
prr.test.ll <- function (formula, data, family=Poisson, subset, weights, na.action, offset,
                                   test.var = NULL, test = "both", n.rep=100, control=fitdist.control(...),
                                   xyz=F, ...)
{
  
  if (is.character(family))
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
     family <- family()
  if (is.null(family$family)) {
     print(family)
     stop("'family' not recognized")
  }
  
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset",
               "na.action", "weights", "offset"),
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  # Extract formula y ~ a + b + c | d + e into ffc = y ~ a + b + c, ffz = ~ d + e
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
    ff <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
    formula[[3]][1] <- call("+")
    mf$formula <- formula
  }
  else {
    ffc <- ff <- formula
    ffz <- NULL
    #ffz[[2]] <- NULL
  }

  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  # Response
  Y <- model.response(mf, "numeric")
  # Weights
  n <- NROW(Y)
  weights <- model.weights(mf)
  if (is.null(weights))
    weights <- 1
  if (length(weights) == 1)
    weights <- rep.int(weights, n)
  weights <- as.vector(weights)
  names(weights) <- rownames(mf)
  # X: count model matrix
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
  if (is.null(offsetx))
    offsetx <- 0
  if (length(offsetx) == 1)
    offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  # Which columns in X are related to test.var?
  X.term.labels  <- attr(mtX, "term.labels")
  X.term.columns <- attr(X, "assign")
  if (attr(mtX, "intercept")) {
    X.term.labels  <- c("(Intercept)", X.term.labels)
    X.term.columns <- X.term.columns + 1
  }

  # Z: zero inflation model matrix
  if(!is.null(ffz)) {
    mtZ <- terms(ffz, data = data)
    mtZ <- terms(update(mtZ, ~.), data = data)
    Z <- model.matrix(mtZ, mf)
    offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
    if (is.null(offsetz))
      offsetz <- 0
    if (length(offsetz) == 1)
      offsetz <- rep.int(offsetz, n)
    offsetz <- as.vector(offsetz)
    # Which columns in X are related to test.var?
    Z.term.labels  <- attr(mtZ, "term.labels")
    Z.term.columns <- attr(Z, "assign")
    if (attr(mtZ, "intercept")) {
      Z.term.labels  <- c("(Intercept)", Z.term.labels)
      Z.term.columns <- Z.term.columns + 1
    }
  } else {
    Z <- NULL
    offsetz <- NULL
  }
  
  X.cols.testvar = X.term.labels[X.term.columns] == test.var
  Z.cols.testvar = if (is.null(ffz)) NULL else Z.term.labels[Z.term.columns] == test.var

  # Check if we should only test the "count" or "zero" component
  if (test == "count")
    Z.cols.testvar = NULL
  if (test == "zero")
    X.cols.testvar = NULL

  # Required count component
  if (!is.null(X.cols.testvar)) {

    # Predict test.var from other variables
    fit.X <- lm.fit(X[,!X.cols.testvar, drop=F], X[,X.cols.testvar, drop=F])
    R.X   <- fit.X$residual
    if (is.null(dim(R.X)) | (length(dim(R.X)) == 1L))
      dim(R.X) <- c(length(R.X), 1)

    # Model matrix without variable of interest
    # and with residuals of variable of interest
    X0 <- X[,!X.cols.testvar, drop=F]
    XR <- cbind(X0, R.X)
  } else {
    XR <- X0 <- X
  }

  # Possible zero-inflation component
  if (!is.null(Z.cols.testvar)) {
    # Predict test.var from other variables
    fit.Z <- lm.fit(Z[,!Z.cols.testvar, drop=F], Z[,Z.cols.testvar, drop=F])
    R.Z   <- fit.Z$residual
    if (is.null(dim(R.Z)) | (length(dim(R.Z)) == 1L))
      dim(R.Z) <- c(length(R.Z), 1)

    # Model matrix without variable of interest
    # and with residuals of variable of interest
    Z0 <- Z[,!Z.cols.testvar, drop=F]
    ZR <- cbind(Z0, R.Z)
  } else {
    ZR <- Z0 <- Z
  }
  
  if (xyz) # return only the model matrix
    return(list(X=X, Z=Z, Y=Y, XR=XR, X0=X0, ZR=ZR, Z0=Z0, offsetx=offsetx, offsetz=offsetz, weights=weights))

  # Calculate observed p-value
  fit.1 <- fitdist(XR, Y, ZR, offsetx, offsetz, weights, family=family, control=control)
  fit.2 <- fitdist(X0, Y, Z0, offsetx, offsetz, weights, family=family, control=control)
  ll.1    <- fit.1$loglik
  ll.2    <- fit.2$loglik
  lr.df <- fit.2$df.residual - fit.1$df.residual
  p.value.obs <- pchisq(2 * (ll.1 - ll.2), lr.df, lower.tail=F)

  # Calculate simulated p-value
  ll.r   <- rep(0, length(n.rep))
  for (i in 1:n.rep) {
    rnd <- sample(1:nrow(X), replace = FALSE)
    if (!is.null(X.cols.testvar))
      XR <- cbind(X0, R.X[rnd,])
    if (!is.null(Z.cols.testvar))
      ZR <- cbind(Z0, R.Z[rnd,])
    fit.r   <- fitdist(XR, Y, ZR, offsetx, offsetz, weights, family=family, control=control)
    ll.r[i]   <- fit.r$loglik
  }
  p.value.sim   <- mean(ll.1 < ll.r)

  results <- list(fit=fit.1, p.value.obs=p.value.obs, p.value.sim=p.value.sim)
  class(results) <- "prr.test.ll"
  return(results)
}
