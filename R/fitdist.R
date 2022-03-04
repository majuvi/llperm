#' Fit an arbitrary distribution to data
#'
#' This is part of the new implementation.
#'
#' @param X model matrix for the count component
#' @param Y response vector
#' @param Z model matrix for the zero component
#' @param offsetx offset for the count component
#' @param offsetz offset for the zero component
#' @param weights weights for each observation
#' @param family specify the distribution to fit
#' @param control extra arguments for the fitting process
#' @return list of likelihood and permutation based p-values
#' @export
fitdist <- function (X, Y, Z=NULL, offsetx=NULL, offsetz=NULL, weights=NULL,
                     family, control = fitdist.control(...), ...)
{

  loglikfun <- family$loglikfun
  gradfun   <- family$gradfun
  zero.inflated  <- family$zero.inflated
  over.dispersed <- family$over.dispersed

  if (control$trace) {
    print(family$family)
  }

  n <- NROW(Y)
  kx <- NCOL(X)
  kz <- ifelse(zero.inflated, NCOL(Z), 0)

  if (is.null(weights))
    weights = rep(1, n)
  if (is.null(offsetx))
    offsetx = rep(0,n)
  if (zero.inflated & is.null(offsetz))
    offsetz = rep(0,n)

  start <- control$start
  if (is.null(start)) {
    if (control$trace)
      cat("generating starting values...")
    start <- family$startfun(X, Y, Z, offsetx, offsetz, weights)
  }

  method <- control$method
  hessian <- control$hessian
  control$method <- control$hessian <- control$start <- NULL
  if (control$trace)
    cat("calling optim() for ML estimation:\n")
  fit <- optim(fn = loglikfun, gr = gradfun, par = c(start$count, start$zero, start$theta),
               X = X, Z = Z, Y = Y, offsetx = offsetx, offsetz = offsetz, weights = weights,
               method = method, hessian = hessian, control = control)
  if (fit$convergence > 0)
    warning("optimization failed to converge")
  coefc <- fit$par[1:kx]
  names(coefc) <- names(start$count) <- colnames(X)
  if (zero.inflated) {
    coefz <- fit$par[(kx + 1):(kx + kz)]
    names(coefz) <- names(start$zero) <- colnames(Z)
  } else {
    coefz <- NULL
  }
  if (over.dispersed)
    theta <- fit$par[(kx + kz + 1)]
  else
    theta <- NULL
  vc <- tryCatch(-solve(as.matrix(fit$hessian)), error = function(e) {
    warning(e$message, call = FALSE)
    k <- nrow(as.matrix(fit$hessian))
    return(matrix(NA, k, k))
  })
  cols.count <- paste("count", colnames(X), sep = "_")
  cols.zero <- if (zero.inflated) paste("zero", colnames(Z), sep = "_") else NULL
  cols.theta <- if (over.dispersed) "theta" else NULL
  colnames(vc) <- rownames(vc) <- c(cols.count, cols.zero, cols.theta)
  mu <- family$count.link$linkinv(X %*% coefc + offsetx)[, 1]
  if (zero.inflated) {
    phi <- family$zero.link$linkinv(Z %*% coefz + offsetz)[, 1]
    Yhat <- (1 - phi) * mu 
  } else {
    Yhat <- mu
  }
  res <- sqrt(weights) * (Y - Yhat)
  nobs <- sum(weights > 0)
  rval <- list(
    coefficients = list(count = coefc, zero = coefz, theta=theta), vcov = vc,
    residuals = res, fitted.values = Yhat, dist = family$family, loglik = fit$value, n = nobs,
    df.null = nobs - 2, df.residual = nobs - (kx + kz + over.dispersed),
    converged = fit$convergence < 1,
    optim = fit, start = start,
    weights = if (identical(weights, rep.int(1, n))) NULL else weights,
    offset = list(count = if (identical(offsetx, rep.int(0, n))) NULL else offsetx,
                  zero = if (identical(offsetz, rep.int(0, n))) NULL else offsetz),
    contrasts = list(count = attr(X, "contrasts"),
                     zero = if(zero.inflated) attr(Z, "contrasts") else NULL))
  return(rval)
}