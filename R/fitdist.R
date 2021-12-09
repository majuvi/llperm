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
#' @param dist specify the distribution
#' @param control extra arguments for the fitting process
#' @param linkobj link object to use for the zero component
#' @return list of likelihood and permutation based p-values
#' @export
fitdist <- function (X, Y, Z=NULL, offsetx=NULL, offsetz=NULL, weights=NULL,
                     dist = c("negbin", "betabin", "zipoisson", "zinegbin", "zibinomial", "zibetabin"),
                     control = fitdist.control(...), 
					 linkobj = make.link("logit"), ...)
{
  dist <- match.arg(dist)
  #TODO: implement gradient for binomial / beta-binomial
  loglikfun <- switch(dist, negbin=NegBin, betabin=BetaBin, zipoisson = ziPoisson, zinegbin = ziNegBin, zibinomial=ziBinomial, zibetabin = ziBetaBin)
  gradfun <- switch(dist, negbin=NULL, betabin=NULL, zipoisson = gradPoisson, zinegbin = gradNegBin, zibinomial=NULL, zibetabin = NULL)

  zero.inflated <- startsWith(dist, "zi")
  over.disperse <- dist %in% c("negbin", "betabin", "zinegbin", "zibetabin")

  if (control$trace) {
    if (zero.inflated)
      cat("Zero-inflated Count Model\n", dist, "\n", sep = "")
    else
      cat("Count Model\n", dist, "\n", sep = "")
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
    if (dist %in% c("negbin", "zipoisson", "zinegbin")) {
      if (zero.inflated)
        model_zero <- glm.fit(Z, as.integer(Y == 0), weights = weights, family = binomial(link = "logit"), offset = offsetz)
      else
        model_zero <- list(coefficients=NULL)
      model_count <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)
      theta <- 1
    } else if (dist %in% c("betabin", "zibinomial", "zibetabin")) {
      y <- as.vector(Y[,1])
      y.c <- as.vector(rowSums(Y))
      phat <- y/y.c
      phat[phat <= 0] <- 1/(2 * y.c[phat == 0])
      phat[phat >= 1] <- 1 - 1/(2 * y.c[phat == 1])
      if (zero.inflated)
        model_zero <- lm.fit(Z, as.numeric(y == 0))
      else
        model_zero <- list(coefficients=NULL)
      model_count <- lm.fit(X, log(phat/(1 - phat)))
      model.p <- mean(y) / mean(y.c)
      phihat <- (1/(mean(y.c) - 1)) * (var(y)/(mean(y.c) * model.p * (1 - model.p)) - 1)
      theta <- phihat/(1 - phihat)
    }
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients,
                  theta = if (over.disperse) log(theta) else NULL)
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
  if (over.disperse)
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
  cols.theta <- if (over.disperse) "theta" else NULL
  colnames(vc) <- rownames(vc) <- c(cols.count, cols.zero, cols.theta)
  mu <- exp(X %*% coefc + offsetx)[, 1]
  if (zero.inflated) {
    phi <- linkobj$linkinv(Z %*% coefz + offsetz)[, 1]
    Yhat <- (1 - phi) * mu
  } else {
    Yhat <- mu
  }
  res <- sqrt(weights) * (Y - Yhat)
  nobs <- sum(weights > 0)
  rval <- list(
    coefficients = list(count = coefc, zero = coefz, theta=theta), vcov = vc,
    residuals = res, fitted.values = Yhat, dist = dist, loglik = fit$value, n = nobs,
    df.null = nobs - 2, df.residual = nobs - (kx + kz + (dist %in% c("negbin", "betabin"))),
    converged = fit$convergence < 1,
    optim = fit, start = start,
    weights = if (identical(weights, rep.int(1, n))) NULL else weights,
    offset = list(count = if (identical(offsetx, rep.int(0, n))) NULL else offsetx,
                  zero = if (identical(offsetz, rep.int(0, n))) NULL else offsetz),
    contrasts = list(count = attr(X, "contrasts"),
                     zero = if(zero.inflated) attr(Z, "contrasts") else NULL))
  #class(rval) <- "zeroinfl"
  return(rval)
}
