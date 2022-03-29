#' Starting values for (ZI)Binomial and (ZI)BetaBinomial
#'
#' This is part of the new implementation.
#'
#' @param X model matrix for the count component
#' @param Y response vector
#' @param Z model matrix for the zero component
#' @param offsetx offset for the count component
#' @param offsetz offset for the zero component
#' @param weights weights for each observation
#' @param zero.inflated calculate starting values for the zero component?
#' @param over.dispersed calculate starting values for the dispersion parameter?
#' @return list of starting values for each component, and possible over-dispersion
#' @export
start_2 <- function(X, Y, Z=NULL, offsetx=0, offsetz=0, weights=1, zero.inflated=FALSE, over.dispersed=FALSE) {
  #y <- as.vector(Y[,1])
  #y.c <- as.vector(rowSums(Y))
  #phat <- y/y.c
  #phat[phat <= 0] <- 1/(2 * y.c[phat == 0])
  #phat[phat >= 1] <- 1 - 1/(2 * y.c[phat == 1])
  #if (zero.inflated)
  #	model_zero <- lm.fit(Z, as.numeric(y == 0))
  #else
  #	model_zero <- list(coefficients=NULL)
  #model_count <- lm.fit(X, log(phat/(1 - phat)))
  #model.p <- mean(y) / mean(y.c)
  #phihat <- (1/(mean(y.c) - 1)) * (var(y)/(mean(y.c) * model.p * (1 - model.p)) - 1)
  #theta <- if (over.dispersed) log(phihat/(1 - phihat)) else NULL
  #list(count = model_count$coefficients, zero = model_zero$coefficients, theta = theta)
  if (zero.inflated)
	model_zero <- glm.fit(Z, as.integer(Y[,1] == 0), weights = weights, family = binomial(link = "logit"), offset = offsetz)
  else
	model_zero <- list(coefficients=NULL)
  model_count <- glm.fit(X, Y, family = binomial(link = "logit"), weights = weights, offset = offsetx)
  y.c <- as.vector(rowSums(Y))
  dispersion <- sum((model_count$weights * model_count$residuals^2)[model_count$weights > 0])/model_count$df.residual
  phihat <- min((dispersion - 1) / (mean(y.c) - 1), (mean(y.c) - 2) / (mean(y.c) - 1))
  theta <- if (over.dispersed) log(phihat/(1 - phihat)) else NULL
  list(count = model_count$coefficients, zero = model_zero$coefficients, theta = theta)
}