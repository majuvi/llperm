#' Starting values for (ZI)Poisson and (ZI)NegativeBinomial 
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
start_1 <- function(X, Y, Z=NULL, offsetx=0, offsetz=0, weights=1, zero.inflated=FALSE, over.dispersed=FALSE) {
  if (zero.inflated)
	model_zero <- glm.fit(Z, as.integer(Y == 0), weights = weights, family = binomial(link = "logit"), offset = offsetz)
  else
	model_zero <- list(coefficients=NULL)
  model_count <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)
  theta <- if (over.dispersed) 0 else NULL
  list(count = model_count$coefficients, zero = model_zero$coefficients, theta = theta)
}