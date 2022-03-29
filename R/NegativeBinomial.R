#' Negative Binomial family
#'
#' This is part of the new implementation.
#'
#' @param count.link link function for the count component
#' @export
NegativeBinomial <- function(count.link="log") {
	count.link <- make.link(count.link)
	list(
		family = "Negative Binomial",
		count.link = count.link, 
		# Log likelihood
		loglikfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  theta <- exp(parms[kx + 1])
			  loglik <- sum(dnbinom(Y, size = theta, mu = mu, log = TRUE) * weights) #(Y * log(mu) - lfactorial(Y) + lgamma(theta + Y) - lgamma(theta) - Y * log(theta + mu) - theta * log(1 + mu/theta))
			  return(loglik)
			},
		# Gradient
		gradfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  mu.d <- count.link$mu.eta(eta)
			  theta <- exp(parms[kx + 1])
			  grad.count = (Y / mu - (Y + theta) / (mu + theta)) * mu.d
			  grad.theta = (digamma(theta + Y) - digamma(theta) - log(1 + mu / theta) + (mu - Y) / (mu + theta)) * theta
			  grad <- colSums(cbind(grad.count * weights * X, grad.theta))
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_1(X, Y, Z, offsetx, offsetz, weights, FALSE, TRUE),
		zero.inflated = FALSE,
		over.dispersed = TRUE
	)	
}