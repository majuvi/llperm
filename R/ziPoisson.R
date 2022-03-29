#' Zero-Inflated Poisson family
#'
#' This is part of the new implementation.
#'
#' @param count.link link function for the count component
#' @param zero.link link function for the zero component
#' @export
ZIPoisson <- function(count.link="log", zero.link="logit") {
	count.link <- make.link(count.link)
	zero.link <- make.link(zero.link)
	list(
		family = "Zero-Inflated Poisson",
		count.link = count.link, 
		zero.link = zero.link,
		# Log likelihood
		loglikfun = function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
			  Y1 <- Y > 0
			  kx <- ncol(X)
			  kz <- ncol(Z)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
			  phi <- zero.link$linkinv(etaz)
			  loglik.0 <- log(phi + (1 - phi) * dpois(0, lambda = mu)) #exp(-mu)
			  loglik.1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE) #Y * log(mu) - mu - lfactorial(Y)
			  loglik <- sum(ifelse(Y1, loglik.1, loglik.0) * weights)
			  return(loglik)
			},
		# Gradient
		gradfun = function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
			  Y1 <- Y > 0
			  kx <- ncol(X)
			  kz <- ncol(Z)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  mu.d <- count.link$mu.eta(eta)
			  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
			  phi <- zero.link$linkinv(etaz)
			  phi.d <- zero.link$mu.eta(etaz)
			  likelihood.0 <- phi + (1 - phi) * dpois(0, lambda = mu) #exp(-mu)
			  grad.count.0 <- -(1 - phi) * exp(-mu) * mu.d / likelihood.0
			  grad.count.1 <- (Y/mu - 1) * mu.d
			  grad.count   <- ifelse(Y1, grad.count.1, grad.count.0)
			  grad.zero.0  <- (1 - exp(-mu)) * phi.d / likelihood.0
			  grad.zero.1  <- -1/(1 - phi) * phi.d
			  grad.zero    <- ifelse(Y1, grad.zero.1, grad.zero.0)
			  grad <- colSums(cbind(grad.count * weights * X, grad.zero * weights * Z))
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_1(X, Y, Z, offsetx, offsetz, weights, TRUE, FALSE),
		zero.inflated = TRUE,
		over.dispersed = FALSE
	)	
}