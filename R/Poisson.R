Poisson <- function(count.link="log") {
	count.link <- make.link(count.link)
	list(
		family = "Poisson",
		count.link = count.link, 
		# Log likelihood
		loglikfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  loglik <- sum(dpois(Y, lambda = mu, log = TRUE) * weights) #(Y * log(mu) - mu - lfactorial(Y))
			  return(loglik)
			},
		# Gradient
		gradfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  mu.d <- count.link$mu.eta(eta)
			  grad = colSums((Y/mu - 1) * mu.d * weights * X)
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_1(X, Y, Z, offsetx, offsetz, weights, FALSE, FALSE),
		zero.inflated = FALSE,
		over.dispersed = FALSE
	)	
}