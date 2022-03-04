Binomial <- function(count.link="logit") {
	count.link <- make.link(count.link)
	list(
		family = "Binomial",
		count.link = count.link, 
		# Log likelihood
		loglikfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  y <- as.vector(Y[,1])
			  y.c <- as.vector(rowSums(Y))
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  loglik <- sum(dbinom(y, y.c, mu, log = TRUE) * weights) #(lchoose(y.c, y) +  y*log(mu) + (y.c - y)*log(1 - mu)
			  return(loglik)
			},
		# Gradient
		gradfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  y <- as.vector(Y[,1])
			  y.c <- as.vector(rowSums(Y))
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  mu.d <- count.link$mu.eta(eta)
			  grad = colSums((y - y.c * mu) / (mu * (1 - mu)) * mu.d * weights * X)
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_2(X, Y, Z, offsetx, offsetz, weights, FALSE, FALSE),
		zero.inflated = FALSE,
		over.dispersed = FALSE
	)	
}