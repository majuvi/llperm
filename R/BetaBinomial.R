BetaBinomial <- function(count.link="logit") {
	count.link <- make.link(count.link)
	llbbinom <- function(y, y.c, a, b) lchoose(y.c, y) +  lbeta(y + a, y.c - y + b) - lbeta(a, b)
	#lgamma(y.c + 1) - lgamma(y + 1) - lgamma(y.c - y + 1) + lgamma(y + a) + lgamma(y.c - y + b) - lgamma(y.c + a + b) + lgamma(a + b) - lgamma(a) - lgamma(b)
	list(
		family = "Beta Binomial",
		count.link = count.link, 
		# Log likelihood
		loglikfun = function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
			  y <- as.vector(Y[,1])
			  y.c <- as.vector(rowSums(Y))
			  kx <- ncol(X)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  theta <- 1/(1+exp(-parms[kx + 1]))
			  a <- mu * (1 - theta) / theta
			  b <- (1 - mu) * (1 - theta) / theta
			  loglik <- sum(llbbinom(y, y.c, a, b) * weights)
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
			  theta <- 1/(1+exp(-parms[kx + 1]))
			  a <- mu * (1 - theta) / theta
			  b <- (1 - mu) * (1 - theta) / theta
			  a.d <- mu.d * (1 - theta) / theta
			  grad.count <- (digamma(y + a) - digamma(y.c - y + b) - digamma(a) + digamma(b)) * a.d
			  grad.theta <- -(1-theta)/theta * (digamma(y + a)*mu + digamma(y.c - y + b)*(1-mu) - digamma(y.c + a + b) + digamma(a + b) - digamma(a)*mu - digamma(b)*(1-mu))
			  grad <- colSums(cbind(grad.count * weights * X, grad.theta))
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_2(X, Y, Z, offsetx, offsetz, weights, FALSE, TRUE),
		zero.inflated = FALSE,
		over.dispersed = TRUE
	)	
}