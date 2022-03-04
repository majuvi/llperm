
ZIBetaBinomial <- function(count.link="logit", zero.link="logit") {
	count.link <- make.link(count.link)
	zero.link <- make.link(zero.link)
	llbbinom <- function(y, y.c, a, b) lchoose(y.c, y) +  lbeta(y + a, y.c - y + b) - lbeta(a, b)
	#lgamma(y.c + 1) - lgamma(y + 1) - lgamma(y.c - y + 1) + lgamma(y + a) + lgamma(y.c - y + b) - lgamma(y.c + a + b) + lgamma(a + b) - lgamma(a) - lgamma(b)
	list(
		family = "Zero-Inflated Negative Binomial",
		count.link = count.link, 
		zero.link = zero.link,
		# Log likelihood
		loglikfun = function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
			  y <- as.vector(Y[,1])
			  y.c <- as.vector(rowSums(Y))
			  Y1 <- y > 0
			  kx <- ncol(X)
			  kz <- ncol(Z)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
			  phi <- zero.link$linkinv(etaz)
			  theta <- 1/(1+exp(-parms[(kx + kz) + 1]))
			  a <- mu * (1 - theta) / theta
			  b <- (1 - mu) * (1 - theta) / theta
			  loglik.0 <- log(phi + (1 - phi) * exp(llbbinom(0, y.c, a, b)))
			  loglik.1 <- log(1 - phi) + llbbinom(y, y.c, a, b)
			  loglik <- sum(ifelse(Y1, loglik.1, loglik.0) * weights)
			  return(loglik)
			},

		# Gradient
		gradfun = function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
			  y <- as.vector(Y[,1])
			  y.c <- as.vector(rowSums(Y))
			  Y1 <- y > 0
			  kx <- ncol(X)
			  kz <- ncol(Z)
			  eta <- as.vector(X %*% parms[1:kx] + offsetx)
			  mu <- count.link$linkinv(eta)
			  mu.d <- count.link$mu.eta(eta)
			  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
			  phi <- zero.link$linkinv(etaz)
			  phi.d <- zero.link$mu.eta(etaz)
			  theta <- 1/(1+exp(-parms[(kx + kz) + 1]))
			  a <- mu * (1 - theta) / theta
			  b <- (1 - mu) * (1 - theta) / theta
			  a.d <- mu.d * (1 - theta) / theta
			  
			  count.0 <- exp(llbbinom(0, y.c, a, b))
			  likelihood.0  <- phi + (1 - phi) * count.0
			  f0.count <- count.0 * (psigamma(y.c + b) - psigamma(b))
			  f0.theta <- count.0 * (psigamma(y.c + b) * (1-mu) + psigamma(a + b) - psigamma(y.c + a + b) - psigamma(b)*(1-mu))
			  grad.count.0 <- - (1 - phi) *f0.count / likelihood.0 * a.d 
			  grad.count.1 <- (digamma(y + a) - digamma(y.c - y + b) - digamma(a) + digamma(b)) * a.d
			  grad.count   <- ifelse(Y1, grad.count.1, grad.count.0)
			  grad.zero.0  <- (1 - count.0)/ likelihood.0  * phi.d 
			  grad.zero.1  <- -1/(1 - phi) * phi.d
			  grad.zero    <- ifelse(Y1, grad.zero.1, grad.zero.0)			  
			  grad.theta.0 <- -(1 - theta)/theta * (1 - phi) * f0.theta / likelihood.0  
			  grad.theta.1 <- -(1 - theta)/theta * (digamma(y + a)*mu + digamma(y.c - y + b)*(1-mu) - digamma(y.c + a + b) + digamma(a + b) - digamma(a)*mu - digamma(b)*(1-mu))
			  grad.theta   <- ifelse(Y1, grad.theta.1, grad.theta.0) 
			  grad <- colSums(cbind(grad.count * weights * X, grad.zero * weights * Z, grad.theta)) 
			  
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_2(X, Y, Z, offsetx, offsetz, weights, TRUE, TRUE),
		zero.inflated = TRUE,
		over.dispersed = TRUE
	)	
}