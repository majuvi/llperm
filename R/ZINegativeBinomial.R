ZINegativeBinomial <- function(count.link="log", zero.link="logit") {
	count.link <- make.link(count.link)
	zero.link <- make.link(zero.link)
	list(
		family = "Zero-Inflated Negative Binomial",
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
			  theta <- exp(parms[(kx + kz) + 1])
			  loglik.0 <- log(phi + (1 - phi) * dnbinom(0, size = theta, mu = mu))
			  loglik.1 <- log(1 - phi) + dnbinom(Y, size = theta, mu = mu, log = TRUE)
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
			  theta <- exp(parms[(kx + kz) + 1])
			  count.0 <- dnbinom(0, size = theta, mu = mu)
			  likelihood.0 <- phi + (1 - phi) * count.0
			  grad.count.0 <- -(1 - phi) / (1 + mu / theta) * count.0 / likelihood.0 * mu.d 
			  grad.count.1 <- (Y / mu - (Y + theta) / (mu + theta)) * mu.d
			  grad.count   <- ifelse(Y1, grad.count.1, grad.count.0)
			  grad.zero.0  <- (1 - count.0)/ likelihood.0  * phi.d 
			  grad.zero.1  <- -1/(1 - phi) * phi.d
			  grad.zero    <- ifelse(Y1, grad.zero.1, grad.zero.0)
			  grad.theta.0 <- ((1 - phi) * (log(theta/(mu + theta)) + mu/(mu + theta)) * count.0 / likelihood.0) * theta 
			  grad.theta.1 <- (digamma(theta + Y) - digamma(theta) - log(1 + mu / theta) + (mu - Y) / (mu + theta)) * theta
			  grad.theta   <- ifelse(Y1, grad.theta.1, grad.theta.0)
			  grad <- colSums(cbind(grad.count * weights * X, grad.zero * weights * Z, grad.theta)) 
			  return(grad)
			    wres_theta <- theta * ifelse(Y1, digamma(Y + theta) -
                                 digamma(theta) + log(theta) - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 -
                                                                                   muz) + clogdens0) * (log(theta) - log(mu + theta) +
                                                                                                          1 - theta/(mu + theta)))
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_1(X, Y, Z, offsetx, offsetz, weights, TRUE, TRUE),
		zero.inflated = TRUE,
		over.dispersed = TRUE
	)	
}