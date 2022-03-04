ZIBinomial <- function(count.link="logit", zero.link="logit") {
	count.link <- make.link(count.link)
	zero.link <- make.link(zero.link)
	list(
		family = "Zero-Inflated Binomial",
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
			  loglik.0 <- log(phi + (1 - phi) * dbinom(0, y.c, mu)) #(1-mu)^y.c
			  loglik.1 <- log(1 - phi) + dbinom(y, y.c, mu, log = TRUE) #lchoose(y.c, y) +  y*log(mu) + (y.c - y)*log(1 - mu)
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
			  count.0 <- dbinom(0, y.c, mu)#(1-mu)^y.c
			  likelihood.0 <- phi + (1 - phi) * count.0
			  grad.count.0 <- -(1 - phi) * y.c * (1 - mu)^(y.c - 1) * mu.d / likelihood.0
			  grad.count.1 <- (y - y.c * mu) / (mu * (1 - mu)) * mu.d
			  grad.count   <- ifelse(Y1, grad.count.1, grad.count.0)
			  grad.zero.0  <- (1 - count.0) * phi.d / likelihood.0
			  grad.zero.1  <- -1/(1 - phi) * phi.d
			  grad.zero    <- ifelse(Y1, grad.zero.1, grad.zero.0)
			  grad <- colSums(cbind(grad.count * weights * X, grad.zero * weights * Z))
			  return(grad)
			},
		startfun = function(X, Y, Z, offsetx, offsetz, weights) start_2(X, Y, Z, offsetx, offsetz, weights, TRUE, FALSE),
		zero.inflated = TRUE,
		over.dispersed = FALSE
	)	
}
