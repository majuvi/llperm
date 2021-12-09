#Negative Binomial likelihood P(N=n)
NegBin <- function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
  kx <- ncol(X)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  theta <- exp(parms[kx + 1])
  loglik <- sum(weights * suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE)))
  return(loglik)
}

