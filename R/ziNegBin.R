#Negative Binomial zero-inflated likelihood P(N=n)
ziNegBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y0 <- Y <= 0
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkobj$linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  theta <- exp(parms[(kx + kz) + 1])
  loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0, size = theta, mu = mu, log = TRUE))))
  loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
  loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
  return(loglik)
}
