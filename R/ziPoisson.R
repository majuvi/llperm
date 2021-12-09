#Poisson zero-inflated likelihood P(N=n)
ziPoisson <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y0 <- Y <= 0
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
  phi <- as.vector(linkobj$linkinv(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  loglik0 <- log(phi + exp(log(1 - phi) - mu))
  loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
  loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
  return(loglik)
}
