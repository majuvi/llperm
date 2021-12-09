# Gradient
gradNegBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  Y1 <- Y > 0
  kx <- ncol(X)
  kz <- ncol(Z)
  eta <- as.vector(X %*% parms[1:kx] + offsetx)
  mu <- exp(eta)
  etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
  muz <- linkobj$linkinv(etaz)
  theta <- exp(parms[(kx + kz) + 1])
  clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
  dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +  clogdens0)
  wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta),
                       -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) -
                              log(mu + theta) + log(mu)))
  wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
                      (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
  wres_theta <- theta * ifelse(Y1, digamma(Y + theta) -
                                 digamma(theta) + log(theta) - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 -
                                                                                   muz) + clogdens0) * (log(theta) - log(mu + theta) +
                                                                                                          1 - theta/(mu + theta)))
  grad <- colSums(cbind(wres_count * weights * X, wres_zero * weights *
                          Z, wres_theta))
  return(grad)
}

