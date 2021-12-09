# Beta-Binomial zero-inflated likelihood P(N=n)
ziBetaBin <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  kz <- NCOL(Z)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  phi <- linkobj$linkinv(as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  theta <- linkobj$linkinv(parms[kx + kz + 1])
  s1 <- as.vector(mu * (1 - theta)/theta)
  s2 <- as.vector((1 - theta)/theta * (1 - mu))
  Y0 <- y == 0
  loglike0 <- log(phi + exp(log(1 - phi) + lbeta(s1, y.c + s2) - lbeta(s1, s2)))
  loglike1 <- log(1 - phi) + lchoose(y.c, y) +  lbeta(y + s1, y.c - y + s2) - lbeta(s1, s2)
  loglike <- sum(weights[Y0] * loglike0[Y0]) + sum(weights[!Y0] * loglike1[!Y0])
  return(loglike)
}
