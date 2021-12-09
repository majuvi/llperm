# Binomial zero-inflated likelihood P(N=n)
ziBinomial <- function(parms, X, Y, Z, offsetx=0, offsetz=0, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  kz <- NCOL(Z)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  phi <- linkobj$linkinv(as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
  Y0 <- y == 0
  loglike0 <- log(phi + exp(log(1 - phi) + y.c*log(1 - mu)))
  loglike1 <- log(1 - phi) + lchoose(y.c, y) +  y*log(mu) + (y.c - y)*log(1 - mu)
  loglike <- sum(weights[Y0] * loglike0[Y0]) + sum(weights[!Y0] * loglike1[!Y0])
  return(loglike)
}

