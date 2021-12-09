# Beta-Binomial likelihood P(N=n)
BetaBin <- function(parms, X, Y, Z=NULL, offsetx=0, offsetz=NULL, weights=1) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  kx <- NCOL(X)
  mu <- linkobj$linkinv(as.vector(X %*% parms[1:kx] + offsetx))
  theta <- linkobj$linkinv(parms[kx + 1])
  s1 <- as.vector(mu * (1 - theta)/theta)
  s2 <- as.vector((1 - theta)/theta * (1 - mu))
  loglike <- sum(weights * (lchoose(y.c, y) +  lbeta(y + s1, y.c - y + s2) - lbeta(s1, s2)))
  return(loglike)
}

