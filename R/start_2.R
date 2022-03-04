# Starting values for (ZI)Binomial and (ZI)BetaBinomial
start_2 <- function(X, Y, Z=NULL, offsetx=0, offsetz=0, weights=1, zero.inflated=FALSE, over.dispersed=FALSE) {
  y <- as.vector(Y[,1])
  y.c <- as.vector(rowSums(Y))
  phat <- y/y.c
  phat[phat <= 0] <- 1/(2 * y.c[phat == 0])
  phat[phat >= 1] <- 1 - 1/(2 * y.c[phat == 1])
  if (zero.inflated)
	model_zero <- lm.fit(Z, as.numeric(y == 0))
  else
	model_zero <- list(coefficients=NULL)
  model_count <- lm.fit(X, log(phat/(1 - phat)))
  model.p <- mean(y) / mean(y.c)
  phihat <- (1/(mean(y.c) - 1)) * (var(y)/(mean(y.c) * model.p * (1 - model.p)) - 1)
  theta <- if (over.dispersed) log(phihat/(1 - phihat)) else NULL
  list(count = model_count$coefficients, zero = model_zero$coefficients, theta = theta)
}