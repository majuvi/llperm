# Starting values for (ZI)Poisson and (ZI)NegativeBinomial TODO: maybe change these to lm.fit also
start_1 <- function(X, Y, Z=NULL, offsetx=0, offsetz=0, weights=1, zero.inflated=FALSE, over.dispersed=FALSE) {
  if (zero.inflated)
	model_zero <- glm.fit(Z, as.integer(Y == 0), weights = weights, family = binomial(link = "logit"), offset = offsetz)
  else
	model_zero <- list(coefficients=NULL)
  model_count <- glm.fit(X, Y, family = poisson(), weights = weights, offset = offsetx)
  theta <- if (over.dispersed) 0 else NULL
  list(count = model_count$coefficients, zero = model_zero$coefficients, theta = theta)
}