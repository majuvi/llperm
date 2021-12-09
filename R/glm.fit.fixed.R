#' Extend the base GLM object with the likelihood attribute
#'
#' This is part of the new implementation.
#'
#' @param y Response vector
#' @param x Model matrix 
#' @return the fitted model
#' @export
glm.fit.fixed <- function (x, y, ...) {
  fit <- glm.fit(x = x, y = y, ...)
  class(fit) <- c("glm", "lm")
  fit$loglik <- logLik(fit)[1]
  return(fit)
}
