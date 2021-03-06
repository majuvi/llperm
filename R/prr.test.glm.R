#' New Permutation of Regression Residuals (PRR) test
#'
#' This is part of the new implementation.
#'
#' @param formula Formula that specifies the model
#' @param data Data frame with the data
#' @param family Specify the distribution to fit
#' @param subset Take a given subset of the observations
#' @param weights Weights to use for the observations
#' @param na.action what to do with NA values
#' @param offset offset to use for each observation
#' @param contrasts contrasts to use for the model matrix
#' @param test.var variable in the formula to permute
#' @param n.rep number of permutation repetitions
#' @param method function that fits the model
#' @return list of likelihood and permutation based p-values
#' @export
		
prr.test.glm <- function (formula, data, family=gaussian, subset, weights, na.action, start=NULL, etastart=NULL, mustart=NULL, control = glm.control(...),
                          offset, contrasts=NULL, test.var=NULL, n.rep=100, seed=12345, model=TRUE, method="glm.fit", x=FALSE, y=TRUE, ...)
{

  if (is.character(family))
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
     family <- family()
  if (is.null(family$family)) {
     print(family)
     stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  #if (identical(method, "glm.fit"))
  #    control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  # Which columns in X are related to test.var?
  term.labels  <- attr(mt, "term.labels")
  term.columns <- attr(X, "assign")
  if (attr(mt, "intercept")) {
    term.labels  <- c("(Intercept)", term.labels)
    term.columns <- term.columns + 1
  }
  cols.testvar <- term.labels[term.columns] == test.var
  # Predict test.var from other variables
  fit <- lm.fit(X[,!cols.testvar, drop=F], X[,cols.testvar, drop=F])
  R   <- fit$residual
  if (is.null(dim(R)) | (length(dim(R)) == 1L))
    dim(R) <- c(length(R), 1)

  # Calculate observed p-value
  X0 <- X[,!cols.testvar, drop=F]
  XR <- cbind(X0, R)
  f     <- call(if (is.function(method)) "method" else method)
  fit.1 <- eval(as.call(c(as.list(f), list(x = XR, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt, "intercept") > 0),
                          substitute(...()))))
  fit.2 <- eval(as.call(c(as.list(f), list(x = X0, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt, "intercept") > 0),
                          substitute(...()))))
  if (method == "glm.fit") { # loglik attribute is missing from base glm.fit
    class(fit.1) <- c("glm", "lm")
    fit.1$loglik <- logLik(fit.1)[1]
    class(fit.2) <- c("glm", "lm")
    fit.2$loglik <- logLik(fit.2)[1]
  }
  ll.1  <- fit.1$loglik
  ll.2  <- fit.2$loglik
  lr.df <- fit.2$df.residual - fit.1$df.residual
  p.value.obs <- pchisq(2 * (ll.1 - ll.2), lr.df, lower.tail=F)

  # Permutations: calculate simulated p-value
  set.seed(seed)
  ll.r <- rep(0, length(n.rep))
  for (i in 1:n.rep) {
    rnd <- sample(1:nrow(X), replace = FALSE)
    XR  <- cbind(X0, R[rnd,])
    fit.r <- eval(as.call(c(as.list(f), list(x = XR, y = Y, weights = weights, start = start,
        etastart = etastart, mustart = mustart, offset = offset, family = family, control = control, intercept = attr(mt, "intercept") > 0),
                            substitute(...()))))
    if (method == "glm.fit") {
      class(fit.r) <- c("glm", "lm")
      fit.r$loglik <- logLik(fit.r)[1]
    }
    ll.r[i] <- fit.r$loglik
  }
  p.value.sim <- mean(ll.1 < ll.r)
  
  # Save model frame, model matrix, and response vector?
  if (model)
    fit.1$model <- mf
  fit.1$na.action <- attr(mf, "na.action")
  if (x)
    fit.1$x <- X
  if (!y)
    fit.1$y <- NULL
  
  results <- list(fit=fit.1, p.value.obs=p.value.obs, p.value.sim=p.value.sim)
  class(results) <- "prr.test.glm"
  return(results)
}