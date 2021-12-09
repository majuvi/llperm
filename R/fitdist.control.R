#' Extra arguments for fitdist
#'
#' This is part of the new implementation.
#'
#' @param method Optimization method
#' @param maxit Maximum number of iterations
#' @param trace Verbose T/F
#' @param start Initialization parameters
#' @return the extra arguments
#' @export
fitdist.control <- function (method = "BFGS", maxit = 10000, trace = FALSE, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, trace = trace, start = start)
  rval <- c(rval, list(...))
  rval$fnscale <- -1
  rval$hessian <- TRUE
  if (is.null(rval$reltol))
    rval$reltol <- .Machine$double.eps^(1/1.6)
  return(rval)
}