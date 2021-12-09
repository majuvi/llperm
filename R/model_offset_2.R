#' Return an offset from model matrix
#'
#' This is part of the new implementation.
#'
#' @param x model matrix
#' @param terms terms
#' @param offset offset
#' @return offset term
model_offset_2 <- function(x, terms = NULL, offset = TRUE)
{
  if(is.null(terms)) terms <- attr(x, "terms")
  offsets <- attr(terms, "offset")
  if(length(offsets) > 0) {
    ans <- if(offset) x$"(offset)" else NULL
    if(is.null(ans)) ans <- 0
    for(i in offsets) ans <- ans + x[[deparse(attr(terms, "variables")[[i + 1]])]]
    ans
  }
  else {
    ans <- if(offset) x$"(offset)" else NULL
  }
  if(!is.null(ans) && !is.numeric(ans)) stop("'offset' must be numeric")
  ans
}
