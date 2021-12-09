#' Cox family object
#'
#' This is part of the implementation in the original glmperm package (GPL-3)
#'
#' @return the family object
Cox <- function() {
        structure(list(family = "Cox"), class = "family")
}
