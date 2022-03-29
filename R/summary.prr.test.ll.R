#' Print PRR test statistics (likelihood)
#'
#' This is part of the new implementation
#'
#' @param object the fitted PRR test object
#' @param digits number of digits to print
#' @export
summary.prr.test.ll <-
function(object, digits = max(3, getOption("digits") - 3),...)
{
    cat("Permutation of Regressor Residuals Test:\n")
    cat("observed likelihood ratio test p-value: ", format(round(object$p.value.obs, digits)),"\n")
    cat("permutation test for simulated p-value: ", format(round(object$p.value.sim, digits)),"\n")
    invisible(object)
}

