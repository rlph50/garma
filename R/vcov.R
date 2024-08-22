#' Covariance matrix
#'
#' Covariance matrix of parameters if available
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) estimated variance-covariance matrix of the parameter estimates
#' @export
vcov.garma_model <- function(object, ...) {
  return(object$var.coef)
}
