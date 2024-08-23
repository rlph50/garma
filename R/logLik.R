#' Log Likelihood
#'
#' The approximate likelihood for the model.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return Object of class "logLik" with values for the (approx) log-likelihood for the model
#' @export
logLik.garma_model <- function(object, ...) {
  # Need to figure out how to indicate these are REML estimates not true LL.
  res <- structure(object$loglik, df = length(object$y) - 1, nobs = length(object$y), class = "logLik")
  return(res)
}
