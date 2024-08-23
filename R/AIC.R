#' AIC for model
#'
#' Approximate AIC for model.
#'
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) Approximate AIC - uses approximation of whichever methoid is used to find model params.
#' @export
AIC.garma_model <- function(object, ...) {
  return(object$aic)
}
