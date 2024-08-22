#' Model Coefficients
#'
#' Model Coefficients/parameters.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) array of parameter value estimates from the fitted model.
#' @export
coef.garma_model <- function(object, ...) {
  return(object$coef[1, ])
}
