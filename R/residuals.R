#' Residuals
#'
#' Response Residuals from the model.
#' @param object The garma_model object
#' @param type (chr) The type of residuals. Must be 'response'.
#' @param h (int) The number of periods ahead for the residuals. Must be 1.
#' @param ... Other parameters. Ignored.
#' @return (double) array of resideuals from the model.
#' @export
residuals.garma_model <- function(object, type = "response", h = 1, ...) {
  .byRef(object) # allows us to update the values of object
  if (!missing(type)) {
    if (type != "response") stop("Only response residuals are available.")
  }
  if (!missing(h)) {
    if (h != 1) stop("Only h=1 response residuals are available.")
  }
  if (!"residuals" %in% names(object)) {
    object <- internal_fitted_values(object)
  }
  return(object$residuals)
}
