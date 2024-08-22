#' Extract fitted values
#'
#' Fitted values are 1-step ahead predictions.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) array of 1-step ahead fitted values for the model.
#' @export
fitted.garma_model <- function(object, ...) {
  .byRef(object) # allows us to update the values of object
  if (!"fitted" %in% names(object)) {
    object <- internal_fitted_values(object)
  }
  return(object$fitted)
}
