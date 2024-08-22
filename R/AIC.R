#' AIC for model
#'
#' AIC for model if available.
#'
#' Note - like the R "arima" function, if a xreg parameter is supplied, the Log Likelihood and the AIC
#' do not include any component of the regression. It relates strictly to the residuals from
#' the regression which are used to fit the GARMA model.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) Approximate AIC - uses approximation of whichever methoid is used to find model params.
#' @export
AIC.garma_model <- function(object, ...) {
  return(object$aic)
}
