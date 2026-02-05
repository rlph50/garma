#' Forecast future values.
#'
#' The forecast function predicts future values of a "garma_model" object, and is exactly the same as the "predict" function with slightly different parameter values.
#' @param object (garma_model) The garma_model from which to forecast the values.
#' @param h (int) The number of time periods to predict ahead. Default: 1
#' @param newdata (real vector or matrix) If the original model was fitted with the 'xreg=' option then this will provide the xreg
#'   values for predictions. If this is a vector then its length should be 'h'; if it is a matrix then it should have 'h' rows.
#'
#'   It should have columns with the same names as the original xreg matrix.
#' @param ... Other parameters passed to the forecast function. For "garma_model" objects, these are ignored.
#' @return - a "ts" object containing the requested forecasts.
#' @examples
#' library(forecast)
#'
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' mdl <- garma(ap, order = c(9, 1, 0), k = 0, method = "CSS", include.mean = FALSE)
#' forecast(mdl, h = 12)
#' @export
forecast.garma_model <- function(object, h = 1, newdata = NULL, ...) {
  temp_pred <- predict.garma_model(object, n.ahead = h, newdata = NULL, ...)

  x_res <- ts(object$y, start = object$y_start, frequency = object$y_freq)
  result <- list(
    mean = temp_pred$pred,
    model = object,
    fitted = fitted.garma_model(object),
    residuals = residuals.garma_model(object),
    x = x_res
  )
  class(result) <- c("garma_forecast", "forecast")
  return(result)
}

#' Return the fitted values for a GARMA forecast.
#' @param object (garma_model) A 'fit' object from the `garma()` function.
#' @param ... other parameters, which are ignored.
#' @export
fitted.garma_forecast <- function(object, ...) {
  result <- object$fitted
  return(result)
}

#' Return the residuals for a GARMA forecast.
#' @param object (garma_model) A 'fit' object from the `garma()` function.
#' @param ... other parameters, which are ignored.
#' @export
residuals.garma_forecast <- function(object, ...) {
  result <- object$residuals
  return(result)
}
