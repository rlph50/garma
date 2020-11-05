#' Diagnostic fit of a garma_model.
#'
#' Produces diagnostic plots of the model fit.
#' This function uses 'tsdiag' of the "stats" package.
#'
#' @param object (garma_model) The garma_model to produce the diagnostic plots for.
#' @param gof.lag (int) the maximum number of lags for a Portmanteau goodness-of-fit test.
#' @param ... further arguments to be passed to particular methods.
#' @return None. Diagnostics are generated.
#' @seealso The stats package tsdiag function: \url{https://stat.ethz.ch/R-manual/R-patched/library/stats/html/tsdiag.html}.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' tsdiag(mdl)
#' @export
tsdiag.garma_model<-function(object, gof.lag = 10, ...) {
  class(object)<-'Arima'
  stats::tsdiag(object,gof.lag=gof.lag, ...)
}
