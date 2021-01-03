#' Diagnostic fit of a garma_model.
#'
#' Produces diagnostic plots of the model fit.
#' This function is copied from stats::tsdiag but modifies the fit_df for the Ljung-Box test for use with garma models.
#'
#' @param object (garma_model) The garma_model to produce the diagnostic plots for.
#' @param gof.lag (int) The number of lags to examine for the Ljung-Box white noise test.
#' @param ... further arguments to be passed to particular methods.
#' @return None. Diagnostics are generated.
#' @seealso The stats package tsdiag function: \url{https://stat.ethz.ch/R-manual/R-patched/library/stats/html/tsdiag.html}.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' tsdiag(mdl)
#' @export
tsdiag.garma_model<-function(object, gof.lag=10, ...) {

  titles <- .generate_default_plot_title(object,h=0)

  ## plot standardized residuals, acf of residuals, Ljung-Box p-values
  oldpar <- par(mfrow = c(2, 1))
  on.exit(par(oldpar))
  rs <- object$residuals
  stdres <- rs/sqrt(object$sigma2)
  plot(stdres, type = "h", main = paste(object$series,"- Standardized Residuals"), ylab = "")
  abline(h = 0)
  acf(as.numeric(object$residuals), plot = TRUE,
      main = paste(object$series," - ACF of Residuals"),
      sub=titles$sub,
      na.action = na.pass)
  # gof(object,gof.lag)
}

