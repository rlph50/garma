#' Diagnostic fit of a garma_model.
#'
#' Produces diagnostic plots of the model fit.
#' This function is copied from stats::tsdiag but modifies the fit_df for the Ljung-Box test for use with garma models.
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
#  class(object)<-'Arima'
#  stats::tsdiag(object,gof.lag=gof.lag, ...)

  ## plot standardized residuals, acf of residuals, Ljung-Box p-values
  oldpar <- par(mfrow = c(3, 1))
  on.exit(par(oldpar))
  rs <- object$residuals
  stdres <- rs/sqrt(object$sigma2)
  plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
  abline(h = 0)
  acf(object$residuals, plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
  nlag <- gof.lag
  pval <- rep(NA,nlag)
  n_param <- sum(object$order)+object$k*2
  if (n_param<1) n_param<-1
  if (n_param>nlag) nlag <- n_param+nlag
  for(i in n_param:nlag) pval[i] <- Box.test(rs, i, type="Ljung-Box", fitdf=n_param)$p.value
  plot(1L:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0,1),
       main = "p values for Ljung-Box statistic")
  abline(h = 0.05, lty = 2, col = "blue")
}
