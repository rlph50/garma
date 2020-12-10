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

  titles <- .generate_default_plot_title(object,h=0)

  ## plot standardized residuals, acf of residuals, Ljung-Box p-values
  oldpar <- par(mfrow = c(3, 1))
  on.exit(par(oldpar))
  rs <- object$residuals
  stdres <- rs/sqrt(object$sigma2)
  plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
  abline(h = 0)
  acf(as.numeric(object$residuals), plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
  #pacf(as.numeric(object$residuals), plot = TRUE, main = "PACF of Residuals", na.action = na.pass)
  nlag <- gof.lag
  pval <- rep(NA,nlag)
  n_param <- object$order[1]+object$order[3]+object$k
  for(i in 1L:nlag) pval[i] <- Box.test(rs, i, type="Ljung-Box", fitdf=ifelse(i>n_param,n_param,0))$p.value
  plot(1L:nlag, pval, xlab = "Lag", ylab = "p value", ylim = c(0,1), xaxp=c(1,nlag,nlag-1),
       main = "p values for Ljung-Box statistic", sub=titles$sub)
  abline(h = 0.05, lty = 2, col = "blue")

  cat('NOTE: The degrees of freedom for the Ljung-Box statistic are determined differently to the standard R "tsdiag".\n')
  cat('      An adjustment is made for the number of model parameters as per Ljung & Box (1978).\n')
}
