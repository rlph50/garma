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
  # nlag <- gof.lag
  # pval <- rep(NA,nlag)
  # n_param <- object$order[1]+object$order[3]+object$k
  # for(i in 1L:nlag) pval[i] <- Box.test(rs, i, type="Ljung-Box", fitdf=ifelse(i>n_param,n_param,0))$p.value
  # plot(1L:nlag, pval, xlab = "Lag", ylab = "p value", ylim = c(0,1), xaxp=c(1,nlag,nlag-1),
  #      main = "p values for Ljung-Box statistic", sub=titles$sub)
  # abline(h = 0.05, lty = 2, col = "blue")
  #
  # cat('NOTE: The degrees of freedom for the Ljung-Box statistic are determined differently\n')
  # cat('      to the standard R "tsdiag". An adjustment is made for the number of model\n')
  # cat('      parameters as per Ljung & Box (1978).\n')
  # cat('      Note also that the Ljung-Box test has not been formally established as valid\n')
  # cat('      for GARMA processes.\n')

  # Fay & Philippe (2002) spectral test for gof.
  ss <- object$spectrum
  n_freq <- length(ss$spec)
  mean_methods <- c('QML','CSS')
  params <- list(y=object$y, orig_y=object$diff_y, ss=ss, p=object$order[1],q=object$order[3],d=object$order[2],k=object$k,
                 include.mean=object$include.mean, include.drift=object$include.drift, drift=object$model$drift,
                 method=object$method,
                 est_mean=ifelse(object$method%in%mean_methods,TRUE,FALSE),
                 scale=sd(object$y),m_trunc=object$m_trunc)

  spec_den_inv <- .calc_spec_den_inv(object$obj_par,params)

  I_f <- ss$spec*spec_den_inv / (2*pi*n_freq)
  s <- log(mean(I_f,na.rm=TRUE)) - mean(log(I_f),na.rm=TRUE) - digamma(1) # (-digamma(1)) == Euler-Mascheroni constant
  s_sigma2 <- pracma::psi(1,1)
  s2 <- s/s_sigma2
  p <- 1-pnorm(s2)

  cat(sprintf('Fay-Philippe Goodness of Fit test.\n\nz-Statistic: %0.4f\np-value:     %0.4f\n',s2,p))
}
