#' garma: A package for estimating and foreasting Gegenbauer time series models.
#'
#' The GARMA package provides the main function "garma" as well as print, summary,
#' predict, forecast and plot/ggplot options.
#'
#' @name garma
#' @docType package
#' @author Richard Hunt
#'
#' @importFrom nloptr cobyla directL
#' @importFrom pracma pinv hessian conv deconv
#' @importFrom signal Arma filter
#' @importFrom zoo zoo index
#' @importFrom lubridate make_date day days_in_month
#' @importFrom forecast forecast ggtsdisplay
#' @importFrom Rsolnp solnp gosolnp
#' @importFrom ggplot2 autoplot ggplot geom_line aes geom_vline theme theme_bw scale_colour_manual geom_text labs element_blank
#' @importFrom graphics abline lines par plot
#' @importFrom stats diffinv end sd start ts tsp var spectrum spec.pgram ks.test ecdf frequency optimise arima Box.test acf pacf tsdiag na.pass optim lm residuals coef fft formula vcov predict
#' @importFrom utils tail packageVersion head globalVariables

#' @aliases garma-package
"_PACKAGE"
