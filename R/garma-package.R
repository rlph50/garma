#' garma: A package for estimating and foreasting Gegenbauer time series models.
#'
#' The GARMA package provides the main function "garma" as well as print, summary,
#' predict, forecast and plot/ggplot options.
#'
#' @name garma
#' @docType package
#' @author Richard Hunt
#'
#' @importFrom nloptr nloptr cobyla directL lbfgs mma auglag
#' @importFrom BB BBoptim
#' @importFrom pso psoptim
#' @importFrom dfoptim hjkb nmkb
#' @importFrom GA de
#' @importFrom pracma pinv hessian conv deconv grad psi
#' @importFrom signal Arma filter
#' @importFrom zoo zoo index
#' @importFrom tswge factor.wge
#' @importFrom lubridate make_date day days_in_month
#' @importFrom forecast forecast ggtsdisplay
#' @importFrom Rsolnp solnp gosolnp
#' @importFrom ggplot2 autoplot ggplot geom_line aes geom_vline theme theme_bw scale_colour_manual geom_text labs element_blank
#' @importFrom graphics abline lines par plot
#' @importFrom stats diffinv end sd start ts tsp var spectrum spec.pgram ks.test ecdf frequency optimise arima Box.test acf pacf tsdiag na.pass optim lm residuals coef fft
#' @importFrom utils tail packageVersion head globalVariables

#' @keywords internal
#' @aliases garma-package
"_PACKAGE"

