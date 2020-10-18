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
#' @importFrom pracma pinv hessian conv deconv
#' @importFrom signal Arma filter
#' @importFrom zoo zoo index
#' @importFrom lubridate make_date day days_in_month
#' @importFrom forecast forecast ggtsdisplay
#' @importFrom Rsolnp solnp gosolnp
#' @importFrom ggplot2 ggplot geom_line aes geom_vline theme theme_bw scale_colour_manual geom_text labs element_blank
#' @importFrom graphics abline lines par plot
#' @importFrom stats diffinv end sd start ts tsp var spectrum frequency optimise arima
#' @importFrom utils tail packageVersion head globalVariables
NULL

