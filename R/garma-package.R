#' garma: A package for estimating and foreasting Gegenbauer time series models.
#'
#' The GARMA package provides the main function "garma" as well as print, summary,
#' predict, forecast and plot/ggplot options.
#'
#' @name garma
#' @docType package
#' @author Richard Hunt
#'
#' @importFrom nloptr cobyla directL lbfgs
#' @importFrom pracma hessian
#' @importFrom signal Arma filter
#' @importFrom zoo zoo
#' @importFrom lubridate make_date day days_in_month
#' @importFrom forecast forecast
#' @importFrom Rsolnp solnp
#' @importFrom pso psoptim
#' @importFrom BB BBoptim
#' @importFrom dfoptim nmkb hjkb
#' @importFrom ggplot2 ggplot geom_line xlab ylab aes geom_vline theme theme_bw scale_color_brewer
NULL

