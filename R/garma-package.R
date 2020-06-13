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
#' @importFrom pracma hessian pinv
#' @importFrom signal Arma filter
#' @importFrom zoo zoo
#' @importFrom lubridate make_date day days_in_month
#' @importFrom forecast forecast ggtsdisplay
#' @importFrom Rsolnp solnp
#' @importFrom ggplot2 ggplot geom_line xlab ylab aes geom_vline theme theme_bw scale_colour_manual annotate geom_text
#' @importFrom graphics abline lines par plot
#' @importFrom stats diffinv end sd start ts var spectrum frequency optimise
#' @importFrom utils tail
NULL

