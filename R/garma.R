#' @docType package
#' @name garma
#' GARMA: A package for estimating Gegenbauer long memory (GARMA) models.
#'
#' The GARMA package provides the main function "garma" as well as print, summary,
#' predict, forecast and plot/ggplot options.
#'
#' @section garma functions:
#' garma: estimate the parameters of a GARMA model, and return a "garma_model" object.
#' print: (summarised) print of the contents of a "garma_model" object.
#' summary: (detailed) print of the contents of a "garma_model" object.
#' predict: predict future values of a "garma_model" object.
#' forecast: predict future values of a "garma_model" object.
#' plot: predicts and plots the future values of a "garma_model" object.
#' ggplot: predicts and plots the future values of a "garma_model" object using the ggplot2 package.
#'
#' @docType package
#' @name garma

library(nloptr)
library(FKF)
library(assertthat)
library(zoo)
library(ggplot2)


