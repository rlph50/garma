#' Extract semiparametric estimates of the Gegenbauer factors.
#'
#' For a Gegenbauer process, use semi-parametric methods to estimate the Gegenbauer frequency and fractional differencing.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param periods (num) This parameter can be used to specify a fixed period or set of periods for the
#' Gegenbauer periodicity. For instance if you have monthly data, then it might be sensible (after an examination of the
#' periodogram) to set `periods = 12`. The default value is NULL. Either `periods` or `k` parameters must be specified
#' but not both - `periods` implies fixed period(s) are to be used and `k` implies that the periods should
#' be estimated.
#' @param k (int) This parameter indicates that the algorithm should estimate the `k` frequencies semi-parametrically,
#' before estimating the degree of fractional differencing at each period.
#'
#' An alternative is the `periods` parameter which can be used to specify exactly which periods should be used by
#' the model.
#' @param alpha (num) Default = 0.8 - This is the bandwidth for the semiparametric estimate, and should be between 0 and 1.
#' Robinson (1994) indicated optimality for a (scaled) version of `alpha` = 0.8, at least for the "lpr" `method`.
#' @param method (char) One of "gsp" or "lpr" - lpr is the log-periodogram-regression technique, "gsp" is the Gaussian
#' semi-parametric technique. "gsp" is the default. Refer Arteche & Robinson (1998).
#' @return An object of class "garma_semipara".
#' @references
#' J Arteche and P Robinson. Semiparametric inference in seasonal and cyclical long memory processes. Journal of Time Series
#' Analysis, 21(1):1–25, 2000. DOI: https://doi.org/10.1111/1467-9892.00170
#'
#' P Robinson. Rates of convergence and optimal spectral bandwidth for long range dependence. Probability Theory and Related
#' Fields, 99:443–473, 1994. DOI: https://doi.org/10.1007/BF01199901.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' sp <- ggbr_semipara(ap)
#' print(sp)
#' @export
ggbr_semipara <- function(x, periods = NULL, k = 1, alpha = 0.8, method = "gsp") {
  x <- as.numeric(x)
  k <- as.integer(k)

  if (!is.numeric(x)) {
    rlang::abort("x should be numeric or a ts object.\n")
  }
  if (any(is.na(x))) {
    rlang::abort("x should not have any missing values.\n")
  }
  if (is.null(periods)) {
    if (is.null(k) | !is.numeric(k) | k <= 0) {
      rlang::abort("The k parameter must be a non-negative integer.\n")
    }
  } else {
    if (!is.numeric(periods) | length(periods) == 0) {
      rlang::abort("The 'periods' parameter must be a numeric vector of at least 1 element.")
    }
    if (any(periods < 0)) {
      rlang::abort("The 'periods' parameter cannot contain negative values.")
    }
  }

  if (!is.character(method)) {
    rlang::abort("The method parameter should be character and one of 'gsp' or 'lpr'.")
  }
  if (length(method) != 1) {
    rlang::abort("The method parameter should be a single value only.")
  }
  if (!method %in% c("gsp", "lpr")) {
    rlang::abort('Invalid method. Should be one of "gsp" or "lpr".')
  }
  if (alpha <= 0 | alpha >= 1) {
    rlang::abort("The alpha parameter should be between 0 and 1.")
  }

  res <- internal_semipara(x, periods, alpha, k, method)

  return(res)
}

#' Print a 'ggbr_factors' object.
#' @param x An object of class ggbr_factors
#' @param ... further parameters for print function
#' @return null
#' @export
print.ggbr_factors <- function(x, ...) {
  printf_9_4 <- function(f) cat(sprintf("%9.4f", f))

  if (length(x) > 0) {
    cat("                      ")
    for (k1 in 1:length(x)) cat(sprintf("%9.9s", paste0("Factor", k1)))
    cat("\nGegenbauer frequency: ")
    for (factor in x) printf_9_4(factor$freq)
    cat("\nGegenbauer Period:    ")
    for (factor in x) printf_9_4(factor$period)
    cat("\nGegenbauer Exponent:  ")
    for (factor in x) printf_9_4(factor$fd)
    cat("\n")
  } else {
    cat("No Gegenbauer factors.\n")
  }
}
