#' Goodness-of-Fit test for a garma_model.
#'
#' Provides a goodness-of-fit test for a GARMA Model, using Bartletts Tp test.
#' This has been justified for long memory and for GARMA models by Delgado, Hidalgo and Velasco (2005).
#'
#' @param object (garma_model) The garma_model to test.
#' @return Invisibly returns the array of p-values from the test.
#' @references
#' M Delgado, J Hidalgo, and C Velasco. Distribution free goodness-of-fit tests for linear processes.
#' The Annals of Statistics, 33(6):2568–2609, 2005. DOI: https://doi.org/10.1214/009053605000000606.
#' @export
gof <- function(object) {
  r <- as.numeric(residuals(object))
  n <- length(r)
  tilde_n <- as.integer(n / 2)

  I_eps <- spec.pgram(r, taper = 0, fast = FALSE, demean = FALSE, detrend = FALSE, plot = FALSE)
  pv <- numeric(0)
  for (i in 2:length(I_eps$spec)) pv <- c(pv, ks.test(I_eps$spec[1:i], ecdf(I_eps$spec))$p.value)
  cat(sprintf(
    "\nBartletts Tp test.\n\nTest for H0: Residuals are White Noise.\nEvaluating %d frequencies, smallest p-value: %0.4f at frequency %.4f period %.4f.\n",
    length(I_eps$spec), min(pv),
    (which.min(pv) + 1) * 2 * pi / n,
    n / ((which.min(pv) + 1))
  ))
  return(invisible(pv))
}
