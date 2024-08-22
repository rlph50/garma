#' summarise a garma_model object.
#'
#' The summary function provides a summary of a "garma_model" object, printed to the output.
#' @param object (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' mdl <- garma(ap, order = c(9, 1, 0), k = 0, method = "CSS", include.mean = FALSE)
#' summary(mdl)
#' @return (null)
#' @export
summary.garma_model <- function(object, ...) {
  internal_print_garma_model(object, verbose = TRUE)
}

#' print a garma_model object.
#'
#' The print function prints a summary of a "garma_model" object, printed to the output.
#' @param x (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @return (null)
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' mdl <- garma(ap, order = c(9, 1, 0), k = 0, method = "CSS", include.mean = FALSE)
#' print(mdl)
#' @export
print.garma_model <- function(x, ...) {
  internal_print_garma_model(x, verbose = FALSE)
}

#' @noRd
internal_print_garma_model <- function(mdl, verbose = TRUE) {
  cat("\nCall:", deparse(mdl$call, width.cutoff = 75L), "", sep = "\n")
  if (!mdl$include.mean) cat("No ")
  cat("Mean term was fitted.\n")
  if (!mdl$include.drift) cat("No ")
  cat("Drift (trend) term was fitted.\n\n")
  if (mdl$method == "Whittle" & mdl$k > 1 & verbose) {
    cat("NOTE: Giraitis, Hidalgo & Robinson (2001) establish consistency and asymptotic Normality only for k=1 processes.\n")
    cat("      Whilst it seems likely that the results also hold for a general k factor process, we are unaware of specific papers\n")
    cat("      which establish this point. The user should be aware therefore that the estimate standard errors etc are provided\n")
    cat("      without the relevant theory being established.\n\n")
  }
  if (verbose) {
    with(
      mdl,
      cat(sprintf(
        "Summary of Gegenbauer Time Series model.\n\nFit using %s method, with Order=(%d,%d,%d) %s\nOptimisation Method(s):  %s\n",
        method, order[1], order[2], order[3],
        ifelse(is.null(mdl$periods), paste0("k=", k), paste0("periods = ", paste0(mdl$periods, collapse = ", "))),
        paste(mdl$opt_method, collapse = ", ")
      ))
    )
    if (mdl$opt_method[[1]] == "best") cat(sprintf("Best optimisation method selected: %s\n", mdl$opt_method_selected))
    cat(sprintf("Convergence Code: %d %s\nOptimal Value found: %0.8f\n",
                mdl$convergence,
                ifelse(!is.null(mdl$conv_message) &!is.na(mdl$conv_message) &nchar(mdl$conv_message) > 1, paste0("(", mdl$conv_message, ")"), ""),
                mdl$obj_value))
  }

  phi_vec <- c(1, -mdl$model$phi)
  theta_vec <- c(1, -mdl$model$theta)
  if (any(Mod(polyroot(phi_vec)) < 1) | any(Mod(polyroot(theta_vec)) < 1)) {
    warning("model estimates are not Stationary! Forecasts may become unbounded.\n")
  }
  cat("\nCoefficients:\n")
  coef <- mdl$coef
  print.default(coef, print.gap = 2, digits = 4)
  cat("\n")

  if (mdl$k > 0) print(mdl$model$ggbr_factors)

  if (mdl$sigma2 > 0) { # make sure we have a valid sigma2 before printing it...
    cat(sprintf("\n\nsigma^2 estimated as %0.4f: log likelihood = %f, aic = %f\n\n", mdl$sigma2, mdl$loglik, mdl$aic))
  }
}
