#' @noRd
internal_fitted_values <- function(object) { # par,params,ggbr_factors,sigma2) { # Generate fitted values and residuals for GARMA process
  y <- object$diff_y
  orig_y <- object$y
  k <- object$k
  include.mean <- object$include.mean
  method <- object$method
  p <- object$order[1]
  id <- object$order[2]
  q <- object$order[3]
  sigma2 <- object$sigma2
  ggbr_factors <- object$model$ggbr_factors
  par <- object$obj_par

  # skip over Gegenbauer parameters
  start <- ifelse(is.null(object$periods), k * 2, length(object$periods))

  n <- length(y)
  phi_vec <- theta_vec <- 1
  if (p > 0) phi_vec <- object$model$phi
  if (q > 0) theta_vec <- object$model$theta

  inf_ar_vec <- c(1, -phi_vec)
  if (k > 0) { # multiply by inverse ggbr expansion polynomial
    for (gf in ggbr_factors) inf_ar_vec <- pracma::conv(inf_ar_vec, internal_ggbr_coef(n, -gf$fd, cos(2*pi*gf$freq)))
    # Divide by theta vec; ignore the remainder. That can be a problem for k=0, q>0 models.
    inf_ar_vec <- pracma::deconv(inf_ar_vec, theta_vec)$q
  }
  # finally get rid of the first "1"
  inf_ar_vec <- c(inf_ar_vec, rep(0, n)) # extra zeros in case this is pure AR
  inf_ar_vec <- inf_ar_vec[2:(n + 1)]

  if (id == 0) y_dash <- y else y_dash <- y - mean(y)
  if (p > 0) start <- p else start <- 1

  fitted <- y_dash[1:start] # Assume we have the first p correct... Thus first p residuals are 0.
  for (i in (start + 1):n) {
    yy_to <- (i - 1)
    yy_from <- ifelse(k == 0, yy_to - p, 1)
    if (yy_from < 1) yy_from <- 1
    yy <- y_dash[yy_from:yy_to]
    vec <- inf_ar_vec[yy_to:yy_from]
    fitted <- c(fitted, -sum(yy * vec))
  }

  if (id > 0) {
    fitted <- diffinv(fitted, differences = id, xi = head(orig_y, id))
    n <- length(orig_y)
    if (length(fitted) > n) fitted <- utils::head(fitted, n)
  } else {
    n <- length(orig_y)
  }

  # Now add on xreg fitted...
  if ("lm_xreg" %in% names(object)) {
    fitted <- fitted + fitted(object$lm_xreg)
  }
  resid <- orig_y - fitted

  object$fitted <- fitted
  object$residuals <- resid
  return(object)
}
