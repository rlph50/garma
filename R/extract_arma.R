#' Extract underlying ARMA process.
#'
#' For a Gegenbauer process, transform to remove Gegenbauer long memory component to get a short memory (ARMA) process.
#' @param x (num) This should be a numeric vector representing the Gegenbauer process.
#' @param ggbr_factors (class) Each element of the list represents a Gegenbauer factor and includes f, u and fd elements.
#' @return An object of same class as x. Any time series attributes of x are copied to the returned object.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' # find semiparametric estimates of the Gegenbauer parameters.
#' sp <- ggbr_semipara(ap)
#' # extract the underlying short-memory ARMA process
#' ap_arma <- extract_arma(ap, sp)
#' summary(arima(ap_arma, order = c(1, 0, 0)))
#' @export
extract_arma <- function(x, ggbr_factors) {
  n <- length(x)
  theta_vec <- 1
  for (gf in ggbr_factors) {
    gc <- internal_ggbr_coef(n, gf$fd, cos(2*pi*gf$freq))
    theta_vec <- pracma::conv(theta_vec, gc)
  }
  if (theta_vec[1] == 1) theta_vec <- theta_vec[2:min(length(theta_vec), n)]
  theta_vec <- rev(theta_vec)
  qk <- length(theta_vec)
  resid <- numeric(0)
  for (i in 1:n) {
    ma_vec <- tail(c(rep(0, qk), resid), qk)
    s <- sum(theta_vec * ma_vec)
    resid <- c(resid, x[i] - s)
  }

  # copy across the ts() attributes if any
  stats::tsp(resid) <- stats::tsp(x)
  return(resid)
}
