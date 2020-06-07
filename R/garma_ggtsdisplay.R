#' For a k=1 Gegenbauer process, use semi-parametric methods to obtain short memory version of the process, then run a ggtsdisplay().
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param ... additional parameters to pass to ggtsdisplay
#' @return A ggplot object.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' garma_ggtsdisplay(ap)
#' @export
garma_ggtsdisplay<-function(x,...) {
  sp <- ggbr_semipara(x)
  arma_process <- extract_arma(x, sp$fd, sp$u)
  ggtsdisplay(arma_process,...)
}
