#' Display the raw periodogram for a time series, not on a log scale.
#' The standard "R" functions display periodograms on a log scale which can make it more difficult to locate high peaks in the spectrum at differning frequencies.
#' This routine will display the peaks on a raw scale.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @return A ggplot object representing the raw periodogram
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' sp <- ggbr_semipara(ap)
#' print(sp)
#' @export
gg_raw_pgram <- function(x) {
  yf <- .yajima_ggbr_freq(x)
  df <- data.frame(Frequency=yf$ssx$freq,Intensity=yf$ssx$spec)
  ggplot2::ggplot(data=df,ggplot2::aes(x=Frequency,y=Intensity)) +
    ggplot2::geom_line() +
    ggplot2::annotate('text',x=yf$ggbr_freq,y=yf$ssx$spec[yf$f_idx],label=sprintf(' Period: %.2f',1.0/yf$ggbr_freq),size=2.5,hjust=0) +
    ggplot2::theme_bw()
}
