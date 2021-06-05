#' Display raw periodogram
#'
#' Display the raw periodogram for a time series, and not on a log scale.
#'
#' The standard "R" functions display periodograms on a log scale which can make it more difficult to locate high peaks in the spectrum at differing frequencies.
#' This routine will display the peaks on a raw scale.
#'
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param k (int) The number of Gegenbauer factors
#' @return A ggplot object representing the raw periodogram
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' sp <- ggbr_semipara(ap)
#' print(sp)
#' @export
gg_raw_pgram <- function(x,k=1) {
  x <- as.numeric(x) # "spectrum" function does weird things if x is a "ts" object

  ssx <- .garma_pgram(as.numeric(x))
  df  <- data.frame(Frequency=ssx$freq,Intensity=ssx$spec)

  sp  <- ggbr_semipara(x,k=k)
  annotate_df <- data.frame(x=numeric(0),y=numeric(0),label=character(0))
  for (factor in sp$ggbr_factors) {
    annotate_df <- rbind(annotate_df,
                         data.frame(x=factor$f,y=ssx$spec[factor$f_idx],label=sprintf(' Period: %.2f',1.0/factor$f)))
  }

  # set up some dummy vars to prevent the RStudio checks from throwing warnings...
  Frequency <- Intensity <- y <- label <- NA
  # now plot it
  ggplot2::ggplot(data=df,ggplot2::aes(x=Frequency,y=Intensity)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste('Raw periogram of',deparse(match.call()$x))) +
    ggplot2::geom_text(data=annotate_df,aes(x=x,y=y,label=label),size=2.5,hjust=0) +
    ggplot2::theme_bw()
}
