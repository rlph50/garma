#' The plot function generates a plot of actuals and predicted values for a "garma_model" object.
#' @param x (garma_model) The garma_model from which to plot the values.
#' @param ... other arguments to be passed to the "plot" function, including h (int) - the number of periods ahgead to forecast.
#' @return An R "plot" object.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' plot(mdl)
#' @export
plot.garma_model<-function(x,...) {
  .plot_garma_model(x,...)
}

#' The ggplot function generates a ggplot of actuals and predicted values for a "garma_model" object.
#' @param mdl (garma_model) The garma_model from which to ggplot the values.
#' @param h (int) The number of time periods to predict ahead. Default: 24
#' @param ... other parameters passed to ggplot.
#' @return A ggplot2 "ggplot" object. Note that the standard ggplot2 "+" notation can be used to enhance the default output.
#' @examples
#' library(ggplot2)
#'
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' ggplot(mdl)
#' @export
ggplot.garma_model<-function(mdl,h=24,...) {
  # plot forecasts from model
  fc <- predict.garma_model(mdl,n.ahead=h)

  if (mdl$y_freq>1) { # then we have actual dates not just an index
    idx <- seq(lubridate::make_date(mdl$y_start[1],mdl$y_start[2],15),by=mdl$y_freq,length.out=(length(mdl$y)+h))
    lubridate::day(idx) <- lubridate::days_in_month(idx)
    cutoff <- lubridate::make_date(mdl$y_end[1],mdl$y_end[2],15)
  } else {
    idx <- (mdl$y_start[1]):(mdl$y_end[1]+h)
    cutoff <- mdl$y_end[1]+1
  }

  df1 <- data.frame(dt=idx,grp='Actuals',value=c(mdl$y,rep(NA,h)))
  df2 <- data.frame(dt=idx,grp='Forecasts',value=c(as.numeric(mdl$fitted),fc$mean))
  df <- rbind(df1,df2)

  # assign some dummy vars to prevent RStudio check from throwing warning msgs
  dt <- value <- grp <- NA

  ggplot2::ggplot(df[!is.na(df$value),],ggplot2::aes(x=dt,y=value,color=grp),...) +
    ggplot2::geom_line() + ggplot2::ylab('') + ggplot2::xlab('') +
    ggplot2::geom_vline(xintercept=cutoff,color='red',linetype=2) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_blank()) +
    #ggplot2::scale_color_brewer(palette="Set1")
    scale_colour_manual(values=c('gray20','mediumblue',rep('gray',10)))
}

.plot_garma_model<-function(mdl,h=24,...) {
  # plot forecasts from model
  actuals <- zoo(stats::ts(mdl$y,start=mdl$y_start,end=mdl$y_end,frequency=mdl$y_freq))
  fitted <- zoo(mdl$fitted)
  fc <- zoo(predict.garma_model(mdl,n.ahead=h)$mean)
  graphics::plot(actuals,col='black',type='l',...)
  graphics::lines(fitted,col='blue')
  graphics::lines(fc,col='blue')
  graphics::abline(v=mdl$y_end,col='red',lty=2)
}

