#' Plot Forecasts from model.
#'
#' The plot function generates a plot of actuals and predicted values for a "garma_model" object.
#' @param x (garma_model) The garma_model from which to plot the values.
#' @param ... other arguments to be passed to the "plot" function, including h (int) - the number of periods ahead to forecast.
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


# This definition allows some vars to be used but not checked by the 'check' function
utils::globalVariables(c('.dt','.value','.grp'))

# default title for plots
.generate_default_plot_title<-function(mdl,h) {
  if (h>0) main <- paste('Forecast for',mdl$series)
  else main <- paste('Actual and Fitted for',mdl$series)
  sub <- sprintf('Model details: order=(%d,%d,%d), k=%d (method: %s)',
                 mdl$order[1],mdl$order[2],mdl$order[3],mdl$k,mdl$method)
  return(list(main=main,sub=sub))
}

#' ggplot of the Forecasts of the model.
#'
#' The ggplot function generates a ggplot of actuals and predicted values for a "garma_model" object.
#' This adds in sensible titles etc as best it can determine.
#'
#' @param mdl (garma_model) The garma_model from which to ggplot the values.
#' @param h (int) The number of time periods to predict ahead. Default: 24
#' @param include_fitted (bool) whether to include the 1-step ahead 'fitted' values in the plot. Default: FALSE
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
ggplot.garma_model<-function(mdl,h=24,include_fitted=FALSE,...) {
  # plot forecasts from model

  if (mdl$y_freq>1) { # then we have actual dates not just an index; set it up on x-axis
    by_str <- 'day'
    if (mdl$y_freq==4) by_str <- 'qtr'
    if (mdl$y_freq==12) by_str<-'month'
    idx <- seq(lubridate::make_date(mdl$y_start[1],mdl$y_start[2],1),by=by_str,length.out=(length(mdl$y)+h))
    lubridate::day(idx) <- lubridate::days_in_month(idx)
    cutoff <- lubridate::make_date(mdl$y_end[1],mdl$y_end[2],1)
  } else {
    idx <- (mdl$y_start[1]):(mdl$y_end[1]+h)
    cutoff <- mdl$y_end[1]+1
  }
  titles <- .generate_default_plot_title(mdl,h)

  if (h>0) {
    fc <- predict.garma_model(mdl,n.ahead=h)
    df1 <- data.frame(.dt=idx,.grp='Actuals',.value=c(as.numeric(mdl$y),rep(NA,h)))
    if (include_fitted) fitted <- as.numeric(mdl$fitted)
    else fitted <- c(rep(NA,length(mdl$fitted)-1),as.numeric(tail(mdl$y,1)))
    df2 <- data.frame(.dt=idx,.grp='Forecasts',.value=c(fitted,as.numeric(fc$pred)))
    df <- rbind(df1,df2)
  } else {
    df1 <- data.frame(.dt=idx,.grp='Actuals',.value=as.numeric(mdl$y))
    df2 <- data.frame(.dt=idx,.grp='Fitted',.value=as.numeric(mdl$fitted))
    df <- rbind(df1,df2)
  }

  ggplot2::ggplot(df[!is.na(df$.value),],ggplot2::aes(x=.dt,y=.value,color=.grp),...) +
    ggplot2::geom_line() + ggplot2::labs(title=titles$main,caption=titles$sub,x='',y='') +
    #ggplot2::ylab('') + ggplot2::xlab('') + ggplot2::ggtitle(title) +
    ggplot2::geom_vline(xintercept=cutoff,color='red',linetype=2) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_blank()) +
    ggplot2::scale_colour_manual(values=c('gray20','mediumblue',rep('gray',10)))
}

.plot_garma_model<-function(mdl,h=24,include_fitted=FALSE,xlab,ylab,main,sub,ylim,...) {
  # plot forecasts from model
  if (missing(xlab)) xlab<-''
  if (missing(ylab)) ylab<-ifelse(is.null(mdl$series),'',mdl$series)
  actuals <- zoo(stats::ts(c(as.numeric(mdl$y),rep(NA,h)),start=mdl$y_start,frequency=mdl$y_freq))
  fitted <- zoo(stats::ts(as.numeric(mdl$fitted),start=mdl$y_start,end=mdl$y_end,frequency=mdl$y_freq))

  # Titles
  titles <- .generate_default_plot_title(mdl,h)
  if (missing(main)) main <- titles$main
  if (missing(sub))  sub  <- titles$sub

  if (missing(ylim)) {
    if (h>0) {
      fc <- zoo(predict.garma_model(mdl,n.ahead=h)$pred)
      # y-limits
      y_min <- min(mdl$y,mdl$fitted,fc)
      y_max <- max(mdl$y,mdl$fitted,fc)
    } else {
      # y-limits
      y_min <- min(mdl$y,mdl$fitted)
      y_max <- max(mdl$y,mdl$fitted)
    }
    # Always include 0
    # if (y_min<0&y_max<0) y_max=0
    # if (y_min>0&y_max>0) y_min=0

    ylim <- c(y_min,y_max)
  }

  graphics::plot(actuals,col='black',type='l',xlab=xlab,ylab=ylab,main=main,sub=sub,ylim=ylim,...)
  if(h==0|include_fitted) graphics::lines(fitted,col='blue')
  if (h>0) { # then draw the predictions.
    fc <- zoo(predict.garma_model(mdl,n.ahead=h)$pred)
    graphics::lines(zoo::index(fc),fc,col='blue')
    graphics::abline(v=mdl$y_end,col='red',lty=2)
  }
}

