#' For a k=1 Gegenbauer process, use semi-parametric methods to estimate the Gegenbauer frequency and fractional differencing.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param alpha (num)
#' @param method (char) One of "gsp" or "lpr" - lpr is the log-periodogram-regression technique, "gsp" is the Gaussian semi-parametric technique. "gsp" is the default. Refer Arteche (1998).
#' @return An object of class "garma_semipara".
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' sp <- ggbr_semipara(ap)
#' print(sp)
#' @export
ggbr_semipara <- function(x,alpha=0.8,method='gsp') {
  if (!method%in%c('gsp','lpr')) stop('Invalid method. Should be one of "gsp" or "lpr".')
  if (alpha<=0|alpha>=1) stop('alpha should be between 0 and 1, but not0 and not 1.')

  if (method=='gsp') res <- .gsp(x,alpha)
  if (method=='lpr') res <- .lpr(x,alpha)

  res$method <- method
  res$alpha  <- alpha
  class(res) <- 'garma_semipara'
  return(res)
}

#' Print a semiparameteric Gegenbauer estimation object.
#' @param x An object of class garma_semipara.
#' @param ... further parameters for print function
#' @export
print.garma_semipara<-function(x,...) {
  cat(sprintf('%s estimation of Gegenbauer process (k=1)\nFrequencies to use: m=%d (alpha=%f)\n\nGegenbauer frequency: %0.4f\nGegenbauer Period:    %0.4f\nFractional Exponent:  %0.4f\n',
              ifelse(x$method=='gsp','Gaussian Semi-Parametric','Log Periodogram Regression'),
              x$m,x$alpha,
              x$f,1/x$f,x$fd))
}

#' For a k=1 Gegenbauer process, transform to remove Gegenbauer long memory component to get a short memory (ARMA) process.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param fd (num) long memory parameter - for stationary process should be in (-1.0,0.5).
#' @param u (num) Cosine of Gegenbauer frequency.
#' @return An object of same class as x.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' # find semiparametric estimates of the Gegenbauer parameters.
#' sp <- ggbr_semipara(ap)
#' # extract the underlying short-memory ARMA process
#' ap_arma <- extract_arma(ap,sp$fd,sp$u)
#' summary(arima(ap_arma,order=c(1,0,0)))
#' @export
extract_arma<-function(x,fd,u) {
  ggbr_filter <- signal::Arma(b=1, a=.ggbr.coef(length(x),fd,u))
  sm          <- signal::filter(ggbr_filter, x)
  return(sm)
}

.yajima_ggbr_freq<-function(x) {
  ssx       <- spectrum(x,plot=F,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)
  f_idx     <- which.max(ssx$spec[1:as.integer(length(x)/2)])
  ggbr_freq <- ssx$freq[f_idx]
  return(list(f_idx=f_idx, ggbr_freq=ggbr_freq, ssx=ssx))
}

.gsp<-function(x,alpha) {
  # as per Arteche 1998 "SEMIPARAMETRIC INFERENCE IN SEASONAL AND CYCLICAL LONG MEMORY PROCESSES"
  # determine "fd"
  c_fcn<-function(fd, omega, spec) {return(mean((omega^(2*fd)) * spec,rm.na=TRUE))}
  r_fcn<-function(fd, f_idx, ssx) {
    omega <- 2*pi*ssx$freq[1:(m-1)]
    spec1 <- ssx$spec[(f_idx+2):(f_idx+m)]
    min_idx <- f_idx - m
    if (m<f_idx+1) spec2 <- ssx$spec[(f_idx-m):(f_idx-2)]
    else spec2 <- c(ssx$spec[(f_idx-2):1], ssx$spec[length(ssx$spec):(length(ssx$spec)-(m-f_idx))])
    return(log(c_fcn(fd, omega, spec1)) + log(c_fcn(fd, omega, spec2)) - 4*fd*mean(log(omega),rm.na=TRUE))
  }
  # first identify the peak - the Gegenbauer frequency
  yf <- .yajima_ggbr_freq(x)
  m  <- as.integer((length(x)/2)^alpha)

  fd <- optimise(r_fcn, f_idx=yf$f_idx, ssx=yf$ssx, lower=-10, upper=10)$minimum / 2
  u  <- cos(2*pi*yf$ggbr_freq)

  return(list(fd=fd,f=yf$ggbr_freq,u=u,m=m))
}

.lpr<-function(x,alpha) {
  # first identify the peak - the Gegenbauer frequency
  yf       <- .yajima_ggbr_freq(x)
  ssx      <- yf$ssx
  f_idx    <- yf$f_idx
  # next, estimate d
  m        <- as.integer((length(x)/2)^alpha)
  v        <- log(1:m) - mean(log(1:m))
  denom    <- 4*sum(v^2)
  spec1    <- ssx$spec[(f_idx+1):(f_idx+m)]
  min_idx  <- f_idx - m
  if (m<f_idx+1) spec2 <- ssx$spec[(f_idx-m):(f_idx-1)]
  else spec2 <- c(ssx$spec[(f_idx-1):1], ssx$spec[length(ssx$spec):(length(ssx$spec)-(m-f_idx))])
  numer    <- sum( v*(log(spec1)+log(spec2)))

  fd       <- (-0.5)*numer/denom
  u        <- cos(2*pi*yf$ggbr_freq)

  return(list(fd=fd,f=yf$ggbr_freq,u=u,m=m))
}
