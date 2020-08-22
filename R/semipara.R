#' For a k=1 Gegenbauer process, use semi-parametric methods to estimate the Gegenbauer frequency and fractional differencing.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param k (int) The number of Gegenabuer frequencies
#' @param alpha (num)
#' @param method (char) One of "gsp" or "lpr" - lpr is the log-periodogram-regression technique, "gsp" is the Gaussian
#' semi-parametric technique. "gsp" is the default. Refer Arteche (1998).
#' @param min_freq (num) The minimum frequency to search through for peaks - default 0.0.
#' @param max_freq (num) The maximum frequency to search through for peaks - default 0.5.
#' @return An object of class "garma_semipara".
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' sp <- ggbr_semipara(ap)
#' print(sp)
#' @export
ggbr_semipara <- function(x,k=1,alpha=0.8,method='gsp',min_freq=0.0,max_freq=0.5) {
  if (!method%in%c('gsp','lpr')) stop('Invalid method. Should be one of "gsp" or "lpr".')
  if (alpha<=0|alpha>=1) stop('alpha should be between 0 and 1, but not0 and not 1.')
  k <- as.integer(k)
  if (k<=0) stop('k should be a small positive integer.')
  peak_idx_list <- c()
  peaks <- c()
  if (method=='gsp')
    for (k1 in 1:k) {
      r <- .gsp(x,alpha,peak_idx_list,min_freq,max_freq)
      peaks <- c(peaks,list(r))
      peak_idx_list <- c(peak_idx_list,r$f_idx)
    }
  if (method=='lpr')
    for (k1 in 1:k) {
      r <- .lpr(x,alpha,peak_idx_list,min_freq,max_freq)
      peaks <- c(peaks,list(r))
      peak_idx_list <- c(peak_idx_list,r$f_idx)
    }
  # look for repeated peaks. Shouldn't happen but sometimes does
  note<-''
  for (i in 1:(length(peaks)-1)) {
    if (i>=length(peaks)) break;
    r1 <- peaks[[i]]
    for (j in length(peaks):(i+1)) {
      if (j>length(peaks)) break;
      r2 <- peaks[[j]]
      if (r1$f_idx==r2$f_idx) {peaks[[j]]<-NULL;note<-'NOTE: Duplicate Peaks were removed.'}
    }
  }

  class(peaks) <- 'ggbr_factors'

  res <- list('ggbr_factors'=peaks, method=method, alpha=alpha, k=k,note=note)
  class(res) <- 'garma_semipara'
  return(res)
}

#' Print a 'ggbr_factors' object.
#' @param x An object of class ggbr_factors
#' @param ... further parameters for print function
#' @export
print.ggbr_factors<-function(x,...) {
  printf_9_4<-function(f) cat(sprintf('%9.4f',f))

  cat('                      ')
  for (k1 in 1:length(x)) cat(sprintf('%9.9s',paste0('Factor',k1)))
  cat('\nGegenbauer frequency: ')
  for (factor in x) printf_9_4(factor$f)
  cat('\nGegenbauer Period:    ')
  for (factor in x) printf_9_4(1/factor$f)
  cat('\nGegenbauer Exponent:  ')
  for (factor in x) printf_9_4(factor$fd)
  cat('\n')
}

#' Print a semiparameteric Gegenbauer estimation object.
#' @param x An object of class garma_semipara.
#' @param ... further parameters for print function
#' @export
print.garma_semipara<-function(x,...) {
  cat(sprintf('%s estimation of Gegenbauer process (k=1)\nFrequencies to use: (alpha=%f)\n\n',
              ifelse(x$method=='gsp','Gaussian Semi-Parametric','Log Periodogram Regression'),
              x$alpha))
  print(x$ggbr_factors)
  if(length(x$note)>0) if (nchar(x$note)>1) cat(paste0('\n',x$note))
}

#' For a k=1 Gegenbauer process, transform to remove Gegenbauer long memory component to get a short memory (ARMA) process.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @param ggbr_factors (list) Each element of the list represents a Gegenbauer factor and includes f, u and fd elements.
#' @return An object of same class as x.
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' # find semiparametric estimates of the Gegenbauer parameters.
#' sp <- ggbr_semipara(ap)
#' # extract the underlying short-memory ARMA process
#' ap_arma <- extract_arma(ap,sp$ggbr_factors)
#' summary(arima(ap_arma,order=c(1,0,0)))
#' @export
extract_arma<-function(x,ggbr_factors) {
  if (!is.null(ggbr_factors)) {
    for (factor in ggbr_factors) {
      ggbr_filter <- signal::Arma(b=1, a=.ggbr.coef(length(x),factor$fd,factor$u))
      sm          <- signal::filter(ggbr_filter, x)
    }
    sm            <- stats::ts(sm,start=stats::start(x),frequency=stats::frequency(x))
  }
  else sm<-x
  return(sm)
}

.garma_pgram<-function(x) {
  return(spectrum(x,plot=F,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE))
}

.yajima_ggbr_freq<-function(x,remove_peaks,min_freq=0.0,max_freq=0.5) {
  ssx       <- .garma_pgram(x)

  if (!is.null(min_freq)&!is.null(max_freq)) {
    from_idx  <- min(which(ssx$freq>=min_freq))
    if(is.null(from_idx)|is.na(from_idx)) from_idx<-0
    to_idx    <- max(which(ssx$freq<=max_freq))
    if(is.null(to_idx)|is.na(to_idx)) to_idx<-length(ssx$freq)
    ssx$spec2 <- ssx$spec[from_idx:to_idx]
  }
  else {
    ssx$spec2 <- ssx$spec
    from_idx <- 1
    to_idx <- length(ssx$spec)
  }

  if (length(remove_peaks)>0) {
    width <- as.integer(length(ssx$spec)/40)
    for (peak in remove_peaks) {
      start_idx <- max(peak-width,1)
      end_idx   <- min(peak+width,length(ssx$spec))
      for (i in start_idx:end_idx) ssx$spec2[i] <- (-1) # 1e-100
    }
  }
  f_idx     <- which.max(ssx$spec2[1:min(as.integer(length(x)/2),length(ssx$spec2))]) + from_idx - 1
  ggbr_freq <- ssx$freq[f_idx]
  return(list(f_idx=f_idx, ggbr_freq=ggbr_freq, ssx=ssx))
}

.gsp<-function(x,alpha,remove_peaks,min_freq=0.0,max_freq=1.0) {
  # as per Arteche 1998 "SEMIPARAMETRIC INFERENCE IN SEASONAL AND CYCLICAL LONG MEMORY PROCESSES"
  # determine "fd"
  c_fcn<-function(fd, omega, spec) {return(mean((omega^(2*fd)) * spec,rm.na=TRUE))}
  r_fcn<-function(fd, f_idx, ssx) {
    omega <- 2*pi*ssx$freq[1:(m-1)]
    spec1 <- ssx$spec[(f_idx+2):(f_idx+m)]
    min_idx <- f_idx - m
    if (m<f_idx+1) spec2 <- ssx$spec[(f_idx-m):(f_idx-2)]
    else {
      if (f_idx>2) spec2 <- c(ssx$spec[(f_idx-2):1], ssx$spec[length(ssx$spec):(length(ssx$spec)-(m-f_idx))])
      else spec2 <- c(ssx$spec[length(ssx$spec):(length(ssx$spec)-(m-f_idx))])
      spec2 <- spec2[1:(m-1)]
    }

    res <- log(c_fcn(fd, omega, spec1)) + log(c_fcn(fd, omega, spec2)) - 4*fd*mean(log(omega),rm.na=TRUE)
    if (is.infinite(res)|is.na(res)|is.null(res)) res<-1e200
    return(res)
  }
  # first identify the peak - the Gegenbauer frequency
  yf <- .yajima_ggbr_freq(x,remove_peaks,min_freq,max_freq)
  m  <- as.integer((length(x)/2)^alpha)

  fd <- stats::optimise(r_fcn, f_idx=yf$f_idx, ssx=yf$ssx, lower=-10, upper=10)$minimum / 2
  u  <- cos(2*pi*yf$ggbr_freq)

  return(list(fd=fd,f=yf$ggbr_freq,u=u,m=m,f_idx=yf$f_idx))
}

.lpr<-function(x,alpha,remove_peaks,min_freq=0.0,max_freq=1.0) {
  # first identify the peak - the Gegenbauer frequency
  yf       <- .yajima_ggbr_freq(x,remove_peaks,min_freq,max_freq)
  ssx      <- yf$ssx
  f_idx    <- yf$f_idx
  # next, estimate d
  m        <- as.integer((length(x)/2)^alpha)
  v        <- log(1:(m-1)) - mean(log(1:(m-1)))
  denom    <- 4*sum(v^2)
  spec1    <- ssx$spec[(f_idx+1):(f_idx+m-1)]
  min_idx  <- f_idx - m
  if (m<f_idx+1) spec2 <- ssx$spec[(f_idx-m):(f_idx-2)]
  else {
    spec2 <- c(ssx$spec[(f_idx-1):1], ssx$spec[length(ssx$spec):(length(ssx$spec)-(m-f_idx))])
    spec2 <- spec2[1:(m-1)]
  }

  numer    <- sum( v*(log(spec1)+log(spec2)))

  fd       <- (-0.5)*numer/denom
  u        <- cos(2*pi*yf$ggbr_freq)

  return(list(fd=fd,f=yf$ggbr_freq,u=u,m=m,f_idx=yf$f_idx))
}
