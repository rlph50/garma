#' For a k=1 Gegenbauer process, use semi-parametric methods to estimate the Gegenbauer frequency and fractional differencing.
#' Then reverse to get a short memory process and return all this.
#' @param x (num) This should be a numeric vector representing the process to estimate.
#' @return An object of class "garma_semipara".
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers,12))
#' sp <- semipara(ap)
#' print(sp)
#' @export
semipara <- function(x) {
  res <- .sm_est(x)
  class(res) <- 'garma_semipara'
  return(res)
}

#' Print a semiparameteric Gegenbauer estimation object.
#' @param x An object of class garma_semipara.
#' @param ... further parameters for print function
#' @export
print.garma_semipara<-function(x,...) {
  cat(sprintf('Semiparametric estimation of Gegenbauer process (k=1)\n\nGegenbauer frequency: %0.4f\nGegenbauer Period:    %0.4f\nFractional Exponent:  %0.4f\n',
              x$f,1/x$f,x$fd))
}

# transform back to short memory
.lm_filter<-function(x,fd,u) {
  n<-length(x)
  xx<-pracma::Toeplitz(.ggbr.coef(n,fd,u),c(1,rep(0,n-1)))
  sm<-solve(xx)%*%x
  return(sm)
}

# estimate fd - fractional differencing at u
.sm_est<-function(x) {
  # Initial estimate for GARMA model with f,u,fd
  # first identify the peak - the Gegenbauer frequency
  ssx      <- spectrum(x,plot=F)
  f_idx    <- which.max(ssx$spec[1:as.integer(length(x)/2)])
  f        <- ssx$freq[f_idx]
  # next, estimate d
  m        <- as.integer(length(x)^0.8)
  v        <- log(1:m)
  v_mean   <- sum(v)/m
  v        <- (v-v_mean)
  denom    <- 2*sum(v^2)
  numer    <- 0
  for (j in 1:m) {
    idx1 <- f_idx+j
    if (idx1>length(ssx$spec)) idx1<-length(ssx$spec)-(idx1-length(ssx$spec))
    idx2 <- f_idx-j
    if (idx2<1) idx2<-abs(idx2)+1
    numer<-numer+v[j]*(log(ssx$spec[idx1])+log(ssx$spec[idx2]))
  }

  fd <- (-0.5)*numer/denom
  u  <- cos(2*pi*f)
  sm <- .lm_filter(x,fd,u)
  return(list(fd=fd,f=(2*pi*f),u=u,sm_x=sm))
}

