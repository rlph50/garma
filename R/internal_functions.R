.ggbr.coef<-function(n,d,u) {
  if (n>0) {
    cf<-c(1,2*d*u,2*d*(d+1)*(u^2)-d)
    if (n>4) for (j in 4:n) cf[j]<-(2*u*(j+d-2)*cf[j-1]-(j-1+2*d-2)*cf[j-2])/(j-1)
    cf <- cf[1:n]
  } else cf <- NA
  return(cf)
}

.a_fcn<-function(a_vec,freq){
  # this is a utility function used to find the short-memory spectral density.
  n_freq <- length(freq)
  cos_sum <- rep(0.0,n_freq)
  a_len<-length(a_vec)
  for (i in 1:a_len) cos_sum <- cos_sum + a_vec[i]*cos(2*i*pi*freq)
  sin_sum <- rep(0.0,n_freq)
  for (i in 1:a_len) sin_sum <- sin_sum + a_vec[i]*sin(2*i*pi*freq)

  return( (1+cos_sum)^2+sin_sum^2)
}

.a_fcn_cos<-function(a_vec,freq){
  # this is a utility function used to find the short-memory spectral density.
  n_freq <- length(freq)
  cos_sum <- rep(0.0,n_freq)
  a_len<-length(a_vec)
  for (i in 1:a_len) cos_sum <- cos_sum + a_vec[i]*cos(2*i*pi*freq)

  return(1+cos_sum)
}
.a_fcn_sin<-function(a_vec,freq){
  # this is a utility function used to find the short-memory spectral density.
  n_freq <- length(freq)
  sin_sum <- rep(0.0,n_freq)
  a_len<-length(a_vec)
  for (i in 1:a_len) sin_sum <- sin_sum + a_vec[i]*sin(2*i*pi*freq)

  return(sin_sum)
}

.getPackageVersion<-function(pkgname) {
  return(paste0(crayon::black("\n\nPackage "),crayon::blue(pkgname),': ',crayon::black(utils::packageVersion('garma')),'\n'))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(.getPackageVersion(pkgname))
}

# next function allows us to update the object
.byRef <- function(..., envir=parent.frame(), inherits=TRUE) {
  cl <- match.call(expand.dots = TRUE)
  cl[c(1, match(c("envir", "inherits"), names(cl), 0L))] <- NULL
  for (x in as.list(cl)) {
    s <- substitute(x)
    sx <- do.call(substitute, list(s), envir=envir)
    dx <- deparse(sx)
    expr <- substitute(assign(dx, s, envir=parent.frame(), inherits=inherits))
    do.call(on.exit, list(expr, add=TRUE), envir=envir)
  }
}

