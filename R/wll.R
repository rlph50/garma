# Estimate a Ggbr model using the Whittle-like-Log (WLL) method.
#
# wll.ggbr.obj - objective function to be minimised to get WLL estimates.
# It is not anticipated that an end-ser would have a need to call this function directly.
# @param par - the parameters to evaluate the function at
# @param params - other parameters - including the p, q, k, and scale parameters and (ss) the spectrum .
# @return The value of the objective at the point par.
.wll.ggbr.obj<-function(par,params) {
  # Objective function to be minimised for the BNP estimates
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k
  scale <- params$scale

  n_freq<-length(ss$freq)
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  # if (k==1) {
  #   u     <- par[start]
  #   fd    <- par[start+1]
  #   start <- start+2
  # }
  u <- c()
  fd <- c()
  if (k>0) for (k1 in 1:k) {
    u     <- c(u,par[start])
    fd    <- c(fd,par[start+1])
    start <- start+2
  }
  if (p>0) phi_vec   <- (-par[start:(start+p-1)])         else phi_vec<-1
  if (q>0) theta_vec <- (-par[(p+start):(length(par)-1)]) else theta_vec<-1
  sigma2 <- par[length(par)] / (scale^2)

  cos_2_pi_f <- cos(2*pi*freq)
  mod_phi <- mod_theta <- 1
  if (p>0) mod_phi   <- .a_fcn(phi_vec,ss$freq)
  if (q>0) mod_theta <- .a_fcn(theta_vec,ss$freq)

  spec_den_inv <- 2.0*pi / sigma2 * mod_phi / mod_theta    # Inverse of spectral density
  if (k>0) for (k1 in 1:k) spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u[k1])^2))^fd[k1]

  spec_den_inv[is.infinite(spec_den_inv)] <- NA #1.0e500
  spec_den_inv[spec_den_inv<=0] <- NA
  #print(spec_den_inv)
  I_f <- spec*spec_den_inv
  res <- sum((log(I_f))^2,na.rm=TRUE)
  #cat(sprintf("ret %.4f\n",res))
  return(res)
}
.wll_d_se<-function(u,ss) {
  x_j <- log(4*((cos(2*pi*ss$freq)-u)^2))
  return (pi^2/(6*sqrt(sum(x_j^2,na.rm=T))))
}
