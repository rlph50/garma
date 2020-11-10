# Evaluate the Whittle (log-) likelihood
# whittle.ggbr.obj - objective function to be minimised to get Whittle estimates.
# called from function "garma"
# @param par - the parameters to evaluate the function at
# @param params - other parameters - including the p, q, k, and scale parameters and (ss) the spectrum .
# @return The value of the objective at the point par.
.whittle.ggbr.obj<-function(par,params) {
  # objective function to be minimised for Whittle estimates
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k

  n_freq<-length(ss$freq)
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  u <- fd <- c()
  for (k1 in seq_len(k)) {
      u <- c(u,par[start])
      fd <- c(fd,par[start+1])
      start<-start+2
    }

  if (p>0) phi_vec <- (-par[start:(start+p-1)])      else phi_vec<-1
  if (q>0) theta_vec <- (par[(p+start):length(par)]) else theta_vec<-1

  cos_2_pi_f <- cos(2*pi*freq)
  mod_phi <- mod_theta <- 1
  if (p>0) mod_phi   <- .a_fcn(phi_vec,freq)
  if (q>0) mod_theta <- .a_fcn(theta_vec,freq)
  spec_den_inv <- mod_phi / mod_theta    # Inverse of spectral density

  for (k1 in seq_len(k)) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }

  # spec_den_inv[spec_den_inv==0] <- NA
  I_f <- spec*spec_den_inv
  res <- sum(I_f,na.rm=TRUE) - sum(log(spec_den_inv[spec_den_inv>0]),na.rm=TRUE)

  return(res)
}

.whittle.ggbr.obj.short<-function(par,params) {
  # Just the first part of the objective function - used for calculating the variance of the process.
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k

  n_freq<-length(ss$freq)
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  u <- fd <- c()
  for (k1 in seq_len(k)) {
    u <- c(u,par[start])
    fd <- c(fd,par[start+1])
    start<-start+2
  }

  if (p>0) phi_vec <- (-par[start:(start+p-1)])      else phi_vec<-1
  if (q>0) theta_vec <- (par[(p+start):length(par)]) else theta_vec<-1

  cos_2_pi_f <- cos(2*pi*freq)
  mod_phi <- mod_theta <- 1
  if (p>0) mod_phi   <- .a_fcn(phi_vec,freq)
  if (q>0) mod_theta <- .a_fcn(theta_vec,freq)
  spec_den_inv <- mod_phi / mod_theta    # Inverse of spectral density

  for (k1 in seq_len(k)) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }

  # spec_den_inv[spec_den_inv==0] <- NA
  I_f <- spec*spec_den_inv
  res <- sum(I_f,na.rm=TRUE)

  return(res)
}

.whittle.ggbr.grad<-function(par,params) {
  # objective function to be minimised for Whittle estimates
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k

  n_freq<-length(ss$freq)
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  u <- fd <- c()
  for (k1 in seq_len(k)) {
    u <- c(u,par[start])
    fd <- c(fd,par[start+1])
    start<-start+2
  }

  if (p>0) phi_vec <- (-par[start:(start+p-1)])      else phi_vec<-1
  if (q>0) theta_vec <- (par[(p+start):length(par)]) else theta_vec<-1

  cos_2_pi_f <- cos(2.0*pi*freq)
  mod_phi <- mod_theta <- 1
  if (p>0) mod_phi   <- .a_fcn(phi_vec,freq)
  if (q>0) mod_theta <- .a_fcn(theta_vec,freq)
  spec_den_inv <- mod_phi / mod_theta    # Inverse of spectral density

  for (k1 in seq_len(k)) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    temp <- (4.0*((cos_2_pi_f-u_k)^2.0))
    spec_den_inv <- spec_den_inv * (temp^fd_k)
  }

  spec_den_inv[spec_den_inv==0] <- NA
  I_f <- spec*spec_den_inv

  grad <- c()
  for (k1 in seq_len(k)) {
    grad_u <- sum(2*fd[k1]/(cos_2_pi_f-u[k1])*(1.0-I_f),na.rm=TRUE)
    grad_fd <- sum(-log((cos_2_pi_f-u[k1])*(cos_2_pi_f-u[k1])*4)*(1.0-I_f),na.rm=TRUE)
    grad <- c(grad,grad_u,grad_fd)
  }
  if (p>0) {
    phi_cos_vec <- .a_fcn_cos(phi_vec,freq)
    phi_sin_vec <- .a_fcn_sin(phi_vec,freq)
    for (p1 in seq_len(p)) {
      grad_phi <- 2.0*sum((phi_cos_vec*cos(p1*2*pi*freq)+phi_sin_vec*sin(p1*2*pi*freq)) / mod_phi * (1.0-I_f), na.rm=TRUE)
      grad <- c(grad,grad_phi)
    }
  }
  if (q>0) {
    theta_cos_vec <- .a_fcn_cos(theta_vec,freq)
    theta_sin_vec <- .a_fcn_sin(theta_vec,freq)
    for (q1 in seq_len(q)) {
      grad_theta <- 2.0*sum((theta_cos_vec*cos(q1*2*pi*freq)+theta_sin_vec*sin(q1*2*pi*freq)) / mod_theta * (I_f+1.0), na.rm=TRUE)
      grad <- c(grad,grad_theta)
    }
  }

  return(grad)
}

# params <- list(y=diff(co2_ts), orig_y=co2_ts, ss=stats::spectrum((diff(co2_ts)-mean(diff(co2_ts))),plot=FALSE,detrend=TRUE,demean=TRUE,method='pgram',taper=0,fast=FALSE),
#                p=0,q=0,d=1,k=1,
#                include.mean=FALSE,
#                method='Whittle',
#                est_mean=FALSE,
#                scale=sd(diff(co2_ts)),m_trunc=50)
# par <- c(cos(0.7), 0.2436360904)
#
# .whittle.ggbr.grad(par,params)
# pracma::grad(.whittle.ggbr.obj,par,params=params)
#
# nloptr::check.derivatives(par,.whittle.ggbr.obj,.whittle.ggbr.grad,params=params)

# params <- list(y=diff(co2_ts), orig_y=co2_ts, ss=stats::spectrum((diff(co2_ts)-mean(diff(co2_ts))),plot=FALSE,detrend=TRUE,demean=TRUE,method='pgram',taper=0,fast=FALSE),
#                p=3,q=0,d=1,k=1,
#                include.mean=FALSE,
#                method='Whittle',
#                est_mean=FALSE,
#                scale=sd(diff(co2_ts)),m_trunc=50)
# par <- c(cos(0.7), 0.2436360904, 0.2, 0.1, 0.3)
#
# nloptr::check.derivatives(par,.whittle.ggbr.obj,.whittle.ggbr.grad,params=params)
