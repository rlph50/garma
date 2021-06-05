# Evaluate the Whittle (log-) likelihood
# whittle.ggbr.obj - objective function to be minimised to get Whittle estimates.
# called from function "garma"
# @param par - the parameters to evaluate the function at
# @param params - other parameters - including the p, q, k, and scale parameters and (ss) the spectrum .
# @return The value of the objective at the point par.
.whittle.garma.obj<-function(par,params) {
  # objective function to be minimised for Whittle estimates
  ss <- params$ss
  spec<-ss$spec
  n_freq <- length(spec)

  spec_den_inv <- .calc_garma_spec_den_inv(par,params)

  ## I have checked the source for for "spectrum" which is used to calculate the "spec"
  ## The I/f calc includes spectrum and the "f" we calc as above.
  ## Source code for "spectrum" however divides by n but not by 2*pi. Thus we make allowances for that here.

  I_f <- spec*spec_den_inv / (2*pi*n_freq)
  res <- sum(I_f,na.rm=TRUE) - sum(log(spec_den_inv[spec_den_inv>0]),na.rm=TRUE)
  if (is.na(res)|is.infinite(res)) res <- 1.0e200

  return(res)
}

.whittle.garma.obj.short<-function(par,params) {
  # Just the first part of the objective function.
  # This is the "S" function of Giraitis, Hidalgo & Robinson 2001.
  ss <- params$ss
  spec<-ss$spec
  n_freq <- length(spec)

  spec_den_inv <- .calc_garma_spec_den_inv(par,params)

  ## I have checked the source for for "spectrum" which is used to calculate the "spec"
  ## The I/f calc includes spectrum and the "f" we calc. f includes the factor 1/(2 pi)
  ## Source code for "spectrum" however divides by n but not by 2*pi. Thus we make allowances for that here.

  I_f <- spec*spec_den_inv / (2*pi*n_freq)
  res <- sum(I_f,na.rm=TRUE)

  return(res)
}

## @export
## calc_garma_spec_den_inv<-function(par,params) return(.calc_garma_spec_den_inv(par,params))

.calc_garma_spec_den_inv<-function(par,params) {
  # This is the "k" function of Giraitis, Hidalgo & Robinson 2001.
  # true spec den is this times sigma^2/(2*pi)
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k

  n_freq<-length(ss$freq)
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]
  # spec <- params$true_spec

  start=1
  u <- fd <- c()
  for (k1 in seq_len(k)) {
    u <- c(u,par[start])
    fd <- c(fd,par[start+1])
    start<-start+2
  }

  if (p>0) phi_vec <- (-par[start:(start+p-1)])      else phi_vec<-1
  if (q>0) theta_vec <- (-par[(p+start):length(par)]) else theta_vec<-1

  cos_2_pi_f <- cos(2*pi*freq)
  mod_phi <- mod_theta <- 1
  if (p>0) mod_phi   <- .a_fcn(phi_vec,freq)
  if (q>0) mod_theta <- .a_fcn(theta_vec,freq)

  ## As per Giraitis, Hidalgo & Robinson 2001, we do not include sigma^2/2pi factor in the spectral density.
  ## The below is the inverse of the spectral density.
  spec_den_inv <- mod_phi / mod_theta

  for (k1 in seq_len(k)) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }

  return(spec_den_inv)
}


.garma_spec_den_0 <-function(par,params) {
  # Just the first part of the objective function - used for calculating the variance of the process.
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k
  freq <- 0

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
  spec_den_inv <- mod_phi / mod_theta   # Inverse of spectral density

  for (k1 in seq_len(k)) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }
  return(1/spec_den_inv)
}

# for vcov calcs...
.whittle_garma_omega<-function(par,params) {
  # This returns the estimate of the "omega" matrix from GHR (2001). Used for determining se's.

  calc_grad_log_k<-function(par,params) {
    calc_log_k<-function(par,params) {
      k <- .calc_garma_spec_den_inv(par,params)
      return(log(k))
    }
    calc_log_k_j<-function(par,params,j) {
      k <- .calc_garma_spec_den_inv(par,params)
      return(log(k[j]))
    }

    #log_k <- calc_log_k(par,params)
    n_freq <- length(params$ss$spec)
    ## Need to be able to evaluate spec den at a frequency.

    res <- array(data=0,dim=c(n_freq,length(par)))
    # find the "q" frequencies
    q_freq<-numeric(0)
    for (k1 in seq_len(params$k)) {
      idx <- (k1-1)*2+1
      u <- par[idx]       # cosine of ggbr freq
      f <- acos(u)/(2*pi) # Ggbr frequency
      ff <- abs(params$ss$freq-f)
      q_freq <- c(q_freq,which.min(ff))
    }

    # Now as per GHR sum over the non-ggbr frequencies
    freq_list <- setdiff(1:n_freq,q_freq)
    for (j in freq_list) {
      x <- pracma::grad(calc_log_k_j,par,params=params,j=j)
      res[j,]<-x
    }

    return(res)
  }

  n_freq <- length(params$ss$spec)
  grad_log_k <- calc_grad_log_k(par,params)
  omega <- array(0,dim=c(length(par),length(par)))  # starting values
  for (j in 1:n_freq) omega <- omega + (grad_log_k[j,]%*%t(grad_log_k[j,]))
  return(omega/n_freq)
}



#
#
# .whittle.ggbr.grad<-function(par,params) {
#   # objective function to be minimised for Whittle estimates
#   ss <- params$ss
#   p  <- params$p
#   q  <- params$q
#   k  <- params$k
#
#   n_freq<-length(ss$freq)
#   freq<-ss$freq[1:n_freq]
#   spec<-ss$spec[1:n_freq]
#   # spec <- params$true_spec
#
#   start=1
#   u <- fd <- c()
#   for (k1 in seq_len(k)) {
#     u <- c(u,par[start])
#     fd <- c(fd,par[start+1])
#     start<-start+2
#   }
#
#   if (p>0) phi_vec <- (-par[start:(start+p-1)])      else phi_vec<-1
#   if (q>0) theta_vec <- (par[(p+start):length(par)]) else theta_vec<-1
#
#   cos_2_pi_f <- cos(2.0*pi*freq)
#   mod_phi <- mod_theta <- 1
#   if (p>0) mod_phi   <- .a_fcn(phi_vec,freq)
#   if (q>0) mod_theta <- .a_fcn(theta_vec,freq)
#   spec_den_inv <- mod_phi / mod_theta    # Inverse of spectral density
#
#   for (k1 in seq_len(k)) {
#     u_k  <- u[k1]
#     fd_k <- fd[k1]
#     temp <- (4.0*((cos_2_pi_f-u_k)^2.0))
#     spec_den_inv <- spec_den_inv * (temp^fd_k)
#   }
#
#   spec_den_inv[spec_den_inv==0] <- NA
#   I_f <- spec*spec_den_inv
#
#   grad <- c()
#   for (k1 in seq_len(k)) {
#     grad_u <- sum(2*fd[k1]/(cos_2_pi_f-u[k1])*(1.0-I_f),na.rm=TRUE)
#     grad_fd <- sum(-log((cos_2_pi_f-u[k1])*(cos_2_pi_f-u[k1])*4)*(1.0-I_f),na.rm=TRUE)
#     grad <- c(grad,grad_u,grad_fd)
#   }
#   if (p>0) {
#     phi_cos_vec <- .a_fcn_cos(phi_vec,freq)
#     phi_sin_vec <- .a_fcn_sin(phi_vec,freq)
#     for (p1 in seq_len(p)) {
#       grad_phi <- 2.0*sum((phi_cos_vec*cos(p1*2*pi*freq)+phi_sin_vec*sin(p1*2*pi*freq)) / mod_phi * (1.0-I_f), na.rm=TRUE)
#       grad <- c(grad,grad_phi)
#     }
#   }
#   if (q>0) {
#     theta_cos_vec <- .a_fcn_cos(theta_vec,freq)
#     theta_sin_vec <- .a_fcn_sin(theta_vec,freq)
#     for (q1 in seq_len(q)) {
#       grad_theta <- 2.0*sum((theta_cos_vec*cos(q1*2*pi*freq)+theta_sin_vec*sin(q1*2*pi*freq)) / mod_theta * (I_f+1.0), na.rm=TRUE)
#       grad <- c(grad,grad_theta)
#     }
#   }
#
#   return(grad)
# }

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
