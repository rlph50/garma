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

  n_freq<-length(ss$freq)#-1
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  u <- fd <- c()
  if (k>0) for (k1 in 1:k) {
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

  if (k>0) for (k1 in 1:k) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }

  # spec_den_inv[spec_den_inv==0] <- NA
  I_f <- spec*spec_den_inv
  res <- sum(I_f,na.rm=TRUE) #- sum(log(spec_den_inv[spec_den_inv>0]),na.rm=TRUE)

  return(res)
}

.whittle.ggbr.obj.short<-function(par,params) {
  # Just the first part of the objective function - used for calculating the variance of the process.
  ss <- params$ss
  p  <- params$p
  q  <- params$q
  k  <- params$k

  n_freq<-length(ss$freq)#-1
  freq<-ss$freq[1:n_freq]
  spec<-ss$spec[1:n_freq]

  start=1
  u <- fd <- c()
  if (k>0) for (k1 in 1:k) {
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

  if (k>0) for (k1 in 1:k) {
    u_k  <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4*((cos_2_pi_f-u_k)^2))^fd_k
  }

  # spec_den_inv[spec_den_inv==0] <- NA
  I_f <- spec*spec_den_inv
  res <- sum(I_f,na.rm=TRUE)

  return(res)
}

