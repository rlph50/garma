#' Evaluate the Whittle (log-) likelihood
#' internal_whittle_garma_obj - objective function to be minimised to get Whittle estimates.
#' called from function "garma"
#' @param par - the parameters to evaluate the function at
#' @param params - other parameters - including the p, q, k, and scale parameters and (ss) the spectrum .
#' @return The value of the objective at the point par.
#' @noRd
internal_whittle_garma_obj <- function(par, params) {
  # objective function to be minimised for Whittle estimates
  ss <- params$ss
  spec <- ss$spec
  n_freq <- length(spec)

  spec_den_inv <- internal_calc_garma_spec_den_inv(par, params)

  ## I have checked the source for for "spectrum" which is used to calculate the "spec"
  ## The I/f calc includes spectrum and the "f" we calc as above.
  ## Source code for "spectrum" however divides by n but not by 2*pi. Thus we make allowances for that here.

  I_f <- spec * spec_den_inv / (2 * pi * n_freq)
  res <- sum(I_f, na.rm = TRUE) - sum(log(spec_den_inv[spec_den_inv > 0]), na.rm = TRUE)
  if (is.na(res) | is.infinite(res)) res <- 1.0e200

  return(res)
}

internal_whittle_garma_obj_short <- function(par, params) {
  # Just the first part of the objective function.
  # This is the "S" function of Giraitis, Hidalgo & Robinson 2001.
  ss <- params$ss
  spec <- ss$spec
  n_freq <- length(spec)

  spec_den_inv <- internal_calc_garma_spec_den_inv(par, params)

  ## I have checked the source for for "spectrum" which is used to calculate the "spec"
  ## The I/f calc includes spectrum and the "f" we calc. f includes the factor 1/(2 pi)
  ## Source code for "spectrum" however divides by n but not by 2*pi. Thus we make allowances for that here.

  I_f <- spec * spec_den_inv / (2 * pi * n_freq)
  res <- sum(I_f, na.rm = TRUE)

  return(res)
}

internal_calc_garma_spec_den_inv <- function(par, params) {
  # This is the "k" function of Giraitis, Hidalgo & Robinson 2001.
  # true spec den is this times sigma^2/(2*pi)
  ss <- params$ss
  p <- params$p
  q <- params$q
  k <- params$k

  n_freq <- length(ss$freq)
  freq <- ss$freq[1:n_freq]
  spec <- ss$spec[1:n_freq]

  start <- 1
  if (is.null(params$periods)) {
    u <- fd <- numeric(0)
    for (k1 in seq_len(k)) {
      u <- c(u, par[start])
      fd <- c(fd, par[start + 1])
      start <- start + 2
    }
  } else {
    k <- length(params$periods)
    u <- cos(2 * pi / params$periods)
    fd <- par[start:(start + k - 1)]
    start <- start + k
  }

  if (p > 0) phi_vec <- (-par[start:(start + p - 1)]) else phi_vec <- 1
  if (q > 0) theta_vec <- (-par[(p + start):length(par)]) else theta_vec <- 1

  cos_2_pi_f <- cos(2 * pi * freq)
  mod_phi <- mod_theta <- 1
  if (p > 0) mod_phi <- internal_a_fcn(phi_vec, freq)
  if (q > 0) mod_theta <- internal_a_fcn(theta_vec, freq)

  ## As per Giraitis, Hidalgo & Robinson 2001, we do not include sigma^2/2pi factor in the spectral density.
  ## The below is the inverse of the spectral density.
  spec_den_inv <- mod_phi / mod_theta

  for (k1 in seq_len(k)) {
    u_k <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4 * ((cos_2_pi_f - u_k)^2))^fd_k
  }

  return(spec_den_inv)
}


internal_garma_spec_den_0 <- function(par, params) {
  # Just the first part of the objective function - used for calculating the variance of the process.
  ss <- params$ss
  p <- params$p
  q <- params$q
  k <- params$k
  freq <- 0

  start <- 1
  if (is.null(params$periods)) {
    u <- fd <- numeric(0)
    for (k1 in seq_len(k)) {
      u <- c(u, par[start])
      fd <- c(fd, par[start + 1])
      start <- start + 2
    }
  } else {
    k <- length(params$periods)
    u <- cos(2 * pi / params$periods)
    fd <- par[start:(start + k - 1)]
    start <- start + k
  }

  if (p > 0) phi_vec <- (-par[start:(start + p - 1)]) else phi_vec <- 1
  if (q > 0) theta_vec <- (par[(p + start):length(par)]) else theta_vec <- 1

  cos_2_pi_f <- cos(2 * pi * freq)
  mod_phi <- mod_theta <- 1
  if (p > 0) mod_phi <- internal_a_fcn(phi_vec, freq)
  if (q > 0) mod_theta <- internal_a_fcn(theta_vec, freq)
  spec_den_inv <- mod_phi / mod_theta # Inverse of spectral density

  for (k1 in seq_len(k)) {
    u_k <- u[k1]
    fd_k <- fd[k1]
    spec_den_inv <- spec_den_inv * (4 * ((cos_2_pi_f - u_k)^2))^fd_k
  }
  return(1 / spec_den_inv)
}

# for vcov calcs...
internal_whittle_garma_omega <- function(par, params) {
  # This returns the estimate of the "omega" matrix from GHR (2001). Used for determining se's.

  calc_grad_log_k <- function(par, params) {
    calc_log_k <- function(par, params) {
      k <- internal_calc_garma_spec_den_inv(par, params)
      return(log(k))
    }
    calc_log_k_j <- function(par, params, j) {
      k <- internal_calc_garma_spec_den_inv(par, params)
      return(log(k[j]))
    }

    # log_k <- calc_log_k(par,params)
    n_freq <- length(params$ss$spec)
    ## Need to be able to evaluate spec den at a frequency.

    res <- array(data = 0, dim = c(n_freq, length(par)))
    # find the "q" frequencies
    q_freq <- numeric(0)
    if (is.null(params$periods)) {
      for (k1 in seq_len(params$k)) {
        idx <- (k1 - 1) * 2 + 1
        u <- par[idx] # cosine of ggbr freq
        f <- acos(u) / (2 * pi) # Ggbr frequency
        ff <- abs(params$ss$freq - f)
        q_freq <- c(q_freq, which.min(ff))
      }
    } else {
      for (period in params$periods) {
        f <- 1/period # Ggbr frequency
        ff <- abs(params$ss$freq - f)
        q_freq <- c(q_freq, which.min(ff))
      }
    }

    # Now as per GHR sum over the non-ggbr frequencies
    freq_list <- setdiff(1:n_freq, q_freq)
    for (j in freq_list) {
      x <- pracma::grad(calc_log_k_j, par, params = params, j = j)
      res[j, ] <- x
    }

    return(res)
  }

  n_freq <- length(params$ss$spec)
  grad_log_k <- calc_grad_log_k(par, params)
  omega <- array(0, dim = c(length(par), length(par))) # starting values
  for (j in 1:n_freq) omega <- omega + (grad_log_k[j, ] %*% t(grad_log_k[j, ]))
  return(omega / n_freq)
}
