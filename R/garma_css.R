
# #' Estimate a Ggbr model using a CSS (Conditional Sum-of-Squares) method.
# #' css.ggbr.obj - objective function to be minimised to get CSS estimates.
# #' called from function "garma"
# #' @param par - the parameters to evaluate the function at
# #' @param params - other parameters - including the p, q, k, and scale parameters and (ss) the spectrum .
# #' @return The value of the objective at the point par.
.css.ggbr.obj<-function(par,params) {
  # Objective function to be minimised for CSS estimates
  y <- params$y
  p <- params$p
  q <- params$q
  k <- params$k
  include.mean <- params$include.mean

  if (include.mean) {
    beta0  <- par[1]
    start  <- 2
  }
  else {
    beta0  <- 0
    start  <- 1
  }

  u <- c()
  fd <- c()
  if (k>0) for (k1 in 1:k) {
    u     <- c(u,par[start])
    fd    <- c(fd,par[start+1])
    start <- start+2
  }

  y_dash <- y-beta0
  if (p>0) phi_vec   <- c(1,-(par[start:(start+p-1)] ))     else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,(par[(p+start):(length(par))])) else theta_vec <- 1

  arma_filter   <- signal::Arma(a = theta_vec, b = phi_vec)
  eps           <- signal::filter(arma_filter, y_dash)
  if (k>0) for (k1 in 1:k) {
    ggbr_filter <- signal::Arma(b = 1, a = .ggbr.coef(length(y_dash),fd[k1],u[k1]))
    eps         <- signal::filter(ggbr_filter, eps)
  }

  ret <- sum(eps^2,na.rm=TRUE)
  if (!is.finite(ret)|is.na(ret)) ret<-(1e500)

  #return(0.5*log(ret))
  return(ret)
}
