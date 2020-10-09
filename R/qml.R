# Estimate a Ggbr model using a QML (Quasi-Maximum-Likelihood) method.
# qml_matricies - generate the statespace matrices needed for QML method.
.qml_matricies<-function(pars,p,q,k,include.mean,m_trunc) {
  # add in ARMA factors
  start1 <- 1
  if (include.mean) start1<-2
  if (k==1) {
    u     <- pars[start1]
    d     <- pars[start1+1]
    start1 <- start1+2
  } else u<-d<-0.0
  sigma2<-1

  if (p>0) phi_vec   <- c(1,-(pars[start1:(start1+p-1)] ))      else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,-(pars[(p+start1):(p+q+start1-1)])) else theta_vec <- 1
  # initial ggbr factor
  if (k==1) si<-as.vector(.ggbr.coef(m_trunc+1,d,u))
  else si<-as.vector(c(1,rep(0,m_trunc)))
  arma_filter   <- signal::Arma(b = theta_vec, a = phi_vec)
  si2   <- sqrt(sigma2) * (signal::filter(arma_filter, si))

  si_sq <- (si2 %*% t(si2))
  Tt    <- rbind(cbind(rep(0,m_trunc),diag(m_trunc)),rep(0,m_trunc+1))
  Zt    <- matrix(c(1, rep(0,m_trunc)), ncol = m_trunc+1)
  ct    <- matrix(0)
  dt    <- matrix(0, nrow = m_trunc+1)
  GGt   <- matrix(0)
  H     <- matrix(c(si2[1:(m_trunc+1)]), nrow = m_trunc+1)
  HHt   <- H %*% t(H)
  a0    <- c(rep(0,m_trunc+1))
  entry <- vector()
  for(j in 1:(m_trunc-1)){
    entry[j] <- sum(diag(si_sq[-1:-j,-(m_trunc+2-j):-(m_trunc+1)]))
  }
  P0 <- stats::toeplitz(c(sum(diag(si_sq)),entry[1:(m_trunc-1)],si_sq[1,m_trunc+1]))

  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt, HHt = HHt))
}

.qml.ggbr.obj <- function(par,params) {
  # objective function to be minimised for QML estimates
  y  <- params$y
  yt <- rbind(y)
  p  <- params$p
  q  <- params$q
  k  <- params$k
  include.mean <- params$include.mean
  m_trunc <- params$m_trunc

  sp  <- .qml_matricies(par,p,q,k,include.mean,m_trunc)
  if (include.mean) beta0 <- par[1] else beta0<-0

  ans <- FKF::fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = (yt-beta0))
  return(-ans$logLik)
}

.qml.ggbr.se2 <- function(pars,params) {
  # objective function to be minimised for QML estimates
  y  <- params$y
  yt <- rbind(y)
  p  <- params$p
  q  <- params$q
  k  <- params$k
  include.mean <- params$include.mean
  m_trunc <- params$m_trunc

  sp  <- .qml_matricies(pars,p,q,k,include.mean,m_trunc)
  if (include.mean) beta0 <- pars[1] else beta0<-0

  ans <- FKF::fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt, Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = (yt-beta0))

  return(sum(ans$vt[1,]^2/ans$Ft[1,1,],na.rm=TRUE))
}
