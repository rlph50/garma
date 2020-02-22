#' Estimate the parmaeters of a GARMA model.
#'
#' The garma function is the main function for the tsggbr package. Depending on the parameters it will
#' calculate the parameter estimates for the GARMA process, and if available the standard errors (se's)
#' for those parameters.
#'
#' The GARMA model is specified as
#' \deqn{\displaystyle{\phi(B)\prod_{i=1}^{k}(1-2u_{i}B+B^{2})^{d_{i}}(Y_{t}-\mu)= \theta(B) \epsilon _{t}}}{\prod(i=1 to k) (1-2u(i)B+B^2)^d(i) \phi(B) Y(t) = \theta(B) \epsilon(t)}
#'
#' where
#' \itemize{
#' \item \eqn{\phi(B)}{\phi(B)} represents the short-memory Autoregressive component of order p,
#' \item \eqn{\theta(B)}{\theta(B)} represents the short-memory Moving Average component of order q,
#' \item \eqn{(1-2u_{i}B+B^{2})^{d_{i}}}{(1-2u(i)B+B^2)^d(i)} represents the long-memory Gegenbauer component (there may in general be k of these),
#' \item \eqn{Y_{t}}{Y(t)} represents the observed process,
#' \item \eqn{\epsilon_{t}}{\epsilon(t)} represents the random component of the model - these are assumed to be uncorrelated but identically distributed variates.
#'       Generally the routines in this package will work best if these have an approximate Gaussian distribution.
#' \item \eqn{B}{B} represents the Backshift operator, defined by \eqn{B Y_{t}=Y_{t-1}}{B Y(t) = Y(t-1)}.
#' }
#' when k=0, then this is just a short memory model as fit by the stats "arima" function.
#'
#' @param py (num) This should be a numeric vector representing the process to estimate. A minimum length of 96 is required.
#' @param order (list) This should be a list (similar to the stats::arima order parameter) which will give the order of the process to fit.
#'     The format should be list(p,d,q) where p, d, and q are all positive integers. p represents the degree of the
#'     autoregressive process to fit, q represents the order of the moving average process to fit and d is the (integer)
#'     differencing to apply prior to any fitting.
#' @param k (int) This is the number of (multiplicative) Gegenbauer terms to fit. Only 0 or 1 are allowed in this version.
#' @param include.mean (bool) A boolean value indicating whether a mean should be fit. Note if you have any differencing, then
#'     it generally does not make sense to fit a mean term.
#' @param method (character) This defines the estimation method for the routine. The valid values are 'CSS', 'Whittle', 'QML' and 'WLL'.
#'     The default (Whittle) will generally return very accurate estimates quite quickly, provided the asumption of a Gaussian
#'     distribution is even approximately correct, and is probably the method of choice for most users. For the theory behind this, refer Giraitis et. al. (2001)
#'     'CSS' is a conditional 'sum-of-squares' technique and can be quite slow. Reference: Chung (1996).
#'     'QML' is a Quasi-Maximum-Likelihood technique, and can also be quite slow. Reference Dissanayake (2016).
#'     'WLL' is a new technique which appears to work well even if the \eqn{\epsilon_{t}}{\epsilon(t)} are highly skewed and/or have heavy tails (skewed and/or lepto-kurtic).
#'     However the asymptotic theory for the WLL method is not complete and so standard errors are not available for most parameters.
#' @param allow_neg_d (bool) A boolean value indicating if a negative value is allowed for the fractional differencing component
#'     of the Gegenbauer term is allowed. This can be set to FALSE to force the routine to find a positive value.
#' @param maxeval (int) the maximum function eveluations to be allowed during each optimisation.
#' @param opt.method (character) This names the optimisation method used to find the parameter estimates. The default is to use the built-in
#'     R algorithm called 'optim'. For some data or some models, however, other methods may work well. Other allowed values are
#'     \itemize{
#'     \item cobyla algorithm in package nloptr
#'     \item directL algorithm in package nloptr
#'     \item BBoptim from package BB
#'     \item psoptim from package pso
#'     \item hjkb from dfoptim package
#'     \item nmkb from dfoptim package
#'     \item best - this option evaluates all the above options in turn and picks the one which finds the lowest value of the objective. This can be quite time consuming to run.
#'     }
#' @param m_trunc Used for the QML estimation method. This defines the AR-truncation point when evaluating the likelihood function. Refer to Dissanayake et. al. (2016) for details.
#' @return An S3 object of class "ggbr_model".
#'
#' \preformatted{
#' References:
#' C Chung. A generalized fractionally integrated autoregressive moving-average process. Journal of Time Series Analysis, 17(2):111–140, 1996.
#' G Dissanayake, S Peiris, and T Proietti. State space modelling of Gegenbauer processes with long memory. Computational Statistics and Data Analysis, 100:115–130, 2016.
#' L Giraitis, J Hidalgo, and P Robinson. Gaussian estimation of parametric spectral density with unknown pole. The Annals of Statistics, 29(4):987–1023, 2001.
#' }

garma<-function(py,
                order=list(0,0,0),
                k=1,
                include.mean=TRUE,
                method='Whittle',
                allow_neg_d=TRUE,
                maxeval=10000,
                opt.method='optim',
                m_trunc=50) {
  if (length(py)<96)
    stop('y should have at least 96 observations.')
  if (is.data.frame(py)) {
    if (ncol(py)>1)
      stop('y should be a numeric vector - not an entire data frame. Please select a single column and try again.')
    else py<-py[,1]
  }
  py<-as.numeric(py)
  if (!is.numeric(py))
    stop('y should be numeric.')
  if (length(order)!=3)
    stop('order parameter must be a 3 integers only.')
  if ((k!=0)&(k!=1))
    stop('Sorry. Only k=0 or k=1 is supported for now.')
  allowed_methods <- c('CSS','Whittle','WLL','QML')
  if (!method%in%allowed_methods)
    stop('Method must be one of CSS, Whittle, QML or WLL.')

  allowed_optimisations <- c('optim','cobyla','directL','BBoptim','psoptim','hjkb','nmkb','best')
  optimisation_packages <- c('optim'='stats','cobyla'='nloptr','directL'='nloptr','BBoptim'='BB','psoptim'='pso','hjkb'='dfoptim','nmkb'='dfoptim','best'='stats')
  if (!opt.method%in%allowed_optimisations)
    stop('\nSorry - supported packages are:\n', paste0(allowed_optimisations,'\n'))

  if (!is.installed(optimisation_packages[[opt.method]]))
    stop(sprintf('package %s needs to be installed to use method ',optimisation_packages[[opt.method]],opt.method))

  #
  # First calc parameter  estimates
  p=as.integer(order[1])
  d=as.integer(order[2])
  q=as.integer(order[3])
  storage.mode(py) <- 'double'
  if (d>0) y<-diff(py,differences=d) else y<-py
  mean_y <- mean(y)
  sd_y   <- sd(y)
  ss<-spectrum((y-mean_y)/sd_y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)

  # Now set up params and lb (lower bounds) and ub (upper bounds)
  n_pars   <- 0
  pars     <- numeric(0)
  lb       <- numeric(0)
  ub       <- numeric(0)
  lb_finite<- numeric(0)
  ub_finite<- numeric(0)

  mean_methods <- c('QML','CSS')
  if (include.mean&method%in%mean_methods) {
    n_pars    <- n_pars+1
    mean_y    <- mean(y)
    pars      <- c(pars,mean_y)
    lb_finite <- c(lb_finite,ifelse(mean_y<0,2*mean_y,-2*mean_y))
    ub_finite <- c(ub_finite,ifelse(mean_y<0,-2*mean_y,2*mean_y))
    lb        <- c(lb,-Inf)
    ub        <- c(ub,Inf)
  }
  if (k==1) {
    n_pars    <- n_pars+2
    start_u   <-cos(2*pi*ss$freq[which.max(ss$spec)])
    if (start_u< (-1)) start_u<-0
    pars      <- c(pars,start_u,0.25)
    lb        <- c(lb,0.0,ifelse(allow_neg_d,-1,0))
    ub        <- c(ub,1.0,1.0)
    lb_finite <- c(lb_finite,0.0,ifelse(allow_neg_d,-1,0))
    ub_finite <- c(ub_finite,1.0,1.0)
  }
  n_pars    <- n_pars + p + q
  methods_to_estimate_var <- c('WLL')
  if (p+q>0) {
    a    <- arima(y,order=c(p,0,q),include.mean=FALSE)
    pars <- c(pars,a$coef)
    if (method%in%methods_to_estimate_var) pars <- c(pars,a$sigma2)
    if (p==1&q==0) { # special limits for AR(1)
      lb<-c(lb,-1)
      if (method%in%methods_to_estimate_var) lb <- c(lb, 1e-10)
      lb_finite <- c(lb_finite,-1)
      if (method%in%methods_to_estimate_var) lb_finite <- c(lb_finite, 1e-10)
      ub<-c(ub,1)
      if (method%in%methods_to_estimate_var) ub<-c(ub,Inf)
      ub_finite <- c(ub_finite,1)
      if (method%in%methods_to_estimate_var) ub_finite <- c(ub_finite,2*var(y))
    } else {
      lb<- c(lb,rep(-Inf,p+q))
      if (method%in%methods_to_estimate_var) lb <- c(lb,1e-10)
      ub<- c(ub,rep(Inf,p+q))
      if (method%in%methods_to_estimate_var) ub <- c(ub,Inf)
      lb_finite <- c(lb_finite,rep(-10,p+q))
      if (method%in%methods_to_estimate_var) lb_finite <- c(lb_finite,1e-10)
      ub_finite <- c(ub_finite,rep(10,p+q))
      if (method%in%methods_to_estimate_var) ub_finite <- c(ub_finite,2*var(y))
    }
  }

  # create a list of all possible params any 'method' might need. The various objective functions can extract the parameters which are relevant to that method.
  params <- list(y=y, ss=ss, p=p,q=q,k=k,include.mean=include.mean,scale=sd_y,m_trunc=m_trunc)
  message <- c()

  # First we make a first pass at optimisation using "optim".
  # If the method chosen is optim then that finishes things; but otherwise the solution found becomes the starting point for the next optimisation.
  fcns <- list('CSS'=css.ggbr.obj,'Whittle'=whittle.ggbr.obj,'QML'=qml.ggbr.obj,'WLL'=wll.ggbr.obj)
  fit <- optim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params,
               hessian=TRUE,method="L-BFGS-B",control=list(maxit=maxeval,factr=1e-25))
  if (fit$convergence==52) {   # Error in line search, then try again
    #cat('1st optimisation using optim - failed in line search.\n2nd optimisation using nloptr.\n')
    fit2<-lbfgs(x0=fit$par, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxeval=maxeval,xtol_rel=1e-8))
    if (fit2$value<fit$value&fit2$convergence>=0) {
      fit<-fit2
      #cat('Using 2nd fit.\n')
    }
  }
  if (fit$convergence>=0) pars <- fit$par
  hh  <- fit$hessian

  if (opt.method=='cobyla') {
    tryCatch(fit <- cobyla(pars,fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxeval=maxeval,xtol_rel=1e-10)),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='directL') {
    tryCatch(fit <- directL(fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxeval=maxeval,xtol_rel=1e-10)),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='BBoptim') {
    pars[1] <- pars[1]-0.1
    tryCatch(fit <- BB::BBoptim(par=pars, fcns[[method]], lower=lb, upper=ub, control=list(trace=FALSE,maxit=maxeval,ftol=1e-15,gtol=1e-8),
                                params=params,quiet=TRUE),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='psoptim') {
    tryCatch(fit <- psoptim(par=pars, fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxit=maxeval)),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='hjkb') {
    tryCatch(fit <- dfoptim::hjkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='nmkb') {
    tryCatch(fit <- dfoptim::nmkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
             error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
    )
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else if (opt.method=='best') {
    message <- c()
    fit.optim<-fit.cobyla<-fit.directL<-fit.bboptim<-fit.psoptim<-fit.hjkb<-fit.nmkb<-list(value=Inf,message=c(message,cond),convergece=999,par=pars) #default to set environment
    tryCatch(
      fit.optim   <- optim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params,
                           hessian=TRUE,method="L-BFGS-B",control=list(maxit=maxeval,factr=1e-25)),
      error=function(cond) {fit.optim<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.cobyla  <- cobyla(pars,fcns[[method]], lower=lb, upper=ub, params=params,control=list(maxeval=maxeval,xtol_rel=1e-10)),
             error=function(cond) {fit.cobyla<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.directL <- directL(fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params,control=list(maxeval=maxeval,xtol_rel=1e-10)),
             error=function(cond) {fit.directL<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.bboptim <- BB::BBoptim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, control=list(trace=FALSE,maxit=maxeval,ftol=1e-15,gtol=1e-8),
                                        params=params,quiet=TRUE),
             error=function(cond) {fit.bboptim<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.psoptim <- psoptim(par=pars, fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxit=maxeval)),
             error=function(cond) {fit.psoptim<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.hjkb    <- dfoptim::hjkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
             error=function(cond) {fit.hjkb<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )
    tryCatch(fit.nmkb    <- dfoptim::nmkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
             error=function(cond) {fit.nmkb<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
    )

    fit_values  <- list('optim'=fit.optim$value,'cobyla'=fit.cobyla$value,'directL'=fit.directL$value,'bboptim'=fit.bboptim$value,'psoptim'=fit.psoptim$value,'hjkb'=fit.hjkb$value,'nmkb'=fit.nmkb$value)
    best_method <- names(fit_values)[which.min(fit_values)]
    fits        <- list('optim'=fit.optim,'cobyla'=fit.cobyla,'directL'=fit.directL,'bboptim'=fit.bboptim,'psoptim'=fit.psoptim,'hjkb'=fit.hjkb,'nmkb'=fit.nmkb)
    fit         <- fits[[best_method]]
    hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
  }

  if (fit$convergence>=0&method!='WLL') {
    # Next, find the c <- alc se's for coefficients
    start<-1
    se<-c()   # default to set this up in the right environment
    if (include.mean) start<-2
    nn <- nrow(hh)
    tryCatch({
      if (abs(det(hh))<1e-50) {
        nn     <- nrow(hh)-1
        hh     <- hh[1:nn,1:nn]
        hh_det <- abs(det(hh))
        if (hh_det>1e-50) se <- c(sqrt(diag(solve(hh*length(y)))*length(y)),NA)
        else {
          nn     <- nrow(hh)-1
          hh     <- hh[2:(nn+1),2:(nn+1)]
          hh_det <- abs(det(hh))
          if (hh_det>1e-50) se <- c(NA,sqrt(diag(solve(hh*length(y)))*length(y)),NA) else se <- rep(NA,length(fit$par))
        }
      } else se <- sqrt(diag(solve(hh*length(y)))*length(y))
    },
    error=function(cond) {se <- rep(NA,length(fit$par));message<-c(message,cond)}
    )

    if (length(se)<length(fit$par)) se<-c(se,NA)
  }
  if (method=='WLL') {
    se<-rep(NA,length(par))
    if (k>0) se[2] <- bnp_d_se(fit$par[1],ss)
  }
  nm<-list()
  if (include.mean&method%in%mean_methods) nm <- c(nm,'Intercept')
  if (k==1) nm<-c(nm,'u','fd')
  if (p>0) nm<-c(nm,paste0('AR',1:p))
  if (q>0) nm<-c(nm,paste0('MA',1:q))

  if (method=='WLL') {
    # adjust sigma2 for bias...
    fit$par[length(fit$par)] <- fit$par[length(fit$par)]/(2*pi) * exp(-digamma(1))
  }

  n_coef <- ifelse(method%in%methods_to_estimate_var,length(fit$par)-1,length(fit$par))  # ignore var on end if it is there
  coef <- t(matrix(round(c(fit$par[1:n_coef],se[1:n_coef]),4),nrow=n_coef))
  colnames(coef) <- nm
  rownames(coef) <- c('coef','se')
  if (method=='WLL')          sigma2 <- fit$par[length(fit$par)]
  else if (method=='QML')     sigma2 <- qml.ggbr.se2(fit$par, params=params)
  #if (method%in%methods_to_estimate_var) sigma2 <- fit$par[length(fit$par)]
  else if (method=='CSS')     sigma2 <- fit$value/length(y)
  else if (method=='Whittle') sigma2 <- 2/length(y) * var(y) * fit$value  # 1997 Ferrara & Geugen eqn 3.7

  if (k>0) {
    fd <- fit$par[which(nm=='fd')]
    u  <- fit$par[which(nm=='u')]
  }

  res<-list('coefficients'=coef,
            'sigma2'=sigma2,
            'fit_value'=fit$value,
            'convergence'=fit$convergence,
            'conv_message'=c(fit$message,message),
            'method'=method,
            'opt.method'=opt.method,
            'maxeval'=maxeval,
            'order'=order,'k'=k,'y'=py,'mean_y'=mean(y),'m_trunc'=m_trunc)
  if (opt.method=='best') res<-c(res,'opt.method.selected'=best_method)
  if (k==1)
    res<-c(res,
           'ggbr_freq'=acos(u)/2/pi,
           'ggbr_period'=2*pi/acos(u),
           'ggbr_d'=fd)
  class(res)<-'ggbr_model'

  return(res)
}

summary.ggbr_model<-function(mdl) {
  with(mdl,
       cat(sprintf('Summary of a Gegenbauer Time Series model.\n\nFit using %s method.\nOrder=(%d,%d,%d) k=%d %s\n\nOptimisation.\nMethod:  %s\nMaxeval: %d\n',
                   method,order[1],order[2],order[3],k,ifelse(mdl$method=='QML',sprintf('QML Truncation at %d',mdl$m_trunc),''),mdl$opt.method,mdl$maxeval))
  )
  if (mdl$opt.method=='best') cat(sprintf('Best optimisation method selected: %s\n',mdl$opt.method.selected))
  cat(sprintf('Optimal Value: %0.4f\n\n',mdl$fit_value))
  if (mdl$convergence<0) cat(sprintf('Model did not converge.\n\n',mdl$conv_message))
  else {
    if (mdl$convergence>0)
      cat(sprintf('WARNING: Only partial convergence achieved!\n%s reports: %s (%d)\nIt is suggested you increase the maxeval parameter, or try an alternative method.\n\n',
                  ifelse(mdl$opt.method=='best',mdl$opt.method.selected,mdl$opt.method),mdl$conv_message,mdl$convergence))
    print(mdl$coef)

    if (mdl$k>0) {
      cat(sprintf('\nGegenbauer parameters:\nGegenbauer Frequency:  %0.4f\nGegenbauer Period:    %7.4f\nFractional Exponent:   %0.4f\n',
                  mdl$ggbr_freq,mdl$ggbr_period,mdl$ggbr_d))
      if (mdl$ggbr_d>0 & mdl$ggbr_d<0.5) cat(sprintf('Fractional Dimension:  %1.4f\n',1.5-mdl$ggbr_d))
      if (mdl$ggbr_d>0.5) cat('WARNING: Fractional Exponent > 0.5 suggesting the process may not be stationary.\n')
    }
    cat(sprintf('\nsigma^2 estimated as %0.4f\n',mdl$sigma2))
  }
}
print.ggbr_model<-function(mdl) {summary(mdl)}

is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}

predict.ggbr_model<-function(mdl,h=1) {
  assert_that(class(mdl)=='ggbr_model')

  ydm <- mdl$y - mdl$mean_y
  n <- length(ydm)

  coef <- mdl$coefficients[1,]

  if (include.mean) {
    beta0  <- coef[1]
    start  <- 2
  }
  else {
    beta0  <- 0
    start  <- 1
  }
  if (k==1) {
    u      <- coef[start]
    d      <- coef[start+1]
    start  <- start+2
  } else u<-d<-0.0

  if (p>0) phi_vec   <- c(1,-(coef[start:(start+p-1)] ))         else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,-(coef[(p+start):(length(coef)-1)])) else theta_vec <- 1

  print(phi_vec)
  print(theta_vec)

  arma_filter <- signal::Arma(a = theta_vec, b = phi_vec)
  if (k>0) ggbr_filter <- signal::Arma(b = 1, a = ggbr.coef(length(y_dash),d,u))

  #y_dash <- y-beta0
  for (i in 1:h) {
    eps <- signal::filter(arma_filter, ydm)
    if (k>0) eps <- signal::filter(ggbr_filter, eps)
    ydm[n+i] <- eps[length(eps)]
  }
  return(ydm[(n+1):length(ydm)]+mdl$mean_y)
}
