#' Estimate the parameters of a GARMA model.
#'
#' The garma function is the main function for the garma package. Depending on the parameters it will
#' calculate the parameter estimates for the GARMA process, and if available the standard errors (se's)
#' for those parameters.
#'
#' The GARMA model is specified as
#' \deqn{\displaystyle{\phi(B)\prod_{i=1}^{k}(1-2u_{i}B+B^{2})^{d_{i}}(X_{t}-\mu)= \theta(B) \epsilon _{t}}}{\prod(i=1 to k) (1-2u(i)B+B^2)^d(i) \phi(B) (X(t) - \mu) = \theta(B) \epsilon(t)}
#'
#' where
#' \itemize{
#' \item \eqn{\phi(B)}{\phi(B)} represents the short-memory Autoregressive component of order p,
#' \item \eqn{\theta(B)}{\theta(B)} represents the short-memory Moving Average component of order q,
#' \item \eqn{(1-2u_{i}B+B^{2})^{d_{i}}}{(1-2u(i)B+B^2)^d(i)} represents the long-memory Gegenbauer component (there may in general be k of these),
#' \item \eqn{X_{t}}{X(t)} represents the observed process,
#' \item \eqn{\epsilon_{t}}{\epsilon(t)} represents the random component of the model - these are assumed to be uncorrelated but identically distributed variates.
#'       Generally the routines in this package will work best if these have an approximate Gaussian distribution.
#' \item \eqn{B}{B} represents the Backshift operator, defined by \eqn{B X_{t}=X_{t-1}}{B X(t) = X(t-1)}.
#' }
#' when k=0, then this is just a short memory model as fit by the stats "arima" function.
#'
#' @param x (num) This should be a numeric vector representing the process to estimate. A minimum length of 96 is required.
#' @param order (list) This should be a list (similar to the stats::arima order parameter) which will give the order of the process to fit.
#'     The format should be list(p,d,q) where p, d, and q are all positive integers. p represents the degree of the
#'     autoregressive process to fit, q represents the order of the moving average process to fit and d is the (integer)
#'     differencing to apply prior to any fitting.
#' @param k (int) This is the number of (multiplicative) Gegenbauer terms to fit.
#' @param include.mean (bool) A boolean value indicating whether a mean should be fit. Note if you have any differencing, then
#'     it generally does not make sense to fit a mean term. Because of this, the default here is to fit a mean time when d (in the "order" parmaeter)
#'     is zero and otherwise not.
#' @param method (character) This defines the estimation method for the routine. The valid values are 'CSS', 'Whittle', 'QML' and 'WLL'.
#'     The default (Whittle) will generally return very accurate estimates quite quickly, provided the asumption of a Gaussian
#'     distribution is even approximately correct, and is probably the method of choice for most users. For the theory behind this, refer Giraitis et. al. (2001)
#'     'CSS' is a conditional 'sum-of-squares' technique and can be quite slow. Reference: Chung (1996).
#'     'QML' is a Quasi-Maximum-Likelihood technique, and can also be quite slow. Reference Dissanayake (2016). (k>1 is not supported for QML)
#'     'WLL' is a new technique which appears to work well even if the \eqn{\epsilon_{t}}{\epsilon(t)} are highly skewed and/or have heavy tails (skewed and/or lepto-kurtic).
#'     However the asymptotic theory for the WLL method is not complete and so standard errors are not available for most parameters.
#' @param allow_neg_d (bool) A boolean value indicating if a negative value is allowed for the fractional differencing component
#'     of the Gegenbauer term is allowed. This can be set to FALSE to force the routine to find a positive value.
#' @param maxeval (int) the maximum function eveluations to be allowed during each optimisation.
#' @param opt_method (character) This names the optimisation method used to find the parameter estimates. The default is to use 'solnp' from package Rsolnp,
#' which has shown to have good performance in a wide range of circumstances. For some data or some models, however, other methods may work well. Supported algorithms include:
#'     \itemize{
#'     \item cobyla algorithm in package nloptr
#'     \item directL algorithm in package nloptr
#'     \item BBoptim from package BB
#'     \item psoptim from package pso
#'     \item hjkb from dfoptim package
#'     \item nmkb from dfoptim package
#'     \item solnp from Rsolnp package
#'     \item best - this option evaluates all the above options in turn and picks the one which finds the lowest value of the objective. This can be quite time consuming to run,
#'     particularly for the 'CSS' method.
#'     }
#' Note further that if you specify a k>1, then inequality constraints are required, and this will further limit the list of supported routines.
#' @param m_trunc Used for the QML estimation method. This defines the AR-truncation point when evaluating the likelihood function. Refer to Dissanayake et. al. (2016) for details.
#' @return An S3 object of class "garma_model".
#'
#' @references
#' C Chung. A generalized fractionally integrated autoregressive moving-average process.
#' *Journal of Time Series Analysis*, **17**(2):111–140, 1996.
#'
#' G Dissanayake, S Peiris, and T Proietti. State space modelling of Gegenbauer processes with long memory.
#' *Computational Statistics and Data Analysis*, **100**:115–130, 2016.
#'
#' L Giraitis, J Hidalgo, and P Robinson. Gaussian estimation of parametric spectral density with unknown pole.
#' *The Annals of Statistics*, **29**(4):987–1023, 2001.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' print(garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE))
#' # Compare with the built-in arima function
#' print(arima(ap,order=c(9,1,0),include.mean=FALSE))
#' @export
garma<-function(x,
                order=list(0,0,0),
                k=1,
                include.mean=(order[2]==0),
                method='Whittle',
                allow_neg_d=TRUE,
                maxeval=10000,
                opt_method='solnp',
                m_trunc=50) {

  inequality_constraints<-function(par,params) {
    # used when k>1 - to ensure a minimum separation between Gegenbauer frequencies
    k <- params$k
    res <- c()
    if (params$include.mean) start<-2 else start<-1
    if (k>1) for (k1 in 1:(k-1))
      if (k1<k) for (k2 in (k1+1):k) {
        u1<-par[start+k1*2-2]
        u2<-par[start+k2*2-2]
        if(abs(u1)>1) u1=1
        if(abs(u2)>1) u2=1
        sep <- (acos(u1)-acos(u2))/(2*pi)
        res <- c(res, abs(sep)-0.01)
      }
    return(res)
  }

  ## Start of "garma" function logic.
  ## 1. Check parameters
  if (length(x)<96)
    stop('y should have at least 96 observations.')
  if (is.data.frame(x)) {
    if (ncol(x)>1)
      stop('x should be a numeric vector - not an entire data frame. Please select a single column and try again.')
    else x<-x[,1]
  }

  # now save the ts elements, if any - we re-apply them to the output values
  x_start <- stats::start(x)
  x_end   <- stats::end(x)
  x_freq  <- stats::frequency(x)

  x<-as.numeric(x)
  if (!is.numeric(x))
    stop('x should be numeric.')
  if (length(order)!=3)
    stop('order parameter must be a 3 integers only.')
  # if ((k!=0)&(k!=1))
  #   stop('Sorry. Only k=0 or k=1 is supported for now.')
  allowed_methods <- c('CSS','Whittle','WLL','QML')
  if (!method%in%allowed_methods)
    stop('Method must be one of CSS, Whittle, QML or WLL.')
  if (method=='QML'&k>1)
    stop('QML method does not support k>1. It is suggested you try either the CSS or Whittle methods.')

  allowed_optimisations <- c('optim','cobyla','directL','BBoptim','psoptim','hjkb','nmkb','solnp', 'best')
  optimisation_packages <- c('optim'='stats','cobyla'='nloptr','directL'='nloptr','BBoptim'='BB','psoptim'='pso','hjkb'='dfoptim','nmkb'='dfoptim','solnp'='Rsolnp','best'='stats')
  if (!opt_method%in%allowed_optimisations)
    stop('\nSorry - supported packages are:\n', paste0(allowed_optimisations,'\n'))

  if (!.is.installed(optimisation_packages[[opt_method]]))
    stop(sprintf('package %s needs to be installed to use method ',optimisation_packages[[opt_method]],opt_method))

  ##
  ## 2. Next calc parameter  estimates
  p=as.integer(order[1])
  d=as.integer(order[2])
  q=as.integer(order[3])
  #if ((d!=0)&(d!=1))
  #  stop('Sorry. Only d=0 or d=1 is supported for now (for the integer portion of d).\nWe suggest you manually difference the series using diff() if you need more than this.')
  storage.mode(x) <- 'double'

  if (d>0) y<-diff(x,differences=d) else y<-x
  mean_y <- mean(y)
  sd_y   <- stats::sd(y)
  ss<-stats::spectrum((y-mean_y)/sd_y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)

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
    lb_finite <- c(lb_finite,ifelse(mean_y<0, 2*mean_y, -2*mean_y))
    ub_finite <- c(ub_finite,ifelse(mean_y<0,-2*mean_y,  2*mean_y))
    lb        <- c(lb,-Inf)
    ub        <- c(ub,Inf)
  }

  # temp_spec and temp_freq used to determine starting values for Gegenbauer params
  temp_spec <- ss$spec
  temp_freq <- ss$freq
  if (k>0) {
    for (k1 in 1:k) {
      n_pars    <- n_pars+2
      max_spec_idx <- which.max(temp_spec)
      start_u   <- cos(2*pi*temp_freq[max_spec_idx])
      min_resolution_idx <- as.integer(length(temp_spec)/40)
      start_idx <- max(max_spec_idx-min_resolution_idx,1)
      end_idx   <- min(max_spec_idx+min_resolution_idx,length(temp_spec))
      for (i in start_idx:end_idx) temp_spec[i] <- 0 # make sure next time through we find the next highest spec, with min separation

      if (start_u< (-1)) start_u<-0
      pars      <- c(pars,start_u,0.25)
      lb        <- c(lb,0.0,ifelse(allow_neg_d,-1,0))
      ub        <- c(ub,1.0,1.0)
      lb_finite <- c(lb_finite,0.0,ifelse(allow_neg_d,-1,0))
      ub_finite <- c(ub_finite,1.0,1.0)
    }
  }

  n_pars    <- n_pars + p + q
  methods_to_estimate_var <- c('WLL')     # WLL method estimates VAR along with other params; For other methods this falls out of the objective value
  if (p+q>0) {
    a    <- stats::arima(y,order=c(p,0,q),include.mean=FALSE)
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
      if (method%in%methods_to_estimate_var) ub_finite <- c(ub_finite,2*stats::var(y))
    } else {
      lb<- c(lb,rep(-Inf,p+q))
      if (method%in%methods_to_estimate_var) lb <- c(lb,1e-10)
      ub<- c(ub,rep(Inf,p+q))
      if (method%in%methods_to_estimate_var) ub <- c(ub,Inf)
      lb_finite <- c(lb_finite,rep(-10,p+q))
      if (method%in%methods_to_estimate_var) lb_finite <- c(lb_finite,1e-10)
      ub_finite <- c(ub_finite,rep(10,p+q))
      if (method%in%methods_to_estimate_var) ub_finite <- c(ub_finite,2*stats::var(y))
    }
  }

  # create a list of all possible params any 'method' might need. The various objective functions can extract the parameters which are relevant to that method.
  params <- list(y=y, orig_y=x, ss=ss, p=p,q=q,d=d,k=k,include.mean=include.mean,est_mean=ifelse(method%in%mean_methods,TRUE,FALSE),scale=sd_y,m_trunc=m_trunc)
  message <- c()

  # First we make a first pass at optimisation using "optim".
  # If the method chosen is optim then that finishes things; but otherwise the solution found becomes the starting point for the next optimisation.
  fcns <- list('CSS'=.css.ggbr.obj,'Whittle'=.whittle.ggbr.obj,'QML'=.qml.ggbr.obj,'WLL'=.wll.ggbr.obj)
  n_constraints <- k*(k+1)/2-k

  if (k>1) { # separate logic since for k>1 we need inequality constraints and not all non-linear optimisers support this.
    fit<-list(value=Inf,par=pars,pars=pars,convergence= -999,value=Inf,message="Error")
    if (opt_method=='cobyla') {
      tryCatch(fit <- nloptr::cobyla(pars,fcns[[method]],
                                     lower=lb,
                                     upper=ub,
                                     hin=inequality_constraints,
                                     params=params,
                                     control=list(maxeval=maxeval,xtol_rel=1e-10)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece= -999,par=pars)}
      )
    } else {
      tryCatch(fit <- Rsolnp::solnp(pars,
                                    fcns[[method]],
                                    LB=lb,
                                    UB=ub,
                                    ineqfun=inequality_constraints,
                                    ineqLB=rep(0,n_constraints),
                                    ineqUB=rep(1,n_constraints), # this should always be true but supplied as solnp requires it.
                                    control=list(tol=1e-12,trace=0),
                                    params=params),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece= -999,par=pars,pars=pars)}
      )
      fit$par <- fit$pars  # copy across for consistency with other optimisation methods
    }
    hh <- pracma::hessian(fcns[[method]], fit$par, params=params)
  } else { # k==0 or k==1
    fit <- stats::optim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params,
                        hessian=TRUE,method="L-BFGS-B",control=list(maxit=maxeval,factr=1e-25))
    if (fit$convergence==52) {   # Error in line search, then try again using nloptr
      fit2<-nloptr::lbfgs(x0=fit$par, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxeval=maxeval,xtol_rel=1e-8))
      if (fit2$value<fit$value&fit2$convergence>=0) fit<-fit2
    }
    if (fit$convergence>=0) {
      pars <- fit$par
      hh   <- fit$hessian
    }

    if (opt_method=='cobyla') {
      tryCatch(fit <- nloptr::cobyla(pars,fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxeval=maxeval,xtol_rel=1e-10)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece= -999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='solnp') {
      tryCatch(fit <- Rsolnp::solnp(pars,fcns[[method]], LB=lb, UB=ub, control=list(tol=1e-12,trace=0), params=params),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece= -999,par=pars)}
      )
      fit$par <- fit$pars  # remap as solnp has slightly different naming convention.
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='directL') {
      tryCatch(fit <- nloptr::directL(fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxeval=maxeval,xtol_rel=1e-10)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece= -999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='BBoptim') {
      pars[1] <- pars[1]-0.1
      tryCatch(fit <- BB::BBoptim(par=pars, fcns[[method]], lower=lb, upper=ub, control=list(trace=FALSE,maxit=maxeval,ftol=1e-15,gtol=1e-8),
                                  params=params,quiet=TRUE),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='psoptim') {
      tryCatch(fit <-pso::psoptim(par=pars, fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxit=maxeval)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='hjkb') {
      tryCatch(fit <- dfoptim::hjkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='nmkb') {
      tryCatch(fit <- dfoptim::nmkb(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxfeval=maxeval)),
               error=function(cond) {fit<-list(value=Inf,message=cond,convergece=999,par=pars)}
      )
      hh  <- pracma::hessian(fcns[[method]], fit$par, params=params)
    } else if (opt_method=='best') {
      message <- c()
      fit.optim<-fit.cobyla<-fit.directL<-fit.bboptim<-fit.psoptim<-fit.hjkb<-fit.nmkb<-list(value=Inf,convergece=999,par=pars) #default to set environment
      tryCatch(
        fit.optim   <- stats::optim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params,
                                    hessian=TRUE,method="L-BFGS-B",control=list(maxit=maxeval,factr=1e-25)),
        error=function(cond) {fit.optim<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
      )
      tryCatch(fit.cobyla  <- nloptr::cobyla(pars,fcns[[method]], lower=lb, upper=ub, params=params,control=list(maxeval=maxeval,xtol_rel=1e-10)),
               error=function(cond) {fit.cobyla<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
      )
      tryCatch(fit.directL <- nloptr::directL(fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params,control=list(maxeval=maxeval,xtol_rel=1e-10)),
               error=function(cond) {fit.directL<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
      )
      tryCatch(fit.bboptim <- BB::BBoptim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, control=list(trace=FALSE,maxit=maxeval,ftol=1e-15,gtol=1e-8),
                                          params=params,quiet=TRUE),
               error=function(cond) {fit.bboptim<-list(value=Inf,message=c(message,cond),convergece=999,par=pars)}
      )
      tryCatch(fit.psoptim <- pso::psoptim(par=pars, fn=fcns[[method]], lower=lb_finite, upper=ub_finite, params=params, control=list(maxit=maxeval)),
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
  }

  # log lik
  loglik <- numeric(0)
  if (method=='CSS')
    loglik <- -0.5 *(length(y)*fit$value[length(fit$value)] + length(y) + length(y) * log(2 * pi))
  if (method=='QML')
    loglik <- -fit$value[length(fit$value)]
  if (method=='Whittle')
    loglik <- .whittle.ggbr.likelihood(fit$par,params)

  # sigma2
  if (method=='WLL') {
    # adjust sigma2 for theoretical bias...
    sigma2 <- fit$par[length(fit$par)] <- fit$par[length(fit$par)]/(2*pi) * exp(-digamma(1))
  }
  else if (method=='QML')     sigma2 <- sqrt(.qml.ggbr.se2(fit$par, params=params))
  else if (method=='CSS')     sigma2 <- exp(2*fit$value[length(fit$value)])/length(y)
  else if (method=='Whittle') sigma2 <- 2/length(y) * var(y) * fit$value[length(fit$value)]  # 1997 Ferrara & Geugen eqn 3.7

  se <- numeric(length(fit$par))
  if (fit$convergence>=0&method!='WLL'&!is.null(hh)) {
    # Next, find the c <- alc se's for coefficients
    start<-1
    se<-c()   # default to set this up in the right environment
    if (include.mean) start<-2
    if (method=='Whittle') se <- sqrt(diag(solve(hh)))
    if (method=='CSS')     se <- sqrt(diag(solve(hh*length(y))))
    if (method=='QML')     {
      se <- sqrt(diag(solve(hh))*length(y))
      if (k==1) {
        if (include.mean) se <- se[1:3]
        else se <- se[1:2]
      }
    }
    if (length(se)<length(fit$par)) se<-c(se,NA)
  }
  if (method=='WLL') {
    se<-rep(NA,length(par))
    if (k==1) se[2] <- .wll_d_se(fit$par[1],ss)  # result only holds for k=1
  }
  nm<-list()
  if (include.mean) nm <- c(nm,'intercept')
  if (k>0) nm<-c(nm,rep(c('u','fd'),k))
  if (p>0) nm<-c(nm,paste0('ar',1:p))
  if (q>0) nm<-c(nm,paste0('ma',1:q))

  n_coef    <- ifelse(method%in%methods_to_estimate_var,length(fit$par)-1,length(fit$par))  # ignore var on end if it is there
  temp_coef <- fit$par[1:n_coef]
  temp_se   <- se[1:n_coef]
  if (include.mean&!method%in%mean_methods) {# then add an Intercept anyway; use Kunsch 1987 result for se
    temp_coef <- c(mean(y),temp_coef)
    # calc se of mean using Kunsch (1987) thm 1; 1-H = 1-(d+.5) = .5-d
    if (k==0) mean_se <- sqrt(sigma2/length(y))
    else mean_se <- sqrt(sigma2/(length(y)^(0.5-fit$par[2])))
    temp_se   <- c(mean_se, temp_se)
    n_coef <- n_coef+1
  }
  coef <- t(matrix(round(c(temp_coef,temp_se),4),nrow=n_coef))
  colnames(coef) <- nm
  rownames(coef) <- c('','s.e.')

  # get fitted values and residuals
  fitted <- .fitted_values(fit$par,params)

  res<-list('call' = match.call(),
            'coef'=coef,
            'sigma2'=sigma2,
            'obj_value'=fit$value[length(fit$value)],
            'loglik'=loglik,
            'convergence'=fit$convergence,
            'conv_message'=c(fit$message,message),
            'method'=method,
            'opt_method'=opt_method,
            'maxeval'=maxeval,
            'order'=order,
            'k'=k,
            'y'=x,
            'y_start'=x_start,
            'y_end'=x_end,
            'y_freq'=x_freq,
            'include.mean'=include.mean,
            'fitted'=stats::ts(fitted$fitted,start=x_start,end=x_end,frequency=x_freq),
            'residuals'=stats::ts(fitted$residuals,start=x_start,end=x_end,frequency=x_freq),
            # 'par'=fit$par,
            # 'params'=params,
            'm_trunc'=m_trunc)
  if (opt_method=='best') res<-c(res,'opt_method.selected'=best_method)
  if (k>0) {
    fd <- u <- c()
    start <- 1
    if (colnames(coef)[1]=='intercept') start <- 2
    for (k1 in 1:k) {
      u  <- c(u,  coef[1,start+k1*2-2])
      fd <- c(fd, coef[1,start+k1*2-1])
    }
    res<-c(res,
           'ggbr_freq'=list(acos(u)/2/pi),
           'ggbr_period'=list(2*pi/acos(u)),
           'ggbr_d'=list(fd))
  }

  class(res)<-'garma_model'

  return(res)
}




.print_garma_model<-function(mdl,verbose=TRUE) {
  cat("\nCall:", deparse(mdl$call, width.cutoff = 75L), "", sep = "\n")
  if (verbose) {
    with(mdl,
         cat(sprintf('Summary of a Gegenbauer Time Series model.\n\nFit using %s method.\nOrder=(%d,%d,%d) k=%d %s\n\nOptimisation.\nMethod:  %s\nMaxeval: %d\n',
                     method,order[1],order[2],order[3],k,ifelse(mdl$method=='QML',sprintf('QML Truncation at %d',mdl$m_trunc),''),mdl$opt_method,mdl$maxeval))
    )
    if (mdl$opt_method=='best') cat(sprintf('Best optimisation method selected: %s\n',mdl$opt_method.selected))
    cat(sprintf('Optimal Value found: %0.8f\n\n',mdl$obj_value))
  }
  if (mdl$convergence<0) cat(sprintf('Model did not converge.\n\n',mdl$conv_message))
  else {
    if (mdl$convergence>0)
      cat(sprintf('WARNING: Only partial convergence achieved!\n%s reports: %s (%d)\n\n',
                  ifelse(mdl$opt_method=='best',mdl$opt_method.selected,mdl$opt_method),mdl$conv_message,mdl$convergence))
    cat('Coefficients:\n')
    print.default(mdl$coef, print.gap = 2)

    if (mdl$k>0) {
      cat('\nGegenbauer parameters:\n                      ')
      for (k1 in 1:(mdl$k)) cat(sprintf('  Factor%d',k1))
      cat('\nGegenbauer Frequency: ')
      for (k1 in 1:(mdl$k)) cat(sprintf('  %7.4f',mdl$ggbr_freq[k1]))
      cat('\nGegenbauer Period:    ')
      for (k1 in 1:(mdl$k)) cat(sprintf('  %7.4f',mdl$ggbr_period[k1]))
      cat('\nFractional Exponent:  ')
      for (k1 in 1:(mdl$k)) cat(sprintf('  %7.4f',mdl$ggbr_d[k1]))
      cat('\nFractional Dimension: ')
      for (k1 in 1:(mdl$k)) cat(sprintf('  %7.4f',1.5-mdl$ggbr_d[k1]))
    }
    cat(sprintf('\n\nsigma^2 estimated as %0.4f',mdl$sigma2))
    if (mdl$method=='CSS') cat(sprintf(':  part log likelihood = %f',mdl$loglik))
    if (mdl$method=='QML') cat(sprintf(':  log likelihood = %f',mdl$loglik))
    if (mdl$method=='Whittle') cat(sprintf(':  log likelihood = %f',mdl$loglik))
    cat('\n')
  }
}


#' The summary function provides a summary of a "garma_model" object.
#' @param object (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' summary(mdl)
#' @export
summary.garma_model<-function(object,...) {
  .print_garma_model(object,verbose=FALSE)
}

#' The print function prints a summary of a "garma_model" object.
#' @param x (garma_model) The garma_model from which to print the values.
#' @param ... Other arguments. Ignored.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' print(mdl)
#' @export
print.garma_model<-function(x,...) {
  .print_garma_model(x,verbose=TRUE)
}

.is.installed <- function(mypkg){
  is.element(mypkg, utils::installed.packages()[,1])
}

#' The predict function predicts future values of a "garma_model" object.
#' @param object (garma_model) The garma_model from which to predict the values.
#' @param n.ahead (int) The number of time periods to predict ahead. Default: 1
#' @param ... Other parameters. Ignored.
#' @return A "ts" object containing the requested forecasts.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' predict(mdl, n.ahead=12)
#' @export
predict.garma_model<-function(object,n.ahead=1,...) {
  coef <- unname(object$coef[1,])
  p<-object$order[1]
  q<-object$order[3]

  if (object$include.mean) {
    beta0  <- coef[1]
    start  <- 2
  }
  else {
    beta0  <- 0
    start  <- 1
  }
  if (object$k==1) {
    u      <- coef[start]
    d      <- coef[start+1]
    start  <- start+2
  } else u<-d<-0.0

  if (p>0) phi_vec   <- c(-(coef[start:(start+p-1)] ))           else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,-(coef[(p+start):(length(coef)-1)])) else theta_vec <- 1

  if (object$order[2]==0)
    ydm <- object$y - beta0
  else if (object$order[2]==1)
    ydm <- diff(object$y) - beta0

  n <- length(ydm)

  # set up filters
  arma_filter <- signal::Arma(a = theta_vec, b = phi_vec)
  if (object$k>0) ggbr_filter <- signal::Arma(b=1, a=.ggbr.coef(n,d,u))

  # generate forecasts
  for (i in 1:n.ahead) {
    eps <- ydm
    if (object$k>0) eps <- signal::filter(ggbr_filter, eps)
    eps <- signal::filter(arma_filter, eps)
    ydm[n+i] <- (-eps[length(eps)])
  }

  # if (integer) differenced then...
  if (object$order[2]>0) {
    ydm2 <- stats::diffinv(ydm,differences=object$order[2])
  }
  else ydm2 <-ydm

  # Now we have the forecasts, we set these up as a "ts" object - as does "predict.arima"
  y_end = object$y_end
  if(length(y_end)>1) {
    if (object$y_freq >= y_end[2]) {
      y_end[1] <- y_end[1]+1
      y_end[2] <- y_end[2]-object$y_freq+1
    }
    else y_end[2] <- y_end[2] + 1
  } else y_end <- y_end +1

  res <- stats::ts(tail(ydm2,n.ahead)+beta0,start=y_end,frequency=object$y_freq)
  return(list(pred=res))
}

#' The forecast function predicts future values of a "garma_model" object, and is exactly the same as the "predict" function with slightly different parameter values.
#' @param mdl (garma_model) The garma_model from which to forecast the values.
#' @param h (int) The number of time periods to predict ahead. Default: 1
#' @return - a "ts" object containing the requested forecasts.
#' @examples
#' library(forecast)
#'
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' forecast(mdl, h=12)
#' @export
forecast.garma_model<-function(mdl,h=1) {return(predict.garma_model(mdl,n.ahead=h))}

.plot_garma_model<-function(mdl,h=24,...) {
  # plot forecasts from model
  actuals <- zoo(stats::ts(mdl$y,start=mdl$y_start,end=mdl$y_end,frequency=mdl$y_freq))
  fitted <- zoo(mdl$fitted_values)
  fc <- zoo(predict.garma_model(mdl,n.ahead=h)$pred)
  graphics::plot(actuals,col='black',type='l',xlim=c(stats::start(actuals)[1],stats::end(fc)[1]),...)
  graphics::lines(fitted,col='blue')
  graphics::lines(fc,col='blue')
  graphics::abline(v=mdl$y_end,col='red',lty=2)
}

#' The plot function generates a plot of actuals and predicted values for a "garma_model" object.
#' @param x (garma_model) The garma_model from which to plot the values.
#' @param ... other arguments to be passed to the "plot" function, including h (int) - the number of periods ahgead to forecast.
#' @return An R "plot" object.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' plot(mdl)
#' @export
plot.garma_model<-function(x,...) {
  .plot_garma_model(x,...)
}

#' The ggplot function generates a ggplot of actuals and predicted values for a "garma_model" object.
#' @param mdl (garma_model) The garma_model from which to ggplot the values.
#' @param h (int) The number of time periods to predict ahead. Default: 24
#' @param ... other parameters passed to ggplot.
#' @return A ggplot2 "ggplot" object. Note that the standard ggplot2 "+" notation can be used to enhance the default output.
#' @examples
#' library(ggplot2)
#'
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' ggplot(mdl)
#' @export
ggplot.garma_model<-function(mdl,h=24,...) {
  # plot forecasts from model
  fc <- predict.garma_model(mdl,n.ahead=h)

  if (mdl$y_freq>1) { # then we have actual dates not just an index
    idx <- seq(lubridate::make_date(mdl$y_start[1],mdl$y_start[2],15),by=mdl$y_freq,length.out=(length(mdl$y)+h))
    lubridate::day(idx) <- lubridate::days_in_month(idx)
    cutoff <- lubridate::make_date(mdl$y_end[1],mdl$y_end[2],15)
  } else {
    idx <- (mdl$y_start[1]):(mdl$y_end[1]+h)
    cutoff <- mdl$y_end[1]+1
  }

  df1 <- data.frame(dt=idx,grp='Actuals',value=c(mdl$y,rep(NA,h)))
  df2 <- data.frame(dt=idx,grp='Forecasts',value=c(as.numeric(mdl$fitted),fc$pred))
  df <- rbind(df1,df2)
  ggplot2::ggplot(df[!is.na(df$value),],ggplot2::aes(x=dt,y=value,color=grp),...) +
    ggplot2::geom_line() + ggplot2::ylab('') + ggplot2::xlab('') +
    ggplot2::geom_vline(xintercept=cutoff,color='red',linetype=2) +
    ggplot2::theme_bw() + ggplot2::theme(legend.title=ggplot2::element_blank()) +
    #ggplot2::scale_color_brewer(palette="Set1")
    scale_colour_manual(values=c('gray20','dodgerblue4',rep('gray',10)))
}


.fitted_values<-function(par,params) { # Generate fitted values and residuals for GARMA process
  y <- params$y
  p <- params$p
  q <- params$q
  id <- params$d
  k <- params$k
  include.mean <- params$include.mean
  est_mean <- params$est_mean

  if (include.mean&est_mean) {
    beta0  <- par[1]
    start  <- 2
  }
  else {
    beta0  <- 0
    start  <- 1
  }
  # if (k==1) {
  #   u      <- par[start]
  #   fd     <- par[start+1]
  #   start  <- start+2
  # } else u<-fd<-0.0
  u <- c()
  fd <- c()
  if (k>0) for (k1 in 1:k) {
    u     <- c(u,par[start])
    fd    <- c(fd,par[start+1])
    start <- start+2
  }


  y_dash <- y-beta0
  if (p>0) phi_vec   <- c(1,-(par[start:(start+p-1)] ))        else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,-(par[(p+start):(length(par)-1)])) else theta_vec <- 1

  arma_filter   <- signal::Arma(a = theta_vec, b = phi_vec)
  eps           <- signal::filter(arma_filter, y_dash)
  if (k>0) for (k1 in 1:k) {
    ggbr_filter <- signal::Arma(b = 1, a = .ggbr.coef(length(y_dash),fd[k1],u[k1]))
    eps         <- signal::filter(ggbr_filter, eps)
  }

  # if (integer) differenced then...
  eps_y <- eps
  fitted <- (params$y-eps)
  if (id>0) {
    eps_y <- stats::diffinv(fitted,differences = id)
    fitted <- (params$orig_y-eps_y)
  }

  return(list(fitted=fitted,residuals=eps_y))
}

#' @export
fitted.garma_model<-function(object,...) {
  return(object$fitted)
}
#' @export
residuals.garma_model<-function(object,...) {
  return(object$residuals)
}
