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
#'     of the Gegenbauer term is allowed. This can be set to FALSE (the default) to force the routine to find a positive value.
#' @param maxeval (int) the maximum function eveluations to be allowed during each optimisation.
#' @param opt_method (character) This names the optimisation method used to find the parameter estimates.
#' This may be a list of methods, in which case the methods are applied in turn,
#' each using the results of the previous one as the starting point for the next. The default is to use c('directL', 'solnp') when k<2 and 'solnp' when k>=2. The
#' directL algorithm is used to perform a global search for the minima, and solnp to refine the values.
#' For some data or some models, however, other methods may work well.
#' Supported algorithms include:
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
#' @param min_freq (num) When searching for Gegenbauer peaks, this is the minimum frequency used. Default 0. Note that when there is an
#' AR(1) component, the peaks corresponding to the AR(1) can be higher than the Gegenbauer peaks. Setting this parameter to 0.05 or above can help.
#' @param max_freq (num) default 0.5. When searching for Gegenbauer peaks, this is the maximum frequency used.
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
                allow_neg_d=FALSE,
                maxeval=10000,
                opt_method=NULL,
                m_trunc=50,
                min_freq=0,
                max_freq=0.5) {

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
  nloptr_ineq_constr<-function(par) {
    # used when k>1 - to ensure a minimum separation between Gegenbauer frequencies
    # this is essentially the same as the above, but passes the "params" from
    # a global var instead of as a parameter - some nloptr routines appear to have a
    # bug where they don't pass the  params over to the ineq. constr function...
    k <- nloptr_params$k
    res <- c()
    if (nloptr_params$include.mean) start<-2 else start<-1
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
    # This is a bit arbitary but I would be concerned about a process of length even 96.
    # But no real evidence for this apart from simulations showing large std errs.
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
    stop('x should be numeric.\n')
  if (length(order)!=3)
    stop('order parameter must be a 3 integers only.\n')
  allowed_methods <- c('CSS','Whittle','WLL','QML')
  if (!method%in%allowed_methods)
    stop('Method must be one of CSS, Whittle, QML or WLL.\n')
  if (method=='QML'&k>1)
    stop('QML method does not support k>1. It is suggested you try either the CSS or Whittle methods.\n')

  if (is.null(opt_method)) {
    if (k>=2) opt_method <- 'solnp'
    else if (order[3]>0) opt_method<-'cobyla'
    else opt_method <- c('directL','solnp')
  }

  for (om in opt_method) {
    if (!om%in%.supported_optim())
      stop(sprintf('\nError: function %s not available.\n\nSupported functions are:\n%s\n', om, .supported_optim()))

    optimisation_packages <- .optim_packages()
    if (!isNamespaceLoaded(optimisation_packages[[om]]))
      stop(sprintf('Package %s is required to use method %s\n',optimisation_packages[[om]],om))

    if (k>1) {
      if (!om%in%.supported_contr_optim())
        stop(sprintf('For k>1 we need to use contrained optimisation, but algorithm %s does not support that.\nPlease try one of %s\n',
                     om,paste(.supported_contr_optim(),collapse=', ')))
    }
  }
  # check min_freq and max_frerq
  if (!is.numeric(min_freq)|!is.numeric(max_freq)|min_freq<0|min_freq>=0.5|max_freq<=0|max_freq>0.5|min_freq>=max_freq)
    stop('min_freq and max_freq must be numeric and between 0 and 0.5 and min_freq<max_freq.\n')
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
  #temp_spec <- ss$spec
  #temp_freq <- ss$freq
  if (k>0) {
    gf <- ggbr_semipara(y,k=k,min_freq=min_freq,max_freq=max_freq)
    for (k1 in 1:k) {
      n_pars    <- n_pars+2
      gf1 <- gf$ggbr_factors[[k1]]
      start_u <- gf1$u
      start_d <- gf1$fd

      max_u <- cos(2*pi*min_freq)
      min_u <- cos(2*pi*max_freq)
      if (start_u< min_u|start_u > max_u) start_u <- (min_u+max_u)/2
      if (start_d>=0.5)  start_d <- 0.49
      if (allow_neg_d) {
        if (start_d<= -0.5) start_d<- -0.49
      } else {
        if (start_d<=0.0) start_d<-0.01
      }
      pars      <- c(pars,start_u,start_d)
      lb        <- c(lb,min_u,ifelse(allow_neg_d,-1,0))
      ub        <- c(ub,max_u,0.5)
      lb_finite <- c(lb_finite,min_u,ifelse(allow_neg_d,-1,0))
      ub_finite <- c(ub_finite,max_u,0.5)
    }
  }

  n_pars    <- n_pars + p + q
  methods_to_estimate_var <- c('WLL')     # WLL method estimates VAR along with other params; For other methods this falls out of the objective value
  if (p+q>0) {
    # if any ARMA params to be estimated, we use the semipara estimates to get ggbf factors then get the underlying ARMA process,
    # and then ask "arima" for estimates. Semi para estimates should make good starting points for optimisation.
    if (k>0) arma_y <- extract_arma(y,gf$ggbr_factors)
    else arma_y <- y
    a <- stats::arima(arma_y,order=c(p,0,q),include.mean=FALSE)
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
  nloptr_params <- params
  message <- c()

  # First we make a first pass at optimisation using "optim".
  # If the method chosen is optim then that finishes things; but otherwise the solution found becomes the starting point for the next optimisation.
  fcns <- list('CSS'=.css.ggbr.obj,'Whittle'=.whittle.ggbr.obj,'QML'=.qml.ggbr.obj,'WLL'=.wll.ggbr.obj)
  n_constraints <- k*(k+1)/2-k

  if (k>1) { # separate logic since for k>1 we need inequality constraints and not all non-linear optimisers support this.
    fit <- .generic_optim_list(opt_method_list=opt_method, initial_pars=pars, fcn=fcns[[method]],lb=lb,ub=ub, ineq_fcn=inequality_constraints,
                               ineq_lb=rep(0,n_constraints), ineq_ub=rep(1,n_constraints), params=params, max_eval=maxeval)
  } else if (opt_method[[1]]=='best') {
    fit <- .best_optim(initial_pars=pars, fcn=fcns[[method]], lb=lb, ub=ub, lb_finite=lb_finite, ub_finite=ub_finite, params=params, max_eval=maxeval)
  } else { # k==0 or k==1
      fit <- .generic_optim_list(opt_method_list=opt_method, initial_pars=pars, fcn=fcns[[method]],
                                 lb=lb, ub=ub, lb_finite=lb_finite, ub_finite=ub_finite, params=params, max_eval=maxeval)
  }
  if (fit$convergence== -999) stop('Failed to converge.')

  hh <- fit$hessian
  for (col in 1:ncol(hh)) hh[any(is.na(hh[,col]))|any(is.infinite(hh[,col])),col] <- 0

  # sigma2
  if (method=='WLL') {
    # adjust sigma2 for theoretical bias...
    sigma2 <- fit$par[length(fit$par)] <- fit$par[length(fit$par)]/(2*pi) * exp(-digamma(1))
  }
  else if (method=='QML')     sigma2 <- sqrt(.qml.ggbr.se2(fit$par, params=params))
  else if (method=='CSS')     sigma2 <- exp(2*fit$value[length(fit$value)])/length(y)
  else if (method=='Whittle') sigma2 <- 2/length(y) * var(y) * fit$value[length(fit$value)]  # 1997 Ferrara & Geugen eqn 3.7

  # log lik
  loglik <- numeric(0)
  if (method=='CSS')
    loglik <- -0.5* ((fit$value/sigma2) + length(y)*(log(2*pi) + log(sigma2)))
  if (method=='QML')
    loglik <- -fit$value[length(fit$value)]
  if (method=='Whittle')
    loglik <- .whittle.ggbr.likelihood(fit$par,params)

  se <- numeric(length(fit$par))
  if (fit$convergence>=0&method!='WLL'&!is.null(hh)) {
    # Next, find the se's for coefficients
    start<-1
    se<-c()   # default to set this up in the right environment
    if (include.mean) start<-2
    if (method=='Whittle') se <- sqrt(diag(pracma::pinv(hh)))
    if (method=='CSS')     se <- sqrt(diag(pracma::pinv(hh*length(y))))
    if (method=='QML')     {
      se <- sqrt(diag(pracma::pinv(hh))*length(y))
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
  if (k>0) nm<-c(nm,unlist(lapply(1:k,function(x) paste0(c('u','fd'),x))))
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

  # build a ggbr_factors object
  gf <- list()
  if (k>0) {
    if (include.mean) start_idx <- 2
    else start_idx <- 1
    for (k1 in 1:k) {
      u <- coef[1,start_idx]
      gf1 <- list(u=u, f=acos(u)/(2*pi), fd=coef[1,start_idx+1], m=NA, f_idx=NA)
      gf <- c(gf,list(gf1))
      start_idx <- start_idx+2
    }
  }
  class(gf) <- 'ggbr_factors'

  res<-list('call' = match.call(),
            'coef'=coef,
            'sigma2'=sigma2,
            'obj_value'=fit$value,
            'loglik'=loglik,
            'aic'=-2*loglik + 2*(n_coef+1),
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
            'm_trunc'=m_trunc)
  if (opt_method[1]=='best') res<-c(res,'opt_method_selected'=fit$best_method)
  if (k>0) res<-c(res, 'ggbr_factors' = list(gf))

  class(res)<-'garma_model'

  return(res)
}




.print_garma_model<-function(mdl,verbose=TRUE) {
  cat("\nCall:", deparse(mdl$call, width.cutoff = 75L), "", sep = "\n")
  if (verbose) {
    with(mdl,
         cat(sprintf('Summary of a Gegenbauer Time Series model.\n\nFit using %s method.\nOrder=(%d,%d,%d) k=%d %s\n\nOptimisation.\nMethod:  %s\nMaxeval: %d\n',
                     method,order[1],order[2],order[3],k,ifelse(mdl$method=='QML',sprintf('QML Truncation at %d',mdl$m_trunc),''),
                     paste(mdl$opt_method,collapse=', '),mdl$maxeval))
    )
    if (mdl$opt_method[[1]]=='best') cat(sprintf('Best optimisation method selected: %s\n',mdl$opt_method_selected))
    cat(sprintf('Convergence Code: %d\nOptimal Value found: %0.8f\n\n',mdl$convergence,mdl$obj_value))
  }
  if (mdl$opt_method[[1]]=='solnp'&mdl$convergence!=0) cat('ERROR: Convergence not acheived. Please try another method.\n')
  if (mdl$convergence<0) cat(sprintf('Model did not converge.\n\n',mdl$conv_message))
  else {
    if (mdl$convergence>0)
      cat(sprintf('WARNING: Only partial convergence achieved!\n%s reports: %s (%d)\n\n',
                  ifelse(mdl$opt_method=='best',mdl$opt_method_selected,mdl$opt_method),mdl$conv_message,mdl$convergence))
    cat('Coefficients:\n')
    print.default(mdl$coef, print.gap = 2)
    cat('\n')

    if (mdl$k>0) print(mdl$ggbr_factors)

    cat(sprintf('\n\nsigma^2 estimated as %0.4f',mdl$sigma2))
    if (mdl$method=='CSS') cat(sprintf(':  part log likelihood = %f',mdl$loglik))
    if (mdl$method=='QML') cat(sprintf(':  log likelihood = %f',mdl$loglik))
    if (mdl$method=='Whittle') cat(sprintf(':  log likelihood = %f, aic = %f',mdl$loglik, mdl$aic))
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
  if (n.ahead<=0) stop('n.ahead must be g.t. 0.')

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
  # jump over the ggbr params
  start <- start + ((object$k)*2)

  # if (p>0) phi_vec   <- c((coef[start:(start+p-1)] ))       else phi_vec   <- c()
  # if (q>0) theta_vec <- c((coef[(p+start):(length(coef))])) else theta_vec <- c()

  if (object$order[2]==0)
    ydm <- object$y - beta0
  else if (object$order[2]>0)
    ydm <- diff(object$y,object$order[2]) - beta0

  # Next section uses eqn 5.3.9 from Brockwell & Davis (1991)
  # to calculate forecasts for the short memory component.
  # n <- length(ydm)
  # resid <- object$residuals
  # ydm <- c(ydm,rep(0,n.ahead))
  # if (p>0) {
  #   for (h in 1:n.ahead)
  #     ydm[n+h] <- ydm[n+h] + ydm[(n+h-1):(n+h-p)] %*% phi_vec
  #   print(tail(ydm,n.ahead))
  # }
  # if(q>0) {
  #   for (h in 1:min(q,n.ahead)) {
  #     ydm[n+h] <- ydm[n+h] + resid[n:(n-q+h)] %*% theta_vec[h:q]
  #     print(tail(ydm,n.ahead))
  #   }
  # }

  if (p>0) phi_vec   <- c(1,-(coef[start:(start+p-1)] ))      else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,(coef[(p+start):(length(coef))])) else theta_vec <- 1

  n <- length(ydm)
  ydm <- c(object$residuals, rep(0,n.ahead))

  # set up filters
  arma_filter <- signal::Arma(b=theta_vec, a=phi_vec)
  ydm <- signal::filter(arma_filter, ydm)
  if (object$k>0) {
    # for each ggbr factor, we set up a filter and add it to the list
    for (k1 in 1:object$k) {
      gf <- object$ggbr_factors[[k1]]
      gc <- .ggbr.coef(n+n.ahead,gf$fd,gf$u)
      ggbr_filter <- signal::Arma(a=1,b=gc)
      ydm <- signal::filter(ggbr_filter, ydm)
    }
  }

  # if (integer) differenced then...
  if (object$order[2]>0) {
    ydm2 <- stats::diffinv(ydm,differences=object$order[2]) + beta0
  }
  else ydm2 <- ydm + beta0

  # Now we have the forecasts, we set these up as a "ts" object - as does "predict.arima"
  y_end = object$y_end
  if(length(y_end)>1) {
    if (object$y_freq >= y_end[2]) {
      y_end[1] <- y_end[1]+1
      y_end[2] <- y_end[2]-object$y_freq+1
    }
    else y_end[2] <- y_end[2] + 1
  } else y_end <- y_end +1

  res <- stats::ts(tail(ydm2,n.ahead),start=y_end,frequency=object$y_freq)
  return(list(mean=res))
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
  if (p>0) phi_vec   <- c(1,-(par[start:(start+p-1)] ))     else phi_vec   <- 1
  if (q>0) theta_vec <- c(1,(par[(p+start):(length(par))])) else theta_vec <- 1

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
