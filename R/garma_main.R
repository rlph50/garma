#' Estimate the parameters of a GARMA model.
#'
#' The garma function is the main function for the garma package. Depending on the parameters it will
#' calculate the parameter estimates for the GARMA process, and if available the standard errors (se's)
#' for those parameters.
#'
#' The GARMA model is specified as
#' \deqn{\displaystyle{\phi(B)\prod_{i=1}^{k}(1-2u_{i}B+B^{2})^{d_{i}}(1-B)^{id} (X_{t}-\mu)= \theta(B) \epsilon _{t}}}{\prod(i=1 to k) (1-2u(i)B+B^2)^d(i) \phi(B)(1-B)^{id} (X(t) - \mu) = \theta(B) \epsilon(t)}
#'
#' where
#' \itemize{
#' \item \eqn{\phi(B)}{\phi(B)} represents the short-memory Autoregressive component of order p,
#' \item \eqn{\theta(B)}{\theta(B)} represents the short-memory Moving Average component of order q,
#' \item \eqn{(1-2u_{i}B+B^{2})^{d_{i}}}{(1-2u(i)B+B^2)^d(i)} represents the long-memory Gegenbauer component (there may in general be k of these),
#' \item \eqn{id} represents the degree of integer differencing.
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
#'     differencing to apply prior to any fitting. WARNING: Currently only d==0 or d==1 are allowed.
#' @param k (int) This is the number of (multiplicative) Gegenbauer terms to fit. Note the 'QML' method only allows k=1.
#' @param include.mean (bool) A boolean value indicating whether a mean should be fit.
#'     Note that no mean term is fit if the series is integer differenced.
#' @param include.drift (bool) A boolean value indicating whether a 'drift' term should be fit to the predictions.
#'     The default is to fit a drift term to the predictions if the process is integer-differenced.
#' @param method (character) This defines the estimation method for the routine. The valid values are 'CSS', 'Whittle', 'QML' and 'WLL'.
#'     The default (Whittle) will generally return very accurate estimates quite quickly, provided the assumption of a Gaussian
#'     distribution is even approximately correct, and is probably the method of choice for most users. For the theory behind this, refer Giraitis et. al. (2001)
#'     'CSS' is a conditional 'sum-of-squares' technique and can be quite slow. Reference: Chung (1996).
##     'QML' is a Quasi-Maximum-Likelihood technique, and can also be quite slow. Reference Dissanayake (2016). (k>1 is not supported for QML)
#'     'WLL' is a new technique which appears to work well even if the \eqn{\epsilon_{t}}{\epsilon(t)} are highly skewed and/or have heavy tails (skewed and/or lepto-kurtic).
#'     However the asymptotic theory for the WLL method is not complete and so standard errors are not available for most parameters.
## @param allow_neg_d (bool) A boolean value indicating if a negative value is allowed for the fractional differencing component
##     of the Gegenbauer term is allowed. This can be set to FALSE (the default) to force the routine to find a positive value.
#' @param d_lim (list) the limits for the d parameter. The default is c(0,0.5), which restricts the model to be stationary.
#'        However sometimes it is desirable to understand what the unrestricted value might be.
#' @param freq_lim (list) the limits for the frequencies to be searched for Gegenbauer factors.
#' When searching for Gegenbauer peaks, when there is an AR(1) component, the peaks corresponding to the AR(1)
#' can be higher than the Gegenbauer peaks. Setting these limits to c(0.05,0.45) or similar can help the routine find the correct optima.
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
#'     \item psoptim from package pso - Particle Swarm algorithm
#'     \item hjkb from dfoptim package
#'     \item nmkb from dfoptim package
#'     \item solnp from Rsolnp package
#'     \item gosolnp from Rsolnp package
#'     \item ga from package GA - a Genetic Algorithm
#'     \item best - this option evaluates all the above options in turn and picks the one which finds the lowest value of the objective. This can be quite time consuming to run,
#'     particularly for the 'CSS' method.
#'     }
#' Note further that if you specify a k>1, then inequality constraints are required, and this will further limit the list of supported routines.
#' @param m_trunc Used for the QML estimation method. This defines the AR-truncation point when evaluating the likelihood function. Refer to Dissanayake et. al. (2016) for details.
## @param min_freq (num) When searching for Gegenbauer peaks, this is the minimum frequency used. Default 0. Note that when there is an
## AR(1) component, the peaks corresponding to the AR(1) can be higher than the Gegenbauer peaks. Setting this parameter to 0.05 or above can help.
## @param max_freq (num) default 0.5. When searching for Gegenbauer peaks, this is the maximum frequency used.
#' @param fitted (bool) indicates whether fitted values should be generated. For longer processes this can take a while and so this option is provided to disable
#' that process. In any event if a call is made to the 'fitted' or 'resid' methods, they will then be generated on demand, but if the practitioner is just exploring
#' different models then this option can be used to speed up the process. Default: TRUE.
#' @param control (list) list of optimisation routine specific values.
#' @return An S3 object of class "garma_model".
#'
#' @references
#' R Hunt, S Peiris, N Weber. Estimation methods for stationary Gegenbauer processes.
#' Statistical Papers 63:1707-1741 2022
#'
#' L Giraitis, J Hidalgo, and P Robinson. Gaussian estimation of parametric spectral density with unknown pole.
#' The Annals of Statistics, 29(4):987–1023, 2001.
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
                include.drift=FALSE,
                method='Whittle',
                d_lim=c(0,0.5),
                freq_lim=c(0,0.5),
                opt_method=NULL,
                m_trunc=50,
                fitted=TRUE,
                control=NULL) {

  ## Start of "garma" function logic.
  ## 1. Check parameters
  if (length(x)<96)
    # This is a bit arbitary but I would be concerned about a process of length even 96.
    # But no real evidence for this restriction apart from simulations showing large std errs.
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
    stop('x should be numeric or a ts object.\n')
  if (any(is.na(x)))
    stop('x should not have any missing values.\n')
  if (length(order)!=3)
    stop('order parameter must be a 3 integers only.\n')
  if (any(order<0))
    stop('order parameter must consist of positive integers.\n')
  if (k<0)
    stop('k parameter must be a non-negative integer.\n')
  if (order[1]+order[3]+k<=0)
    stop('At least one of p, q or k must be positive.\n')
  if (order[2]>0&include.mean) {
    warning('"include.mean" is ignored since integer differencing is specified.\n')
    include.mean <- FALSE
  }

  if (method=='QML')  # Major problems with this method. It is not giving accurate results.
    stop('The QML method is not currently supported due to likely bugs.\n')

  allowed_methods <- c('CSS','Whittle','WLL','QML')
  if (!method%in%allowed_methods)
    stop('Method must be one of CSS, Whittle, QML or WLL.\n')
  if (method=='QML'&k>1)
    stop('QML method does not support k>1. It is suggested you try either the CSS or Whittle methods.\n')

  if (is.null(opt_method)) {
    if (method=='WLL') opt_method <- c('optim','solnp')
        else opt_method <- c('directL','solnp')
  }

  for (om in opt_method) {
    if (!om%in%.supported_optim())
      stop(sprintf('\nError: function %s not available.\n\nSupported functions are:\n%s\n', om, .supported_optim()))

    optimisation_packages <- .optim_packages()
    if (!isNamespaceLoaded(optimisation_packages[[om]]))
      stop(sprintf('Package %s is required to use method %s\n',optimisation_packages[[om]],om))

  }
  # check d_lim
  if (!is.numeric(d_lim)|length(d_lim)!=2|d_lim[1]>d_lim[2]) {
    stop(cat('parameter d_lim should be a list of 2 numerics, Eg c(0,0.5) and the minimum should be <= maximum.\n'))
  }
  # check min_freq and max_freq
  if (!is.numeric(freq_lim)|length(freq_lim)!=2|freq_lim[1]<0|freq_lim[2]>0.5|freq_lim[1]>freq_lim[2])
    stop('freq_lim must be numeric and between 0 and 0.5 and the minimum should be <= maximum.\n')

  if (missing(control)|is.null(control)) control <- list(tol=1e-15,maxeval=10000,max_eval=10000,maxit=10000,trace=0,delta=1e-10)

  ##
  ## 2. Next calc parameter  estimates
  p=as.integer(order[1])
  d=as.integer(order[2])
  q=as.integer(order[3])
  storage.mode(x) <- 'double'

  # drift term
  if (include.drift) {
    drift_term <- 1:length(x)
    lm1 <- stats::lm(as.numeric(x)~1+drift_term)
    lm_coef <- stats::coef(lm1)
    drift_const <- lm_coef[1]
    drift<-lm_coef[2]
    drift_se <- summary(lm1)$coef[2,2]
    drift_const_se <- summary(lm1)$coef[1,2]
    x1<-as.numeric(stats::residuals(lm1))+lm_coef[1]
  } else {
    x1 <- x
    drift <- 0
    drift_const <- 0
    drift_se<-0
    drift_const_se <- 0
  }

  if (d>0) y <- diff(x1,differences=d) else y <- x1
  mean_y <- mean(y)
  sd_y   <- stats::sd(y)
  ss<-stats::spectrum(as.numeric(y-mean_y),plot=FALSE,detrend=TRUE,demean=TRUE,method='pgram',taper=0,fast=FALSE)

  # Now set up params and lb (lower bounds) and ub (upper bounds)
  n_pars   <- 0
  pars     <- numeric(0)
  lb       <- numeric(0)
  ub       <- numeric(0)

  mean_methods <- c('QML','CSS')
  if (include.mean&method%in%mean_methods) {
    n_pars    <- n_pars+1
    mean_y    <- mean(y)
    pars      <- c(pars,mean_y)
    lb <- c(lb,ifelse(mean_y<0, 2*mean_y, -2*mean_y))
    ub <- c(ub,ifelse(mean_y<0,-2*mean_y,  2*mean_y))
  }

  if (k>0) {# initial parameter estimates for Gegenbauer factors
    gf <- ggbr_semipara(y,k=k,min_freq=freq_lim[1],max_freq=freq_lim[2])
    for (k1 in seq_len(k)) {
      n_pars    <- n_pars+2
      gf1 <- gf$ggbr_factors[[k1]]
      start_u <- gf1$u
      start_d <- gf1$fd

      max_u <- cos(2*pi*freq_lim[1])
      min_u <- cos(2*pi*freq_lim[2])
      #params
      if (start_u< min_u|start_u > max_u) start_u <- (min_u+max_u)/2
      if (start_d>=d_lim[2])  start_d <- d_lim[2]-0.01
      if (start_d<=d_lim[1])  start_d <- d_lim[1]+0.01

      pars      <- c(pars,start_u,start_d)
      lb        <- c(lb,min_u,d_lim[1])
      ub        <- c(ub,max_u,d_lim[2])
    }
  }

  n_pars    <- n_pars + p + q
  methods_to_estimate_var <- c('WLL')     # WLL method estimates VAR along with other params; For other methods this falls out of the objective value
  if (p+q>0) {
    # if any ARMA params to be estimated, we use the semipara estimates to get ggbf factors then get the underlying ARMA process,
    # and then ask "arima" for estimates. Semi para estimates should make good starting points for optimisation.
    # if (k>0) arma_y <- extract_arma(y,gf$ggbr_factors)
    # else arma_y <- y
    # a <- stats::arima(arma_y,order=c(p,0,q),include.mean=FALSE)
    # pars <- c(pars,a$coef)
    pars <- c(pars,rep(0,p+q))
    if (method%in%methods_to_estimate_var) {
      pars <- c(pars,sd(y))
    }
    if (p==1&q==0) { # special limits for AR(1)
      lb<-c(lb,-1)
      if (method%in%methods_to_estimate_var) lb <- c(lb, 1e-10)
      ub<-c(ub,1)
      if (method%in%methods_to_estimate_var) ub<-c(ub,2*stats::var(y))
    } else {
      lb <- c(lb,rep(-10,p+q))
      if (method%in%methods_to_estimate_var) lb <- c(lb,1e-10)
      ub <- c(ub,rep(10,p+q))
      if (method%in%methods_to_estimate_var) ub <- c(ub,2*stats::var(y))
    }
  }

  # create a list of all possible params any 'method' might need.
  # The various objective functions can extract the parameters which are relevant to that method.
  params <- list(y=y, orig_y=x, ss=ss, p=p,q=q,d=d,k=k,
                 include.mean=include.mean, include.drift=include.drift, drift=drift,
                 method=method,
                 est_mean=ifelse(method%in%mean_methods,TRUE,FALSE),
                 scale=sd_y,m_trunc=m_trunc)
  message <- c()

  # Optimisation functions for each method
  fcns  <- list('CSS'=.css.ggbr.obj,'Whittle'=.whittle.garma.obj.short,'QML'=.qml.ggbr.obj,'WLL'=.wll.ggbr.obj)
  if (opt_method[[1]]=='best') {
    fit <- .best_optim(initial_pars=pars, fcn=fcns[[method]],
                       lb=lb, ub=ub, params=params, control=control)
  } else {
      fit <- .generic_optim_list(opt_method_list=opt_method, initial_pars=pars, fcn=fcns[[method]],
                                 lb=lb, ub=ub, params=params, control=control)
  }
  if (fit$convergence== -999) stop('ERROR: Failed to converge.\n')

  hh <- fit$hessian

  # sigma2
  if (method=='WLL') {
    # adjust sigma2 for theoretical bias...
    sigma2 <- fit$par[length(fit$par)] <- fit$par[length(fit$par)]/(2*pi) * exp(-digamma(1))
  }
  # else if (method=='QML')     sigma2 <- .qml.ggbr.se2(fit$par, params=params)
  else if (method=='CSS')     sigma2 <- fit$value[length(fit$value)]/length(y)  # Chung (1996)
  else if (method=='Whittle') sigma2 <- .whittle.garma.obj.short(fit$par,params) # GHR 2001.

  # log lik
  loglik <- numeric(0)
  if (method=='CSS')
    loglik <- -0.5* (fit$value[length(fit$value)]/sigma2 + length(y)*(log(2*pi) + log(sigma2)))
  # if (method=='QML')
  #   loglik <- -fit$value[length(fit$value)]
  if (method=='Whittle')
    loglik <- -0.5*(2*length(y)*log(2*pi)+ .whittle.garma.obj(fit$par,params))  #refer GHR 2001, Whittle (1953) thm 6.

  se <- numeric(length(fit$par))

  # check convergence. Unfortunately "solnp" routine uses positive values to indicate an error.
  conv_ok <- TRUE
  if (opt_method[[1]]=='solnp') conv_ok <- (fit$convergence==0)
  else conv_ok <- (fit$convergence>=0)

  vcov1 <- matrix(nrow=length(fit$par),ncol=length(fit$par))  # default
  if (conv_ok&method!='WLL'&!is.null(hh)) {
    # Next, find the se's for coefficients
    start<-1
    se<-c()   # default to set this up in the right environment
    # next check the hessian
    h_inv_diag <- diag(inv_hessian <- pracma::pinv(hh))
    if (method=='Whittle') {
      # Next line from GHR (2001) thm 3.2
      omega <- .whittle_garma_omega(fit$par,params)
      vcov1 <- pracma::pinv(omega)/length(y)
      se <- suppressWarnings(sqrt(diag(vcov1)))
    }
    if (method=='CSS') {
      se <- suppressWarnings(sqrt(h_inv_diag*sigma2*2))
      vcov1 <- inv_hessian*2*sigma2
    }
    # if (method=='QML') {
    #   se <- suppressWarnings(sqrt(h_inv_diag*length(y)))
    #   vcov1 <- inv_hessian*length(y)
    # }
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
  vcov_nm <- nm

  n_coef    <- ifelse(method%in%methods_to_estimate_var,length(fit$par)-1,length(fit$par))  # ignore var on end if it is there
  temp_coef <- fit$par[1:n_coef]
  temp_se   <- se[1:n_coef]
  if (include.mean&!method%in%mean_methods) {# then add an Intercept anyway; use Kunsch 1987 result for se
    temp_coef <- c(mean(y),temp_coef)
    mean_se <- sigma2 * .garma_spec_den_0(fit$par,params)/(2*pi) # Chung (1996)
    temp_se   <- c(mean_se, temp_se)
    n_coef <- n_coef+1
  }
  if (include.drift) {
    temp_coef <- c(temp_coef,drift)
    temp_se   <- c(temp_se,drift_se)
    nm <- c(nm,'drift')
    n_coef <- n_coef+1
  }
  coef <- t(matrix(c(temp_coef,temp_se),nrow=n_coef))
  colnames(coef) <- nm
  rownames(coef) <- c('','s.e.')
  colnames(vcov1) <- rownames(vcov1) <- tail(vcov_nm,nrow(vcov1))

  # build a ggbr_factors object
  gf <- list()
  if (k>0) {
    if (include.mean) start_idx <- 2
    else start_idx <- 1
    for (k1 in seq_len(k)) {
      u <- coef[1,start_idx]
      gf1 <- list(u=u, f=acos(u)/(2*pi), fd=coef[1,start_idx+1], m=NA, f_idx=NA)
      gf <- c(gf,list(gf1))
      start_idx <- start_idx+2
    }
  }
  class(gf) <- 'ggbr_factors'

  # set up the 'model' list
  mean_y<-beta0<-0
  if (include.mean&order[2]==0) {
    mean_y <- mean(y)
    if (colnames(coef)[1]=='intercept') beta0  <- coef[1]
    else beta0 <- mean_y
  }

  model <- list('phi'=coef[1,substr(colnames(coef),1,2)=='ar'],
                'theta'=coef[1,substr(colnames(coef),1,2)=='ma'],
                'Delta'=order[2],
                'beta0'=beta0,
                'drift'=drift,
                'drift_se'=drift_se,
                'drift_const'=drift_const,
                'drift_const_se'=drift_const_se)
  if (k>0) model <-c(model, 'ggbr_factors' = list(gf))  # get fitted values and residuals

  garma_version <- as.character(utils::packageVersion('garma'))

  res<-list('call' = match.call(),
            'series' = deparse(match.call()$x),
            'coef'=coef,
            'var.coef'=vcov1,
            'sigma2'=sigma2,
            'obj_value'=fit$value,
            'obj_par' = fit$par,
            'loglik'=loglik,
            'aic'=-2*loglik + 2*(n_coef+1),
            'model'=model,
            'convergence'=fit$convergence,
            'conv_message'=c(fit$message,message),
            'method'=method,
            'opt_method'=opt_method,
            'control'=control,
            'order'=order,
            'k'=k,
            'y'=x,
            'diff_y'=y,
            'y_start'=x_start,
            'y_end'=x_end,
            'y_freq'=x_freq,
            'include.mean'=include.mean,
            'include.drift'=include.drift,
            'd_lim'=d_lim,
            'freq_lim'=freq_lim,
            'm_trunc'=m_trunc,
            'fitted_avail'=FALSE,
            'garma_version'=garma_version)
  if (opt_method[1]=='best') res<-c(res,'opt_method_selected'=fit$best_method)

  class(res) <- 'garma_model'

  if (fitted) { # add fitted values
    res <- .fitted_values(res)

  } else {
    fitted_values<-list()
    resid_values<-list()
  }

  return(res)
}




#' Predict future values.
#'
#' Predict ahead using algorithm of (2009) Godet, F
#' "Linear prediction of long-range dependent time series", ESAIM: PS 13 115-134.
#' DOI: 10.1051/ps:2008015
#'
#' @param object (garma_model) The garma_model from which to predict the values.
#' @param n.ahead (int) The number of time periods to predict ahead. Default: 1
#' @param max_wgts (int) The number of past values to use when forecasting ahead.
#' By default, all available data is used.
#' @param ggbr_scale (logical) - whether or not to scale the Gegenbauer weights
#' to add up to 1. By default this is FALSE.
#' @param ... Other parameters. Ignored.
#' @return A "ts" object containing the requested forecasts.
#' @examples
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' predict(mdl, n.ahead=12)
#' @export
predict.garma_model<-function(object,n.ahead=1,max_wgts=length(object$diff_y),ggbr_scale=FALSE,...) {
  ## Start of Function "predict"

  if (n.ahead<=0|!is.numeric(n.ahead)|is.na(n.ahead)) {
    message('n.ahead must be an integer g.t. 0.')
    return(NA)
  }

  coef <- unname(object$coef[1,])
  p  <- object$order[1]
  id <- object$order[2]
  q  <- object$order[3]
  k  <- object$k
  y  <- as.numeric(object$diff_y)
  orig_y <- as.numeric(object$y)
  n <- length(orig_y)
  resid  <- as.numeric(object$resid)

  beta0 <- object$model$beta0
  mean_y <- mean(y)
  phi_vec <- c(1,-object$model$phi)
  theta_vec <- c(1,-object$model$theta)
  if (any(Mod(polyroot(phi_vec))<1)|any(Mod(polyroot(theta_vec))<1))
    warning('model estimates are not Stationary! Forecasts may become unbounded.\n')


  if (k>0) {
    # for (gf in object$model$ggbr_factors)
    #     ggbr_inv_vec <- pracma::conv(ggbr_inv_vec,.ggbr.coef(n+n.ahead+2*k,-gf$fd,gf$u))
    build_weights<-function(len) {
      wgts <-  1
      for (gf in object$model$ggbr_factors) {
        fctr <- .ggbr.coef(len+2,-gf$fd,gf$u)
        wgts <- pracma::conv(wgts,fctr)
      }
      # Next line multiplies and divides the various polynomials to get psi = theta * delta * ggbr / phi
      # pracma::conv gives polynomial multiplication, and pracma::deconv gives polynomial division.
      # we don't bother with the remainder. For non-ggbr models this may be a mistake.
      wgts <- pracma::conv(phi_vec,wgts)
      if (length(theta_vec)>1) wgts <- pracma::deconv(wgts,theta_vec)$q

      return( wgts[(len+1):2] )
    }

    wgts <- build_weights(length(y)+n.ahead+id)
    y_dash <- y - beta0 # if differenced or just include.mean=FALSE then beta0 is zero.
    gf <- object$model$ggbr_factors[[1]]
    totsum <- (1-2.0*gf$u+gf$u^2)^(gf$fd)
    for (h in 1:(n.ahead)) {
      yy <- y_dash
      wgts1 <- tail(wgts,max_wgts)
      yy <- tail(yy,max_wgts)
      next_forecast <- (-sum(yy*wgts1))
      if (ggbr_scale) next_forecast <- next_forecast / abs(sum(wgts1)) * abs(totsum)
      y_dash <- c(y_dash, next_forecast)
    }

    pred<-y_dash[(length(y)+1):length(y_dash)]
    if (id>0) {
      if (object$include.drift) pred <- pred + object$model$drift
      pred<-diffinv(pred+mean_y,differences=id,xi=tail(orig_y,id))
      if (length(pred)>n.ahead) pred <- tail(pred,n.ahead)
    } else {
      pred <- pred + beta0
      if (length(pred)>n.ahead) pred <- tail(pred,n.ahead)
      n <- length(orig_y)
      if (object$include.drift) pred <- pred + ((n+1):(n+length(pred)))*object$model$drift
    }
  } else { # ARIMA forecasting only
    y_dash <- y-beta0
    phi_vec <- rev(-phi_vec[2:length(phi_vec)])
    if (length(theta_vec)>1) theta_vec <- rev(-theta_vec[2:length(theta_vec)])
    else theta_vec<-numeric(0)  # length will be zero. thus not used.
    pp <- length(phi_vec)
    qq <- length(theta_vec)
    pred <- rep(beta0,n.ahead)

    for (i in 1:n.ahead) {
      if (p>0) {
        if (i>1) ar_vec <- tail(c(rep(0,pp),y_dash,pred[1:(i-1)]),pp)
        else ar_vec <- tail(c(rep(0,pp),y_dash),pp)
        pred[i] <- pred[i] + sum(phi_vec*ar_vec)
      }
      if (q>0) {
        if (i>qq) ma_vec <- rep(0,qq)
        else if (i>1) ma_vec <- tail(c(rep(0,qq),resid,rep(0,i-1)),qq)
        else ma_vec <- rep(0,qq)
        pred[i] <- pred[i] + sum(theta_vec*ma_vec)
      }
    }
    if (id>0) {
      if (object$include.drift) pred <- pred + object$model$drift  # mean(y)
      pred<-diffinv(pred,differences=id,xi=tail(orig_y,id))
      if (length(pred)>n.ahead) pred <- tail(pred,n.ahead)
    } else {
      pred <- pred + beta0
      n <- length(orig_y)
      if (object$include.drift) pred <- pred + ((n+1):(n+length(pred)))*object$model$drift
    }
  }

  # Now we have the forecasts, we set these up as a "ts" object
  y_end = object$y_end
  if(length(y_end)>1) {
    if (object$y_freq >= y_end[2]) {
      y_end[1] <- y_end[1]+1
      y_end[2] <- y_end[2]-object$y_freq+1
    }
    else y_end[2] <- y_end[2] + 1
  } else y_end <- y_end +1

  res <- stats::ts(pred,start=y_end,frequency=object$y_freq)
  return(list(pred=res))
}

#' Predict2 future values.
#'
#' Predict ahead using algorithm of (2009) Godet, F
#' "Linear prediction of long-range dependent time series", ESAIM: PS 13 115-134.
#' DOI: 10.1051/ps:2008015
#'
#' @param object (garma_model) The garma_model from which to predict the values.
#' @param n.ahead (int) The number of time periods to predict ahead. Default: 1
#' @return A "ts" object containing the requested forecasts.
#' @export
predict2<-function(object, n.ahead=1) {
  ## Start of Function "predict"

  if (n.ahead<=0|!is.numeric(n.ahead)|is.na(n.ahead)) {
    message('n.ahead must be an integer g.t. 0.')
    return(NA)
  }

  coef <- unname(object$coef[1,])
  p  <- object$order[1]
  id <- object$order[2]
  q  <- object$order[3]
  k  <- object$k
  y  <- as.numeric(object$diff_y)
  orig_y <- as.numeric(object$y)
  n <- length(orig_y)
  resid  <- as.numeric(object$resid)

  beta0 <- object$model$beta0
  mean_y <- mean(y)
  phi_vec <- c(1,-object$model$phi)
  theta_vec <- c(1,-object$model$theta)
  if (any(Mod(polyroot(phi_vec))<1)|any(Mod(polyroot(theta_vec))<1))
    warning('model estimates are not Stationary! Forecasts may become unbounded.\n')

  gf<-object$model$ggbr_factors[[1]]
  #  fctr <- .ggbr.coef(len+2,-gf$fd,gf$u)

  legendre_array<-function(n,u,fd) {
    legendre_initial_1<-function(u,fd) {
      res1 <- ( (1.0+u)/(1.0-u) )^(fd-0.25) / gamma(1.5-2.0*fd)
      res2 <- Re(hypergeo::hypergeo(0.5,0.5,1.5-2.0*fd,(1.0-u)/2))
      return(res1*res2)
    }
    legendre_initial_2<-function(u,fd) {
      res1 <- ( (1.0+u)/(1.0-u) )^(fd-0.25) / gamma(1.5-2.0*fd)
      res3 <- Re(hypergeo::hypergeo(-0.5,1.5,1.5-2.0*fd,(1.0-u)/2))
      return(res1*res3)
    }

    res <- rep(NA,n+3)
    res[1] <- legendre_initial_1(u,fd)
    res[2] <- legendre_initial_2(u,fd)

    b <- 2.0*fd-0.5
    for (j in 1:(n+1)) {
      a <- j-0.5
      res[j+2] <- (2.0*a-1.0)/(a-b)*u*res[j+1] - (a+b-1.0)/(a-b)*res[j]
    }

    return(res)
  }
  calc_acf<-function(n,u,fd) {
    pos_leg_fcn <- legendre_array(n,u,fd)
    neg_leg_fcn <- legendre_array(n,-u,fd)
    nu <- acos(u)
    acf_array <- rep(0,n+3) # don't store acf[0] since this will always be 1.0
    for (j in 0:(n+2)) {
      if (j%%2==0) acf_array[j+1] <- pos_leg_fcn[j+1]+neg_leg_fcn[j+1]
      else acf_array[j+1] <- pos_leg_fcn[j+1]-neg_leg_fcn[j+1]
    }

    sign_sin_nu <- sign(sin(nu))
    acf_array <- acf_array * gamma(1.0-2.0*fd)/(2.0*sqrt(pi)) * abs(2.0*sin(nu))^(0.5-2.0*fd) * sign_sin_nu
    acf_array <- acf_array / acf_array[1]
    return(acf_array)
  }

  print(gf$u)
  print(gf$fd)
  g <- calc_acf(n+n.ahead,gf$u,gf$fd)
  print(g)
  pred <- ltsa::TrenchForecast(y, g, mean_y, n, n.ahead)

  # Now we have the forecasts, we set these up as a "ts" object
  y_end = object$y_end
  if(length(y_end)>1) {
    if (object$y_freq >= y_end[2]) {
      y_end[1] <- y_end[1]+1
      y_end[2] <- y_end[2]-object$y_freq+1
    }
    else y_end[2] <- y_end[2] + 1
  } else y_end <- y_end +1

  res <- stats::ts(pred,start=y_end,frequency=object$y_freq)
  return(list(pred=res))
}

#' Forecast future values.
#'
#' The forecast function predicts future values of a "garma_model" object, and is exactly the same as the "predict" function with slightly different parameter values.
#' @param object (garma_model) The garma_model from which to forecast the values.
#' @param h (int) The number of time periods to predict ahead. Default: 1
#' @param ... Other parameters passed to the forecast function. For "garma_model" objects, these are ignored.
#' @return - a "ts" object containing the requested forecasts.
#' @examples
#' library(forecast)
#'
#' data(AirPassengers)
#' ap  <- as.numeric(diff(AirPassengers,12))
#' mdl <- garma(ap,order=c(9,1,0),k=0,method='CSS',include.mean=FALSE)
#' forecast(mdl, h=12)
#' @export
forecast.garma_model<-function(object,h=1,...) {
  res <- predict.garma_model(object,n.ahead=h)
  return(list(mean=res$pred))
}

.printf<-function(val) {
  if (class(val)[1] == 'integer') fmtstr <- '%s: %d\n'
  else fmtstr <- '%s: %f\n'
  cat(sprintf(fmtstr,as.character(substitute(val)),val))
}

.fitted_values<-function(object) {#par,params,ggbr_factors,sigma2) { # Generate fitted values and residuals for GARMA process
  # y <- as.numeric(params$y)
  # orig_y <- as.numeric(params$orig_y)
  # p <- params$p
  # q <- params$q
  # id <- params$d
  # k <- params$k
  # include.mean <- params$include.mean
  # method <- params$method
  y <- object$diff_y
  orig_y <- object$y
  k <- object$k
  include.mean <- object$include.mean
  method <- object$method
  p <- object$order[1]
  id <- object$order[2]
  q <- object$order[3]
  sigma2 <- object$sigma2
  ggbr_factors <- object$model$ggbr_factors
  par <- object$obj_par

  beta0  <- 0
  start  <- 1
  if (include.mean) {
    if (method%in%c('CSS','QML')) {
      beta0  <- par[1]
      start <- 2
      } else beta0 <- mean(y)
  }

  # skip over Gegenbauer parameters
  start <- start + k*2

  n       <- length(y)
  phi_vec <- theta_vec <- 1
  if (p>0) phi_vec   <- object$model$phi #par[start:(start+p-1)]
  if (q>0) theta_vec <- object$model$theta #par[(p+start):(length(par))]
  # testing
  # if (q==2) theta_vec <- c(0.2357262, -0.2934927)
  # if (p==2) phi_vec <- c(0.2396479, -0.1674436)

  inf_ar_vec <- c(1, -phi_vec)
  if (k>0) { # multiply by inverse ggbr expansion polynomial
    for (gf in ggbr_factors) inf_ar_vec <- pracma::conv(inf_ar_vec,.ggbr.coef(n,-gf$fd,gf$u))
    # Divide by theta vec; ignore the remainder. That can be a problem for k=0, q>0 models.
    inf_ar_vec <- pracma::deconv(inf_ar_vec,theta_vec)$q
  }
  # finally get rid of the first "1"
  inf_ar_vec <- c(inf_ar_vec, rep(0,n))  # extra zeros in case this is pure AR
  inf_ar_vec <- inf_ar_vec[2:(n+1)]

  if (id==0) y_dash <- y-beta0 else y_dash <- y-mean(y)
  if (p>0) start<-p else start<-1

  fitted <- y_dash[1:start] # Assume we have the first p correct... Thus first p residuals are 0.
  for (i in (start+1):n) {
    yy_to   <- (i-1)
    yy_from <- ifelse(k==0,yy_to-p,1)
    if (yy_from<1) yy_from <- 1
    yy <- y_dash[yy_from:yy_to]
    vec <- inf_ar_vec[yy_to:yy_from]
    fitted <- c(fitted, -sum(yy*vec))
  }

  if (id>0) {
    fitted <- diffinv(fitted,differences=id,xi=head(orig_y,id))
    n <- length(orig_y)
    if (length(fitted)>n) fitted <- utils::head(fitted,n)
  } else {
    fitted <- fitted + beta0
    n <- length(orig_y)
  }
  resid=orig_y-fitted

  object$fitted <- fitted
  object$residuals  <- resid
  object$fitted_avail <- TRUE
  return(object)
}

#' Fitted values
#'
#' Fitted values are 1-step ahead predictions.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) array of 1-step ahead fitted values for the model.
#' @export
fitted.garma_model<-function(object,...) {
  .byRef(object)  # allows us to update the values of object
  if (!object$fitted_avail) {
    object <- .fitted_values(object)
  }
  return(object$fitted)
}

#' Residuals
#'
#' Response Residuals from the model.
#' @param object The garma_model object
#' @param type (chr) The type of residuals. Must be 'response'.
#' @param h (int) The number of periods ahead for the residuals. Must be 1.
#' @param ... Other parameters. Ignored.
#' @return (double) array of resideuals from the model.
#' @export
residuals.garma_model<-function(object,type='response',h=1,...) {
  .byRef(object)  # allows us to update the values of object
  if (!missing(type)) {
    if (type!='response') stop('Only response residuals are available.')
  }
  if (!missing(h))
    if (h!=1) stop('Only h=1 response residuals are available.')
  if (!object$fitted_avail) {
    object <- .fitted_values(object)
  }
  return(object$residuals)
}

#' Model Coefficients
#'
#' Model Coefficients/parameters.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) array of parameter value estimates from the fitted model.
#' @export
coef.garma_model<-function(object,...) {
  return(object$coef[1,])
}

#' AIC for model
#'
#' AIC for model if available.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) Approximate AIC - uses approximation of whichever methoid is used to find model params.
#' @export
AIC.garma_model<-function(object,...) {
  return(object$aic)
}

#' Covariance matrix
#'
#' Covariance matrix of parameters if available
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return (double) estimated variance-covariance matrix of the parameter estimates
#' @export
vcov.garma_model<-function(object,...) {
  return(object$var.coef)
}

#' Log Likelihood
#'
#' Log Likelihood, or approximate likelihood or part likelihood, depending on the method.
#' @param object The garma_model object
#' @param ... Other parameters. Ignored.
#' @return Object of class "logLik" with values for the (approx) log-likelihood for the model
#' @export
logLik.garma_model<-function(object,...) {
  # Need to figure out how to indicate these are REML estimates not true LL.
  res <-  structure(object$loglik, df=length(object$y)-1, nobs=length(object$y), class="logLik")
  return(res)
}

#' garma package version
#'
#' The version function returns the garma package version.
#' @return The package version.
#' @examples
#' library(garma)
#' garma::version()
#' @export
version<-function() {message(.getPackageVersion('garma'))}

