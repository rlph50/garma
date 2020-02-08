garma<-function(py,
                order=list(0,0,0),
                k=1,
                include.mean=TRUE,
                method='CSS',
                allow_neg_d=TRUE,
                maxeval=10000,
                opt.method='optim',
                m_trunc=50) {
  if (length(py)<12)
    stop('y should have at least 12 observations.')
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
  allowed_methods <- c('CSS','Whittle','BNP','QML')
  if (!method%in%allowed_methods)
    stop('Method must be one of CSS, Whittle, QML or BNP.')

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
  methods_to_estimate_var <- c('BNP')#,'QML')
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
  params <- list(y=y, ss=ss, p=p,q=q,k=k,include.mean=include.mean,scale=sd_y)
  message <- c()

  # First we make a first pass at optimisation using "optim".
  # If the method chosen is optim then that finishes things; but otherwise the solution found becomes the starting point for the next optimisation.
  fcns <- list('CSS'=css.ggbr.obj,'Whittle'=whittle.ggbr.obj,'QML'=qml.ggbr.obj,'BNP'=bnp.ggbr.obj)
  fit <- optim(par=pars, fn=fcns[[method]], lower=lb, upper=ub, params=params,
               hessian=TRUE,method="L-BFGS-B",control=list(maxit=maxeval,factr=1e-25))
  if (fit$convergence==52) {   # Error in line search, then try again
    cat('1st optimisation using optim - failed in line search.\n2nd optimisation using nloptr.\n')
    fit2<-lbfgs(x0=fit$par, fn=fcns[[method]], lower=lb, upper=ub, params=params, control=list(maxeval=maxeval,xtol_rel=1e-8))
    if (fit2$value<fit$value&fit2$convergence>=0) {
      fit<-fit2
      cat('Using 2nd fit.\n')
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

  if (fit$convergence>=0&method!='BNP') {
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
  if (method=='BNP') {
    se<-rep(NA,length(par))
    if (k>0) se[2] <- bnp_d_se(fit$par[1],ss)
  }
  nm<-list()
  if (include.mean&method%in%mean_methods) nm <- c(nm,'Intercept')
  if (k==1) nm<-c(nm,'u','fd')
  if (p>0) nm<-c(nm,paste0('AR',1:p))
  if (q>0) nm<-c(nm,paste0('MA',1:q))

  if (method=='BNP') {
    # adjust sigma2 for bias...
    fit$par[length(fit$par)] <- fit$par[length(fit$par)]/(2*pi) * exp(-digamma(1))
  }

  n_coef <- ifelse(method%in%methods_to_estimate_var,length(fit$par)-1,length(fit$par))  # ignore var on end if it is there
  coef <- t(matrix(round(c(fit$par[1:n_coef],se[1:n_coef]),4),nrow=n_coef))
  colnames(coef) <- nm
  rownames(coef) <- c('coef','se')
  if (method=='BNP')          sigma2 <- fit$par[length(fit$par)]
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
            "method"=method,
            'opt.method'=opt.method,
            'maxeval'=maxeval,
            'order'=order,'k'=k,'y'=py,'m_trunc'=m_trunc)
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
