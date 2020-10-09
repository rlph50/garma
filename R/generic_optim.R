#
.generic_optim<-function(opt_method, initial_pars, fcn, lb, ub, params, ineq_fcn=NULL, ineq_lb=NULL, ineq_ub=NULL, control) {
  nlopt_control <- list(maxeval=control[['max_eval']],ftol_rel=control[['tol']],xtol_rel=0)
  fit <- save_fit <- list(value=Inf,convergence= -999,par=initial_pars,pars=initial_pars) # default value
  if (opt_method=='solnp'&!is.null(ineq_fcn)) {
    tryCatch(fit <- Rsolnp::solnp(pars=initial_pars,
                                  fun=fcn,
                                  LB=lb,
                                  UB=ub,
                                  ineqfun=ineq_fcn,
                                  ineqLB=ineq_lb,
                                  ineqUB=ineq_ub,
                                  control=control,
                                  params=params),
             warning=function(cond) {fit$message<-cond},
             error=function(cond) {fit<-save_fit;fit$message<-cond}
    )
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }
  if (opt_method=='gosolnp'&!is.null(ineq_fcn)) {
    fit <- Rsolnp::gosolnp(fun=fcn, LB=lb, UB=ub, ineqfun=ineq_fcn, ineqLB=ineq_lb, ineqUB=ineq_ub, n.restarts = 100,
                           n.sim = 10000, control=list(outer.iter = 100, trace = 0), params=params)
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }
  if (opt_method=='solnp'&is.null(ineq_fcn)) {
    fit <- Rsolnp::solnp(pars=initial_pars, fun=fcn, LB=lb, UB=ub, control=control, params=params)
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }
  if (opt_method=='gosolnp'&is.null(ineq_fcn)) {
    fit <- Rsolnp::gosolnp(fun=fcn, LB=lb, UB=ub, n.restarts = 100, n.sim = 10000, control=list(outer.iter = 100, trace = 0), params=params)
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }

  if (opt_method=='slsqp'&!is.null(ineq_fcn)) {
    tryCatch(fit <- nloptr::slsqp(x0=initial_pars, fn=fcn, lower=lb, upper=ub, hin=ineq_fcn, params=params, control=nlopt_control),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='slsqp'&is.null(ineq_fcn)) {
    tryCatch(fit <- nloptr::slsqp(x0=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=nlopt_control),
             error=function(cond) {fit<-save_fit}
    )
  }

  if (opt_method=='cobyla'&!is.null(ineq_fcn)) {
    tryCatch(fit <- nloptr::cobyla(x0=initial_pars, fn=fcn, lower=lb, upper=ub, hin=ineq_fcn, params=params, control=nlopt_control),
           error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='cobyla'&is.null(ineq_fcn)) {
    tryCatch(fit <- nloptr::cobyla(x0=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=nlopt_control),
             error=function(cond) {fit<-save_fit}
    )
  }

  if (opt_method=='directL') {
    tryCatch(fit <- nloptr::directL(fn=fcn, lower=lb, upper=ub, params=params, control=nlopt_control),
             error=function(cond) {fit<-save_fit}
    )
  }

  if (opt_method=='BBoptim') {
    tryCatch(fit <- BB::BBoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, control=control, params=params,quiet=TRUE),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='psoptim') {
    tryCatch(fit <-pso::psoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=control),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='hjkb') {
    tryCatch(fit <- dfoptim::hjkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=control),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='nmkb') {
    tryCatch(fit <- dfoptim::nmkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=control),
             error=function(cond) {fit<-save_fit}
    )
  }

  if (!is.null(fit$hessian)) det_hessian <- det(fit$hessian) else det_hessian <- NA
  if (length(fit$par>0)) {
    if(is.na(det_hessian)|det_hessian==1) fit$hessian <- pracma::hessian(fcn, fit$par, params=params)
    det_hessian <- det(fit$hessian)
    if (is.null(fit$hessian)|is.na(det_hessian)) {
      if (!any(is.na(fit$par))) fit$hessian <- pracma::hessian(fcn, fit$par, params=params)
    }
  }
  if (length(fit$value)>1)
    fit$value <- fit$value[length(fit$value)]

  return(fit)
}

.best_optim <- function(initial_pars, fcn, lb, ub, lb_finite, ub_finite, params, control) {
  best_val    <- Inf
  best_method <- 'none'
  best_fit    <- list()
  for (opt_method in .supported_optim()) {
    if (.is_finite_bounds(opt_method)) {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn, lb=lb_finite, ub=ub_finite, params=params, control=control)
    } else {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn, lb=lb, ub=ub, params=params, control=control)
    }
    if (fit$value<best_val) {
      best_val    <- fit$value
      best_fit    <- fit
      best_method <- opt_method
    }
  }
  best_fit$best_method <- best_method
  return(best_fit)
}

.generic_optim_list <- function(opt_method_list, initial_pars, fcn,
                                lb, ub, lb_finite, ub_finite,
                                ineq_fcn=NULL, ineq_lb=NULL, ineq_ub=NULL,
                                params, control) {
  # chain through the opt methods using the optimal value from the last one as the initial value for the next one.

  for (opt_method in opt_method_list) {
    if (.is_finite_bounds(opt_method)) {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb_finite, ub=ub_finite,
                            ineq_fcn=ineq_fcn, ineq_lb=ineq_lb, ineq_ub=ineq_ub,
                            params=params, control=control)
    } else {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb, ub=ub,
                            ineq_fcn=ineq_fcn, ineq_lb=ineq_lb, ineq_ub=ineq_ub,
                            params=params, control=control)
    }

    initial_pars <- fit$par
  }
  return(fit)
}

.optim_packages<-function() {
  # return a list giving the package required for each method
  return(c('optim'='stats','cobyla'='nloptr','directL'='nloptr', 'slsqp'='nloptr',
           'BBoptim'='BB','psoptim'='pso','hjkb'='dfoptim','nmkb'='dfoptim','solnp'='Rsolnp','gosolnp'='Rsolnp','best'='stats'))
}
.supported_optim<-function() {
  # which methods are supported?
  return(c('best', 'cobyla','directL','slsqp', 'BBoptim','psoptim','hjkb','nmkb','solnp','gosolnp'))
}
.supported_contr_optim<-function() {
  # which methods are supported for inequality contraints?
  return(c('gosolnp','solnp','cobyla','slsqp'))
}
.finite_bounds_methods<-function() {
  # return a lit showing which methods require finite box-bounds
  return(c('directL','psoptim','gosolnp'))
}
.is_finite_bounds<-function(opt_method) {
  # return TRUE or FALSE as to whether the method requires finite box-bounds
  return(opt_method%in%.finite_bounds_methods())
}
