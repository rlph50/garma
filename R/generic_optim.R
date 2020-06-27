#
.generic_optim<-function(opt_method, initial_pars, fcn, lb, ub, params, ineq_fcn=NULL, ineq_lb=NULL, ineq_ub=NULL, max_eval=10000, tol=1e-15) {
  nlopt_control <- list(maxeval=max_eval,ftol_rel=tol)
  fit <- save_fit <- list(value=Inf,convergence= -999,par=initial_pars,pars=initial_pars) # default value
  if (opt_method=='solnp'&!is.null(ineq_fcn)) {
    tryCatch(fit <- Rsolnp::solnp(pars=initial_pars,
                                  fun=fcn,
                                  LB=lb,
                                  UB=ub,
                                  ineqfun=ineq_fcn,
                                  ineqLB=ineq_lb,
                                  ineqUB=ineq_ub,
                                  control=list(trace=0,tol=tol),
                                  params=params),
             error=function(cond) {fit<-save_fit}
    )
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }
  if (opt_method=='solnp'&is.null(ineq_fcn)) {
    tryCatch(fit <- Rsolnp::solnp(pars=initial_pars, fun=fcn, LB=lb, UB=ub, control=list(trace=0,tol=tol), params=params),
             error=function(cond) {fit<-save_fit}
    )
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
    tryCatch(fit <- BB::BBoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, control=list(trace=FALSE,maxit=max_eval,ftol=tol,gtol=tol), params=params,quiet=TRUE),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='psoptim') {
    tryCatch(fit <-pso::psoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=list(maxit=max_eval)),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='hjkb') {
    tryCatch(fit <- dfoptim::hjkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=list(maxfeval=max_eval)),
             error=function(cond) {fit<-save_fit}
    )
  }
  if (opt_method=='nmkb') {
    tryCatch(fit <- dfoptim::nmkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=list(maxfeval=max_eval)),
             error=function(cond) {fit<-save_fit}
    )
  }

  if (is.null(fit$hessian)|diag(fit$hessian)[2]==1)
    fit$hessian <- pracma::hessian(fcn, fit$par, params=params)
  if (length(fit$value)>1)
    fit$value <- fit$value[length(fit$value)]

  return(fit)
}

.best_optim <- function(initial_pars, fcn, lb, ub, lb_finite, ub_finite, params, max_eval=10000, tol=1e-15) {
  best_val    <- Inf
  best_method <- 'none'
  best_fit    <- list()
  for (opt_method in .supported_optim()) {
    if (.is_finite_bounds(opt_method)) {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn, lb=lb_finite, ub=ub_finite, params=params, max_eval=max_eval, tol=tol)
    } else {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn, lb=lb, ub=ub, params=params, max_eval=max_eval, tol=tol)
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
                                params, max_eval=10000, tol=1e-15) {
  # chain through the opt methods using the optimal value from the last one as the initial value for the next one.
  for (opt_method in opt_method_list) {
    if (.is_finite_bounds(opt_method)) {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb_finite, ub=ub_finite,
                            ineq_fcn=ineq_fcn, ineq_lb=ineq_lb, ineq_ub=ineq_ub,
                            params=params, max_eval=max_eval, tol=tol)
    } else {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb, ub=ub,
                            ineq_fcn=ineq_fcn, ineq_lb=ineq_lb, ineq_ub=ineq_ub,
                            params=params, max_eval=max_eval, tol=tol)
    }

    initial_pars <- fit$par
  }
  return(fit)
}

.optim_packages<-function() {
  # return a list giving the package required for each method
  return(c('optim'='stats','cobyla'='nloptr','directL'='nloptr', 'slsqp'='nloptr',
           'BBoptim'='BB','psoptim'='pso','hjkb'='dfoptim','nmkb'='dfoptim','solnp'='Rsolnp','best'='stats'))
}
.supported_optim<-function() {
  # which methods are supported?
  return(c('best', 'cobyla','directL','slsqp', 'BBoptim','psoptim','hjkb','nmkb','solnp'))
}
.supported_contr_optim<-function() {
  # which methods are supported for inequality contraints?
  return(c('solnp','cobyla','slsqp'))
}
.finite_bounds_methods<-function() {
  # return a lit showing which methods require finite box-bounds
  return(c('directL','psoptim'))
}
.is_finite_bounds<-function(opt_method) {
  # return TRUE or FALSE as to whether the method requires finite box-bounds
  return(opt_method%in%.finite_bounds_methods())
}
