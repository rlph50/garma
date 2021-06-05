#
.generic_optim<-function(opt_method, initial_pars, fcn, lb, ub, params, control) {
  error_handling <- function(cond) {
    save_fit$message=cond
    return(save_fit)
  }
  nlopt_control <- list(maxeval=control[['max_eval']],ftol_rel=control[['tol']],xtol_rel=0)#,check_derivatives=TRUE)
  fit <- save_fit <- list(value=Inf,convergence= -999,par=initial_pars,pars=initial_pars) # default value
  if (opt_method=='solnp') {
    fit <- Rsolnp::solnp(pars=initial_pars, fun=fcn, LB=lb, UB=ub, control=control, params=params)
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }
  if (opt_method=='gosolnp') {
    fit <- Rsolnp::gosolnp(fun=fcn, LB=lb, UB=ub, n.restarts = 100, n.sim = 10000,
                           control=list(outer.iter = 100, trace = 0), params=params)
    fit$par <- fit$pars  # copy across for consistency with other optimisation methods
  }

  if (opt_method=='optim') {
    fit <- tryCatch(stats::optim(initial_pars, fn=fcn, method="L-BFGS-B", hessian = TRUE,
                                  lower=lb, upper=ub, params=params),
                    error=error_handling
    )
  }
  if (opt_method=='lbfgs') {
    fit <- tryCatch(nloptr::lbfgs(x0=initial_pars, fn=fcn,
                                  lower=lb, upper=ub, params=params, control=nlopt_control),
                    error=error_handling
    )
  }
  if (opt_method=='slsqp') {
    fit <- tryCatch(nloptr::slsqp(x0=initial_pars, fn=fcn, lower=lb, upper=ub,
                                  params=params, control=nlopt_control),
                    error=error_handling
    )
  }

  if (opt_method=='cobyla') {
    fit <- tryCatch(nloptr::cobyla(x0=initial_pars, fn=fcn, lower=lb, upper=ub,
                                   params=params, control=nlopt_control),
                    error=error_handling
    )
  }

  if (opt_method=='directL') {
    fit <- tryCatch(nloptr::directL(fn=fcn, lower=lb, upper=ub, params=params, control=nlopt_control),
                    error=error_handling
    )
  }

  if (opt_method=='BBoptim') {
    bb_control <- list(maxit=control[['max_eval']],ftol=control[['tol']],trace=control[['trace']],eps=control[['tol']])
    fit <- tryCatch(BB::BBoptim(par=initial_pars, fn=fcn,lower=lb, upper=ub,
                                control=bb_control, params=params, quiet=(control[['trace']]>0)),
                    error=error_handling
    )
    # fit <- BB::BBoptim(par=initial_pars, fn=fcn, gr=grad, lower=lb, upper=ub,
    #                    control=bb_control, params=params, quiet=FALSE)
  }
  if (opt_method=='psoptim') {
    ps_control<-list(maxf=control[['max_eval']],maxit=control[['max_eval']],reltol=control[['tol']],trace=control$trace)#,hybrid=TRUE
    fit <- tryCatch(pso::psoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=ps_control),
                    error=error_handling
    )
    # fit <- pso::psoptim(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=ps_control)

  }
  if (opt_method=='ga') { # Evolutionary algorithm
    # TODO: enable parallel=TRUE
    fitness_fcn <- function(...) return(-fcn(...,params=params))
    res <- GA::de(fitness=fitness_fcn,lower=lb,upper=ub,monitor=FALSE)
    fit <- list(value=res@fitnessValue,
                convergence=0,
                par=res@solution[1,],
                hessian=pracma::hessian(fcn, res@solution[1,], params=params))
  }
  if (opt_method=='hjkb') {
    hjkb_control <- list(maxfeval=control[['max_eval']],tol=control[['tol']],info=(control[['trace']]>0))
    fit <- tryCatch(dfoptim::hjkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=hjkb_control),
                    error=error_handling
    )
  }
  if (opt_method=='nmkb') {
    nmkb_control <- list(maxfeval=control[['max_eval']],tol=control[['tol']],trace=control[['trace']],restarts.max=10)
    fit <- tryCatch(dfoptim::nmkb(par=initial_pars, fn=fcn, lower=lb, upper=ub, params=params, control=nmkb_control),
                    error=error_handling
    )
  }

  if (!is.null(fit$hessian)) det_hessian <- det(fit$hessian) else det_hessian <- NA
  if (length(fit$par>0)) {
    if(is.na(det_hessian)|det_hessian==1|is.null(fit$hessian)) fit$hessian <- pracma::hessian(fcn, fit$par, params=params)
    det_hessian <- det(fit$hessian)
    if (is.null(fit$hessian)|is.na(det_hessian)) {
      if (!any(is.na(fit$par))) fit$hessian <- pracma::hessian(fcn, fit$par, params=params)
    }
  }
  if (length(fit$value)>1)
    fit$value <- fit$value[length(fit$value)]

  return(fit)
}

.best_optim <- function(initial_pars, fcn, lb, ub, params, control) {
  best_val    <- Inf
  best_method <- 'none'
  best_fit    <- list()
  for (opt_method in .supported_optim()) {
    fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                          lb=lb, ub=ub, params=params, control=control)
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
                                lb, ub, params, control) {
  # chain through the opt methods using the optimal value from the last one as the initial value for the next one.

  for (opt_method in opt_method_list) {
    if (.is_finite_bounds(opt_method)) {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb, ub=ub, params=params, control=control)
    } else {
      fit <- .generic_optim(opt_method=opt_method, initial_pars=initial_pars, fcn=fcn,
                            lb=lb, ub=ub, params=params, control=control)
    }

    initial_pars <- fit$par
  }
  return(fit)
}

.optim_packages<-function() {
  # return a list giving the package required for each method
  return(c('optim'='stats','lbfgs'='nloptr', 'cobyla'='nloptr','directL'='nloptr', 'slsqp'='nloptr',
           'BBoptim'='BB','psoptim'='pso','hjkb'='dfoptim','nmkb'='dfoptim',
           'solnp'='Rsolnp','gosolnp'='Rsolnp','ga'='GA','best'='stats'))
}
.supported_optim<-function() {
  # which methods are supported?
  return(c('best', 'optim', 'lbfgs', 'cobyla','directL','slsqp', 'BBoptim','psoptim','hjkb','nmkb','solnp','gosolnp','ga'))
}
.supported_contr_optim<-function() {
  # which methods are supported for inequality contraints?
  return(c('gosolnp','solnp','cobyla','slsqp'))
}
.finite_bounds_methods<-function() {
  # return a lit showing which methods require finite box-bounds
  return(c('directL','psoptim','gosolnp','ga'))
}
.is_finite_bounds<-function(opt_method) {
  # return TRUE or FALSE as to whether the method requires finite box-bounds
  return(opt_method%in%.finite_bounds_methods())
}
