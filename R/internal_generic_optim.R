#' @noRd
internal_generic_optim <- function(opt_method, initial_pars, fcn, lb, ub, params, control) {
  error_handling <- function(cond) {
    save_fit$message <- cond
    return(save_fit)
  }

  nlopt_control <- list(maxeval = control[["max_eval"]], ftol_rel = control[["tol"]], xtol_rel = 0) # ,check_derivatives=TRUE)
  fit <- save_fit <- list(value = Inf, convergence = -999, par = initial_pars, pars = initial_pars) # default value
  if (opt_method == "solnp") {
    fit <- Rsolnp::solnp(pars = initial_pars, fun = fcn, LB = lb, UB = ub, control = control, params = params)
    nm <- names(fit)
    nm[which(nm == "pars")] <- "par"
    names(fit) <- nm
  }
  if (opt_method == "gosolnp") {
    fit <- Rsolnp::gosolnp(
      fun = fcn, LB = lb, UB = ub, n.restarts = 100, n.sim = 10000,
      control = list(outer.iter = 100, trace = 0), params = params
    )
    nm <- names(fit)
    nm[which(nm == "pars")] <- "par"
    names(fit) <- nm
  }

  if (opt_method == "cobyla") {
    fit <- tryCatch(
      nloptr::cobyla(
        x0 = initial_pars, fn = fcn, lower = lb, upper = ub,
        params = params, control = nlopt_control
      ),
      error = error_handling
    )
  }

  if (opt_method == "directL") {
    fit <- tryCatch(nloptr::directL(fn = fcn, lower = lb, upper = ub, params = params, control = nlopt_control),
      error = error_handling
    )
  }

  if (opt_method=='lbfgs') {
    fit <- tryCatch(nloptr::lbfgs(x0=initial_pars, fn=fcn,
                                  lower=lb, upper=ub, params=params, control=nlopt_control),
                    error=error_handling
    )
  }


  if (!is.null(fit$hessian)) fit$hessian <- pracma::hessian(fcn, fit$par, params = params)

  return(fit)
}


internal_generic_optim_list <- function(opt_method_list, initial_pars, fcn,
                                        lb, ub, params, control) {
  # chain through the opt methods using the optimal value from the last one as the initial value for the next one.

  results <- list()
  for (opt_method in opt_method_list) {
    if (internal_is_finite_bounds(opt_method)) {
      fit <- internal_generic_optim(
        opt_method = opt_method, initial_pars = initial_pars, fcn = fcn,
        lb = lb, ub = ub, params = params, control = control
      )
    } else {
      fit <- internal_generic_optim(
        opt_method = opt_method, initial_pars = initial_pars, fcn = fcn,
        lb = lb, ub = ub, params = params, control = control
      )
    }

    initial_pars <- fit$par
    results <- c(results, list(fit))
  }

  return(results)
}

internal_optim_packages <- function() {
  # return a list giving the package required for each method
  return(c("cobyla" = "nloptr", "directL" = "nloptr", "lbfgs" = "nloptr", "solnp" = "Rsolnp", "gosolnp" = "Rsolnp", "best" = "stats"))
}
internal_supported_optim <- function() {
  # which methods are supported?
  return(c("best", "cobyla", "directL", "solnp", "gosolnp", "lbfgs"))
}
internal_supported_contr_optim <- function() {
  # which methods are supported for inequality contraints?
  return(c("gosolnp", "solnp", "cobyla"))
}
internal_finite_bounds_methods <- function() {
  # return a lit showing which methods require finite box-bounds
  return(c("directL", "gosolnp"))
}
internal_is_finite_bounds <- function(opt_method) {
  # return TRUE or FALSE as to whether the method requires finite box-bounds
  return(opt_method %in% internal_finite_bounds_methods())
}
