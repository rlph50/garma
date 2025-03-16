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
#' \item \eqn{(1-2u_{i}B+B^{2})^{d_{i}}}{(1-2u(i)B+B^2)^d(i)} represents the long-memory Gegenbauer component (there may in
#' general be k of these),
#' \item \eqn{id} represents the degree of integer differencing, where as \eqn{d_i} represents the degree of fractional
#' differencing. Note that \eqn{id} is a value supplied by the user (the second number on the `order=` parameter - similarly
#' to the way that the base R `arima` function works) whereas \eqn{d_i} is estimated by this function.
#' \item \eqn{X_{t}}{X(t)} represents the observed process,
#' \item \eqn{\epsilon_{t}}{\epsilon(t)} represents the random component of the model - these are assumed to be uncorrelated but
#'       identically distributed variates. Generally the routines in this package will work best if these have an approximate
#'       Gaussian distribution.
#' \item \eqn{B}{B} represents the Backshift operator, defined by \eqn{B X_{t}=X_{t-1}}{B X(t) = X(t-1)}.
#' }
#' when k=0, then this is just a short memory model as fit by the stats "arima" function.
#'
#' @param x (num) This should be a numeric vector representing the process to estimate. A minimum length of 96 is required.
#' @param order (numeric vector) This should be a vector (similar to the stats::arima order parameter) which will give the order
#'   of the process to fit. The format should be list(p,d,q) where p, d, and q are all positive integers. p represents the degree
#'   of the autoregressive process to fit, q represents the order of the moving average process to fit and d is the (integer)
#'   differencing to apply prior to any fitting. WARNING: Currently only d==0 or d==1 are allowed.
#' @param periods (num) This parameter can be used to specify a fixed period or set of periods for the
#' Gegenbauer periodicity. For instance if you have monthly data, then it might be sensible (after an examination of the
#' periodogram) to set periods = 12. The default value is NULL. Either `periods` or `k` parameters must be specified
#' but not both - `periods` implies fixed period(s) are to be used and `k` implies that the periods should
#' be estimated.
#' @param k (int) This parameter indicates that the algorithm should estimate the `k` frequencies as a part of the model.
#' An alternative is the `periods` parameter which can be used to specify exactly which periods should be used by
#' the model.
#'
#' This parameter can also be interpreted as specifying the number of (multiplicative) Gegenbauer terms to fit in the model.
#' @param include.mean (bool) A boolean value indicating whether a mean should be fit.
#'     Note that no mean term is fit if the series is integer differenced.
#' @param include.drift (bool) A boolean value indicating whether a 'drift' term should be fit to the predictions.
#'     The default is to fit a drift term to the predictions if the process is integer-differenced.
#' @param xreg (numeric matrix) A numerical vector or matrix of external regressors, which must have the same number of rows as x.
#'   It should not have any NA values. It should not be a data frame. The default value is NULL.
#'
#'    Note that the algorithm used here is that if any `xreg` is supplied, then a linear regression model is fit first, and the
#'    GARMA model is then based on the residuals from that regression model.
#' @param method (character) This defines the estimation method for the routine. The valid values are 'CSS', 'Whittle', and
#'   'WLL'. The default ('Whittle') method will generally return very accurate estimates quite quickly, provided the assumption
#'   of a Gaussian distribution is even approximately correct, and is probably the method of choice for most users. For the
#'   theory behind this, refer Giraitis et. al. (2001).
#'
#'   The 'CSS' method is a conditional 'sum-of-squares' technique and can be quite slow.
#'   Reference: Robinson (2006), Chung (1996). Note that the paper of Chung (1996) was partially critisised by Giraitis et.
#'   al. (2001), however still contains useful results.
#'
#'   'WLL' is a new technique, originally developed by the author of this package and which appears to work well even if the
#'   \eqn{\epsilon_{t}}{\epsilon(t)} are highly skewed and/or have heavy tails (skewed and/or lepto-kurtic). However the
#'   asymptotic theory for the WLL method is not complete and so standard errors are not available for most parameters.
#'   Refer Hunt et. al. (2021).
#' @param d_lim (list) the limits for the d parameter. The default is `c(0,0.5)`, which restricts the model to be stationary.
#'        However sometimes it is desirable to understand what the unrestricted value might be.
#' @param opt_method (character) This names the optimisation method used to find the parameter estimates.
#' This may be a list of methods, in which case the methods are applied in turn,
#' each using the results of the previous one as the starting point for the next. The default is to use c('solnp', 'cobyla').
#' For some data or some models, however, other methods may work well.
#'
#' Supported algorithms include:
#'  \itemize{
#'  \item 'cobyla' algorithm in package nloptr
#'  \item 'directL' algorithm in package nloptr
#'  \item 'solnp' from Rsolnp package
#'  \item 'gosolnp' from Rsolnp package.
#'  }
#'
#' Note that the algorithms are selected to be those which do not require derivatives, even numerically calculated
#' derivatives. The function being optimised by `garma()` has a point of discontinuity at the minimum value - the point
#' we are trying to find. This means that standard algorithms like BFGS et al. perform very poorly here.
#'
#' Note further that if you specify a value of `k` > 1, then inequality constraints are required, and this will further limit
#' the list of supported routines.
#' @param control (list) list of optimisation routine specific values.
#' @return An S3 object of class "garma_model".
#'
#' @references
#' C Chung. A generalized fractionally integrated autoregressive moving-average process.
#' Journal of Time Series Analysis, 17(2):111-140, 1996. DOI: https://doi.org/10.1111/j.1467-9892.1996.tb00268.x
#'
#' L Giraitis, J Hidalgo, and P Robinson. Gaussian estimation of parametric spectral density with unknown pole.
#' The Annals of Statistics, 29(4):987–1023, 2001. DOI: https://doi.org/10.1214/AOS/1013699989
#'
#' R Hunt, S Peiris, N Webe. A General Frequency Domain Estimation Method for Gegenbauer Processes.
#' Journal of Time Series Econometrics, 13(2):119-144, 2021. DOI: https://doi.org/10.1515/jtse-2019-0031
#'
#' R Hunt, S Peiris, N Weber. Estimation methods for stationary Gegenbauer processes.
#' Statistical Papers 63:1707-1741, 2022. DOI: https://doi.org/10.1007/s00362-022-01290-3
#'
#' P. Robinson. Conditional-sum-of-squares estimation of models for stationary time series with long memory.
#' IMS Lecture Notes Monograph Series, Time Series and Related Topics, 52:130-137, 2006.
#' DOI: https://doi.org/10.1214/074921706000000996.
#'
#' @examples
#' data(AirPassengers)
#' ap <- as.numeric(diff(AirPassengers, 12))
#' print(garma(ap, order = c(9, 1, 0), k = 0, method = "CSS", include.mean = FALSE))
#' # Compare with the built-in arima function
#' print(arima(ap, order = c(9, 1, 0), include.mean = FALSE))
#' @export
garma <- function(x,
                  order = c(0L, 0L, 0L),
                  periods = NULL,
                  k = 1,
                  include.mean = (order[2] == 0L),
                  include.drift = FALSE,
                  xreg = NULL,
                  method = "Whittle",
                  d_lim = c(0, 0.5),
                  opt_method = c("cobyla", "solnp"),
                  control = NULL) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check parameters
  # Check x
  if (is.data.frame(x)) {
    if (ncol(x) > 1) {
      rlang::abort("x should be a numeric vector - not an entire data frame. Please select a single column and try again.")
    } else {
      x <- x[, 1]
    }
  }
  if (length(x) < 96) {
    # This is a bit arbitary but I would be concerned about fitting to a process of length even 96.
    # But no real evidence for this restriction apart from simulations showing large std errs.
    rlang::abort("x should have at least 96 observations.")
  }

  # now save the ts elements, if any - we re-apply them to the output values
  x_start <- stats::start(x)
  x_end <- stats::end(x)
  x_freq <- stats::frequency(x)

  x <- as.numeric(x)
  if (!is.numeric(x)) {
    rlang::abort("x should be numeric or a ts object.\n")
  }
  if (any(is.na(x))) {
    rlang::abort("x should not have any missing values.\n")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check periods
  if (is.null(periods)) {
    if (is.null(k) | !is.numeric(k) | k < 0) {
      rlang::abort("The k parameter must be a non-negative integer.\n")
    }
  } else {
    if (!is.numeric(periods) | length(periods) == 0) {
      rlang::abort("The 'periods' parameter must be a numeric vector of at least 1 element.")
    }
    if (any(periods < 0)) {
      rlang::abort("The 'periods' parameter cannot contain negative values.")
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check 'order'
  if (length(order) != 3) {
    rlang::abort("The 'order' parameter must be a 3 integers only.\n")
  }
  if (any(order < 0)) {
    rlang::abort("The 'order' parameter must consist of positive integers.\n")
  }
  if (order[1] + order[3] + k + length(periods) <= 0) {
    rlang::abort("At least one of p, q or k (or periods) must be positive.\n")
  }
  if (order[2] > 0 & include.mean) {
    rlang::warn("The parameter 'include.mean' is ignored since integer differencing is specified in the order= parameter.")
    include.mean <- FALSE
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check xreg
  if (!is.null(xreg)) {
    if (!is.numeric(xreg)) {
      rlang::abort("The parameter 'xreg' should be a numeric vector or matrix.")
    }
    if (!is.vector(xreg) & ! is.matrix(xreg)) {
      rlang::abort("The parameter 'xreg' shoud be a vector or a matrix. Dataframes are specifically not supported.")
    }
    if (is.vector(xreg)) {
      if (length(xreg) != length(x)) {
        rlang::abort("The parameter 'xreg' should be the same length as 'x'.")
      }
      # convert vector to matrix
      xreg_name <- as.character(deparse(substitute(xreg)))
      xreg <- as.matrix(xreg)
      colnames(xreg) <- internal_make_column_name(xreg_name)

    } else if (is.matrix(xreg)) {
      if (nrow(xreg) != length(x)) {
        rlang::abort("The parameter 'xreg' should have the same number of rows as the length of 'x'.")
      }
    }
    if (any(is.na(xreg))) {
      rlang::abort("The parameter 'xreg' should not have any NA values.")
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check estimation method
  allowed_methods <- c("CSS", "Whittle", "WLL")
  if (!method %in% allowed_methods) {
    rlang::abort("The parameter 'method' must be one of CSS, Whittle or WLL.\n")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Check optimisation method
  unrecognised_opt_method <- setdiff(opt_method, internal_supported_optim())
  if (length(unrecognised_opt_method) > 0) {
      rlang::abort(
        paste0("The following optimisation routines are not supported: ",
               paste(unrecognised_opt_method, collapse = ", "), ".")
      )
  }
  rm(unrecognised_opt_method)

  optimisation_packages <- internal_optim_packages()
  required_opt_packages <- sapply(opt_method, function(x) optimisation_packages[[x]])
  required_opt_packages <- unique(required_opt_packages)
  packages_loaded <- sapply(required_opt_packages, isNamespaceLoaded)
  if (any(!packages_loaded)) {
    idx <- which(!packages_loaded)
    rlang::abort(paste0("The following packages are required but are not available: ", paste(names(packages_loaded)[idx], collapse = ", "), "."))
  }
  rm(optimisation_packages, required_opt_packages, packages_loaded)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # check d_lim
  if (any(is.na(d_lim))) {
    rlang::abort("The parameter 'd_lim' should not have any NA values.")
  }
  if (!is.numeric(d_lim) | length(d_lim) != 2 | d_lim[1] > d_lim[2]) {
    rlang::abort("The parameter 'd_lim' should be a list of 2 numerics, Eg c(0,0.5) and the minimum should be < maximum.\n")
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set up a default 'control'
  if (missing(control) | is.null(control)) {
    control <- list(tol = 1e-8, maxeval = 10000, max_eval = 10000, maxit = 10000, trace = 0, delta = 1e-8,
                    inner.iter = 100, outer.iter = 100)
  }


  result <- list(
    "call" = match.call(),
    "series" = deparse(match.call()$x),
    "method" = method,
    "opt_method" = opt_method,
    "control" = control,
    "order" = order,
    "k" = k,
    "periods" = periods,
    "y" = x,
    "y_start" = x_start,
    "y_end" = x_end,
    "y_freq" = x_freq,
    "include.mean" = include.mean,
    "include.drift" = include.drift,
    "xreg" = xreg,
    "d_lim" = d_lim,
    "garma_version" = as.character(utils::packageVersion("garma"))
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set up data for calculating parameter estimates
  p <- as.integer(order[1])
  d <- as.integer(order[2])
  q <- as.integer(order[3])
  storage.mode(x) <- "double"
  if (!is.null(xreg)) {
    # ensure we have column names for xreg
    if (is.null(colnames(xreg))) colnames(xreg) <- paste0("x", 1:ncol(xreg))
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # mean term
  if (include.mean) {
    mean_matrix <- matrix(data=rep(1,length(x)), ncol = 1)
    colnames(mean_matrix) <- "intercept"
    if (!is.null(xreg)) {
      xreg <- cbind(mean_matrix, xreg)
    } else {
      xreg <- mean_matrix
    }
    rm(mean_matrix)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # drift term
  result$lm_xreg <- lm1 <- NULL

  if (include.drift) {
    drift_matrix = matrix(data = 1:length(x), ncol = 1)
    colnames(drift_matrix) <- "drift"
    if (!is.null(xreg)) {
      xreg <- cbind(drift_matrix, xreg)
    } else {
      xreg <- drift_matrix
    }
    rm(drift_matrix)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # xreg estimate
  if (!is.null(xreg)) {
    df <- as.data.frame(xreg)
    df$temp_y_ <- x
    f <- stats::formula(paste0("temp_y_ ~ ", paste(colnames(xreg), collapse = "+"), "- 1"))
    result$lm_xreg <- lm1 <- stats::lm(f, data = df, model = FALSE)
    rm(df, f)

    # now we use the residuals for the time series model.
    y <- as.numeric(lm1$residuals)
    coef <- lm1$coefficients
    se <- sqrt(diag(vcov(lm1)))
  } else {
    y <- x
    coef <- numeric(0)
    se <- numeric(0)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # First adjust for integer differencing
  if (d > 0) {
    y <- diff(y, differences = d)
  } else {
    y <- y
  }
  # save this away in the result object
  result$diff_y <- y

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Now set up params and lb (lower bounds) and ub (upper bounds)
  pars <- numeric(0)
  lb <- numeric(0)
  ub <- numeric(0)
  sd_y <- stats::sd(y)
  ss <- internal_garma_pgram(y) # periodogram

  # First, Gegenbauer pars
  if (is.null(periods)) {
    # need to estimate both periods and fractional differencing
    gf <- internal_semipara(y, k = k, pgram = ss)
    for (k1 in seq_len(k)) {
      gf1 <- gf[[k1]]
      start_u <- cos(2*pi*gf1$freq)
      start_d <- gf1$fd

      max_u <- 1 - 1e-15
      min_u <- 1e-15

      if (is.null(start_u)) start_u <- (min_u + max_u) / 2
      if (is.na(start_u)) start_u <- (min_u + max_u) / 2
      if (start_u < min_u | start_u > max_u) start_u <- (min_u + max_u) / 2

      if (is.null(start_d)) start_d <- mean(d_lim)
      if (is.na(start_d)) start_d <- mean(d_lim)
      if (start_d >= d_lim[2]) start_d <- d_lim[2] - 0.01
      if (start_d <= d_lim[1]) start_d <- d_lim[1] + 0.01

      pars <- c(pars, start_u, start_d)
      lb <- c(lb, min_u, d_lim[1])
      ub <- c(ub, max_u, d_lim[2])
    }
  } else {
    # periods is not null
    gf <- internal_semipara(y, periods = periods, pgram = ss)
    for (gf1 in gf) {
      start_d <- gf1$fd
      if (is.null(start_d)) start_d <- mean(d_lim)
      if (is.na(start_d)) start_d <- mean(d_lim)
      if (start_d >= d_lim[2]) start_d <- d_lim[2] - 0.01
      if (start_d <= d_lim[1]) start_d <- d_lim[1] + 0.01

      pars <- c(pars, start_d)
      lb <- c(lb, d_lim[1])
      ub <- c(ub, d_lim[2])
    }
  } # end of: if is null(periods) ... else ...

  if (p + q > 0) {
    # ARMA params - experiment suggests starting with 0 is as effective as any other method.
    pars <- c(pars, rep(0, p + q))
    lb <- c(lb, rep(-1, p + q))
    ub <- c(ub, rep(1, p + q))
  }

  if (method == "WLL") {
    # The WLL method needs to estimate the VAR along with other parameters.
    pars <- c(pars, stats::sd(y))
    lb <- c(lb, 1e-10)
    ub <- c(ub, 2 * stats::var(y))
  }

  # create a list of all possible params any 'method' might need.
  # The various objective functions can extract the parameters which are relevant to that method.
  params <- list(
    y = y,
    ss = ss,
    p = p,
    q = q,
    d = d,
    k = k,
    periods = periods,
    method = method
  )
  message <- character(0)

  # Optimisation functions for each method
  fcns <- list(
    "CSS" = internal_css_ggbr_obj,
    "Whittle" = internal_whittle_garma_obj_short,
    "WLL" = internal_wll_ggbr_obj
  )

  # apply a list of optimisation methods in sequence
  # result will be a list of the "fit" results for each algorithm
  opt_result <- internal_generic_optim_list(
    opt_method_list = opt_method,
    initial_pars = pars,
    fcn = fcns[[method]],
    lb = lb,
    ub = ub,
    params = params,
    control = control
  )
  fit <- opt_result[[length(opt_result)]]
  result$opt_result <-list(opt_result)
  result$convergence <- fit$convergence
  result$conv_message <- fit$message

  hh <- fit$hessian

  # sigma2
  if (method == "WLL") {
    # adjust sigma2 for theoretical bias...
    # Note: digamma() function is part of base R.
    sigma2 <- fit$par[length(fit$par)] <- fit$par[length(fit$par)] / (2 * pi) * exp(-digamma(1))
  } else if (method == "CSS") {
    sigma2 <- fit$value[length(fit$value)] / length(y)
  } # Chung (1996)
  else if (method == "Whittle") sigma2 <- internal_whittle_garma_obj_short(fit$par, params) # GHR 2001.

  result$sigma2 <- sigma2

  if (!is.null(lm1)) {
    coef <- c(lm1$coefficients, fit$par)
    se <- sqrt(diag(vcov(lm1)))
  } else {
    coef <- fit$par
    se <- numeric(0)
  }
  se1 <- rep(NA_real_, length(fit$par)) # the se values for the GARMA part of the model - default

  if (method == "WLL") {
    vcov1 <- matrix(nrow = length(fit$par) - 1, ncol = length(fit$par) - 1) # default
    # "-1" above since we don't want to include the se calc in the vcov matrix
  } else {
    vcov1 <- matrix(nrow = length(fit$par), ncol = length(fit$par))
  }
  if (method != "WLL" & !is.null(hh)) {
    # Next, find the se's for coefficients
    start <- 1
    se1 <- numeric(0) # default to set this up in the right environment
    # next check the hessian
    h_inv_diag <- diag(inv_hessian <- pracma::pinv(hh))
    if (method == "Whittle") {
      # Next line from GHR (2001) thm 3.2
      omega <- internal_whittle_garma_omega(fit$par, params)
      vcov1 <- pracma::pinv(omega) / length(y)
      se1 <- suppressWarnings(sqrt(diag(vcov1)))
    }
    if (method == "CSS") {
      se1 <- suppressWarnings(sqrt(h_inv_diag * sigma2 * 2))
      vcov1 <- inv_hessian * 2 * sigma2
    }
  }
  se <- c(se, se1)
  rm(se1)

  # set up names - nm

  nm <- character(0)
  if (!is.null(xreg)) {
    nm <- colnames(xreg)
  }

  if (is.null(periods)) {
    if (k > 0) nm <- c(nm, as.vector(unlist(lapply(1:k, function(x) paste0(c("u", "fd"), x)))))
  } else {
    nm <- c(nm, paste0("fd", 1:length(periods)))
  }
  if (p > 0) nm <- c(nm, paste0("ar", 1:p))
  if (q > 0) nm <- c(nm, paste0("ma", 1:q))
  if (method == "WLL") {
    nm <- c(nm, "Var")
  }

  coef <- matrix(data=c(coef, se), nrow = 2, byrow = TRUE)
  colnames(coef) <- nm
  rownames(coef) <- c("", "s.e.")
  result$coef <- coef
  result$var.coef <- vcov1
  colnames(vcov1) <- rownames(vcov1) <- tail(nm, nrow(vcov1))

  # build a ggbr_factors object
  gf <- list()
  if (is.null(xreg)) {
    start_idx <- 1
  } else {
    start_idx <- ncol(xreg) + 1
  }

  if (is.null(periods)) {
    for (k1 in seq_len(k)) {
      u <- coef[1, start_idx]
      gf1 <- list(freq = acos(u) / (2 * pi), period = (2 * pi) / acos(u), fd = coef[1, start_idx + 1])
      gf <- c(gf, list(gf1))
      start_idx <- start_idx + 2
    }
  } else {
    for (period in periods) {
      gf1 <- list(freq = 1/period, period = period, fd = coef[1, start_idx])
      gf <- c(gf, list(gf1))
      start_idx <- start_idx + 1
    }
  }
  class(gf) <- "ggbr_factors"

  model <- list(
    "phi" = coef[1, substr(colnames(coef), 1, 2) == "ar"],
    "theta" = coef[1, substr(colnames(coef), 1, 2) == "ma"],
    "ggbr_factors" = gf
  )
  if (!is.null(xreg)) {
    model$xreg = coef[1, 1:ncol(xreg)]
  }

  result$model <- model

  # fitted values and residuals
  result <- internal_fitted_values(result)

  # log lik
  result$loglik <-
    (-length(y) / 2 * log(2 * pi)
     - length(y) / 2 * log(sigma2)
     - 1 / (2 * sigma2) * sum(result$residuals^2))
  result$aic <- -2 * result$loglik + 2 * (p + q + ifelse(is.null(periods), 2*k, length(periods)) + 1)

  class(result) <- "garma_model"
  return(result)
}

# TO DO:
# 3. predict to include xreg factors
