internal_ggbr_coef <- function(n, d, u) {
  if (n > 0) {
    cf <- c(1, 2 * d * u, 2 * d * (d + 1) * (u^2) - d)
    if (n > 4) for (j in 4:n) cf[j] <- (2 * u * (j + d - 2) * cf[j - 1] - (j - 1 + 2 * d - 2) * cf[j - 2]) / (j - 1)
    cf <- cf[1:n]
  } else {
    cf <- NA
  }
  return(cf)
}

internal_a_fcn <- function(a_vec, freq) {
  # this is a utility function used to find the short-memory spectral density.
  n_freq <- length(freq)
  cos_sum <- rep(0.0, n_freq)
  a_len <- length(a_vec)
  for (i in 1:a_len) cos_sum <- cos_sum + a_vec[i] * cos(2 * i * pi * freq)
  sin_sum <- rep(0.0, n_freq)
  for (i in 1:a_len) sin_sum <- sin_sum + a_vec[i] * sin(2 * i * pi * freq)

  return((1 + cos_sum)^2 + sin_sum^2)
}

internal_garma_pgram <- function(x) {
  x <- as.numeric(x) # If we don't do this, the spectrum function returns 'different' freq.
  return(spectrum(x, plot = FALSE, detrend = TRUE, demean = TRUE, method = "pgram", taper = 0, fast = FALSE))
}

internal_yajima_ggbr_freq <- function(x, k=1, pgram = NULL, bandwidth = 1L) {
  x <- as.numeric(x) # If we don't do this, the spectrum function returns 'different' frequencies.
  if (is.null(pgram)) {
    pgram <- internal_garma_pgram(x)
  }

  spec <- pgram$spec # create a copy
  peak_indexes <- integer(0)
  peak_freq <- numeric(0)
  for (k1 in 1:k) {
    idx <- which.max(spec)
    peak_indexes <- c(peak_indexes, idx)
    peak_freq <- c(peak_freq, pgram$freq[idx])
    # Now zero out the spec around this point so we can find another peak.
    idxm1 <- idx - bandwidth
    idxp1 <- idx + bandwidth
    if (idxm1 < 1) idxm1 <- 1
    if (idxp1 > length(spec)) idxp1 <- length(spec)
    spec[idxm1:idxp1] <- 0
  }

  return(list(peak_indexes = peak_indexes, peak_freq = peak_freq, peak_periods = 1/peak_freq))
}

internal_semipara <- function(x, periods = NULL, alpha = 0.8, k = 1, method = "gsp", pgram = NULL) {
  # as per Arteche 1998 "SEMIPARAMETRIC INFERENCE IN SEASONAL AND CYCLICAL LONG MEMORY PROCESSES"
  # If we need to estimate the cycles in the data then (a) periods will be NULL and (b) k>=1
  # Then for each peak we will estimate the fd value.

  x <- as.numeric(x) # If we don't do this, the spectrum function returns 'different' frequencies.

  if (is.null(pgram)) {
    pgram <- internal_garma_pgram(x)
  }

  yf <- NULL
  if (is.null(periods)) {
    # we will have to estimate the positions of the periods
    yf <- internal_yajima_ggbr_freq(x, k = k, pgram = pgram)
    periods <- yf$peak_periods
  } else {
    k <- length(periods)
  }

  # Now estimate the fd for each period
  if (method == "gsp") fd_list <- internal_gsp_estimate(x, alpha, periods, pgram = pgram)
  if (method == "lpr") fd_list <- internal_lpr_estimate(x, alpha, periods, pgram = pgram)

  result <- list()
  for (k1 in 1:k) {
    res <- list(
      freq = ifelse(is.null(yf), 1/periods[k1], yf$peak_freq[k1]),
      periods = periods[k1],
      fd = fd_list[k1]
    )
    result <- c(result, list(res))
  }
  class(result) <- "ggbr_factors"

  return(result)
}

# For a given vector of "periods", find the gsp (Gaussian Semi-Parametric) estimate of the
# degree of fractional differencing
internal_gsp_estimate <- function(x, alpha, periods, pgram = NULL) {
  # as per Arteche 1998 "SEMIPARAMETRIC INFERENCE IN SEASONAL AND CYCLICAL LONG MEMORY PROCESSES"
  # determine "fd" for a set of fixed periods.
  c_fcn <- function(fd, omega, spec) {
    return(mean((omega^(2 * fd)) * spec, na.rm = TRUE))
  }
  r1_fcn <- function(fd, f_idx, pgram, m) {
    omega <- 2 * pi * pgram$freq[1:m] # Frequencies

    # Spec to use, as offset from ggbr_freq. These are specs above ggbr_freq.
    spec_2pi <- c(pgram$spec, rev(pgram$spec))
    spec1 <- spec_2pi[f_idx:(f_idx + m)]
    spec1 <- spec1[1:m]

    res <- log(c_fcn(fd, omega, spec1)) - 2 * fd * mean(log(omega), na.rm = TRUE)
    if (is.infinite(res) | is.na(res) | is.null(res)) res <- 1e200
    return(res)
  }

  x <- as.numeric(x) # If we don't do this, the spectrum function returns frequencies based on the ts() settings.
  m <- as.integer((length(x) / 2L)^alpha)
  if (is.null(pgram)) {
    pgram <- internal_garma_pgram(x)
  }

  fd_list <- numeric(0)
  for (period in periods) {
    f <- 1/period
    idx <- which.min(abs(pgram$freq - f)) # Find closest index to period.
    fd <- stats::optimise(r1_fcn, f_idx = idx, pgram = pgram, m = m, lower = -10, upper = 10)$minimum
    fd_list <- c(fd_list, fd)
  }

  return(fd_list)
}

# For a given vector of "periods", find the lpr (Least-squares Periodogram Regression) estimate of the
# degree of fractional differencing
internal_lpr_estimate <- function(x, alpha, periods, pgram = NULL) {
  # as per Arteche 1998 "SEMIPARAMETRIC INFERENCE IN SEASONAL AND CYCLICAL LONG MEMORY PROCESSES"
  # determine "fd" for a set of fixed periods.
  x <- as.numeric(x) # If we don't do this, the spectrum function returns frequencies based on the ts() settings.
  m <- as.integer((length(x) / 2L)^alpha)
  if (is.null(pgram)) {
    pgram <- internal_garma_pgram(x)
  }

  v <- log(1:(m - 1)) - mean(log(1:(m - 1)))
  denom <- 4 * sum(v^2)

  # we want to be able to easily grab the Spectrum for negative indicies if needed, so
  # here we form a larger vector of spectrum values.
  all_spec <- c(rev(pgram$spec[2:length(pgram$spec)]), pgram$spec)
  start_index <- length(pgram$spec)

  fd_list <- numeric(0)
  for (period in periods) {
    f <- 1/period
    idx <- which.min(abs(pgram$freq - f)) # Find closest index to period.
    spec1 <- all_spec[(start_index + idx):(start_index + idx + m - 2)]
    spec2 <- all_spec[(start_index + idx):(start_index + idx - m + 2)]

    numer <- sum(v * (log(spec1) + log(spec2)))
    fd <- (-0.5) * numer / denom

    fd_list <- c(fd_list, fd)
  }

  return(fd_list)
}


.getPackageVersion <- function(pkgname) {
  return(paste0(crayon::black("\n\nPackage "), crayon::blue(pkgname), ": ", crayon::black(utils::packageVersion("garma")), "\n"))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(.getPackageVersion(pkgname))
}

# next function allows us to update the object
.byRef <- function(..., envir = parent.frame(), inherits = TRUE) {
  cl <- match.call(expand.dots = TRUE)
  cl[c(1, match(c("envir", "inherits"), names(cl), 0L))] <- NULL
  for (x in as.list(cl)) {
    s <- substitute(x)
    sx <- do.call(substitute, list(s), envir = envir)
    dx <- deparse(sx)
    expr <- substitute(assign(dx, s, envir = parent.frame(), inherits = inherits))
    do.call(on.exit, list(expr, add = TRUE), envir = envir)
  }
}

.printf <- function(val) {
  if (class(val)[1] == "integer") {
    fmtstr <- "%s: %d\n"
  } else {
    fmtstr <- "%s: %f\n"
  }
  cat(sprintf(fmtstr, as.character(substitute(val)), val))
}

internal_make_column_name <- function(name) {
  # if we have a vector supplied as xreg= or newdata= then we can use
  # deparse(substitute) to get the name for the column when we put it in a  matrix or dataframe.
  # but if this came from a dataframe Eg deparse(substitute()) is 'mydf$col'
  # so we remove everything before the '$' (why? maybe the user has a df called train and another called test...)

  name <- gsub(".*\\$", "", name)

  return(name)
}
