# Test coefficients
testthat::test_that("Test 1-COEF. Short Memory AR model coef match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
    sqrt(mean((coef(garma(dap, order = c(2, 0, 0), k = 0, method = "CSS", include.mean = F)) -
      coef(arima(dap, order = c(2, 0, 0), include.mean = F, method = "CSS")))^2)) < 0.005
  )
})

testthat::test_that("Test 2-COEF. Short Memory MA model coef match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
    sqrt(mean((coef(garma(dap, order = c(0, 0, 2), k = 0, method = "CSS", include.mean = F)) -
      coef(arima(dap, order = c(0, 0, 2), include.mean = F, method = "CSS")))^2)) < 0.005
  )
})


# RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
testthat::test_that("Test 3-COEF. Short Memory ARMA model coef match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  g_c <- coef(garma(dap, order = c(2, 0, 2), k = 0, method = "CSS", include.mean = T))
  g_c2 <- c(g_c[2:5], g_c[1])
  a_c <- coef(arima(dap, order = c(2, 0, 2), include.mean = T, method = "CSS"))

  testthat::expect_true(sqrt(mean((g_c2 - a_c)^2)) < 0.05)
})


# Test residuals
testthat::test_that("Test 1-RESID. Short Memory AR model resid match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((residuals(garma(dap, order = c(2, 0, 0), k = 0, method = "CSS", include.mean = F)) -
      residuals(arima(dap, order = c(2, 0, 0), method = "CSS", include.mean = F)))^2)) < 0.03
  )
})

testthat::test_that("Test 2-RESID. Short Memory AR model resid match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap, order = c(2, 1, 0), k = 0, method = "CSS", include.mean = F)) -
      residuals(arima(dap, order = c(2, 1, 0), method = "CSS", include.mean = F)))^2)) < 0.10
  )
})

testthat::test_that("Test 3-RESID. Short Memory MA model resid match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap, order = c(0, 0, 2), k = 0, method = "CSS", include.mean = F)) -
      residuals(arima(dap, order = c(0, 0, 2), method = "CSS", include.mean = F)))^2)) < 0.05
  )
})

testthat::test_that("Test 4-RESID. Short Memory ARMA model resid match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  gmdl <- garma(dap, order = c(2, 0, 2), k = 0, method = "Whittle", include.mean = F)
  amdl <- arima(dap, order = c(2, 0, 2), method = "CSS", include.mean = F)
  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((resid(gmdl) - resid(amdl))^2)) < 0.05
  )
})

testthat::test_that("Test 5-RESID. Short Memory ARMA model with intercept resid match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap, order = c(2, 0, 2), k = 0, method = "CSS", include.mean = T)) -
      residuals(arima(dap, order = c(2, 0, 2), method = "CSS", include.mean = T)))^2)) < 0.07
  )
})

# Test predictions
testthat::test_that("Test 1-PRED. Short Memory AR model. Check that pred match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(dap, order = c(2, 0, 0), k = 0, method = "CSS", include.mean = F), n.ahead = 12)$pred -
      predict(arima(dap, order = c(2, 0, 0), method = "CSS", include.mean = F), n.ahead = 12)$pred)^2)) < 0.001
  )
})

testthat::test_that("Test 2-PRED. Short Memory MA model. Check that pred match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(dap, order = c(0, 0, 2), k = 0, method = "CSS", include.mean = F), n.ahead = 12)$pred -
      predict(arima(dap, order = c(0, 0, 2), method = "CSS", include.mean = F), n.ahead = 12)$pred)^2)) < 0.03
  )
})

testthat::test_that("Test 3-PRED. Short Memory AR model with diff. Check that pred match 'arima'", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(ap, order = c(2, 1, 0), k = 0, method = "CSS", include.mean = F, include.drift = FALSE), n.ahead = 12)$pred -
      predict(arima(ap, order = c(2, 1, 0), method = "CSS", include.mean = F), n.ahead = 12)$pred)^2)) < 0.001
  )
})

testthat::test_that("Parameter checks on garma() function", {
  df <- data.frame(x=runif(120), y =runif(120))
  testthat::expect_error(
    garma(df),
    regexp = "x should be a numeric vector - not an entire data frame. Please select a single column and try again."
  )

  x <- runif(120)
  x[30] <- NA_real_
  testthat::expect_error(
    garma(x),
    regexp = "x should not have any missing values"
  )

  x <- runif(120)
  testthat::expect_error(
    garma(x, k=(-1)),
    regexp = "The k parameter must be a non-negative integer"
  )

  testthat::expect_error(
    garma(x, periods = "A"),
    regexp = "The 'periods' parameter must be a numeric vector of at least 1 element"
  )

  testthat::expect_error(
    garma(x, periods = numeric(0)),
    regexp = "The 'periods' parameter must be a numeric vector of at least 1 element"
  )

  testthat::expect_error(
    garma(x, periods = (-5)),
    regexp = "The 'periods' parameter cannot contain negative values"
  )

  testthat::expect_error(
    garma(x, order = c(1, 1)),
    regexp = "The 'order' parameter must be a 3 integers only"
  )

  testthat::expect_error(
    garma(x, order = c(-1, -1, -1)),
    regexp = "The 'order' parameter must consist of positive integers"
  )

  testthat::expect_error(
    garma(x, k = 0),
    regexp = "At least one of p, q or k \\(or periods\\) must be positive"
  )

  testthat::expect_error(
    garma(x, xreg = df),
    regexp = "The parameter 'xreg' should be a numeric vector or matrix"
  )

  testthat::expect_error(
    garma(x, xreg = runif(150)),
    regexp = "The parameter 'xreg' should be the same length as 'x'"
  )

  testthat::expect_error(
    garma(x, xreg = c(runif(119), NA_real_)),
    regexp = "The parameter 'xreg' should not have any NA values"
  )

  testthat::expect_error(
    garma(x, method = "invalid"),
    regexp = "The parameter 'method' must be one of CSS, Whittle or WLL"
  )

  testthat::expect_error(
    garma(x, opt_method = "invalid"),
    regexp = "The following optimisation routines are not supported: invalid"
  )

  testthat::expect_error(
    garma(x, d_lim = 0),
    regexp = "The parameter 'd_lim' should be a list of 2 numerics, Eg c\\(0,0.5\\) and the minimum should be < maximum"
  )

  testthat::expect_error(
    garma(x, d_lim = c(0.4, 0.1)),
    regexp = "The parameter 'd_lim' should be a list of 2 numerics, Eg c\\(0,0.5\\) and the minimum should be < maximum"
  )

  testthat::expect_error(
    garma(x, d_lim = c(NA, 1)),
    regexp = "The parameter 'd_lim' should not have any NA values"
  )

})

testthat::test_that("garma xreg", {
  # Set up some synthetic data and compare the results

  set.seed(31415926L)
  x <- runif(200)
  m <- matrix(runif(400), ncol = 2)
  m.ahead <- matrix(runif(20), ncol = 2)
  colnames(m) <- colnames(m.ahead) <- c("R1", "R2")

  fit <- garma(x, xreg = m)
  testthat::expect_snapshot({
    print(fit)
    predict(fit, n.ahead = 10, newdata = m.ahead)
  })
})
