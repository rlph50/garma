context("GARMA")
library(garma)

data(AirPassengers)
ap <- log(AirPassengers)
dap<-diff(ap)

# Test coefficients
test_that("Test 1-COEF. Short Memory AR model coef match 'arima'", {
  expect_true(
    # RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
    sqrt(mean((coef(garma(dap,order=c(2,0,0),k=0,method='CSS',include.mean=F)) -
                 coef(arima(dap,order=c(2,0,0),include.mean=F,method='CSS')))^2)) < 0.005
  )
})
test_that("Test 2-COEF. Short Memory MA model coef match 'arima'", {
  expect_true(
    # RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
    sqrt(mean((coef(garma(dap,order=c(0,0,2),k=0,method='CSS',include.mean=F)) -
                 coef(arima(dap,order=c(0,0,2),include.mean=F,method='CSS')))^2)) < 0.005
  )
})

# calc next bits in advance, and re-arrange for the test.
g_c <- coef(garma(dap,order=c(2,0,2),k=0,method='CSS',include.mean=T))
g_c2 <- c(g_c[2:5],g_c[1])
a_c <- coef(arima(dap,order=c(2,0,2),include.mean=T,method='CSS'))

# RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
test_that("Test 3-COEF. Short Memory ARMA model coef match 'arima'", {
  expect_true(sqrt(mean( (g_c2-a_c)^2)) < 0.05)
})


# Test residuals
test_that("Test 1-RESID. Short Memory AR model resid match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((residuals(garma(dap,order=c(2,0,0),k=0,method='CSS',include.mean=F)) -
                 residuals(arima(dap,order=c(2,0,0),method='CSS',include.mean=F)))^2)) < 0.001
  )
})

test_that("Test 2-RESID. Short Memory AR model resid match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap,order=c(2,1,0),k=0,method='CSS',include.mean=F)) -
                 residuals(arima(dap,order=c(2,1,0),method='CSS',include.mean=F)))^2)) < 0.08
  )
})

test_that("Test 3-RESID. Short Memory MA model resid match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap,order=c(0,0,2),k=0,method='CSS',include.mean=F)) -
                 residuals(arima(dap,order=c(0,0,2),method='CSS',include.mean=F)))^2)) < 0.15
  )
})

test_that("Test 4-RESID. Short Memory ARMA model resid match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap,order=c(2,0,2),k=0,method='CSS',include.mean=F)) -
                 residuals(arima(dap,order=c(2,0,2),method='CSS',include.mean=F)))^2)) < 0.07
  )
})

test_that("Test 5-RESID. Short Memory ARMA model with intercept resid match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((residuals(garma(dap,order=c(2,0,2),k=0,method='CSS',include.mean=T)) -
                 residuals(arima(dap,order=c(2,0,2),method='CSS',include.mean=T)))^2)) < 0.15
  )
})

# Test predictions
test_that("Test 1-PRED. Short Memory AR model. Check that pred match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(dap,order=c(2,0,0),k=0,method='CSS',include.mean=F),n.ahead=12)$pred -
                 predict(arima(dap,order=c(2,0,0),method='CSS',include.mean=F),n.ahead=12)$pred)^2)) < 0.001
  )
})

test_that("Test 2-PRED. Short Memory MA model. Check that pred match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(dap,order=c(0,0,2),k=0,method='CSS',include.mean=F),n.ahead=12)$pred -
                 predict(arima(dap,order=c(0,0,2),method='CSS',include.mean=F),n.ahead=12)$pred)^2)) < 0.03
  )
})

test_that("Test 3-PRED. Short Memory AR model with diff. Check that pred match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small
    sqrt(mean((predict(garma(ap,order=c(2,1,0),k=0,method='CSS',include.mean=F,include.drift=FALSE),n.ahead=12)$pred -
                 predict(arima(ap,order=c(2,1,0),method='CSS',include.mean=F),n.ahead=12)$pred)^2)) < 0.001
  )
})

test_that("Test 4-PRED. Short Memory ARIMA model with diff. Check that pred match 'arima'", {
  expect_true(
    # RMSE for difference in residuals between GARMA and ARIMA is reasonably small with differencing.
    sqrt(mean((predict(garma(ap,order=c(2,1,2),k=0,method='CSS',include.mean=F,include.drift=FALSE),n.ahead=6)$pred -
                 predict(arima(ap,order=c(2,1,2),method='CSS',include.mean=F),n.ahead=6)$pred)^2)) < 0.3
  )
})
