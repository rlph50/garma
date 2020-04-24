context("GARMA")
library(garma)

data(AirPassengers)
ap  <- as.numeric(diff(AirPassengers,12))

test_that("Short Memory model coef match 'arima'", {
  expect_true(
    # RMSE for difference in coefficients between GARMA and ARIMA is reasonably small
    sqrt(mean((garma(ap,order=c(9,1,0),k=0,method='Whittle',include.mean=F)$coef[1,] - arima(ap,order=c(9,1,0))$coef)^2)) < 0.05
  )
})
