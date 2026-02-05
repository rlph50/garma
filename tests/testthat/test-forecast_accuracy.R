testthat::test_that("accuracy function works", {

  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)
  train_dap <- dap[1:131]
  test_dap <- dap[132:143]

  fit <- garma(train_dap, order = c(9, 0, 0), k = 1)
  fc <- forecast(fit, 12)
  result <- forecast::accuracy(fc, x=test_dap)

  expected_result <- matrix(
    data = c(5.01801128694943e-05, 0.0911768449068468, 0.0703926772247195, NaN, Inf, 0.648959489052162, -0.0301943944340089,
             5.00625748787946e-03, 0.0867980763982614, 0.0759770295082921, 44.9578939522492, 89.1661901082397, 0.700442378288856, NA),
    byrow = TRUE,
    nrow = 2,
    dimnames = list(c("Training set", "Test set"),
                    c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1"))
  )

  testthat::expect_equal(result, expected_result)
})
