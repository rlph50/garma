testthat::test_that("accuracy function works", {

  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)
  train_dap <- dap[1:131]
  test_dap <- dap[132:143]

  fit <- garma(train_dap, order = c(9, 0, 0), k = 1)
  fc <- forecast(fit, 12)
  result <- forecast::accuracy(fc, x=test_dap)

  # this accuracy() function has varying results depending on whether the
  # Intel MKL library is being used with R. To work around this we test only
  # specific elements of the 'result' matrix which are similar between the
  # two versions.
  testthat::expect_true(abs(result[1,1]) < 1e-4)
  testthat::expect_true(abs(result[2,1]) < 1e-2)
  testthat::expect_true(abs(round(result[1,2], 2)-0.0910) < 1e-2)
  testthat::expect_true(abs(round(result[2,2], 3)-0.0867) < 1e-3)
  testthat::expect_true(abs(round(result[1,3], 2)-0.0700) < 1e-2)
  testthat::expect_true(abs(round(result[2,3], 3)-0.0759) < 1e-3)
  testthat::expect_true(is.nan(result[1,4]))
  testthat::expect_true(abs(round(result[2,4], 0)-45) < 1.0)
  testthat::expect_true(is.infinite(result[1,5]))
  testthat::expect_true(abs(round(result[2,5], 0)-88) < 2.0)
  testthat::expect_true(abs(round(result[1,6], 2)-0.645) < 0.006)
  testthat::expect_true(abs(round(result[2,6], 2)-0.7) < 1e-8)
  testthat::expect_true((abs(round(result[1,7], 2))-0.03) < 1e-8)
  testthat::expect_true(is.na(result[2,7]))

  # For this test, the results are markedly different depending on whether
  # the Intel MKL library is being used. Here we check for that and set up
  # the expected results accordingly.
  # if (sessionInfo()$matprod == "default") {
  #   expected_result <- matrix(
  #     data = c(5.01801128694943e-05, 0.0911768449068468, 0.0703926772247195, NaN, Inf, 0.648959489052162, -0.0301943944340089,
  #              5.00625748787946e-03, 0.0867980763982614, 0.0759770295082921, 44.9578939522492, 89.1661901082397, 0.700442378288856, NA),
  #     byrow = TRUE,
  #     nrow = 2,
  #     dimnames = list(c("Training set", "Test set"),
  #                     c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1"))
  #   )
  # } else  {
  #   expected_result <- matrix(
  #     data = c(7.332368e-05, 0.09061627, 0.06983208, NaN, Inf, 0.6437913, -0.03172361,
  #              5.139126e-03, 0.08673775, 0.07553035, 44.74868, 87.70533, 0.6963244, NA),
  #     byrow = TRUE,
  #     nrow = 2,
  #     dimnames = list(c("Training set", "Test set"),
  #                     c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1"))
  #   )
  # }
  #
  # testthat::expect_equal(result, expected_result)
})
