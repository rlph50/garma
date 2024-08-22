testthat::test_that("Check Output", {
  data(AirPassengers)
  ap <- log(AirPassengers)
  dap <- diff(ap)

  testthat::expect_snapshot({
    fit <- garma(dap, order = c(9, 0, 0), k = 1)
    print(fit)
  })
})
