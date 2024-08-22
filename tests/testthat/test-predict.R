testthat::test_that("predict", {
  fit <- garma(runif(120))
  testthat::expect_error(
    predict(fit, n.ahead = (-1)),
    regexp = "The parameter 'n.ahead' must be an integer g.t. 0"
  )

  testthat::expect_error(
    predict(fit, n.ahead = 5, newdata = runif(5)),
    regexp = "You have provided the 'newdata' parameter but the original 'fit' did not include an xreg"
  )

  fit <- garma(runif(120), xreg = runif(120))
  testthat::expect_error(
    predict(fit, n.ahead = 5, newdata = runif(4)),
    regexp = "The length of newdata must be 5"
  )

})
