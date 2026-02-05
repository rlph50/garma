# Extracted from test-garma2.R:9

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "garma", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
data(AirPassengers)
ap <- log(AirPassengers)
dap <- diff(ap)
testthat::expect_snapshot({
    fit <- garma(dap, order = c(9, 0, 0), k = 1)
    print(fit)
  })
