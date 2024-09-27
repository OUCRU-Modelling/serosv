library(testthat)

test_that("sir_static_model returns expected results", {
  expected = list(
    time=40,
    s=0.01829732,
    i=5.026737e-05,
    r=0.9816524)
  state <- c(s=0.999,i=0.001,r=0)
  parameters <- c(lambda = 0.1, nu=36.5)
  times <- seq(0, 40, by=0.01)

  output <- sir_static_model(times, state, parameters)$output
  actual <- as.list(tail(output, 1))

  expect_equal(actual, expected, tolerance = 0.000001)
})
