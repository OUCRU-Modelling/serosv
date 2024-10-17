library(testthat)

test_that("sir_basic_model returns expected results", {
  S_equilibrium = 2026.67
  I_equilibrium =   39.12
  state <- c(S=4999, I=1, R=0)
  parameters <- c(
    mu=1/75,
    alpha=0,
    beta=0.0005,
    nu=1,
    p=0
  )
  times <- seq(0, 1000, by=0.1)

  output <- sir_basic_model(times, state, parameters)$output

  expect_equal(round(tail(output, 1)$S, 2), S_equilibrium)
  expect_equal(round(tail(output, 1)$I, 2), I_equilibrium)

  # make sure utilities work
  expect_no_error(plot(output))
})
