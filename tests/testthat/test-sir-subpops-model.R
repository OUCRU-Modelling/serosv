test_that("sir_subpops_model returns expected results", {
  expected <- list(
    time=10000,
    s1=0.6869925,
    s2=0.6869925,
    i1=0.009141676,
    i2=0.009141676,
    r1=0.3038659,
    r2=0.3038659
  )

  k <- 2
  state <- c(
    s = c(0.8, 0.8),
    i = c(0.2, 0.2),
    r = c(  0,   0)
  )
  beta_matrix <- c(
    c(0.05, 0.00),
    c(0.00, 0.05)
  )
  parameters <- list(
    beta = matrix(beta_matrix, nrow=k, ncol=k, byrow=TRUE),
    nu = c(1/30, 1/30),
    mu = 0.001,
    k = k
  )
  times<-seq(0,10000,by=10)

  output <- sir_subpops_model(times, state, parameters)$output
  actual <- as.list(tail(output, 1))

  expect_equal(actual, expected, tolerance = 0.000001)
  # make sure utilities work
  expect_no_error(plot(output))
})
