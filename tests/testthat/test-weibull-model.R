test_that("weibull_model returns same result as in the book", {
  expected_coefs <- c(-0.27596492, 0.38073667) # page 97
  expected_beta_0_hat <- 0.759

  df <- hcv_be_2006[order(hcv_be_2006$dur), ]
  colnames(df) <- c("t", "status")

  model <- weibull_model(df)

  actual_coefs <- unname(c(
    coef(model$info)[1], # intercept
    coef(model$info)[2]
  ))
  actual_beta_0_hat <- exp(unname(coef(model$info)[1]))

  expect_equal(actual_coefs, expected_coefs, tolerance=0.000001)
  expect_equal(actual_beta_0_hat, expected_beta_0_hat, tolerance=0.001)

  # make sure helper works fine
  expect_no_error(compute_ci.weibull_model(model))
  expect_no_error(plot(model))
})
