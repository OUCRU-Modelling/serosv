test_that("bivariate dale model returns the same result as original code", {
  # --- Expected values obtained using the code for the paper
  expected <- c(
    # coefficients
    -0.79898001, -0.28972390, 0.47942086, 0.13138027, 0.15787313, 0.01700757,
    # deviance and degrees of freedom
    131.2903, 117.3236
  )

  # --- Get result from package
  data <- rubella_mumps_uk
  model <- bivariate_dale_model(age = data$age, y = data[, c("NN", "NP", "PN", "PP")], monotonized=TRUE)

  actual <- unname(coef(model$info))
  actual <- append(actual, deviance(model$info))
  actual <- append(actual, df.residual(model$info))

  expect_equal(expected, actual, tolerance = 0.00001)
})
