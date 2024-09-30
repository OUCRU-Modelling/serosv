library(testthat)

test_that("mixture_model returns expected result", {
  expected_params_dist1 <- c(0.1088, 2.349, 0.6804)
  expected_params_dist2 <- c(0.8912, 6.439, 0.9437)
  expected_se_dist1 <- c(0.006238, 0.04522, 0.03423)
  expected_se_dist2 <- c(0.006238, 0.02019, 0.01538)

  df <- vzv_be_2001_2003[vzv_be_2001_2003$age < 40.5,]
  data <- df$VZVmIUml[order(df$age)]
  model <- mixture_model(data)

  actual_params_dist1 <- as.double(model$info$parameters[1,])
  actual_params_dist2 <- as.double(model$info$parameters[2,])
  actual_se_dist1 <- as.double(model$info$se[1,])
  actual_se_dist2 <- as.double(model$info$se[2,])

  expect_equal(expected_params_dist1, actual_params_dist1, tolerance = 0.0001)
  expect_equal(expected_params_dist2, actual_params_dist2, tolerance = 0.0001)
  expect_equal(expected_se_dist1, actual_se_dist1, tolerance = 0.0001)
  expect_equal(expected_se_dist2, actual_se_dist2, tolerance = 0.001)
})

# TODO: update expected return value
test_that("estimate_from_mixture returns expected result", {
  df <- vzv_be_2001_2003[vzv_be_2001_2003$age < 40.5,]
  data <- df$VZVmIUml[order(df$age)]
  model <- mixture_model(data)

  # est_mixture <- estimate_from_mixture(df$age[order(df$age)], data, mixture_model = model)
  expect_no_error(
    estimate_from_mixture(df$age[order(df$age)], data, mixture_model = model)
  )

})
