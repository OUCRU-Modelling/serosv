test_that("no error while using hierarchical bayesian model", {
  df <- mumps_uk_1986_1987

  # making sure all models work without error
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(df, type="far2", warmup = 500, iter=2000)
    )
  )
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(df, type="far3", warmup = 500, iter=2000)
    )
  )
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(df, type="log_logistic", warmup = 500, iter=2000)
    )
  )

})

test_that("test utility functions for hierarchical bayesian model", {
  model <- hierarchical_bayesian_model(mumps_uk_1986_1987, type="far2", warmup = 500, iter=2000)

  pred <- predict(model, data.frame(age = mumps_uk_1986_1987$age))
  expect_equal(
    model$sp,
    pred$y,
    tolerance = 0.01
  )
  expect_no_error(print(model))
  expect_no_error(compute_ci(model))
  expect_no_error(plot(model, foi_ci=TRUE))
})
