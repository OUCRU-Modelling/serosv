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
