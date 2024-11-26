test_that("no error while using hierarchical bayesian model", {
  df <- mumps_uk_1986_1987

  # making sure all models work without error
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="far2")
    )
  )
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="far3")
    )
  )
  expect_no_error(
    suppressWarnings(
      hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="log_logistic")
    )
  )

})
