test_that("no error while using hierarchical bayesian model", {
  df <- mumps_uk_1986_1987

  # making sure all types work
  expect_no_error(hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="far2"))
  expect_no_error(hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="far3"))
  expect_no_error(hierarchical_bayesian_model(age = df$age, pos = df$pos, tot = df$tot, type="log_logistic"))
})
