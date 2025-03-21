test_that("penalized_spline_model returns correct format", {
  df <- vzv_parvo_be
  subset <- (df$age > 0.5) &
    (df$age < 76) &
    (!is.na(df$age)) & !is.na(df$parvo_res)
  df <- df[subset, c("age", "parvo_res")]
  colnames(df) <- c("age", "status")

  glmm_model <- penalized_spline_model(df, s = "tp", framework = "glmm")
  pl_model <- penalized_spline_model(df, s = "tp", framework = "pl")

  expect_equal(class(glmm_model$info)[1], "gamm")
  expect_equal(class(pl_model$info)[1], "gam")
})

test_that("penalized_spline_model works without error (pl framework)", {
  df <- vzv_parvo_be
  subset <-
    (df$age > 0.5) &
    (df$age < 76) &
    (!is.na(df$age)) & !is.na(df$parvo_res)
  df <- df[subset, c("age", "parvo_res")]
  colnames(df) <- c("age", "status")

  expect_no_error(penalized_spline_model(df, s = "tp", framework = "pl"))
  model <- penalized_spline_model(df, s = "tp", framework = "pl")

  # making sure compute ci and plot function works without error
  expect_no_error(compute_ci.penalized_spline_model(model))
  expect_no_error(plot(model))
})

test_that("penalized_spline_model works without error (glmm framework)", {
  # expected_coef <- c(0.7080086, 1.4665966, 5.0837472, 0.2509096, -2.3617926, 0.1543315, -3.1849006,
  #              0.9776792, 8.0273085, 3.5284318)

  df <- vzv_parvo_be
  subset <-
    (df$age > 0.5) &
    (df$age < 76) &
    (!is.na(df$age)) & !is.na(df$parvo_res)
  df <- df[subset, c("age", "parvo_res")]
  colnames(df) <- c("age", "status")

  model <- penalized_spline_model(df, s = "tp", framework = "glmm")
  actual_coef <- unname(model$info$gam$coefficients)

  # expect_equal(actual_coef,expected_coef, tolerance=0.00001)

  # making sure compute ci and plot function works without error
  expect_no_error(compute_ci.penalized_spline_model(model))
  expect_no_error(plot(model))
})

test_that("penalized_spline_model works with aggregated data", {
  df <- hav_bg_1964[order(hav_bg_1964$age), ]

  expect_no_error(suppressWarnings(penalized_spline_model(df)))
})
