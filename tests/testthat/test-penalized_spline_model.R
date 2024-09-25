# TODO: add test
test_that("penalized_spline_model returns same result as in the book (ps framework)", {
  df <- vzv_parvo_be
  subset <-
    (df$age > 0.5) &
    (df$age < 76) &
    (!is.na(df$age)) & !is.na(df$parvo_res)
  df <- df[subset, ]

  expect_no_error(penalized_spline_model(df$age, status = df$parvo_res, s = "tp", framework = "pl"))
})

test_that("penalized_spline_model returns same result as in the book (glmm framework)", {
  expected_coef <- c(0.7080086, 1.4665966, 5.0837472, 0.2509096, -2.3617926, 0.1543315, -3.1849006,
               0.9776792, 8.0273085, 3.5284318)

  df <- vzv_parvo_be
  subset <-
    (df$age > 0.5) &
    (df$age < 76) &
    (!is.na(df$age)) & !is.na(df$parvo_res)
  df <- df[subset, ]
  model <- penalized_spline_model(df$age, status = df$parvo_res, s = "tp", framework = "glmm")
  actual_coef <- unname(model$info$gam$coefficients)

  expect_equal(expected_coef, actual_coef, tolerance=0.00001)
})

test_that("penalized_spline_model works with aggregated data", {
  df <- hav_bg_1964[order(hav_bg_1964$age), ]

  expect_no_error(suppressWarnings(penalized_spline_model(df$age, pos = df$pos, tot = df$tot)))
})
