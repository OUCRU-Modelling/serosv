library(testthat)

test_that("fp_model returns same result as in the book (Hepatitis A (BG))", {
  expected_coefs <- c(-1.09473686, 0.02622843, -0.01613128)
  expected_D <- 77.748963

  df <- hav_bg_1964
  model <- fp_model(
    df$age, pos = df$pos, tot = df$tot,
    p=c(1.9, 2.0), link="logit"
  )
  actual_coefs <- unname(c(
    coef(model$info)[1], # intercept
    coef(model$info)[2],
    coef(model$info)[3]
  ))
  actual_D <- model$info$deviance

  expect_equal(actual_coefs, expected_coefs)
  expect_equal(actual_D, expected_D)
})


test_that("find_best_fp_powers returns same result as in the book (non-monotone)", {
  expected_p <- c(1.9, 2.0)

  df <- hav_be_1993_1994
  output <- suppressWarnings(find_best_fp_powers(
    df$age, pos = df$pos, tot = df$tot,
    p=seq(-2,3,0.1), mc=F, degree=2, link="logit"
  ))

  actual_p <- output$p

  expect_equal(actual_p, expected_p)
})

test_that("find_best_fp_powers returns same result as in the book (monotone)", {
  expected_p <- c(1.0, 1.6)

  df <- hav_be_1993_1994
  output <- suppressWarnings(find_best_fp_powers(
    df$age, pos = df$pos, tot = df$tot,
    p=seq(-2,3,0.1), mc=T, degree=2, link="logit"
  ))

  actual_p <- output$p

  expect_equal(actual_p, expected_p)
})

