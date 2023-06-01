library(testthat)
library(stats4)

test_that("farrington_model returns same result as in the book", {
  expected <- c(
    alpha=0.07034904,
    beta=0.20243950,
    gamma=0.03665599
  )

  df <- rubella_uk_1986_1987
  model <- farrington_model(
      df$age, df$pos, df$tot,
      start=list(alpha=0.07,beta=0.1,gamma=0.03)
    )
  actual <- c(
    coef(model$info)[1],
    coef(model$info)[2],
    coef(model$info)[3]
    )

  expected == actual
  expect_equal(actual, expected, tolerance=0.000001)
})
