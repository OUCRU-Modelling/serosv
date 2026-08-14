library(testthat)
library(stats4)

test_that("farrington_model returns same result as in the book", {
  expected <- c(
    alpha=0.07034904,
    beta=0.20243950,
    gamma=0.03665599
  )

  model <- suppressWarnings(farrington_model(
    rubella_uk_1986_1987,
    start=list(alpha=0.07,beta=0.1,gamma=0.03)
  ))
  actual <- c(
    model$info@coef[1],
    model$info@coef[2],
    model$info@coef[3]
  )

  expect_equal(actual, expected, tolerance=0.001)

  # make sure utilities work
  expect_no_error(plot(model))
})

test_that("farrington_model works with linelisting data", {
  df <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age),]
  df$status <- df$seropositive

  expect_no_error(
    suppressWarnings(farrington_model(
      df,
      start=list(alpha=0.07,beta=0.1,gamma=0.03)
    ))
  )
})

test_that("farrington_model utils function work", {
  df <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age),]

  model <- suppressWarnings(farrington_model(
    df,
    status_col = "seropositive",
    start=list(alpha=0.07,beta=0.1,gamma=0.03)
  ))

  # test plot function
  expect_no_error(plot(model, foi_ci=TRUE))

  # test print function
  expect_no_error(
    capture.output(print(model))
  )

  # test predict function
  expect_equal(
    predict(model, newdata = data.frame(
      age = model$df$age
    )),
    model$sp
  )
})


