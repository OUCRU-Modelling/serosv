library(testthat)
library(stats4)

test_that("polynomial_model works with line listing data", {

  data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age),]
  expect_no_error(polynomial_model(data, k=1,
                                   status_col = "seropositive"))

})

test_that("polynomial_model returns same result as in the book (Muench)", {
  expected <- c(-0.0505004)

  model <- polynomial_model(
    hav_bg_1964,
    k = 1
  )
  actual <- unname(c(
    coef(model$info)[1]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})

test_that("polynomial_model returns same result as in the book (Griffiths)", {
  expected <- c(-0.0442615740, -0.0001888796)

  model <- polynomial_model(
    hav_bg_1964,
    k = 2
  )
  actual <- unname(c(
    coef(model$info)[1],
    coef(model$info)[2]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})


test_that("polynomial_model returns same result as in the book (Grenfell & Anderson)", {
  expected <- c(-5.325918e-02, 5.065095e-04, -1.018736e-05)

  model <- polynomial_model(
    hav_bg_1964,
    k = 3
  )
  actual <- unname(c(
    coef(model$info)[1],
    coef(model$info)[2],
    coef(model$info)[3]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})

test_that("test utility functions for polynomial_model", {
  model <- polynomial_model(
    hav_bg_1964,
    k = 3
  )

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
