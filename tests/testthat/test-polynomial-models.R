library(testthat)
library(stats4)

test_that("farrington_model returns same result as in the book", {
  expected <- c(
    alpha=0.07034904,
    beta=0.20243950,
    gamma=0.03665599
  )

  df <- rubella_uk_1986_1987
  model <- suppressWarnings(farrington_model(
      df$age, pos = df$pos, tot = df$tot,
      start=list(alpha=0.07,beta=0.1,gamma=0.03)
    ))
  actual <- c(
    model$info@coef[1],
    model$info@coef[2],
    model$info@coef[3]
    )

  expect_equal(actual, expected, tolerance=0.000001)

  # make sure utilities work
  expect_no_error(plot(model))
})

test_that("farrington_model works with linelisting data", {
  df <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age),]

  expect_no_error(
    suppressWarnings(farrington_model(
      df$age, status = df$seropositive,
      start=list(alpha=0.07,beta=0.1,gamma=0.03)
    ))
  )
})


test_that("polynomial_model works with line listing data", {

  data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age),]
  expect_no_error(polynomial_model(data$age, status = data$seropositive, type = "Muench"))


})

test_that("polynomial_model returns same result as in the book (Muench)", {
  expected <- c(-0.0505004)

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
    type= "Muench"
  )
  actual <- unname(c(
    coef(model$info)[1]
  ))

  expect_equal(actual, expected, tolerance=0.000001)

  # make sure utilities work
  expect_no_error(compute_ci(model))
  expect_no_error(plot(model))
})

test_that("polynomial_model returns same result as in the book (Muench)", {
  expected <- c(-0.0505004)

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
    k = 1
  )
  actual <- unname(c(
    coef(model$info)[1]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})

test_that("polynomial_model returns same result as in the book (Griffiths)", {
  expected <- c(-0.0442615740, -0.0001888796)

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
    type = "Griffith"
  )
  actual <- unname(c(
    coef(model$info)[1],
    coef(model$info)[2]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})



test_that("polynomial_model returns same result as in the book (Griffiths)", {
  expected <- c(-0.0442615740, -0.0001888796)

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
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

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
    type = "Grenfell"
  )
  actual <- unname(c(
    coef(model$info)[1],
    coef(model$info)[2],
    coef(model$info)[3]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})

test_that("polynomial_model returns same result as in the book (Grenfell & Anderson)", {
  expected <- c(-5.325918e-02, 5.065095e-04, -1.018736e-05)

  df <- hav_bg_1964
  model <- polynomial_model(
    df$age, pos = df$pos, tot = df$tot,
    k = 3
  )
  actual <- unname(c(
    coef(model$info)[1],
    coef(model$info)[2],
    coef(model$info)[3]
  ))

  expect_equal(actual, expected, tolerance=0.000001)
})
