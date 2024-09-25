library(testthat)

test_that("transform_data returns a data frame", {
  df <- transform_data(c(1, 2, 3), c(1, 0, 1))
  expect_type(df, "list")
})

test_that("transform_data returns the correct number of rows", {
  df <- transform_data(c(1, 2, 3), c(1, 0, 1))
  expect_equal(nrow(df), 3)
})

test_that("transform_data returns the correct column names", {
  df <- transform_data(c(1, 2, 3), c(1, 0, 1))
  expect_equal(colnames(df), c("t", "pos", "tot"))
})

test_that("transform_data returns the correct values", {
  df <- transform_data(c(1, 1, 2, 2, 3), c(1, 0, 1, 1, 0))
  expect_equal(df$pos, c(1, 2, 0))
  expect_equal(df$tot, c(2, 2, 1))
})

test_that("est_foi returns a numeric vector", {
  t <- c(1, 2, 3, 4)
  sp <- c(0.1, 0.2, 0.3, 0.4)
  expect_type(est_foi(t, sp), "double")
})

test_that("est_foi returns the correct output", {
  t <- c(1, 2, 3, 4)
  sp <- c(0.1, 0.2, 0.3, 0.4)
  expect_equal(est_foi(t, sp), c(0.125, 0.14285714))
})
