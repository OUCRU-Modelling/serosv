library(testthat)
library(dplyr)

test_that("transform returns a data frame", {
  df <- transform(c(1, 2, 3), c(1, 0, 1))
  expect_type(df, "list")
})

test_that("transform returns the correct number of rows", {
  df <- transform(c(1, 2, 3), c(1, 0, 1))
  expect_equal(nrow(df), 3)
})

test_that("transform returns the correct column names", {
  df <- transform(c(1, 2, 3), c(1, 0, 1))
  expect_equal(colnames(df), c("t", "pos", "tot"))
})

test_that("transform returns the correct values", {
  df <- transform(c(1, 1, 2, 2, 3), c(1, 0, 1, 1, 0))
  expect_equal(df$pos, c(1, 2, 0))
  expect_equal(df$tot, c(2, 2, 1))
})
