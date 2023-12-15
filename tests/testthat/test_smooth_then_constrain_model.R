library(testthat)
test_that("smooth_then_constrain_model returns expected results", {
  expected_foi_summary <- c(0, 0.05875577, 0.09195305, 0.9057742)
  expected_sp_summary <- c(0.1736335, 0.7792654, 0.6995438, 0.9999992)
  df <- hav_bg_1964
  model <- stc(df$age,df$pos,df$tot,
               alpha = 0.35,family = "binomial")
  actual_foi_summary <- c(
    min(model$foi), median(model$foi), mean(model$foi), max(model$foi)
  )
  actual_sp_summary  <- c(
    min(model$sp), median(model$sp), mean(model$sp), max(model$sp)
  )
  expect_equal(actual_foi_summary, expected_foi_summary, tolerance=0.000001)
  expect_equal(actual_sp_summary, expected_sp_summary, tolerance=0.000001)
})
