library(testthat)
test_that("keiding_model returns expected results", {
  expected_foi_summary <- c(0.02330453, 0.05276397, 0.05504527, 0.08295355)
  expected_sp_summary <- c(0.06664893, 0.9155623, 0.7704096, 0.9892495)
  df <- hav_bg_1964
  model <- keiding_model(
    df$age,df$pos,df$tot, 
    kernel = "normal", bw = 30)
  actual_foi_summary <- c(
    min(model$foi), median(model$foi), mean(model$foi), max(model$foi)
  )
  actual_sp_summary  <- c(
    min(model$sp), median(model$sp), mean(model$sp), max(model$sp)
  )
  expect_equal(actual_foi_summary, expected_foi_summary, tolerance=0.000001)
  expect_equal(actual_sp_summary, expected_sp_summary, tolerance=0.000001)
})
