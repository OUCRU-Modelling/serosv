library(testthat)
library(locfit)

test_that("lp_model works with linelisting data", {
  df <- hav_be_2002[order(hav_be_2002$age),]

  expect_no_error(lp_model(df$age, status = df$seropositive))
})

test_that("lp_model returns expected results", {
  expected_foi_summary <- c(0.0018361, 0.0856287, 0.107952, 0.294777)
  expected_sp_summary <- c(0.165227, 0.968281, 0.884160, 0.985884)

  df <- mumps_uk_1986_1987
  model <- lp_model(
      df$age, pos = df$pos, tot = df$tot,
      nn=0.7, kern="tcub"
    )

  actual_foi_summary <- c(
    min(model$foi), median(model$foi), mean(model$foi), max(model$foi)
    )
  actual_sp_summary  <- c(
    min(model$sp), median(model$sp), mean(model$sp), max(model$sp)
  )

  expect_equal(actual_foi_summary, expected_foi_summary, tolerance=0.000001)
  expect_equal(actual_sp_summary, expected_sp_summary, tolerance=0.000001)

  # make sure utilities work
  expect_no_error(compute_ci.lp_model(model))
  expect_no_error(plot(model))
})
