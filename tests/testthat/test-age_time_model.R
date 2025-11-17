library(dplyr)
library(magrittr)

test_that("test age time model", {
  # preparing data
  tb_nl <- tb_nl_1966_1973 %>%
    mutate(
      survey_year = age + birthyr,
      survey_time = as.Date(paste0(survey_year, "-01-01"))
    ) %>% select(-birthyr) %>%
    group_by(age, survey_year, survey_time) %>%
    summarize(pos = sum(pos), tot = sum(tot), .groups = "drop")

  # try running model with PAVA for monotonization
  expect_no_error(
    out <- suppressWarnings(age_time_model(tb_nl, time_col = "survey_time", grouping_col = "survey_year", monotonize_method = "pava"))
  )

  # try running model with scam model for monotonization
  expect_no_error(
    out <- suppressWarnings(age_time_model(tb_nl, time_col = "survey_time", grouping_col = "survey_year", monotonize_method = "scam"))
  )

  # try running model with facet settings
  expect_no_error(plot(out))
  expect_no_error(plot(out, facet = FALSE))
})
