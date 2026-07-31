# tests/testthat/test-compare_models.R
#
# Tests for compare_models(), aic_bic(), and cv()

library(testthat)
library(dplyr)

# ===== Set up ========
# A minimal mock "model" constructor + S3 methods
mock_model <- function(dat, datatype = "aggregated") {
  structure(
    list(dat = dat, datatype = datatype, info = NULL),
    class = "mock_model"
  )
}

predict.mock_model <- function(object, newdata, ...) {
  rep(0.5, nrow(newdata))
}

# A dummy metric function with extra argument for testing
dummy_metric <- function(dat, mod_func, extra = 0) {
  fit <- mod_func(dat)
  data.frame(
    type = class(fit),
    score = 1 + extra
  )
}

# A metric function returning something that is NOT a data.frame, to test
# the assert_that() guard on the return type.
bad_metric <- function(dat, mod_func) {
  list(not = "a data.frame")
}

# ====== Test compare_models function =======
test_that("compare_models works with a custom metric function", {
  res <- compare_models(
    data = hbv_ru_1999,
    method = dummy_metric,
    model_a = ~ mock_model(.x),
    model_b = ~ mock_model(.x)
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("label", "type", "score") %in% names(res)))
  expect_equal(res$label, c("model_a", "model_b"))
  expect_true(all(res$score == 1))
})

test_that("compare_models forwards method_args to the metric function", {
  res <- compare_models(
    data = hbv_ru_1999,
    method = dummy_metric,
    method_args = list(extra = 5),
    model_a = ~ mock_model(.x)
  )

  expect_equal(res$score, 6)
})


# ====== Test built-in metric functions =======
test_that("test aic_bic function behavior", {
  # data
  dat <- hav_be_2002 |> rename(status = seropositive)

  # try metrics w/ different model
  metric_1 <- aic_bic(dat, \(dat) polynomial_model(dat, k=3))
  metric_2 <- aic_bic(dat, lp_model)

  # check returned metrics
  expect_true(all(c("AIC", "BIC", "logLik") %in% names(metric_1)))
  expect_true(all(c("AIC", "BIC", "logLik") %in% names(metric_2)))

  # check datatype
  expect_s3_class(metric_1, "data.frame")
  expect_s3_class(metric_2, "data.frame")

  # AIC/BIC should exist for polynomial
  expect_true(all(
    !is.na(metric_1)
  ))

  # AIC/BIC expected to not available for local polynomial
  expect_true(all(
    is.na(metric_2[, c("AIC", "BIC", "logLik")])
  ))
})

test_that("test cv function behavior", {
  # data
  linelist_dat <- hav_be_2002 |> rename(status = seropositive)
  aggregated_dat <- transform_data(linelist_dat)

  # metrics for different data formats
  linelist_metric <- cv(linelist_dat, \(dat) polynomial_model(dat, k=3))
  aggr_metric <- cv(aggregated_dat, \(dat) polynomial_model(dat, k=3))

  # expect return logloss and auc for serostatus data
  expect_true(all(c("logloss", "auc") %in% names(linelist_metric)))

  # expect return logloss and mse for aggregated data
  expect_true(all(c("logloss", "mse") %in% names(aggr_metric)))

  # expect function return data.frame
  expect_s3_class(linelist_metric, "data.frame")
  expect_s3_class(aggr_metric, "data.frame")

  # expect no NAs
  expect_true(all(!is.na(linelist_metric)))
  expect_true(all(!is.na(aggr_metric)))
})

# ====== Test compare_models example use cases =======
test_that("test compare_model with aic_bic", {
  data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ] %>%
    rename(status = seropositive)

  aggregated <- transform_data(data, stratum_col = "age", status_col="status")

  expect_no_error(
    compare_models(
      data = data,
      method = "AIC/BIC",
      griffith = ~polynomial_model(.x, k=2),
      penalized_spline = penalized_spline_model,
      farrington = ~farrington_model(.x,
                                     start=list(beta=0.12,gamma=0.05),
                                     fixed=list(alpha=0.1)),
      local_polynomial = lp_model # expect to not return any values
    ) %>% suppressWarnings()
  )

  expect_no_error(
    compare_models(
      data = aggregated,
      method = "AIC/BIC",
      griffith = ~polynomial_model(.x, k=2),
      penalized_spline = penalized_spline_model,
      farrington = ~farrington_model(.x,
                                     start=list(beta=0.1,gamma=0.03),
                                     fixed=list(alpha=0.07)),
      local_polynomial = lp_model # expect to not return any values
    ) %>% suppressWarnings()
  )
})

test_that("test compare_model with CV", {
  data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ] %>%
    rename(status = seropositive)

  aggregated <- transform_data(data, stratum_col = "age", status_col="status")

  expect_no_error(
    compare_models(
      data,
      method = "CV",
      griffith = ~polynomial_model(.x, k=2),
      penalized_spline = penalized_spline_model,
      farrington = ~farrington_model(.x,
                                     start=list(beta=0.12,gamma=0.05),
                                     fixed=list(alpha=0.1)),
      local_polynomial = lp_model,
      method_args = list(
        k = 4 # adjust the fold k for the CV function
      )
    ) %>% suppressWarnings()
  )

  expect_no_error(
    compare_models(
      data = aggregated,
      method = "CV",
      griffith = ~polynomial_model(.x, k=2),
      penalized_spline = penalized_spline_model,
      farrington = ~farrington_model(.x,
                                     start=list(beta=0.1,gamma=0.03),
                                     fixed=list(alpha=0.07)),
      local_polynomial = lp_model,
      method_args = list(
        k = 4 # adjust the fold k for the CV function
      )
    ) %>% suppressWarnings()
  )
})
