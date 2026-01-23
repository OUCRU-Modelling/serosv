#' Compare models
#'
#' @param data input data to fit into the models
#' @param method method to compare models. Can be one of the built-in methods or a function to compute the returned metrics (see Details).
#' @param ... models to be compared. Must be models created by serosv. If models' names are not provided, indices will be used instead for the `model` column in the returned data.frame.
#'
#'
#' @return
#' a data.frame of 4 columns
#'   \item{label}{name or index of the model}
#'   \item{type}{model type of the given model (a serosv model name)}
#'   \item{AIC}{AIC value for the model (lower value indicates better fit)}
#'   \item{BIC}{BIC value for the model (lower value indicates better fit)}
#'
#' @details
#' Built-in comparison methods include:
#' - computing AIC and BIC, which returns AIC, BIC values of the model if available
#' - cross validation, which reutns
#'
#' @importFrom magrittr %>%
#' @importFrom purrr imap_dfr as_mapper
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#'
#' @export
compare_models <- function(data, method="AIC/BIC",...){
  list(...) %>%
    imap_dfr(~ {
      # return error if input contains non-serosv models
      # if(!all(str_detect(class(.x), "_model"))) {
      #   stop("Inputs must be serosv models")
      # }

      # get function to compute comparison metrics
      metric_func <- switch(
        method,
        "AIC/BIC" = aic_bic,
        "CV" = cv,
        method
      )

      assert_that(is.function(metric_func),
                  msg = "Function to compute the metrics must be provided")

      # TODO: apply cross validation to get comparison metrics instead
      out <- metric_func(data, as_mapper(.x))

      assert_that("data.frame" %in% class(out),
                  msg = "Function to compute the metrics must return a data.frame")

      out %>% mutate(
        label = .y,
        .before = 1
      )

    })
}

# function returning goodness-of-fit metrics such as AIC/BIC, likelihood (with degree-of-freedom)
#' @importFrom stats4 logLik AIC BIC
aic_bic <- function(dat, mod_func){
  out <- mod_func(dat)

  aic <- tryCatch(stats4::AIC(out$info),
                  error = \(e){NULL})
  bic <- tryCatch(stats4::BIC(out$info),
                  error = \(e){NULL})
  ll <- tryCatch(stats4::logLik(out$info),
                     error = \(e){NULL})

  tibble(
    type = class(out),
    AIC = aic,
    BIC = bic,
    logLik = if (!is.null(ll)) as.numeric(ll) else NA,
    df = if (!is.null(ll) && !is.null(attr(ll, "df"))) attr(ll, "df") else NA,
    mod_out = list(out)
  )
}


# function to compute metrics from cross validation
# assess the generalization/prediction of the model
#' @importFrom stats4 logLik AIC BIC
#' @importFrom stats predict.glm
cv <- function(dat, mod_func){
  # making sure samples are the same across models
  set.seed(123)

  idx <- sample(nrow(dat), size = floor(0.2*nrow(dat)), replace=FALSE)
  fit_dat <- dat[-idx, ]
  test_dat <- dat[idx, ]

  out <- mod_func(fit_dat)

  # generate prediction
  pred <- predict(out, data.frame(age=test_dat[,1]), type="response")

  metrics <- list()
  metrics$type <- class(out)

  if(out$datatype == "aggregated"){
    # if data is aggregated
    seroprev_obs <- test_dat$pos/test_dat$tot

    # MAE - agnostic to sample size
    metrics$mae <- sum(pred - seroprev_obs)**2/nrow(test_dat)
    # compute logloss
    metrics$logloss <- sum(dbinom(test_dat$pos, test_dat$tot, prob=pred, log=TRUE))

  }else{
    # if data is linelisting
    # compute log-loss
    metrics$logloss <- sum(test_dat$status*log(pred) + (1 - test_dat$status)*log(1-pred))
    # and estimate auc
    metrics$auc <- as.numeric(pROC::auc(test_dat$status, pred))
  }

  metrics$mod_out <- list(out)

  as_tibble(metrics)
}


