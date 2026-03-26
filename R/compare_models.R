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
#' - cross validation, which returns MSE and logloss (negative )
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
      metric_func <- if(is.character(method)){
        switch(
          method,
          "AIC/BIC" = aic_bic,
          "CV" = cv,
          method
        )
      }else{
        method
      }

      assert_that(is.function(metric_func),
                  msg = "Function to compute the metrics must be provided")

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
    mod_out = list(out),
    plots = list(plot(out)+ggtitle(paste("Fitted model using", class(out))))
  )
}


# function to compute metrics from cross validation
# assess the generalization/prediction of the model
#' @importFrom stats4 logLik AIC BIC
#' @importFrom stats predict.glm
#' @import tidyr dplyr pROC
cv <- function(dat, mod_func, k=4){
  # assign each row of data to each fold
  idx_fold <- sort(rep(1:k, length.out=nrow(dat)))

  metrics <- lapply(1:k, \(fold){
    curr_metric <- list()

    # split data
    fit_dat <- dat[idx_fold != fold, ]
    test_dat <- dat[idx_fold == fold, ]

    # get model info
    out <- mod_func(fit_dat)
    curr_metric$type <- class(out)
    # generate prediction
    pred <- predict(out, data.frame(age=test_dat[,1]), type="response")

    if(out$datatype == "aggregated"){
      # if data is aggregated
      seroprev_obs <- test_dat$pos/test_dat$tot

      # MSE
      curr_metric$mse <- sum((pred - seroprev_obs)**2)/nrow(test_dat)
      # compute logloss (negative binomial loglikelihood)
      curr_metric$logloss <- -sum(
        dbinom(test_dat$pos, test_dat$tot, prob=pred, log=TRUE),
        na.rm = TRUE)

    }else{
      # if data is linelisting
      # make sure pred is slightly higher than 0 and lower than 1 to avoid log(0)
      eps <- .Machine$double.eps  # smallest positive floating point number ~ 2.2e-16
      pred <- pmax(eps, pmin(1 - eps, pred))

      # compute logloss (negative bernoulli loglikelihood)
      curr_metric$logloss <- -sum(test_dat$status*log(pred) + (1 - test_dat$status)*log(1-pred),
                                  na.rm = TRUE)
      # and estimate auc
      curr_metric$auc <- as.numeric(pROC::auc(test_dat$status, pred))
    }

    curr_metric
  }) %>%
  bind_rows() %>%
  summarise(
    # compute average of the metrics
    across(where(is.numeric), mean),
    # for type, simply get the first one
    type = first(type))

  # also return the model when it is fitted using the whole data
  out <- mod_func(dat)

  metrics %>%
    mutate(
      mod_out = list(out),
      plots = list(plot(out)+ggtitle(paste("Fitted model using", class(out))))
    )
}


