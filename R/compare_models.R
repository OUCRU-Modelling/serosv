#' Generate table of metrics for model comparison
#'
#' @param data input data to fit into the models
#' @param method method to compare models. Can be one of the built-in methods or a function to compute the returned metrics (see Details).
#' @param method_args additional arguments to be passed to the method function.
#' @param ... models to be compared. Must be models created by serosv. If models' names are not provided, indices will be used instead for the `model` column in the returned data.frame.
#'
#'
#' @return
#' a data.frame with the following columns
#'   \item{label}{name or index of the model}
#'   \item{type}{model type of the given model (a serosv model name)}
#'   \item{mod_out}{the fitted models}
#'   \item{plots}{the plots for each of the fitted model}
#'   \item{metrics columns}{the columns for metrics of comparison, the number of which depends on the function that generate these metrics}
#'
#' @details
#' Built-in comparison methods include:
#' \itemize{
#' \item{computing AIC and BIC, which returns AIC, BIC values of the model if available}
#' \item{cross validation (perform k-fold validation), which returns MSE and logloss (negative log Binomial likelihood)
#'  for aggregated data, or AUC and logloss (negative log Bernoulli likelihood)
#'  for linelisting data}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom purrr imap_dfr as_mapper
#' @importFrom stringr str_detect
#' @importFrom assertthat assert_that
#'
#' @examples
#' comparison_table <- suppressWarnings(
#'   compare_models(
#'     data = hav_bg_1964,
#'     method = "CV",
#'     polynomial_mod = ~polynomial_model(.x, k=1),
#'     penalized_spline = penalized_spline_model,
#'     farrington = ~farrington_model(.x, start=list(alpha=0.3,beta=0.1,gamma=0.03))
#'   )
#' )
#' # view table of metrics
#' comparison_table
#' # view the model fitted with the whole dataset
#' comparison_table$plots
#' @export
compare_models <- function(data, method="AIC/BIC", method_args=list(), ...){
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

      # out <- metric_func(data, as_mapper(.x))
      out <- do.call(
        metric_func,
        c(
          list(dat = data, mod_func = as_mapper(.x)),
          method_args
        )
      )

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
    AIC = if (!is.null(aic)) as.numeric(aic) else NA,
    BIC = if (!is.null(bic)) as.numeric(bic) else NA,
    logLik = if (!is.null(ll)) as.numeric(ll) else NA,
    df = if (!is.null(ll) && !is.null(attr(ll, "df"))) attr(ll, "df") else NA,
    mod_out = list(out),
    plots = list(plot(out)+ggtitle(paste("Fitted model using", class(out))))
  )
}


# function to compute metrics from cross validation
# assess the generalization/prediction of the model
#' @importFrom stats4 logLik AIC BIC
#' @importFrom stats predict.glm dbinom
#' @import tidyr dplyr pROC
cv <- function(dat, mod_func, k=4){
  # resolve no visible binding NOTE during check()
  type <- NULL

  # assign each row of data to each fold
  idx_fold <- sample(rep(1:k, length.out=nrow(dat)))

  metrics <- lapply(1:k, \(fold){
    curr_metric <- list()

    # split data
    fit_dat <- dat[idx_fold != fold, ]
    test_dat <- dat[idx_fold == fold, ]

    # get model info
    out <- mod_func(fit_dat)
    curr_metric$type <- class(out)
    # generate prediction
    pred <- predict(out, data.frame(age=test_dat[,1]))

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
      # compute logloss (negative bernoulli loglikelihood)
      curr_metric$logloss <- -sum(dbinom(test_dat$status, 1, prob=pred, log=TRUE),
                                  na.rm = TRUE)
      # and estimate auc
      curr_metric$auc <- as.numeric(pROC::auc(test_dat$status, pred, quiet=TRUE))
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


