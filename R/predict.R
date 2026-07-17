# ====== Predict function for serosv models ======
#' Prediction for serosv polynomial model
#'
#' A wrapper of predict.glm for direct prediction from polynomial_model object
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @import dplyr
#'
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.polynomial_model <- function(object, newdata=NULL, ...){
  # return seroprevalence
  1 - predict.glm(object$info, newdata, type="response", ...)
}

#' Prediction for serosv fractional polynomial model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.fp_model <- function(object, newdata=NULL, ...){
  # return seroprevalence
  predict.glm(object$info,newdata=newdata, type="response", ...)
}

#' Prediction for serosv Weibull model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.weibull_model <- function(object, newdata=NULL, ...){
  predict.glm(object$info,data.frame("age" = newdata$age),type="response", ...)
}


#' Prediction for serosv local polynomial model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @return prediction output
#' @export
predict.lp_model <- function(object, newdata=NULL,...){
  predict(object$info, data.frame(age = newdata$age), type="response", ...)
}


#' Prediction for serosv penalized spline model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @importFrom mgcv predict.gam
#' @return prediction output
#' @seealso
#' [mgcv::predict.gam()] for more information on the predict function
#' @export
predict.penalized_spline_model <- function(object, newdata=NULL,...){

  # handle different output for different frameworks
  if(object$framework == "pl"){
    gam_obj <- object$info
  }else{
    gam_obj <- object$info$gam
  }

  predict.gam(gam_obj, newdata, type="response", ...)
}

#' Prediction for serosv Farrington model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @return prediction output
#' @export
predict.farrington_model <- function(object, newdata=NULL,...){
  alpha <- object$info@fullcoef[1]
  beta  <- object$info@fullcoef[2]
  gamma <- object$info@fullcoef[3]

  # 1-exp(
  #   (alpha/beta)*newdata[[1]]*exp(-beta*newdata[[1]])
  #   +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*newdata[[1]])-1)
  #   -gamma*newdata[[1]])
  object$sp_mod(newdata[[1]], alpha, beta, gamma)
}

#' Predict from an hierarchical bayesian model
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary arguments
#' @import dplyr
#'
#' @return list of confidence interval for seroprevalence and foi. Each confidence interval dataframe with 4 variables, x and y for the fitted values and ymin and ymax for the confidence interval
#' @export
predict.hierarchical_bayesian_model <- function(object, newdata=NULL, ...){
  out_x <- object$df$age
  out.DF <- NULL

  if (object$type == "far3"){
    alpha1 <- summary(object$info)$summary["alpha1", "50%"]
    alpha2 <- summary(object$info)$summary["alpha2", "50%"]
    alpha3 <- summary(object$info)$summary["alpha3", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = object$sp_func(out_x, alpha1, alpha2, alpha3),
    )
  }else if(object$type == "far2"){
    alpha1 <- summary(object$info)$summary["alpha1", "50%"]
    alpha2 <- summary(object$info)$summary["alpha2", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = object$sp_func(out_x, alpha1, alpha2),
    )

  }else if(object$type == "log_logistic"){
    alpha1 <- summary(object$info)$summary["alpha1", "50%"]
    alpha2 <- summary(object$info)$summary["alpha2", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = object$sp_func(out_x, alpha1, alpha2),
    )
  }else{
    warning('Expect model type to be one of the following: "far3", "far2", "log_logistic"')
  }

  out.DF
}

#' Predict from the age_time_mdoel
#'
#' @param object serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param modtype either "monotonized" (to predict using monotonized model) or "non-monotonized"
#' @param ... arbitrary argument
#'
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @return confidence interval dataframe with n_group x 3 cols, the columns are `group`, `sp_df`, `foi_df`
#' @export
predict.age_time_model <- function(object, newdata, modtype="monotonized", ...){
  # resolve no visible binding note
  df <- monotonized_info <- monotonized_ci_mod <- age <- info <- fit <- se.fit <- sp_df <- foi_df <- NULL
  data <- age_df <- NULL

  # check which type of model user wants to predict
  modtype <- if (is.null(list(...)[["modtype"]])) "monotonized" else list(...)$modtype
  assert_that(
    modtype == "monotonized" | modtype == "non-monotonized",
    msg = "modtype argument must be eithers 'monotonized' or 'non-monotonized'"
  )

  p <- (1 - object$ci) / 2

  # check whether newdata match the requirement
  if(!all(c(object$grouping_col, "age") %in% colnames(newdata)) ){
    stop(paste0(
      "Data must have `",
      object$grouping_col,
      "`, `age` columns"
    ))
  }


  # generate the newdata by survey time for prediction
  out <- newdata %>%
    group_by(.data[[object$grouping_col]]) %>%
    nest() %>%
    rename(age_df = data) %>%
    left_join(
      # can only predict for the survey time fitted to the model
      object$out,
      join_by(!!sym(object$grouping_col))
    )

  # --- use the monotonized model for prediction and ci-----
  if(modtype == "monotonized"){
    out <- out %>%
      mutate(
        sp_df = pmap(list(monotonized_info, monotonized_ci_mod, age_df), \(mod, ci_mod, grid){
          data.frame(
            x = grid$age,
            y = predict(mod, grid, type = "response"),
            ymin = predict(ci_mod$ymin, grid, type = "response"),
            ymax = predict(ci_mod$ymax, grid, type = "response")
          )
        })
      )
  }else{
    # --- if user specify non-monotonized then simply compute CI from gam model-----
    out <- out %>%
      mutate(
        sp_df = map2(info, age_df, \(mod, grid){
          link_inv <- mod$family$linkinv
          dataset <- mod$model[,1:2]
          n <- nrow(dataset) - length(mod$coefficients)

          predict(mod, grid, se.fit = TRUE)  %>%
            as_tibble()  %>%
            select(fit, se.fit) %>%
            mutate(
              x = grid,
              ymin = link_inv(fit + qt(    p, n) * se.fit),
              ymax = link_inv(fit + qt(1 - p, n) * se.fit),
              y = link_inv(fit)
            )  %>%
            select(- se.fit)
        })
      )
  }

  # --- finally, compute FOI -----
  out <- out %>%
    mutate(
      foi_df = map2(age_df, sp_df, \(grid, sp){
        foi_x <- sort(unique(grid$age))
        foi_x <- foi_x[c(-1, -length(foi_x) )]

        data.frame(
          x = foi_x,
          y = est_foi(grid$age, sp$y)
        )
      })
    ) %>%
    select(!!sym(object$grouping_col), sp_df, foi_df)

  out
}

