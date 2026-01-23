# ====== Predict function for serosv models ======
#' Prediction for serosv polynomial model
#'
#' A wrapper of predict.glm for direct prediction from polynomial_model object
#'
#' @param x serosv models
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @import dplyr
#'
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.polynomial_model <- function(x, newdata=NULL, ...){
  predict.glm(x$info, newdata, ...)
}

#' Prediction for serosv fractional polynomial model
#'
#' @param x serosv models
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.fp_model <- function(x, newdata=NULL, ...){

  predict.glm(x$info,newdata=newdata, ...)
}

#' Prediction for serosv Weibull model
#'
#' @param x serosv models
#' @param ... arbitrary argument
#'
#' @importFrom stats predict.glm
#' @return prediction output
#' @seealso
#' [stats::predict.glm()] for more information on the predict function
#' @export
predict.weibull_model <- function(x, newdata=NULL, ...){
  predict.glm(x$info,data.frame("log(t)" = newdata$`log(t)`), ...)
}


#' Prediction for serosv local polynomial model
#'
#' @param x serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @return prediction output
#' @export
predict.lp_model <- function(x, newdata=NULL,...){
  predict(x$info, data.frame(age = newdata[[1]]), ...)
}


#' Prediction for serosv penalized spline model
#'
#' @param x serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @importFrom mgcv predict.gam
#' @return prediction output
#' @seealso
#' [mgcv::predict.gam()] for more information on the predict function
#' @export
predict.penalized_spline_model <- function(x, newdata=NULL,...){

  # handle different output for different frameworks
  if(x$framework == "pl"){
    gam_obj <- x$info
  }else{
    gam_obj <- x$info$gam
  }

  predict.gam(gam_obj, newdata, ...)
}

#' Prediction for serosv Farrington model
#'
#' @param x serosv models
#' @param newdata data.frame with age column to generate prediction
#' @param ... arbitrary argument
#'
#' @return prediction output
#' @export
predict.farrington_model <- function(x, newdata=NULL,...){
  alpha <- x$info@coef[1]
  beta  <- x$info@coef[2]
  gamma <- x$info@coef[3]

  1-exp(
    (alpha/beta)*newdata[[1]]*exp(-beta*newdata[[1]])
    +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*newdata[[1]])-1)
    -gamma*newdata[[1]])
}

#' Predict from an hierarchical bayesian model
#'
#' @param x serosv models
#' @param ... arbitrary arguments
#' @import dplyr
#'
#' @return list of confidence interval for seroprevalence and foi. Each confidence interval dataframe with 4 variables, x and y for the fitted values and ymin and ymax for the confidence interval
#' @export
predict.hierarchical_bayesian_model <- function(x, newdata=NULL, ...){
  out_x <- x$df$age
  out.DF <- NULL

  if (x$type == "far3"){
    alpha1 <- x$info["alpha1", "50%"]
    alpha2 <- x$info["alpha2", "50%"]
    alpha3 <- x$info["alpha3", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = x$sp_func(out_x, alpha1, alpha2, alpha3),
    )
  }else if(x$type == "far2"){
    alpha1 <- x$info["alpha1", "50%"]
    alpha2 <- x$info["alpha2", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = x$sp_func(out_x, alpha1, alpha2),
    )

  }else if(x$type == "log_logistic"){
    alpha1 <- x$info["alpha1", "50%"]
    alpha2 <- x$info["alpha2", "50%"]

    out.DF <- data.frame(
      x = out_x,
      y = x$sp_func(out_x, alpha1, alpha2),
    )
  }else{
    warning('Expect model type to be one of the following: "far3", "far2", "log_logistic"')
  }

  out.DF
}

#' Predict from the age_time_mdoel
#'
#' @param x serosv models
#' @param ... arbitrary argument
#'
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @return confidence interval dataframe with n_group x 3 cols, the columns are `group`, `sp_df`, `foi_df`
#' @export
predict.age_time_model <- function(x, ci=0.95, le = 100, ...){
  # resolve no visible binding note
  df <- monotonized_info <- monotonized_ci_mod <- age <- info <- fit <- se.fit <- sp_df <- foi_df <- NULL

  # check which type of model user wants to visualize
  modtype <- if (is.null(list(...)[["modtype"]])) "monotonized" else list(...)$modtype
  assert_that(
    modtype == "monotonized" | modtype == "non-monotonized",
    msg = "modtype argument must be eithers 'monotonized' or 'non-monotonized'"
  )

  p <- (1 - ci) / 2

  # use model to generate seroprev (with CI) and FOI on a finer grid for plotting
  age_range <- range(bind_rows(x$out$df)$age)
  out <- x$out %>%
    mutate(
      age = map(df, \(dat){
        seq(age_range[1], age_range[2], length.out = le)
      })
    )

  # --- use the monotonized model for prediction and ci-----
  if(modtype == "monotonized"){
    out <- out %>%
      mutate(
        sp_df = pmap(list(monotonized_info, monotonized_ci_mod, age), \(mod, ci_mod, grid){
          data.frame(
            x = grid,
            y = predict(mod, list(age = grid), type = "response"),
            ymin = predict(ci_mod$ymin, list(age = grid), type = "response"),
            ymax = predict(ci_mod$ymax, list(age = grid), type = "response")
          )
        })
      )
  }else{
    # --- if user specify non-monotonized then simply compute CI from gam model-----
    out <- out %>%
      mutate(
        sp_df = map2(info, age, \(mod, grid){
          link_inv <- mod$family$linkinv
          dataset <- mod$model[,1:2]
          n <- nrow(dataset) - length(mod$coefficients)

          predict(mod, data.frame(age = grid), se.fit = TRUE)  %>%
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
      foi_df = map2(age, sp_df, \(grid, sp){
        foi_x <- sort(unique(grid))
        foi_x <- foi_x[c(-1, -length(foi_x) )]

        tibble(
          x = foi_x,
          y = est_foi(grid, sp$y)
        )
      })
    ) %>%
    select(!!sym(x$grouping_col), sp_df, foi_df)
}

