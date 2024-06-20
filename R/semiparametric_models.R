#' Fitting penalized spline model with Generalized Linear Mixture Model framework
#'
#' @param age the age vector.
#' @param res the serostatus vector.
#' @param s smoothing basis to use.
#' @param link link function to use.
#' @import mgcv
#' @importFrom stats binomial
#'
#' @export
#'
#' @examples
#' data <- parvob19_be_2001_2003
#' mod <- glmm_ps_model(data$age, data$seropositive)
#' mod$gam$info
glmm_ps_model <- function(age, res, s = "bs", link = "logit"){
  model <- list()

  model$info <- mgcv::gamm(res ~ s(age, bs = s), family = binomial(link = link))

  model$sp <- model$info$gam$fitted.values
  model$foi <- est_foi(age, model$sp)

  # aggregate data after fitting
  data <- transform(age, res)

  model$df <- data.frame(age=data$t, pos = data$pos, tot = data$tot)
  class(model) <- "glmm_ps_model"
  model
}

# pl_ps_model <- function(age, res, s = "bs", link = "logit"){
#   model <- list()
#
#   model$info <- gam::gam(res ~ s(age, bs = s), family = binomial(link = link))
#
#   model$sp <- model$info$fitted.values
#   model$sp <- est_foi(age, model$sp)
#
#   model$df <- data.frame(age=data$t, pos = data$pos, tot = data$tot)
#   class(model) <- "glmm_ps_model"
#   model
# }


# Merged
#' Penalized Spline model
#'
#' @param age the age vector
#' @param res the serostatus vector
#' @param s smoothing basis to use
#' @param link link function to use
#' @param framework which approach to fit the model ("ps" for penalized spline framework, "glmm" for generalized linear mixed model framework)
#'
#' @import mgcv
#' @importFrom stats binomial
#' @importFrom glue glue
#'
#' @return a penalized_spline_model object
#' @export
#'
#' @examples
#' data <- parvob19_be_2001_2003
#' mod <- penalized_spline_model(data$age, data$seropositive, framework="glmm")
#' mod$gam$info
penalized_spline_model <- function(age, res, s = "bs", link = "logit", framework = "ps"){
  model <- list()

  if (framework == "ps"){
    model$info <- mgcv::gam(res ~ s(age, bs = s), family = binomial(link = link))
    model$sp <- model$info$fitted.values
  }else if(framework == "glmm"){
    model$info <- mgcv::gamm(res ~ s(age, bs = s), family = binomial(link = link))
    model$sp <- model$info$gam$fitted.values
  }else{
    stop(glue('Invalid value for framework. Expected "ps" or "glmm", got "{framework}"'))
  }

  model$foi <- est_foi(age, model$sp)

  # aggregate data after fitting for plotting
  data <- transform(age, res)

  model$df <- data.frame(age=data$t, pos = data$pos, tot = data$tot)
  model$foi <- est_foi(age, model$sp)
  model$framework <- framework

  class(model) <- "penalized_spline_model"
  model
}
