#' Fitting penalized spline model with Generalized Linear Mixture Model framework
#'
#' @param age the age vector.
#' @param res the serostatus vector.
#' @param s smoothing basis to use.
#' @param link link function to use.
#' @import mgcv
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
