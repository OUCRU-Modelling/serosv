#' Penalized Spline model
#'
#' @param age the age vector
#' @param spos the seropositive vector
#' @param s smoothing basis to use
#' @param link link function to use
#' @param framework which approach to fit the model ("pl" for penalized likelihood framework, "glmm" for generalized linear mixed model framework)
#'
#' @import mgcv
#' @importFrom stats binomial
#'
#' @return a penalized_spline_model object
#' @export
#'
#' @examples
#' data <- parvob19_be_2001_2003
#' mod <- penalized_spline_model(data$age, status = data$seropositive, framework="glmm")
#' mod$gam$info
penalized_spline_model <- function(age, pos=NULL,tot=NULL,status=NULL, s = "bs", link = "logit", framework = "pl", sp = NULL){
  stopifnot("Values for either `pos & tot` or `status` must be provided" = !is.null(pos) & !is.null(tot) | !is.null(status) )
  model <- list()
  age <- as.numeric(age)

  # check input whether it is line-listing or aggregated data
  if (!is.null(pos) & !is.null(tot)){
    pos <- as.numeric(pos)
    tot <- as.numeric(tot)
    model$datatype <- "aggregated"
  }else{
    pos <- as.numeric(status)
    tot <- rep(1, length(pos))
    model$datatype <- "linelisting"
  }

  # s <- mgcv:::s
  spos <- pos/tot

  if (framework == "pl"){
    model$info <- mgcv::gam(spos ~ s(age, bs = s, sp=sp), family = binomial(link = link))
    model$sp <- model$info$fitted.values
  }else if(framework == "glmm"){
    model$info <- mgcv::gamm(spos ~ s(age, bs = s, sp=sp), family = binomial(link = link))
    model$sp <- model$info$gam$fitted.values
  }else{
    stop(paste0('Invalid value for framework. Expected "pl" or "glmm", got ', framework))
  }

  # aggregate data after fitting for plotting
  model$df <- data.frame(age=age, pos = pos, tot = tot)
  model$foi <- est_foi(age, model$sp)
  model$framework <- framework

  class(model) <- "penalized_spline_model"
  model
}
