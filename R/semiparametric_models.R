#' Penalized Spline model
#'
#' @param age the age vector
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
#' @param s smoothing basis to use
#' @param sp smoothing parameter
#' @param link link function to use
#' @param framework which approach to fit the model ("pl" for penalized likelihood framework, "glmm" for generalized linear mixed model framework)
#'
#' @import mgcv
#' @importFrom stats binomial
#'
#' @return a list of class penalized_spline_model with 6 attributes
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{framework}{either pl or glmm}
#'   \item{info}{fitted "gam" model when framework is pl or "gamm" model when framework is glmm}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'
#' @seealso [mgcv::gam()], [mgcv::gamm()] for more information the fitted gam and gamm model
#'
#' @export
#'
#' @examples
#' data <- parvob19_be_2001_2003
#' model <- penalized_spline_model(data$age, status = data$seropositive, framework="glmm")
#' model$gam$info
#' plot(model)
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
