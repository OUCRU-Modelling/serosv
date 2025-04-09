#' Penalized Spline model
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` column for aggregated data OR `age`, `status` for linelisting data
#' @param s smoothing basis to use
#' @param sp smoothing parameter
#' @param link link function to use
#' @param framework which approach to fit the model ("pl" for penalized likelihood framework, "glmm" for generalized linear mixed model framework)
#'
#' @importFrom mgcv gam gamm
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
#' data$status <- data$seropositive
#' model <- penalized_spline_model(data, framework="glmm")
#' model$info$gam
#' plot(model)
penalized_spline_model <- function(data, s = "bs", link = "logit", framework = "pl", sp = NULL){
  model <- list()

  data <- check_input(data)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

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
