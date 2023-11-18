#' A local polynomial model.
#'
#' Refers to section 7.1. and 7.2.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#' @param nn Nearest neighbor component of the smoothing parameter.
#' Default value is 0.7, unless either h is provided, in which case the default is 0.
#' @param h The constant component of the smoothing parameter. Default: 0.
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' df <- mumps_uk_1986_1987
#' model <- lp_model(
#'   df$age, df$pos, df$tot,
#'   nn=0.7, kern="tcub"
#'   )
#' plot(model)
#'
#' @importFrom locfit locfit
#' @importFrom locfit lp
#' @importFrom graphics par
#'
#' @export
# library(locfit)
lp_model <- function(age, pos, tot, kern="tcub", nn=0, h=0, deg=2) {
  if (missing(nn) & missing(h)) {
    nn <- 0.7 # default nn from lp()
  }

  model <- list()

  y <- pos/tot
  estimator <- lp(age, deg=deg, nn=nn, h=h)
  model$pi  <- locfit(y~estimator, family="binomial", kern=kern)
  model$eta <- locfit(y~estimator, family="binomial", kern=kern, deriv=1)
  model$sp  <- fitted(model$pi)
  model$foi <- fitted(model$eta)*fitted(model$pi) # λ(a)=η′(a)π(a)
  model$df  <- list(age=age, pos=pos, tot=tot)

  class(model) <- "lp_model"
  model
}
