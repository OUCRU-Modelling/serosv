#' A local polynomial model.
#'
#' Refers to section 7.1. and 7.2.
#'
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#'
#' @param nn Nearest neighbor component of the smoothing parameter.
#' Default value is 0.7, unless either h is provided, in which case the default is 0.
#'
#' @param h The constant component of the smoothing parameter. Default: 0.
#'
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' model <- lp_model()
#' mumps_uk_1986_1987 %>%
#'   model$fit() %>%
#'   plot_p_foi_wrt_age()
#'
#' @importFrom locfit locfit
#'
#' @export
# library(locfit)
lp_model <- function(kern="tcub", nn=0, h=0, deg=2) {
  if (is_missing(nn) & is_missing(h)) {
    nn <- 0.7 # default nn from lp()
  }
  model <- list()
  model$parameters <- list(kern, nn, h, deg)
  model$fit <- function(df) {
    with(c(df, model$parameters), {
      y <- pos/tot
      predictor <- lp(age, deg=deg, nn=nn, h=h)
      lpfit   <- locfit(y~predictor, family="binomial", kern=kern)
      lpfitd1 <- locfit(y~predictor, family="binomial", kern=kern, deriv=1)
      model$sp <- fitted(lpfit)
      model$foi <- fitted(lpfitd1) * fitted(lpfit)
      model$df <- df
      model
    })
  }
  model
}
