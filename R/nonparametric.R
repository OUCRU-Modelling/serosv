#' A local polynomial model.
#'
#' Refers to section 7.1. and 7.2.
#'
#' @param age the age vector.
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
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
#'   df$age, pos = df$pos, tot = df$tot,
#'   nn=0.7, kern="tcub"
#'   )
#' plot(model)
#'
#' @importFrom locfit locfit
#' @importFrom locfit lp
#' @importFrom graphics par
#' @importFrom stats fitted
#'
#' @return a list of class lp_model with 6 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{pi}{fitted locfit object for pi}
#'   \item{eta}{fitted locfit object for eta}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [locfit::locfit()] for more information on the fitted locfit object
#'
#' @export
lp_model <- function(age, pos=NULL, tot=NULL, status=NULL, kern="tcub", nn=0, h=0, deg=2) {
  if (missing(nn) & missing(h)) {
    nn <- 0.7 # default nn from lp()
  }

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
