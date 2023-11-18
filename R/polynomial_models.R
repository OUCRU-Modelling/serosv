X <- function(t, degree) {
  X_matrix <- matrix(rep(1, length(t)), ncol = 1)
  if (degree > 1) {
    for (i in 2:degree) {
      X_matrix <- cbind(X_matrix, i * t^(i-1))
    }
  }
  -X_matrix
}

predictor <- function(degree) {
  formula <- "cbind(tot-pos, pos)~-1"
  for (i in 1:degree) {
    formula <- paste0(formula, "+I(age^", i, ")")
  }
  formula
}

#' A polynomial model.
#'
#' Refers to section 6.1.1.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.
#' @param deg the degree of the model.
#'
#' @examples
#' df <- hav_bg_1964
#' model <- polynomial_model(
#'   df$age, df$pos, df$tot,
#'   deg=2
#'   )
#' plot(model)
#'
#' @export
polynomial_model <- function(age, pos, tot, deg=1) {
  model <- list()

  f <- predictor(deg)
  model$info <- glm(
    as.formula(f),
    family=binomial(link="log")
  )
  X <- X(age, deg)
  model$sp <- 1 - model$info$fitted.values
  model$foi <- X%*%model$info$coefficients
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "polynomial_model"
  model
}

#' The Farrington (1990) model.
#'
#' Refers to section 6.1.2.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.
#' @param start Named list of vectors or single vector.
#' Initial values for optimizer.
#' @param fixed Named list of vectors or single vector.
#' Parameter values to keep fixed during optimization.
#'
#' @examples
#' df <- rubella_uk_1986_1987
#' model <- farrington_model(
#'   df$age, df$pos, df$tot,
#'   start=list(alpha=0.07,beta=0.1,gamma=0.03)
#'   )
#' plot(model)
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(age, pos, tot, start, fixed=list())
{
  farrington <- function(alpha,beta,gamma) {
    p=1-exp((alpha/beta)*age*exp(-beta*age)
            +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)-gamma*age)
    ll=pos*log(p)+(tot-pos)*log(1-p)
    return(-sum(ll))
  }

  model <- list()

  model$info <- mle(farrington, fixed=fixed, start=start)
  alpha <- model$info@coef[1]
  beta  <- model$info@coef[2]
  gamma <- model$info@coef[3]
  model$sp <- 1-exp(
    (alpha/beta)*age*exp(-beta*age)
    +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
    -gamma*age)
  model$foi <- (alpha*age-gamma)*exp(-beta*age)+gamma
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "farrington_model"
  model
}

