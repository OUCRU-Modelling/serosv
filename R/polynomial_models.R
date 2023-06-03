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
#' @param degree the degree of the model.
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

plot.polynomial_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#' The Farrington (1990) model.
#'
#' Refers to section 6.1.2.
#'
#' @param parameters the parameters of the model.
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
  alpha <- coef(model$info)[1]
  beta  <- coef(model$info)[2]
  gamma <- coef(model$info)[3]
  model$sp <- 1-exp(
    (alpha/beta)*age*exp(-beta*age)
    +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
    -gamma*age)
  model$foi <- (alpha*age-gamma)*exp(-beta*age)+gamma
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "farrington_model"
  model
}

plot.farrington_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}
