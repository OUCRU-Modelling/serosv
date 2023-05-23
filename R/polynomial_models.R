#' Plot the prevalence and force of infection against age.
#'
#' Refers to section 6.1.
#'
#' @param model a model returned from this package's function.
#'
#' @examples
#' model <- polynomial_model(degree=3)
#' hav_be_1993_1994 %>%
#'   model$fit() %>%
#'   plot_p_foi_wrt_age()
#'
#' @export
plot_p_foi_wrt_age <- function(model)
{
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(model$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, model$sp, lwd=2)
    lines(age, model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(model$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#' A polynomial model.
#'
#' Refers to section 6.1.1.
#'
#' @param degree the degree of the model.
#'
#' @examples
#' muench_model <- polynomial_model(degree=1)
#' griffiths_model <- polynomial_model(degree=2)
#' grenfell_anderson_model <- polynomial_model(degree=3)
#'
#' @export
polynomial_model <- function(degree) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        as.formula(linear_predictor_str()),
        family=binomial(link="log")
      )
      X <- generate_X(age)
      model$sp <- 1 - model$info$fitted.values
      model$foi <- X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  linear_predictor_str <- function() {
    formula <- "cbind(tot-pos, pos)~-1"
    for (i in 1:degree) {
      formula <- paste0(formula, "+I(age^", i, ")")
    }
    formula
  }
  generate_X <- function(age) {
    X <- matrix(rep(1, length(age)), ncol = 1)
    if (degree > 1) {
      for (i in 2:degree) {
        X <- cbind(X, i * age^(i-1))
      }
    }
    -X
  }
  model
}

#' The Farrington (1990) model.
#'
#' Refers to section 6.1.2.
#'
#' @param parameters the parameters of the model.
#'
#' @examples
#' parameters <- c(alpha=0.07,beta=0.1,gamma=0.03)
#' model <- farrington_model(parameters)
#' model
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(parameters)
{
  model <- list()
  model$parameters <- parameters
  model$fit <- function(df) {
    with(c(df, parameters), {
      farrington <- function(alpha,beta,gamma) {
        p=1-exp((alpha/beta)*age*exp(-beta*age)
                +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)-gamma*age)
        ll=pos*log(p)+(tot-pos)*log(1-p)
        return(-sum(ll))
      }

      model$info <- mle(farrington, start=as.list(parameters))
      alpha <- coef(model$info)[1]
      beta  <- coef(model$info)[2]
      gamma <- coef(model$info)[3]
      model$sp <- 1-exp(
        (alpha/beta)*age*exp(-beta*age)
        +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
        -gamma*age)
      model$foi <- (alpha*age-gamma)*exp(-beta*age)+gamma
      model$df <- df
      model
    })
  }
  model
}
