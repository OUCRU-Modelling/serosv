#' The Weibull model.
#'
#' Refers to section 6.1.2.
#'
#' @param t the time vector.
#' @param spos the seropositive vector.
#'
#' @examples
#' df <- hcv_be_2006[order(hcv_be_2006$dur), ]
#' model <- weibull_model(
#'   t=df$dur,
#'   spos=df$seropositive
#'   )
#' plot(model)
#'
#' @export
weibull_model <- function(t, spos)
{
  model <- list()

  model$info <- glm(
    spos~log(t),
    family=binomial(link="cloglog")
    )
  b0 <- coef(model$info)[1]
  b1 <- coef(model$info)[2]
  model$foi <- exp(b0)*b1*exp(log(t))^(b1-1)
  model$sp <- 1-exp(-exp(b0)*t^b1)
  model$df <- data.frame(t=t, spos=spos)

  class(model) <- "weibull_model"
  model
}


