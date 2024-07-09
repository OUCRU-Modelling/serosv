#' The Weibull model.
#'
#' Refers to section 6.1.2.
#'
#' @param t the time vector.
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
#'
#' @importFrom stats coef
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
weibull_model <- function(t, status=NULL, pos=NULL, tot=NULL)
{
  stopifnot("Values for either `pos & tot` or `status` must be provided" = !is.null(pos) & !is.null(tot) | !is.null(status) )
  model <- list()
  t <- as.numeric(t)

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

  spos <- pos/tot
  model$info <- glm(
    spos~log(t),
    family=binomial(link="cloglog")
    )
  b0 <- coef(model$info)[1]
  b1 <- coef(model$info)[2]
  model$foi <- exp(b0)*b1*exp(log(t))^(b1-1)
  model$sp <- 1-exp(-exp(b0)*t^b1)
  model$df <- data.frame(age=t, pos=pos, tot=tot)

  class(model) <- "weibull_model"
  model
}


