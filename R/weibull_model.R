#' The Weibull model.
#'
#' Refers to section 6.1.2.
#'
#' @param data the input data frame, must either have `t`, `pos`, `tot` column for aggregated data OR `t`, `status` for linelisting data
#'
#' @importFrom stats coef
#'
#' @examples
#' df <- hcv_be_2006[order(hcv_be_2006$dur), ]
#' df$t <- df$dur
#' df$status <- df$seropositive
#' model <- weibull_model(df)
#' plot(model)
#'
#' @return list of class weibull_model with the following items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "glm" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'
#' @seealso [stats::glm()] for more information on the fitted "glm" object
#'
#' @export
weibull_model <- function(data)
{
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data, heterogeneity_col = "t")
  t <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

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


