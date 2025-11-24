#' The Weibull model.
#'
#' @description Model seroprevalence as a function of duration since vaccination using the Weibull
#' model, where the force of infection is assumed to vary monotonically with duration.
#'
#' @details
#' For a Weibull model, the prevalence is given by
#' \deqn{
#'  \pi (d) = 1 - e^{ - \beta_0 d ^ {\beta_1}}
#' }
#' Where \eqn{d} is exposure time (difference between age of vaccination and age at test)
#'
#' Which implies the force of infection to be the monotonic function
#' \deqn{
#'  \lambda(d) = \beta_0 \beta_1 d^{\beta_1 - 1}
#' }
#'
#' Refer to section 6.1.2. of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
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
  data <- check_input(data, stratum_col = "t")
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


