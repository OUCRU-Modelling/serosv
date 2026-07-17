#' The Weibull model.
#'
#' @description Model seroprevalence as a function of age using the Weibull
#' model, where the force of infection is assumed to vary monotonically with age.
#'
#' @details
#' For a Weibull model, the prevalence is given by
#' \deqn{
#'  \pi (a) = 1 - e^{ - \beta_0 a ^ {\beta_1}}
#' }
#' Where \eqn{a} is the age, which may refer to biological age or a
#' time scale of interest (e.g., time since vaccination).
#'
#' Which implies the force of infection to be the monotonic function
#' \deqn{
#'  \lambda(a) = \beta_0 \beta_1 a^{\beta_1 - 1}
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
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR
#'  `age`, `status` (for linelisting data)
#' @param age_col name of the `age` column (default age_col="age")
#' @param pos_col name of the `pos` column (default pos_col="pos")
#' @param tot_col name of the `tot` column (default tot_col="tot")
#' @param status_col name of the `status` column (default status_col="status")
#' @param ... additional arguments to be passed to `glm()` function that fits the model
#'
#' @importFrom stats coef
#'
#' @examples
#' df <- hcv_be_2006[order(hcv_be_2006$dur), ]
#' model <- weibull_model(df, age_col="dur", status_col="seropositive")
#' plot(model)
#'
#' @return list of class weibull_model with the following items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "glm" object}
#'   \item{sp}{estimated seroprevalence}
#'   \item{foi}{estimated force of infection}
#'   \item{foi_mod}{function to generate FoI given age, and parameter values}
#'
#' @seealso [stats::glm()] for more information on the fitted "glm" object
#'
#' @export
weibull_model <- function(data,
                          age_col="age",pos_col="pos", tot_col="tot", status_col="status",
                          ...)
{
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data, stratum_col = age_col, pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  spos <- pos/tot
  model$info <- glm(
    spos~log(age),
    family=binomial(link="cloglog"),
    ...
  )
  b0 <- coef(model$info)[1]
  b1 <- coef(model$info)[2]

  model$sp <- model$info$fitted.values
  model$foi_mod <- function(age, b0, b1){ exp(b0)*b1*exp(log(age))^(b1-1) }
  model$foi <- model$foi_mod(age, b0, b1)
  model$df <- data.frame(age=age, pos=pos, tot=tot)

  class(model) <- "weibull_model"
  model
}


