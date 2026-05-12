#' The Farrington (1990) model.
#'
#' @description Fit age-stratified seroprevalence data using the Farrington (1990) model, which assumes the
#' force of infection increases linearly with age and subsequently decreases exponentially.
#'
#' @details
#' The force of infection is defined as followed
#'
#' \deqn{
#' \lambda(a) = (\alpha a - \gamma)e^{-\beta a} + \gamma
#' }
#' Where \eqn{\gamma} is called the long term residual for FOI,
#' as \eqn{a \rightarrow \infty} , \eqn{\lambda (a) \rightarrow \gamma}
#'
#' The seroprevalence can thus be estimated using the non-linear model
#' \deqn{
#'  \pi(a) = 1 - exp\{ \frac{\alpha}{\beta}ae^{-\beta a} +
#'  \frac{1}{\beta}(\frac{\alpha}{\beta} -
#'  \gamma)(e^{-\beta a} - 1) -\gamma a \}
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
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR `age`, `status` (for linelisting data)
#' @param start Named list of vectors or single vector.
#' Initial values for optimizer.
#' @param fixed Named list of vectors or single vector.
#' Parameter values to keep fixed during optimization.
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#' @return a list of class farrington_model with 5 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "mle" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [stats4::mle()] for more information on the fitted mle object
#'
#' @examples
#' df <- rubella_uk_1986_1987
#' model <- farrington_model(
#'   df,
#'   start=list(alpha=0.07,beta=0.1,gamma=0.03)
#'   )
#' plot(model)
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(data, start, fixed=list(),
                             age_col="age",pos_col="pos", tot_col="tot", status_col="status")
{
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  # model for seroprevalence in terms of age, alpha, beta and gamma
  seroprev_mod <- function(age, alpha, beta, gamma){
    1-exp(
      (alpha/beta)*age*exp(-beta*age)
      +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
      -gamma*age)
  }

  farrington <- function(alpha,beta,gamma) {
    p <- seroprev_mod(age, alpha, beta, gamma)
    # ll=pos*log(p)+(tot-pos)*log(1-p)
    # compute loglikelihood with dbinom instead
    ll <- dbinom(pos, size = tot, prob = p, log = TRUE)
    return(-sum(ll))
  }

  model$info <- mle(farrington, fixed=fixed, start=start)
  alpha <- model$info@coef[1]
  beta  <- model$info@coef[2]
  gamma <- model$info@coef[3]
  # functions to estimate seroprev and foi given age and parameters
  model$sp_mod <- seroprev_mod
  model$foi_mod <- function(age, alpha, beta, gamma){
    (alpha*age-gamma)*exp(-beta*age)+gamma
  }
  model$sp <- seroprev_mod(age, alpha, beta, gamma)
  model$foi <- model$foi_mod(age, alpha, beta, gamma)
  model$df <- data.frame(age=age, pos=pos, tot=tot)

  class(model) <- "farrington_model"
  model
}
