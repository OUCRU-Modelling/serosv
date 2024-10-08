#' MSEIR model
#'
#' Refers to section 3.4.
#'
#' @param a age sequence
#' @param gamma time in maternal class.
#' @param lambda time in susceptible class.
#' @param sigma time in latent class.
#' @param nu time in infected class.
#'
#' @examples
#' model <- mseir_model(
#'   a=seq(from=1,to=20,length=500), # age range from 0 -> 20 yo
#'   gamma=1/0.5, # 6 months in the maternal antibodies
#'   lambda=0.2,  # 5 years in the susceptible class
#'   sigma=26.07, # 14 days in the latent class
#'   nu=36.5      # 10 days in the infected class
#' )
#' model
#'
#' @return list of class mseir_model with the following parameters
#'   \item{parameters}{list of parameters used for fitting the model}
#'   \item{output}{matrix of proportion for each compartment over time}
#'
#' @export
mseir_model <- function(a, gamma, lambda, sigma, nu)
{


  ma  <- exp(-gamma*a)
  sa  <- (gamma/(gamma-lambda))*(exp(-lambda*a)-exp(-gamma*a))
  ea  <- ((lambda*gamma)/(gamma-lambda))*
    (
      ((exp(-sigma*a)-exp(-lambda*a))/(lambda-sigma))
      -((exp(-sigma*a)-exp(-gamma*a))/(gamma-sigma))
    )
  ia  <- (sigma*lambda*gamma)*
    (
      ((exp(-nu*a)-exp(-sigma*a))/((lambda-sigma)*(gamma-sigma)*(sigma-nu)))
      +((exp(-nu*a)-exp(-lambda*a))/((lambda-gamma)*(lambda-sigma)*(lambda-nu)))
      +((exp(-nu*a)-exp(-gamma*a))/((gamma-lambda)*(gamma-sigma)*(gamma-nu)))
    )


  model <- list()
  model$parameters <- list(gamma = gamma, lambda = lambda, sigma = sigma, nu = nu)
  model$output <- data.frame(
    a =  c(0, a),
    m = c(1, ma),
    s = c(0, sa),
    e = c(0, ea),
    i = c(0, ia),
    r = c(0, 1 - ma - sa - ea - ia)
  )

  class(model) <- "mseir_model"
  model
}
