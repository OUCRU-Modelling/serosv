sir_basic <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
     N  <- S + I + R   # total population
     lambda <- beta*I  # force of infection
     dS <- N*mu*(1-p) - (lambda + mu)*S
     dI <- lambda*S - (nu+alpha+mu)*I
     dR <- N*mu*p + nu*I - mu*R
     list(c(dS, dI, dR))
  })
}

#' Basic SIR model
#'
#' Refers to section 3.1.3.
#'
#' @details
#' In \code{state}:
#'
#'  - \code{S}: number of susceptible
#'
#'  - \code{I}: number of infected
#'
#'  - \code{R}: number of recovered
#'
#' In \code{parameters}:
#'
#'  - \code{alpha}: disease-related death rate
#'
#'  - \code{mu}: natural death rate (= 1/life expectancy)
#'
#'  - \code{beta}: transmission rate
#'
#'  - \code{nu}: recovery rate
#'
#'  - \code{p}: percent of population vaccinated at birth
#'
#' @param times time sequence.
#'
#' @param state the initial state of the model.
#'
#' @param parameters the parameters of the model.
#'
#' @examples
#' state <- c(S=4999, I=1, R=0)
#' parameters <- c(
#'   mu=1/75, # 1 divided by life expectancy (75 years old)
#'   alpha=0, # no disease-related death
#'   beta=0.0005, # transmission rate
#'   nu=1, # 1 year for infected to recover
#'   p=0 # no vaccination at birth
#' )
#' times <- seq(0, 250, by=0.1)
#' model <- sir_basic_model(times, state, parameters)
#' model
#'
#' @importFrom deSolve ode
#' @return list of class sir_basic_model with the following items
#'   \item{parameters}{list of parameters used for fitting the model}
#'   \item{output}{matrix of population for each compartment over time}
#'
#'
#' @export
sir_basic_model <- function(times, state, parameters)
{
  model <- list()
  model$parameters <- parameters
  model$output <- as.data.frame(ode(y=state,times=times,func=sir_basic,parms=parameters))

  class(model) <- "sir_basic_model"
  model
}
