sir_static <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
    ds <- -lambda*s
    di <- lambda*s - nu*i
    dr <- nu*i
    list(c(ds, di, dr))
  })
}

#' SIR static model (age-heterogeneous, endemic equilibrium)
#'
#' Refers to section 3.2.2.
#'
#' @details
#' In \code{state}:
#'
#'   - \code{s}: proportion susceptible
#'
#'   - \code{i}: proportion infected
#'
#'   - \code{r}: proportion recovered
#'
#' In \code{parameters}:
#'
#'   - \code{lambda}: natural death rate
#'
#'   - \code{nu}: recovery rate
#'
#' @param a age sequence.
#'
#' @param state the initial state of the system.
#'
#' @param parameters the model's parameter.
#'
#' @examples
#' state <- c(s=0.99,i=0.01,r=0)
#' parameters <- c(
#'   lambda = 0.05,
#'   nu=1/(14/365) # 2 weeks to recover
#' )
#' ages<-seq(0, 90, by=0.01)
#' model = sir_static_model(ages, state, parameters)
#' model
#'
#' @return list of class sir_static_model with the following items
#'   \item{parameters}{list of parameters used for fitting the model}
#'   \item{output}{matrix of proportion for each compartment over time}
#'
#'
#' @export
sir_static_model <- function(a, state, parameters)
{

  model <- list()
  model$parameters <- parameters
  model$output <- as.data.frame(ode(y=state,times=a,func=sir_static,parms=parameters))

  class(model) <- "sir_static_model"
  model
}
