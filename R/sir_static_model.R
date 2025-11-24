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
#' @description Simulate transmission model with constant force of infection and consists of 3 compartments: susceptible (S), infected (I), recovered (R)
#'
#' @details
#' Follow the SIR model at endemic state described in the book by Hens et al. (section 3.2.2.)
#'
#' Assumptions:
#'
#' - Time homogeneity
#'
#' - Age heterogeneity
#'
#' - Constant force of infection
#'
#' The model is described by a system of 3 differential equations
#'
#' \deqn{
#' \begin{cases}
#' \frac{ds(a)}{da} = -\lambda s(a) \\
#' \frac{di(a)}{da} = \lambda s(a) - \nu i(a)  \\
#' \frac{dr(a)}{da} =  \nu i(a)
#' \end{cases}
#' }
#'
#' Where:
#'
#' -   \eqn{s(a), i(a), r(a)} are proportion of susceptible, infected, recovered population of age group \eqn{a} respectively
#'
#' -   \eqn{\lambda} is the force of infection
#'
#' -   \eqn{\nu} is the recovery rate
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' @param a age sequence.
#'
#' @param state the initial state of the system. A named vector with the following
#'   \itemize{
#'     \item \code{s}: initial susceptible proportion
#'     \item \code{i}: initial infected proportion
#'     \item \code{r}: initial recovered proportion
#'   }
#'
#' @param parameters the model's parameter. A named vector with the following
#'  \itemize{
#'     \item \code{lambda}: natural death rate
#'     \item \code{nu}: recovery rate
#'  }
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
