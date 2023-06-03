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
#' @param times time sequence.
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
#' @export
sir_static_model <- function(times, state, parameters)
{
  as.data.frame(ode(y=state,times=times,func=sir_static,parms=parameters))
}
