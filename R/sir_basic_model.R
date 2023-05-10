# install.packages("deSolve")
require(deSolve)

sir_basic <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
     dS <- N*mu*(1-p) - beta*I*S + mu*S
     dI <- beta*I*S - (nu+alpha+mu)*I
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
#'  - \code{beta}: transmission rate (= lambda/I)
#'
#'  - \code{nu}: recovery rate
#'
#'  - \code{N}: total population
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
#'   mu=1/75,
#'   alpha=0,
#'   beta=0.0005,
#'   nu=1,
#'   N=5000,
#'   p=0
#' )
#' times <- seq(0, 250, by=0.1)
#' model <- sir_basic_model(times, state, parameters)
#' model
#'
#' @export
sir_basic_model <- function(times, state, parameters)
{
  as.data.frame(ode(y=state,times=times,func=sir_basic,parms=parameters))
}

#
# plot_SIR_basic <- function(df)
# {
#   par(mfrow=c(1,3),cex.axis=1.3,cex.lab=1.5,lwd=3,las=1,mgp=c(3, 0.5, 0))
#
#   plot(
#     df[,1],
#     df[,2]/(df[,2]+df[,3]+df[,4]),
#     xlab="time", ylab="proportion susceptible", type="l", lwd=2
#   )
#   plot(
#     df[,1],
#     df[,3]/(df[,2]+df[,3]+df[,4]),
#     xlab="time", ylab="proportion infected", type="l", lwd=2
#   )
#   plot(
#     df[,2],
#     df[,3],
#     xlab="number susceptible", ylab="number infected", type="l", lwd=2
#   )
# }
#
# state <- c(S=4999, I=1, R=0)
# parameters <- c(
#   mu=1/75,
#   alpha=0,
#   beta=0.0005,
#   nu=1,
#   N=5000,
#   p=0
# )
# times <- seq(0, 250, by=0.1)
# model <- get_SIR_basic(times, state, parameters)
# plot_SIR_basic(model)
