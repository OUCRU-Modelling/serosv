ds <- function(state, parameters, i)
{
  with(as.list(c(state, parameters)), {
    sum_beta_i <- 0
    for (j in 1:k) {
      sum_beta_i <- sum_beta_i + beta[i,j]*get(paste0("i", j))
    }
    -sum_beta_i*get(paste0("s", i)) + mu - mu*get(paste0("s", i))
  })
}

di <- function(state, parameters, i)
{
  with(as.list(c(state, parameters)), {
    sum_beta_i <- 0
    for (j in 1:k) {
      sum_beta_i <- sum_beta_i + beta[i,j]*get(paste0("i", j))
    }
    sum_beta_i*get(paste0("s", i)) - nu[i]*get(paste0("i", i)) - mu*get(paste0("i", i))
  })
}

dr <- function(state, parameters, i)
{
  with(as.list(c(state, parameters)), {
    nu[i]*get(paste0("i", i)) - mu*get(paste0("r", i))
  })
}

sir_subpop <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    s_states <- c()
    i_states <- c()
    r_states <- c()

    for (i in 1:k) {
      s_states <- c(s_states, ds(state, parameters, i))
      i_states <- c(i_states, di(state, parameters, i))
      r_states <- c(r_states, dr(state, parameters, i))
    }

    list(c(s_states, i_states, r_states))
  })
}

#' SIR Model with Interacting Subpopulations
#'
#' Refers to section 3.5.1.
#'
#' @details
#' In \code{state}:
#'
#' - \code{s}: Percent susceptible
#'
#' - \code{i}: Percent infected
#'
#' - \code{r}: Percent recovered
#'
#' In \code{parameters}:
#'
#' - \code{mu}: natural death rate (1/L).
#'
#' - \code{beta}: transmission rate w.r.t population (beta tilde)
#'
#' - \code{nu}: recovery rate
#'
#' - \code{k}: number of subpopulations
#'
#' @param times time sequence.
#'
#' @param state the initial state of the model.
#'
#' @param parameters the parameters of the model.
#'
#' @examples
#' \donttest{
#' k <- 2
#' state <- c(
#'   s = c(0.8, 0.8),
#'   i = c(0.2, 0.2),
#'   r = c(  0,   0)
#' )
#' beta_matrix <- c(
#'   c(0.05, 0.00),
#'   c(0.00, 0.05)
#' )
#' parameters <- list(
#'   beta = matrix(beta_matrix, nrow=k, ncol=k, byrow=TRUE),
#'   nu = c(1/30, 1/30),
#'   mu = 0.001,
#'   k = k
#' )
#' times<-seq(0,10000,by=0.5)
#' model <- sir_subpops_model(times, state, parameters)
#' model
#' }
#'
#' @return list of class sir_subpops_model with the following items
#'
#'   \item{parameters}{list of parameters used for fitting the model}
#'   \item{output}{matrix of proportion for each compartment over time}
#'
#'
#' @export
sir_subpops_model <- function(times, state, parameters) {
  model <- list()
  model$parameters <- parameters
  model$output <- as.data.frame(
    ode(y=state,times=times,func=sir_subpop,parms=parameters)
  )

  class(model) <- "sir_subpops_model"
  model
}
