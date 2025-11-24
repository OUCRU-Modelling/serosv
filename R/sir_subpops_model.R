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
#' @description An extension of the basic SIR model that incorporates interaction between sub-populations
#'
#' @details
#' Follow the SIR model with sub populations described in the book by Hens et al. (section 3.5.1.)
#'
#' With K subpopulations, the WAIFW matrix or mixing matrix is given by
#'
#' \deqn{
#' C = \begin{bmatrix}
#' \beta_{11} & \beta_{12}  & ... & \beta_{1K} \\
#' \beta_{21} & \beta_{22}  & ... & \beta_{2K} \\
#' \vdots & \vdots  & ... & \vdots \\
#' \beta_{K1} & \beta_{K2}  & ... & \beta_{KK} \\
#' \end{bmatrix}
#' }
#'
#' And the \eqn{i^{th}} sub population is described by the following system of differential equations
#' \deqn{
#' \begin{cases}
#' \frac{dS_i(t)}{dt} = -(\sum^K_{j=1}\beta_{ij}I_j(t)) S_i(t) + N_i\mu_i - \mu_i S_i(t) \\
#' \frac{dI_i(t)}{dt} = (\sum^K_{j=1}\beta_{ij}I_j(t)) S_i(t)  - (\nu_i + \mu_i) I_i(t)  \\
#' \frac{dR_i(t)}{dt} = \nu_i I_i(t)  - \mu_i R_i(t)
#' \end{cases}
#' }
#'
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' @param times time sequence.
#'
#' @param state the initial state of the model. A named vector with the following
#' \itemize{
#'    \item \code{s}: initial susceptible proportion
#'    \item \code{i}: initial infected proportion
#'    \item \code{r}: initial recovered proportion
#'  }
#'
#' @param parameters the parameters of the model. A named vector with the following
#' \itemize{
#'    \item \code{mu}: natural death rate (1/L).
#'    \item \code{beta}: the WAIFW matrix, with dimensions \code{[K, K]}.
#'    \item \code{nu}: recovery rate
#'  }
#'
#' @examples
#' \donttest{
#' state <- c(
#'   s = c(0.8, 0.8),
#'   i = c(0.2, 0.2),
#'   r = c(  0,   0)
#' )
#' beta_matrix <- matrix(
#'   c(0.05, 0.00,
#'   0.00, 0.05),
#'   2
#' )
#' parameters <- list(
#'   beta = beta_matrix,
#'   nu = c(1/30, 1/30),
#'   mu = 0.001
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
  # Infer k here, beta should be in the dimension [k, k]
  parameters$k <- dim(parameters$beta)[1]

  model$parameters <- parameters
  model$output <- as.data.frame(
    ode(y=state,times=times,func=sir_subpop,parms=parameters)
  )

  class(model) <- "sir_subpops_model"
  model
}
