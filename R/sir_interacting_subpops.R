# dY1 <- -(beta11*Y2+beta12*Y5)*Y1+mu-mu*Y1
# dY4 <- -(beta21*Y2+beta22*Y5)*Y4+mu-mu*Y4
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

# dY2 <- (beta11*Y2+beta12*Y5)*Y1-v1*Y2-mu*Y2
# dY5 <- (beta21*Y2+beta22*Y5)*Y4-v2*Y5-mu*Y5
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

# dY3 <- v1*Y2 - mu*Y3
# dY6 <- v2*Y5 - mu*Y6
dr <- function(state, parameters, i)
{
  with(as.list(c(state, parameters, i)), {
    nu[i]*get(paste0("i", i)) - mu*get(paste0("r", i))
  })
}

SIR_interacting_subpops <- function(t, state, parameters) {
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

k <- 2
state <- c(
  s = c(0.8, 0.8),
  i = c(0.2, 0.2),
  r = c(  0,   0)
)
beta_matrix <- c(
  c(0.05, 0.00),
  c(0.00, 0.05)
)
parameters <- list(
  beta = matrix(beta_matrix, nrow=k, ncol=k, byrow=TRUE),
  nu = c(1/30, 1/30),
  mu = 0.001,
  k = k
)
times<-seq(0,10000,by=0.5)

sir_interacting_subpops <- as.data.frame(ode(y=state,times=times,func=SIR_interacting_subpops,parms=parameters))
head(sir_interacting_subpops)
summary(sir_interacting_subpops)

