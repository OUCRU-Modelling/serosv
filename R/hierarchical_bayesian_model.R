#' Hierarchical Bayesian Model
#'
#' Refers to section 10.3
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param type type of model ("far2", "far3" or "log_logistic")
#' @param chains number of Markov chains
#' @param warmup number of warmup runs
#' @param iter number of iterations
#'
#' @importFrom rstan sampling summary
#' @importFrom boot inv.logit
#'
#' @return a list of class hierarchical_bayesian_model with 6 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{type}{type of bayesian model far2, far3 or log_logistic}
#'   \item{info}{parameters for the fitted model}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @export
#'
#' @examples
#' \donttest{
#' df <- mumps_uk_1986_1987
#' model <- hierarchical_bayesian_model(df, type="far3")
#' model$info
#' plot(model)
#' }
hierarchical_bayesian_model <- function(data,
                            type="far3",chains = 1,warmup = 1500,iter = 5000){
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data)
  model$datatype <- data$type
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  # prepare data for stan model
  data <- list(age = age,
               posi = pos,
               ni = tot,
               Nage = as.numeric(length(age)))

  if (type == "far3"){
    # file <-  file.path(getwd(), "R", "stan_code", "fra_3.stan")
    fit <- rstan::sampling(stanmodels$fra_3, data=data, chains=chains,warmup=warmup, iter=iter)
  }
  else if (type == "far2"){
    # file <-  file.path(getwd(), "R", "stan_code", "fra_3.stan")
    fit <- rstan::sampling(stanmodels$fra_2, data=data, chains=chains,warmup=warmup, iter=iter)
  }
  else if (type == "log_logistic"){
    # file <-  file.path(getwd(), "R", "stan_code", "fra_3.stan")
    fit <- rstan::sampling(stanmodels$log_logistic, data=data, chains=chains,warmup=warmup, iter=iter)
  }
  else {
    stop('Model is not defined. Please choose "far3", "far2" or "log_logistic"')
  }

  model$info <- summary(fit)$summary

  theta <- data.frame(
    sp = rep (0,length(age)),
    foi = rep (0,length(age))
    )

  if (type == "far3"){
    alpha1 <- model$info["alpha1",c("mean")]
    alpha2 <- model$info["alpha2",c("mean")]
    alpha3 <- model$info["alpha3",c("mean")]

    for (i in 1:data$Nage){
      theta$sp[i] = 1-exp((alpha1/alpha2)*data$age[i]*exp(-alpha2*data$age[i])+
                            (1/alpha2)*((alpha1/alpha2)-alpha3)*(exp(-alpha2*data$age[i])-1)-
                            alpha3*data$age[i])
      theta$foi[i] = (alpha1*data$age[i]-alpha3)*exp(-alpha2*data$age[i])+alpha3
    }

  }

  if (type == "far2"){
    alpha1 <- model$info["alpha1",c("mean")]
    alpha2 <- model$info["alpha2",c("mean")]

    for (i in 1:data$Nage){
      theta$sp[i] = 1-exp((alpha1 / alpha2) * data$age[i] * exp(-alpha2 * data$age[i]) +
                            (1 / alpha2) * ((alpha1 / alpha2)) * (exp(-alpha2 * data$age[i]) - 1));
      theta$foi[i] = (alpha1*data$age[i])*exp(-alpha2*data$age[i])
    }

  }

  if (type == "log_logistic"){
    alpha1 <- model$info["alpha1",c("mean")]
    alpha2 <- model$info["alpha2",c("mean")]

    for (i in 1:data$Nage){
      theta$sp[i] = inv.logit(alpha2+alpha1*log(data$age[i]))
      theta$foi[i] = alpha1*exp(alpha2)*(data$age[i]^(alpha1-1))*(1-theta$sp[i])
    }
  }

  model$sp <- theta$sp
  model$foi <- theta$foi
  model$df <- data.frame(age = age, tot = tot, pos = pos)
  model$type <- type

  class(model) <- "hierarchical_bayesian_model"

  model
}
