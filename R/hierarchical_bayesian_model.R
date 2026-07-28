#' Hierarchical Bayesian Model
#'
#'
#' @description Fit age-stratified seroprevalence to parametric hierarchical Bayesian models.
#' Supported models including Farrington model (2 and 3 parameters variants)
#' and Log Logistic model
#'
#' @details Consider a model for prevalence that has a parametric form
#' \eqn{\pi(a_i, \alpha)} where \eqn{\alpha} is a parameter vector
#'
#' Under a Bayesian framework, we can constraint the parameter space of the prior distribution \eqn{P(\alpha)}
#' to achieve monotonicity of the posterior distribution \eqn{P(\pi_1, \pi_2, ..., \pi_m|y,n)}
#'
#' Where:
#'
#' - \eqn{n = (n_1, n_2, ..., n_m)} and \eqn{n_i} is the sample size at age \eqn{a_i}
#'
#' - \eqn{y = (y_1, y_2, ..., y_m)} and \eqn{y_i} is the number of infected individual from the \eqn{n_i} sampled subjects
#'
#' For \bold{Farrington} model with 3 parameters, prevalence is formulated as follow
#'
#' \deqn{
#' \pi (a) = 1 - exp\{ \frac{\alpha_1}{\alpha_2}ae^{-\alpha_2 a} +
#' \frac{1}{\alpha_2}(\frac{\alpha_1}{\alpha_2} - \alpha_3)(e^{-\alpha_2 a} - 1) -\alpha_3 a \}
#' }
#'
#' The likelihood model is defined as \eqn{y_i \sim Bin(n_i, \pi_i), \text{  for } i = 1,2,3,...m}
#'
#' The constraint on the parameter space can be incorporated by assuming
#'  truncated normal distribution for the components of \eqn{\alpha},
#' \eqn{\alpha = (\alpha_1, \alpha_2, \alpha_3)} in \eqn{\pi_i = \pi(a_i,\alpha)}
#'
#' The flat hyperpriors are defined as \eqn{\mu_j \sim \mathcal{N}(0, 10000)} and
#' \eqn{\tau^{-2}_j \sim \Gamma(100,100)}
#'
#' For \bold{Farrington} model with 2 parameters, it is equivalent to the previous model with \eqn{\alpha_3 = 0}
#'
#' For \bold{Log logistic model}, seroprevalence is instead defined as
#'
#' \deqn{\pi(a) = \frac{\beta a^\alpha}{1 + \beta a^\alpha}, \text{ } \alpha, \beta > 0}
#'
#' The likelihood is similarly defined as \eqn{y_i \sim Bin(n_i, \pi_i))}
#'
#' The prior model of \eqn{\alpha_1} is specified as \eqn{\alpha_1 \sim \text{truncated  } \mathcal{N}(\mu_1, \tau_1)}
#' with flat hyperpriors as in Farrington model
#'
#' \eqn{\beta} is constrained to be positive by specifying \eqn{\alpha_2 \sim \mathcal{N}(\mu_2, \tau_2)}
#'
#' Refer to section Chapter 10.3 of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR `age`, `status` (for linelisting data)
#' @param type type of model ("far2", "far3" or "log_logistic")
#' @param chains number of Markov chains
#' @param warmup number of warmup runs
#' @param iter number of iterations
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#' @importFrom rstan sampling summary
#' @importFrom boot inv.logit
#'
#' @return a list of class hierarchical_bayesian_model with the following items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{type}{type of bayesian model "far2", "far3" or "log_logistic"}
#'   \item{info}{a stanfit object for the fitted result}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'   \item{sp_func}{function to compute seroprevalence given age and model parameters}
#'   \item{foi_func}{function to compute force of infection given age and model parameters}
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
                            age_col="age",pos_col="pos", tot_col="tot", status_col="status",
                            type="far3",chains = 1,warmup = 1500,iter = 5000){
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
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

  model$info <- fit

  theta <- list()

  if (type == "far3"){
    alpha1 <- summary(model$info)$summary["alpha1",c("mean")]
    alpha2 <- summary(model$info)$summary["alpha2",c("mean")]
    alpha3 <- summary(model$info)$summary["alpha3",c("mean")]

    theta$sp_func <- \(age, alpha1, alpha2, alpha3){
      1 - exp((alpha1/alpha2)*age*exp(-alpha2*age)+
              (1/alpha2)*((alpha1/alpha2)-alpha3)*(exp(-alpha2*age)-1)-
              alpha3*age)
    }

    theta$foi_func <- \(age, alpha1, alpha2, alpha3){
      (alpha1*age-alpha3)*exp(-alpha2*age)+alpha3
    }

    theta$sp <-  theta$sp_func(data$age, alpha1, alpha2, alpha3)
    theta$foi <-  theta$foi_func(data$age, alpha1, alpha2, alpha3)

  }

  if (type == "far2"){
    alpha1 <- summary(model$info)$summary["alpha1",c("mean")]
    alpha2 <- summary(model$info)$summary["alpha2",c("mean")]

    theta$sp_func <- \(age, alpha1, alpha2){
      1-exp((alpha1 / alpha2) * age * exp(-alpha2 * age) +
              (1 / alpha2) * ((alpha1 / alpha2)) * (exp(-alpha2 * age) - 1))
    }
    theta$foi_func <- \(age, alpha1, alpha2){
      (alpha1*age)*exp(-alpha2*age)
    }

    theta$sp <- theta$sp_func(data$age, alpha1, alpha2)
    theta$foi <- theta$foi_func(data$age, alpha1, alpha2)

  }

  if (type == "log_logistic"){
    alpha1 <- summary(model$info)$summary["alpha1",c("mean")]
    alpha2 <- summary(model$info)$summary["alpha2",c("mean")]

    theta$sp_func <- \(age, alpha1, alpha2){
      inv.logit(alpha2+alpha1*log(age))
    }
    theta$foi_func <- \(age, sp, alpha1, alpha2){
      alpha1*exp(alpha2)*(age^(alpha1-1))*(1-sp)
    }

    theta$sp <- theta$sp_func(data$age, alpha1, alpha2)
    theta$foi <- theta$foi_func(data$age, theta$sp, alpha1, alpha2)
  }

  model$sp <- theta$sp
  model$foi <- theta$foi
  model$df <- data.frame(age = age, tot = tot, pos = pos)
  model$type <- type
  model$sp_func <- theta$sp_func
  model$foi_func <- theta$foi_func

  class(model) <- "hierarchical_bayesian_model"

  model
}
