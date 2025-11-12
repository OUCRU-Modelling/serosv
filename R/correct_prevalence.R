#' Estimate the true sero prevalence using Frequentist/Bayesian estimation
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param bayesian whether to adjust sero-prevalence using the Bayesian or frequentist approach. If set to `TRUE`, true sero-prevalence is estimated using MCMC.
#' @param init_se sensitivity of the serological test
#' @param init_sp specificity of the serological test
#' @param study_size_se (applicable when `bayesian=TRUE`) study size for sensitivity validation study (i.e., number of confirmed infected patients in the study)
#' @param study_size_sp (applicable when `bayesian=TRUE`) study size for specificity validation study (i.e., number of confirmed non-infected patients in the study)
#' @param chains (applicable when `bayesian=TRUE`) number of Markov chains
#' @param warmup (applicable when `bayesian=TRUE`) number of warm up runs
#' @param iter (applicable when `bayesian=TRUE`) number of iterations
#'
#' @importFrom rstan sampling summary
#' @importFrom dplyr mutate
#' @importFrom stats prop.test
#' @importFrom magrittr %>%
#'
#' @return a list of 3 items
#'   \item{info}{estimated parameters (when `bayesian = TRUE`) or formula to compute corrected prevalence (when `bayesian = FALSE`)}
#'   \item{df}{data.frame of input data (in aggregated form)}
#'   \item{corrected_sero}{data.frame containing age, the corresponding estimated seroprevalance with 95\% confidence/credible interval, and adjusted tot and pos}
#' @export
#'
#' @examples
#' data <- rubella_uk_1986_1987
#' correct_prevalence(data)
correct_prevalence <- function(data, bayesian=TRUE,
                         init_se = 0.95, init_sp = 0.8, study_size_se = 1000, study_size_sp = 1000,
                         chains = 1, warmup = 1000, iter = 2000){
  # resolve no visible binding note
  lower <- upper <- prop.test <- NULL

  output <- list()

  data <- check_input(data)
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  if (data$type == "linelisting"){
    transform_df <- transform_data(age, pos)
    age <- transform_df$t
    pos <- transform_df$pos
    tot <- transform_df$tot
  }

  if (bayesian){

    # format data
    data <- list(
      posi = pos,
      ni = tot,
      init_se = init_se,
      init_sp = init_sp,
      study_size_se = study_size_se,
      study_size_sp = study_size_sp,
      Nage = length(age)
    )

    fit <- rstan::sampling(stanmodels$prevalence_correction, data=data, chains=chains,warmup=warmup, iter=iter)

    output$info <- summary(fit)$summary

    # subsetting theta estimates, also get the credible interval
    thetas <- data.frame(
      lower = summary(fit)$summary[3:(length(age) + 2), "2.5%"],
      fit = summary(fit)$summary[3:(length(age) + 2), "50%"],
      upper = summary(fit)$summary[3:(length(age) + 2), "97.5%"]
    )
  }else{
    thetas <- data.frame(
      # use prop.test to get confidence interval
      lower = mapply(\(pos, tot){prop.test(pos, tot, conf.level = 0.95)$conf.int[1]}, pos, tot),
      fit = pos/tot,
      upper = mapply(\(pos, tot){prop.test(pos, tot, conf.level = 0.95)$conf.int[2]}, pos, tot)
    ) %>%
    mutate(
      # estimate true prevalence and lower, upper bound
      lower = pmax(0, (lower - 1 + init_sp)/(init_se + init_sp - 1)),
      fit =  pmax(0, (fit + init_sp - 1) / (init_se + init_sp - 1)) %>% pmin(1),
      upper = pmin(1, (upper - 1 + init_sp)/(init_se + init_sp - 1))
    )

    # info for frequentist approach
    output$info <- "Formula: real_sero = (observed_sero + sp - 1) / (se + sp -1)"
  }

  output$df <- data.frame(
    age = age,
    pos = pos,
    tot = tot
  )

  output$corrected_se <- data.frame(
    age = age,
    sero = thetas$fit,
    pos = thetas$fit*tot, #adjusted pos with estimated sero
    tot = tot,
    sero_lwr = thetas$lower,
    sero_upr = thetas$upper
  )

  output$method <- if(bayesian) "bayesian" else "frequentist"

  output
}
