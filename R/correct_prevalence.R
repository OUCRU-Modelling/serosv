#' Estimate the true sero prevalence using Bayesian estimation
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
#'
#' @return
#' a list of 2 items
#'   \item{info}{estimated parameters}
#'   \item{corrected_sero}{data.frame containing age, the corresponding estimated seroprevalance, adjusted tot and pos}
#' @export
#'
#' @examples
#' data <- rubella_uk_1986_1987
#' correct_prevalence(data)
correct_prevalence <- function(data, bayesian=TRUE,
                         init_se = 0.95, init_sp = 0.8, study_size_se = 1000, study_size_sp = 1000,
                         chains = 1, warmup = 1000, iter = 2000){
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

    # subsetting theta estimates (use mean value as corrected seroprevalence)
    thetas <- summary(fit)$summary[3:(length(age) + 2), "mean"]
  }else{
    # simply compute thetas using the formula
    thetas <- (pos/tot + init_sp - 1) / (init_se + init_sp - 1)
    thetas <- pmin(thetas, 1) # keep upper bound as 1
    thetas <- pmax(thetas, 0) # keep lower bound as 0
    # info for frequentist approach
    output$info <- "Formula: real_sero = (apparent_sero + sp - 1) / (se + sp -1)"
  }

  output$corrected_se <- data.frame(
    age = age,
    sero = thetas,
    pos = thetas*tot, #adjusted pos with estimated sero
    tot = tot
  )

  output
}
