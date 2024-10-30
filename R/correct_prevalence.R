#' Estimate the true sero prevalence using Bayesian estimation
#'
#' @param age the age vector
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
#' @param init_se sensitivity of the serological test
#' @param init_sp specificity of the serological test
#' @param study_size_se study size for sensitivity validation study (i.e., number of confirmed infected patients in the study)
#' @param study_size_sp study size for specificity validation study (i.e., number of confirmed non-infected patients in the study)
#' @param chains number of Markov chains
#' @param warmup number of warm up runs
#' @param iter number of iterations
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
#' correct_prevalence(data$age, pos = data$pos, tot = data$tot)
correct_prevalence <- function(age, pos=NULL, tot=NULL, status=NULL,
                         init_se = 0.95, init_sp = 0.8, study_size_se = 1000, study_size_sp = 1000,
                         chains = 1, warmup = 1000, iter = 2000){
  output <- list()
  if (!is.null(pos) & !is.null(tot)){
    pos <- as.numeric(pos)
    tot <- as.numeric(tot)
  }else{
    # automatically transform to aggregated data if line listing data is given
    transform_df <- transform_data(age, status)
    age <- transform_df$t
    pos <- as.numeric(transform_df$pos)
    tot <- as.numeric(transform_df$tot)
  }

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

  # subsetting theta estimates
  thetas <- summary(fit)$summary[3:(length(age) + 2), "mean"]
  output$corrected_se <- data.frame(
    age = age,
    sero = thetas,
    pos = thetas*tot, #adjusted pos with estimated sero
    tot = tot
  )

  output
}
