#' Fit a mixture model to classify serostatus
#'
#' @param antibody_level - vector of the corresponding raw antibody level
#' @param breaks - number of intervals which the antibody_level are grouped into
#' @param pi - proportion of susceptible, infected
#' @param mu - a vector of means of component distributions (vector of 2 numbers in ascending order)
#' @param sigma-  a vector of standard deviations of component distributions (vector of 2 number)
#'
#' @importFrom mixdist mix mixgroup mixparam
#'
#' @export
#'
#' @examples
#' df <- vzv_be_2001_2003[vzv_be_2001_2003$age < 40.5,]
#' data <- df$VZVmIUml[order(df$age)]
#' model <- mixture_model(antibody_level = data)
#' model$info
mixture_model <- function (antibody_level, breaks=40, pi=c(0.2, 0.8), mu=c(2,6), sigma=c(0.5, 1)) {
  model <- list()

  # add 1 to avoid 0 when computing logs
  log_antibody <- log(antibody_level + 1)
  data <- mixgroup(log_antibody, breaks = breaks)
  starting_values <- mixparam(pi=pi,mu=mu,sigma=sigma)

  model$info <- mix(data,starting_values,dist="norm")
  model$df <- data.frame(antibody_level = data$X, count = data$count)
  model$susceptible <- fitted(model$info)$joint[,1]
  model$infected <- fitted(model$info)$joint[,2]

  class(model) <- "mixture_model"
  model
}
