compute_ci <- function(x, ci = 0.95, ...){
  UseMethod("compute_ci")
}

#' Compute confidence interval
#'
#' @param x - serosv models
#' @param ci - confidence interval
#' @param le - number of data for computing confidence interval
#' @param ... - arbitrary argument
#'
#' @export
compute_ci <- function(x, ci = 0.95, le = 100, ...){
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- x$info$data
  n <- nrow(dataset) - length(x$info$coefficients)
  age_range <- range(dataset$Age)
  ages <- seq(age_range[1], age_range[2], le = le)

  mod1 <- predict.glm(x$info,data.frame(Age = ages), se.fit = TRUE)
  n1 <- mod1 |> extract(c("fit", "se.fit")) %>%
    c(age = list(ages), .) %>%
    as_tibble() |>
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)

  out.DF <- data.frame(x = n1$age, y = 1- n1$fit, ymin= 1-  n1$lwr, ymax=1- n1$upr)
}

#' Compute confidence interval for fractional polynomial model
#'
#' @param x - serosv models
#' @param ci - confidence interval
#' @param le - number of data for computing confidence interval
#' @param ... - arbitrary argument
#'
#' @export
compute_ci.fp_model <- function(x, ci = 0.95, le = 100, ...){
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- data.frame(x$df)
  n <- nrow(dataset) - length(x$info$coefficients)
  age_range <- range(dataset$age)
  ages <- seq(age_range[1], age_range[2], le = le)

  mod1 <- predict.glm(x$info,data.frame(age = ages), se.fit = TRUE)
  n1 <- data.frame(mod1)[,-3] %>%
    c(age = list(ages), .) %>%
    as_tibble() |>
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)
  out.DF <- data.frame(x = n1$age, y = n1$fit,
                       ymin= n1$lwr, ymax= n1$upr)
}

#' Compute confidence interval for Weibull model
#'
#' @param x - serosv models
#' @param ci - confidence interval
#' @param ... - arbitrary argument
#'
#' @export
compute_ci.weibull_model <- function(x, ci = 0.95, ...){
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- x$info$model
  n <- nrow(dataset) - length(x$info$coefficients)
  age_range <- range(dataset$`log(t)`)
  exposure_time <- dataset$`log(t)`

  mod1 <- predict.glm(x$info,data.frame("log(t)" = exposure_time), se.fit = TRUE)
  n1 <- mod1 |> extract(c("fit", "se.fit")) %>%
    c(exposure = list(exposure_time), .) %>%
    as_tibble() |>
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)

  out.DF <- data.frame(x = x$df$t, y = n1$fit,
                       ymin= n1$lwr, ymax= n1$upr)
}


#' Compute confidence interval for local polynomial model
#'
#' @param x - serosv models
#' @param ci - confidence interval
#' @param ... - arbitrary arguments
#' @export
compute_ci.lp_model <- function(x,ci = 0.95, ...){
  ages <- x$df$age
  crit<- crit(x$pi,cov = ci)$crit.val
  mod1 <- predict(x$pi, data.frame(a = ages),se.fit = TRUE)
  out.DF <- data.frame(x = ages, y = mod1$fit,ymin= mod1$fit-crit*(mod1$se.fit/100),
                       ymax= mod1$fit+crit*(mod1$se.fit/100))
  out.DF
}


#' Compute confidence interval for GLMM
#'
#' @param x - serosv models
#' @param ci - confidence interval
#' @param ... - arbitrary arguments
#' @export
compute_ci.glmm_ps_model <- function(x,ci = 0.95, ...){
  m <- 1

  p <- (1 - ci) / 2
  link_inv <- x$info$gam$family$linkinv
  dataset <- x$info$gam$model[,1:2]
  n <- nrow(dataset) - length(x$info$gam$coefficients)
  ages <- dataset[2]

  mod <- predict.gam(x$info$gam, data.frame(a = ages), se.fit = TRUE) |>
    extract(c("fit", "se.fit")) %>%
    c(age = list(ages), .) |>
    as_tibble() |>
    mutate(lwr = m * link_inv(fit + qt(    p, n) * se.fit),
           upr = m * link_inv(fit + qt(1 - p, n) * se.fit),
           fit = m * link_inv(fit)) |>
    select(- se.fit)

  out.DF <- data.frame(x = ages, y = mod$fit,
                       ymin= mod$lwr, ymax = mod$upr)
  out.FOI <- data.frame(x = x$df$age[c(-1, -length(x$df$age) )],
                        y = est_foi(ages[[1]], mod$fit),
                        ymin= est_foi(ages[[1]],mod$lwr),
                        ymax = est_foi(ages[[1]],mod$upr)
                        )

  return(list(out.DF, out.FOI))
}

#' Compute confidence interval for mixture model
#'
#' @param x - serosv mixture_model object
#' @param ci - confidence interval
#' @param ... - arbitrary arguments
#' @export
compute_ci.mixture_model <- function(x,ci = 0.95, ...){

  susceptible <- x$info$parameters[1, ]
  infected <- x$info$parameters[2, ]

  lower_q <- (1 - ci)/2
  upper_q <- 1 - lower_q
  susceptible <- list(lower_bound = qnorm(lower_q, mean = susceptible$mu, sd = susceptible$sigma),
                      upper_bound = qnorm(upper_q, mean = susceptible$mu, sd = susceptible$sigma))

  infected <- list(lower_bound = qnorm(lower_q, mean = infected$mu, sd = infected$sigma),
                      upper_bound = qnorm(upper_q, mean = infected$mu, sd = infected$sigma))

  return(list(susceptible= susceptible, infected=infected))
}



