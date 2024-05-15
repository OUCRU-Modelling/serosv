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



