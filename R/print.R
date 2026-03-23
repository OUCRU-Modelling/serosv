#' @export
print.polynomial_model <- function(x, ...){
  cat("Polynomial model\n\n")
  cat("Input type: ", x$datatype, "\n")
  cat("Degree (k): ",x$k, "\n")
  print(x$info)
}

#' @export
print.fp_model <- function(x, ...){
  cat("Fractional polynomial model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat("Powers: ", paste(x$p, collapse=", "), "\n")
  print(x$info)
}

#' @export
print.age_time_model <- function(x, ...){
  cat("Age-time varying seroprevalence model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat("Grouping variable: ", x$grouping_col,"\n")
  cat("Monotonization method: ", x$monotonize_method, "\n")
  cat("Monotonize across: ", if(x$age_correct) "birth cohort" else "age group", "\n")
  print(x$out)
}

#' @export
print.farrington_model <- function(x, ...){
  cat("Farrington model \n\n")
  cat("Input type: ", x$datatype, "\n")
  # cat("Fitted parameters: ",
  #     paste(c("alpha", "beta", "gamma"), sprintf("%.4g", x$info@coef[1:3]), sep="=", collapse=", "),
  #     "\n\n")
  print(x$info)
}

#' @importFrom purrr compact
#' @export
print.hierarchical_bayesian_model <- function(x, ...){
  pars <- intersect(c("alpha1", "alpha2", "alpha3"), rownames(x$info))
  fitted_pars <- x$info[pars,  "mean"]
  # get the CrI
  lower_pars <- x$info[pars,  "2.5%"]
  upper_pars <- x$info[pars,  "97.5%"]
  # get the sd
  sd_pars <- purrr::compact(x$info[pars,  "sd"])

  cat("Hierarchical Bayesian model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat("Model: ",
      switch(x$type,
             far2 = "Farrington model with 2 parameters",
             far3 = "Farrington model with 3 parameters",
             log_logistic = "Log-logistic model"),
      "\n\n")
  cat("Fitted parameters:\n",
      paste(names(fitted_pars),
            # print fitted parameter in the format mean (95% CrI=lower-upper, sd=)
            paste(
              sprintf("%.4g", fitted_pars),
              " (95% CrI [",
              sprintf("%.4g", lower_pars),
              ", ",
              sprintf("%.4g", upper_pars),
              "], sd = ",
              sprintf("%.4g", sd_pars),
              ")",
              sep = ""
            ),
            sep=" = ", collapse="\n "),
      "\n")
}

#' @export
print.mixture_model <- function(x, ...){
  cat("Mixture model \n\n")
  cat("Estimated proportion:\n",
      paste(
        c("Susceptible", "Infected"),
        sprintf("%.4g", x$info$parameters$pi),
        sep="=", collapse = ", "
      ),
      "\n\n"
      )
  cat("Estimated mean Log(Antibody):\n",
      paste(
        c("Susceptible", "Infected"),
        sprintf("%.4g", x$info$parameters$mu),
        sep="=", collapse = ", "
      ),
      "\n")
}

#' @export
print.estimate_from_mixture <- function(x, ...){
  cat("Age-varying seroprevalence estimated from mixture model \n\n")
  cat("Monotonized seroprevalence: ", x$monotonize)
  print(x$info)
}

#' @export
print.lp_model <- function(x, ...){
  cat("Local polynomial model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat(
    "Configs: ",
    paste(
      c("nn", "bandwidth(h)", "degree", "kernel"),
      c(as.character(round(
        as.numeric(x[c("nn", "h", "deg")]), digits = 4
      )), x$kern),
      sep = "=",
      collapse = ", "
    ),
    "\n\n"
  )
  print(x$info)
}

#' @export
print.penalized_spline_model <- function(x, ...){
  cat("Penalized spline model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat("Framework: ", if(x$framework=="pl") "Penalized likelihood" else "Mixed model", "\n")
  print(x$info)
}

#' @export
print.weibull_model <- function(x, ...){
  cat("Weibull model \n\n")
  cat("Input type: ", x$datatype, "\n")
  cat(
    paste(
      c("b0", "b1"),
      sprintf("%.4g", coef(x$info)[1:2]),
      sep="=", collapse=", "
    ),
    "\n"
  )

  print(x$info)
}

