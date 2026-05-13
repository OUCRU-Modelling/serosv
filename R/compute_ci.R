compute_ci <- function(x, ci = 0.95, le = 100, ...){
  UseMethod("compute_ci")
}

# ======= Quantify CI helper functions =======
#' @param foi_func a function that takes coefficents, newdat, and return estimated FOI
#' @param newdat new age-range to generate FOI
#' @param coef estimated coefficients
#' @param vcov variance-covariance matrix of coefficients
#' @param nb number of bootstrap iterations
#' @param alpha significance level
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats quantile setNames
#' @importFrom purrr map_dfc
parametric_bootstrapping <- function(foi_func,
                                     newdat,
                                     coef, vcov,
                                     nb=9999, alpha=.025){
  rowsplit <- function(df) split(df, 1:nrow(df))

  # Sample parameters value
  sampling_out <- rmvnorm(
    nb,
    mean = coef,
    sigma = vcov
  ) %>%
    as.data.frame() %>%
    rowsplit() %>%
    map(as.numeric)

  # Get the bounds for FOI
  sampling_out %>%
    map_dfc(~foi_func(newdat, .x)) %>%
    apply(1, quantile, c(alpha, 1 - alpha)) %>%
    t() %>% as.data.frame() %>%
    cbind(foi_func(newdat, coef)) %>%
    setNames(c("ymin", "ymax","y")) %>%
    cbind(newdat)
}

# TODO: General function for parametric delta method

#' @param mod fitted serosv model
#' @param refit_func function to refit the model
#' @param newdat new age-range for estimating FoI
#' @param nb number of bootstrap iterations
#' @param ci confidence level for the interval
#'
#' @importFrom boot boot boot.ci
#' @importFrom stats setNames
#' @importFrom map map_dfr
nonparametric_bootstrapping <- function(mod, refit_func,
                                        newdat,
                                        nb=200, ci=.95){
  message("Running nonparametric bootstrap for FoI confidence intervals, this may take a while")
  dat <- mod$df

  # dat is the data returned by ran.gen (i.e., resampled data)
  stat_func <- function(resampled_dat, indices=NULL,refit_func, newdat){
    # refit model
    refit_mod <- refit_func(resampled_dat)

    # return foi estim over the age range we're interested in
    sp_pred <- predict(refit_mod, newdat)
    foi_pred <- est_foi(newdat[[1]], sp_pred)

    foi_pred
  }

  resample_func <- function(dat, mod){
    resampled_dat <- dat
    # check fitted datatype
    if(mod$datatype == "aggregated"){
      # aggregated -> resample success count
      resampled_dat$pos <- rbinom(nrow(dat),
                                 size=resampled_dat$tot, prob=resampled_dat$pos/resampled_dat$tot)
    }else{
      # linelisting -> resample rows
      indices <- sample(1:nrow(dat), nrow(dat), replace = TRUE)
      resampled_dat <- resampled_dat[indices, ]
    }

    resampled_dat
  }

  boot_out <- boot::boot(
    dat, statistic = stat_func, R = nb,
    # specify custom ran.gen function to implement unsupported nonparametric
    # bootstrapping scheme
    sim = "parametric", ran.gen = resample_func,
    parallel = "multicore", ncpus = parallel::detectCores(),
    # argument for resample_func
    mle = mod,
    # arguments for stat_func
    refit_func = refit_func,
    newdat = newdat
  )

  setNames(
    map_dfr(
      1:ncol(boot_out$t),
      \(i) {
        ci_out <- boot::boot.ci(boot_out, type = "perc", index = i, conf=ci)
        as.data.frame(ci_out$percent)[,4:5]
      }
    ),
    c("ymin", "ymax")
  ) |>
    bind_cols(
      # adjust length after numerical differentiation (for FoI)
      newdat[c(-1, -nrow(newdat)), ,drop=FALSE]
    )
}

#' Compute confidence interval for a model of serosv
#'
#' Computes CI for Seroprevalence from model standard errors, and (optionally)
#' for Force of Infection via parametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param le number of data for computing confidence interval
#' @param foi_ci whether to compute CI for FoI
#' @param ... arbitrary argument
#'
#' @importFrom stats qt predict.glm coef vcov
#' @import dplyr
#'
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#'
#' @export
compute_ci.default <- function(x, foi_ci=TRUE, ci = 0.95, le = 100, ...){
  # resolve no visible binding issue with CRAN check
  fit <- se.fit <- NULL

  # Set up
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- x$info$data
  n <- nrow(dataset) - length(x$info$coefficients)
  age_range <- range(dataset$age)
  ages <- seq(age_range[1], age_range[2], le = le)

  # Quantify CI of seroprevalence estimate
  mod1 <- predict.glm(x$info,data.frame(age = ages), se.fit = TRUE)
  n1 <- mod1 %>% as_tibble() %>%  select(fit, se.fit) %>%
    mutate(age = ages ) %>%
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)
  out.DF <- data.frame(x = n1$age, y = 1- n1$fit, ymin= 1-  n1$lwr, ymax=1- n1$upr)

  # Quantify CI of FoI if specified
  out.FOI <- if(foi_ci){
    parametric_bootstrapping(
      # make sure foi_func match the expected function signature
      foi_func = \(newdat, coef){
        x$foi_mod(newdat$x, coef)
      },
      newdat = data.frame(x = ages),
      coef = coef(x$info), vcov = vcov(x$info),
      alpha = p
    )
  }else{
    data.frame(
      x = ages,
      y = x$foi_mod(ages, x$info$coefficients)
    )
  }

  list(out.DF, out.FOI)
}



# ======= Parametric model ===========
#' Compute confidence interval for fractional polynomial model
#'
#' Computes CI for Seroprevalence from model standard errors, and (optionally)
#' for Force of Infection via nonparametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param foi_ci whether to compute CI for FoI
#' @param le number of data for computing confidence interval
#' @param ... arbitrary argument
#'
#' @import dplyr
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#'
#' @export
compute_ci.fp_model <- function(x, ci = 0.95, le = 100, foi_ci=FALSE, ...){
  # resolve no visible binding issue with CRAN check
  fit <- se.fit <- NULL

  # Set up
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- data.frame(x$df)
  n <- nrow(dataset) - length(x$info$coefficients)
  age_range <- range(dataset$age)
  ages <- seq(age_range[1], age_range[2], le = le)

  # Quantify CI of seroprevalence estimation
  mod1 <- predict.glm(x$info,data.frame(age = ages), se.fit = TRUE)
  n1 <- data.frame(mod1)[,-3] %>%
    mutate(age = ages) %>%
    as_tibble() %>%
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)
  out.DF <- data.frame(x = n1$age, y = n1$fit,
                       ymin= n1$lwr, ymax= n1$upr)

  # Quantify CI of FOI if specified
  foi_x <- sort(unique(ages))
  foi_x <- foi_x[c(-1, -length(foi_x) )]

  out.FOI <- data.frame(
    x = foi_x,
    y = est_foi(ages, out.DF$y)
  )

  out.FOI <- if(foi_ci){
    bootstrap_out <- nonparametric_bootstrapping(
      mod = x, newdat = data.frame(age = ages),
      refit_func = \(df){
        do.call(
          fp_model,
          c(
            list(
              data = df,
              p = x$p
            ),
            x$pars
          )
        )
      },
      ci=ci, ...
    )
    cbind(out.FOI, bootstrap_out)
  }else{
    out.FOI
  }

  list(out.DF, out.FOI)
}

#' Compute confidence interval for Weibull model
#'
#' Computes CI for Seroprevalence from model standard errors, and (optionally)
#' for Force of Infection via parametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param foi_ci whether to compute CI for FoI
#' @param ... arbitrary argument
#'
#' @importFrom stats vcov coef
#' @import dplyr
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#' @export
compute_ci.weibull_model <- function(x, ci = 0.95, foi_ci=TRUE, ...){
  # resolve no visible binding issue with CRAN check
  fit <- se.fit <- NULL

  # set up
  p <- (1 - ci) / 2
  link_inv <- x$info$family$linkinv
  dataset <- data.frame(x$df)
  n <- nrow(dataset) - length(x$info$coefficients)

  mod1 <- predict.glm(x$info,data.frame(t = dataset$age), se.fit = TRUE)
  n1 <- mod1 %>% as_tibble() %>%
    select(fit, se.fit) %>%
    mutate(t = dataset$age) %>%
    mutate(lwr = link_inv(fit + qt(    p, n) * se.fit),
           upr = link_inv(fit + qt(1 - p, n) * se.fit),
           fit = link_inv(fit)) %>%
    select(-se.fit)

  out.DF <- data.frame(x = dataset$age, y = n1$fit,
                       ymin= n1$lwr, ymax= n1$upr)

  # estimate FOI CI if specified
  out.FOI <- if(foi_ci){
    parametric_bootstrapping(
      # make sure foi_func match the expected function signature
      foi_func = \(newdat, coef){
        x$foi_mod(newdat$x, coef[1], coef[2])
      },
      newdat = data.frame(x = dataset$age),
      coef = coef(x$info), vcov = vcov(x$info),
      alpha = p
    )
  }else{
    data.frame(
      x = ages,
      y = x$foi_mod(ages, coef(x$info)[1], coef(x$info)[2])
    )
  }

  list(out.DF, out.FOI)
}

#' Compute confidence interval for Farrington model
#'
#' Computes CI for Seroprevalence and (optionally)
#' for Force of Infection via parametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param nb number of samples for parametric bootstrapping
#' @param foi_ci whether to compute CI for FoI
#' @param ... arbitrary argument
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom purrr map_dfc
#' @importFrom stats setNames vcov coef quantile formula
#'
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#'
#' @export
compute_ci.farrington_model <- function(x, ci = 0.95, nb=9999, foi_ci=TRUE,...){
  rowsplit <- function(df) split(df, 1:nrow(df))

  mod <- x$info
  age <- unique(x$df$age)

  alpha <- (1-ci)/2

  # CIs for seroprevalence and FoI are quantified using parametric bootstrapping
  sampling_out <-  rmvnorm(
    nb,
    mean = mod@coef,
    sigma = mod@vcov
  ) %>%
    as.data.frame() %>%
    rowsplit() %>%
    map(as.list) %>%
    map(~ c(.x, data.frame(age = age)))

  # ----- Estimate CI for seroprevalence
  out.DF <- sampling_out %>%
    map_dfc(~do.call(x$sp_mod, .x)) %>%
    apply(1, quantile, c(alpha, 1 - alpha)) %>%
    t() %>% as.data.frame() %>%
    setNames(c("ymin", "ymax")) %>%
    cbind(
      data.frame(
        x = age,
        # use the estimated parameter to compute estimated seroprev
        y = do.call(x$sp_mod, c(
          list(age=age),
          mod@coef
        ))
      )
    )

  # ----- Estimate CI for FOI
  out.FOI <- data.frame(
    x = age,
    # use the estimated parameter to compute estimated FOI
    y = do.call(x$foi_mod, c(
      list(age=age),
      mod@coef
    ))
  )
  # compute CI if specified
  out.FOI <- if(foi_ci){
    sampling_out %>%
      map_dfc(~do.call(x$foi_mod, .x)) %>%
      apply(1, quantile, c(alpha, 1 - alpha)) %>%
      t() %>% as.data.frame() %>%
      setNames(c("ymin", "ymax")) %>%
      cbind(out.FOI)
  }else{
    out.FOI
  }

  list(out.DF, out.FOI)
}

#' Compute 95\% credible interval for hierarchical Bayesian model
#'
#' Return CrI for Seroprevalence and Force of Infection via parameters' posterior distributions.
#'
#' @param x serosv models
#' @param ... arbitrary arguments
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper credible interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), \code{ymin} and
#'       \code{ymax} (lower and upper credible interval bounds)
#' @export
compute_ci.hierarchical_bayesian_model <- function(x, ...){
  out_x <- x$df$age
  out.DF <- NULL
  out.FOI <- NULL

  if (x$type == "far3"){
    alpha1 <- x$info["alpha1",c("2.5%","50%", "97.5%")]
    alpha2 <- x$info["alpha2",c("2.5%","50%", "97.5%")]
    alpha3 <- x$info["alpha3",c("2.5%","50%", "97.5%")]

    out.DF <- data.frame(
      x = out_x,
      ymin = x$sp_func(out_x, alpha1[1], alpha2[1], alpha3[1]),
      y = x$sp_func(out_x, alpha1[2], alpha2[2], alpha3[2]),
      ymax = x$sp_func(out_x, alpha1[3], alpha2[3], alpha3[3])
    )
    out.FOI <- data.frame(
      x = out_x,
      ymin = x$foi_func(out_x, alpha1[1], alpha2[1], alpha3[1]),
      y = x$foi_func(out_x, alpha1[2], alpha2[2], alpha3[2]),
      ymax = x$foi_func(out_x, alpha1[3], alpha2[3], alpha3[3])
    )
  }else if(x$type == "far2"){
    alpha1 <- x$info["alpha1",c("2.5%","50%", "97.5%")]
    alpha2 <- x$info["alpha2",c("2.5%","50%", "97.5%")]

    out.DF <- data.frame(
      x = out_x,
      ymin = x$sp_func(out_x, alpha1[1], alpha2[1]),
      y = x$sp_func(out_x, alpha1[2], alpha2[2]),
      ymax = x$sp_func(out_x, alpha1[3], alpha2[3])
    )
    out.FOI <- data.frame(
      x = out_x,
      ymin = x$foi_func(out_x, alpha1[1], alpha2[1]),
      y = x$foi_func(out_x, alpha1[2], alpha2[2]),
      ymax = x$foi_func(out_x, alpha1[3], alpha2[3])
    )
  }else if(x$type == "log_logistic"){
    alpha1 <- x$info["alpha1",c("2.5%","50%", "97.5%")]
    alpha2 <- x$info["alpha2",c("2.5%","50%", "97.5%")]

    out.DF <- data.frame(
      x = out_x,
      ymin = x$sp_func(out_x, alpha1[1], alpha2[1]),
      y = x$sp_func(out_x, alpha1[2], alpha2[2]),
      ymax = x$sp_func(out_x, alpha1[3], alpha2[3])
    )
    out.FOI <- data.frame(
      x = out_x,
      ymin = x$foi_func(out_x, out.DF$ymin, alpha1[1], alpha2[1]),
      y = x$foi_func(out_x, out.DF$y, alpha1[2], alpha2[2]),
      ymax = x$foi_func(out_x, out.DF$ymax, alpha1[3], alpha2[3])
    )
  }else{
    warning('Expect model type to be one of the following: "far3", "far2", "log_logistic"')
  }

  list(out.DF, out.FOI)
}

# =========== Nonparametric =============
#' Compute confidence interval for local polynomial model
#'
#' Computes CI for Seroprevalence from model standard errors, and (optionally)
#' for Force of Infection via nonparametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param foi_ci whether to compute CI for FoI (default to FALSE)
#' @param ... arbitrary arguments
#'
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#'
#' @export
compute_ci.lp_model <- function(x,ci = 0.95, foi_ci=FALSE, ...){
  ages <- x$df$age
  crit<- crit(x$info,cov = ci)$crit.val
  mod1 <- predict(x$info, data.frame(a = ages),
                  se.fit = TRUE, band="local",
                  what="coef")

  # get the fit in predictor scale (before inverse link) to work with SE
  pred_raw <- log(mod1$fit/(1-mod1$fit))

  out.DF <- data.frame(
    x = ages,
    y = mod1$fit,
    # quantify CI
    ymin = x$info$trans(pred_raw - crit * mod1$se.fit),
    ymax = x$info$trans(pred_raw + crit * mod1$se.fit)
  )

  foi_x <- sort(unique(ages))
  foi_x <- foi_x[c(-1, -length(foi_x) )]

  out.FOI <- data.frame(
    x = foi_x,
    y = est_foi(ages, out.DF$y)
  )

  out.FOI <- if(foi_ci){
    bootstrap_res <- nonparametric_bootstrapping(
      mod = x, newdat = data.frame(age = ages),
      refit_func = \(df){
        do.call(
          lp_model,
          list(
            data = df,
            nn = x$nn, h = x$h, deg = x$deg, kern = x$kern
          )
        )
      },
      ci=ci, ...
    )

    cbind(out.FOI, bootstrap_res)
  }else{
    out.FOI
  }

  list(out.DF, out.FOI)
}

# ========== Semi-parametric ==========
#' Compute confidence interval for penalized_spline_model
#'
#' Computes CI for Seroprevalence from model standard errors, and (optionally)
#' for Force of Infection via nonparametric bootstrap.
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param foi_ci whether to compute CI for FoI (default to FALSE)
#' @param ... arbitrary arguments
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @return a list of 2 data frames:
#'   \itemize{
#'     \item seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#'     \item FoI estimates with columns: \code{x} (age), \code{y}
#'       (fitted FoI), and if \code{foi_ci = TRUE}, \code{ymin} and
#'       \code{ymax} (lower and upper confidence interval bounds)
#'
#' @export
compute_ci.penalized_spline_model <- function(x,ci = 0.95, foi_ci=FALSE, ...){
  # resolve no visible binding issue with CRAN check
  fit <- se.fit <- NULL

  m <- 1
  p <- (1 - ci) / 2

  # handle different output for different frameworks
  if(x$framework == "pl"){
    link_inv <- x$info$family$linkinv
    dataset <- x$info$model
    n <- nrow(dataset) - length(x$info$coefficients)
    gam_obj <- x$info
  }else{
    link_inv <- x$info$gam$family$linkinv
    dataset <- x$info$gam$model
    n <- nrow(dataset) - length(x$info$gam$coefficients)
    gam_obj <- x$info$gam
  }

  ages <- unique(x$df$age)
  # print(head(ages))

  mod <- predict.gam(gam_obj, data.frame(age = ages), se.fit = TRUE)  %>%
    as_tibble()  %>%
    select(fit, se.fit) %>%
    mutate(age = ages)  %>%
    mutate(lwr = m * link_inv(fit + qt(    p, n) * se.fit),
           upr = m * link_inv(fit + qt(1 - p, n) * se.fit),
           fit = m * link_inv(fit))  %>%
    select(- se.fit)

  out.DF <- data.frame(x = ages, y = mod$fit,
                       ymin= mod$lwr, ymax = mod$upr)

  # print(ages)

  foi_x <- sort(ages)
  foi_x <- foi_x[c(-1, -length(foi_x) )]
  out.FOI <- data.frame(x = foi_x, y = est_foi(ages, mod$fit))

  out.FOI <- if(foi_ci){
    bootstrap_res <- nonparametric_bootstrapping(
      mod = x, newdat = data.frame(age = ages),
      refit_func = \(df){
        do.call(
          penalized_spline_model,
          c(
            list(
              data = df,
              framework = x$framework
            ),
            x$pars
          )
        )
      },
      ci=ci, ...
    )

    cbind(out.FOI, bootstrap_res)
  }else{
    out.FOI
  }

  return(list(out.DF, out.FOI))
}

#' Compute confidence interval for time age model
#'
#' @param x serosv models
#' @param ci confidence level for the interval
#' @param le number of data for computing confidence interval
#' @param ... arbitrary argument
#'
#' @importFrom mgcv predict.gam
#' @import dplyr
#'
#' @return confidence interval dataframe with n_group x 3 cols, the columns are `group`, `sp_df`, `foi_df`
#' @export
compute_ci.age_time_model <- function(x, ci=0.95, le = 100, ...){
  # resolve no visible binding note
  df <- monotonized_info <- monotonized_ci_mod <- age <- info <- fit <- se.fit <- sp_df <- foi_df <- NULL

  # check which type of model user wants to visualize
  modtype <- if (is.null(list(...)[["modtype"]])) "monotonized" else list(...)$modtype
  assert_that(
    modtype == "monotonized" | modtype == "non-monotonized",
    msg = "modtype argument must be eithers 'monotonized' or 'non-monotonized'"
  )

  p <- (1 - ci) / 2

  # use model to generate seroprev (with CI) and FOI on a finer grid for plotting
  age_range <- range(bind_rows(x$out$df)$age)
  out <- x$out %>%
    mutate(
      age = map(df, \(dat){
        seq(age_range[1], age_range[2], length.out = le)
      })
    )

  # --- use the monotonized model for prediction and ci
  if(modtype == "monotonized"){
    out <- out %>%
      mutate(
        sp_df = pmap(list(monotonized_info, monotonized_ci_mod, age), \(mod, ci_mod, grid){
          data.frame(
            x = grid,
            y = predict(mod, list(age = grid), type = "response"),
            ymin = predict(ci_mod$ymin, list(age = grid), type = "response"),
            ymax = predict(ci_mod$ymax, list(age = grid), type = "response")
          )
        })
      )
  }else{
    # --- if user specify non-monotonized then simply compute CI from gam model
    out <- out %>%
      mutate(
        sp_df = map2(info, age, \(mod, grid){
          link_inv <- mod$family$linkinv
          dataset <- mod$model[,1:2]
          n <- nrow(dataset) - length(mod$coefficients)

          predict(mod, data.frame(age = grid), se.fit = TRUE)  %>%
            as_tibble()  %>%
            select(fit, se.fit) %>%
            mutate(
              x = grid,
              ymin = link_inv(fit + qt(    p, n) * se.fit),
              ymax = link_inv(fit + qt(1 - p, n) * se.fit),
              y = link_inv(fit)
            )  %>%
            select(- se.fit)
        })
      )
  }

  # ------- finally, compute FOI
  out <- out %>%
    mutate(
      foi_df = map2(age, sp_df, \(grid, sp){
        foi_x <- sort(unique(grid))
        foi_x <- foi_x[c(-1, -length(foi_x) )]

        tibble(
          x = foi_x,
          y = est_foi(grid, sp$y)
        )
      })
    ) %>%
    select(!!sym(x$grouping_col), sp_df, foi_df)

  out
}


# ========== Mixture model ========
#' Compute confidence interval for mixture model
#'
#' @param x serosv mixture_model object
#' @param ci confidence level for the interval
#' @param ... arbitrary arguments
#' @importFrom stats qnorm
#'
#' @return list of confidence interval for susceptible and infected. Each confidence interval is a list with 2 items for lower and upper bound of the interval.
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

#' Compute confidence interval for the prevalence estimate from mixture model
#'
#' CI for prevalence is the transform CI for mu(a) estimator, and not accounting for
#' the uncertainty in mu_I and mu_S estimates
#'
#' @param x serosv mixture_model object
#' @param ci confidence level for the interval
#' @param ... arbitrary arguments
#' @importFrom stats qnorm
#' @importFrom dplyr mutate
#'
#' @return a data frames of seroprevalence estimates with columns: \code{x} (age),
#'       \code{y} (fitted seroprevalence), \code{ymin} and \code{ymax}
#'       (lower and upper confidence interval bounds)
#' @export
compute_ci.estimate_from_mixture <- function(x, ci=.95, ...){
  # resolve no visible binding issue with CRAN check
  fit <- se.fit <- NULL

  # set up
  p <- (1 - ci) / 2

  link_inv <- x$info$family$linkinv
  dataset <- x$info$model[,1:2]
  n <- nrow(dataset) - length(x$info$coefficients)
  ages <- sort(unique(x$df$age))
  mu_s <- x$mu_s
  mu_i <- x$mu_i

  # estimate the CI of mu(a)
  mu_a <- predict.gam(x$info, data.frame(age = ages), se.fit = TRUE)  %>%
    as.data.frame()  %>%
    select(fit, se.fit) %>%
    mutate(x = ages)  %>%
    mutate(ymin = link_inv(fit + qt(    p, n) * se.fit),
           ymax = link_inv(fit + qt(1 - p, n) * se.fit),
           y = link_inv(fit))  %>%
    select(- se.fit)

  # transform mu(a) estimate to get CI for prevalence
  out.DF <- mu_a %>%
    mutate(
      ymin = (ymin - mu_s)/(mu_i - mu_s),
      ymax = (ymax - mu_s)/(mu_i - mu_s),
      y = (y - mu_s)/(mu_i - mu_s)
    )

  # if monotonized, apply pava to estimates
  if(x$monotonize){
    out.DF <- out.DF %>%
      mutate(
        ymin = pava(lwr)$pai2,
        ymax = pava(upr)$pai2,
        y = pava(fit)$pai2
      )
  }

  out.DF
}



