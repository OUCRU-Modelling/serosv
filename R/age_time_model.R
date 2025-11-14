#' Age-time varying seroprevalence
#'
#'
#' @description Fit age-stratified seroprevalence across multiple time points. Also try to monotonize age (or birth cohort) - specific seroprevalence.
#'
#' @param data - input data, must have`age`, `status`, time, group columns, where group column determines how data is aggregated
#' @param time_col - name of the column for time (default to `date`)
#' @param grouping_col - name of the column for time (default to `group`)
#' @param age_correct - a boolean, if `TRUE`, monotonize age-specific prevalence. Monotonize birth cohort-specific seroprevalence otherwise.
#' @param le - number of bins to generate age grid, used when monotonizing data
#' @param ci - confidence interval for smoothing
#' @param monotonize_method - either "pava" or "scam"
#' @import mgcv scam assertthat
#'
#' @return a list of class time_age_model with 4 items
#'   \item{out}{a data.frame with dimension n_group x 9, where columns `info`, `sp`, `foi` store output for non-monotonized
#' data and `monotonized_info`, `monotonized_sp`,  `monotonized_foi`,  `monotonized_ci_mod` store output for monotonized data}
#'   \item{grouping_col}{name of the column for grouping}
#'   \item{age_correct}{a boolean indicating whether the data is monotonized across age or cohort}
#'   \item{datatype}{whether the input data is aggregated or line-listing data}
#' @export
age_time_model <- function(data, time_col="date", grouping_col="group", age_correct=F, le=512, ci = 0.95, monotonize_method = "pava"){
  # work around to resolve no visible binding note NOTE during check()
  x <- label <- family <- fit <- se.fit <- ymin <- ymax <- y <- mean_time <- prevalence <- sim_data <- NULL
  age <- ys <- shift_no <- cohort <- col_time <- monotonized_mod <- df <- info <- sp <- monotonized_info <- monotonized_sp <- NULL


  assert_that(
    (monotonize_method == "pava") | (monotonize_method == "scam"),
    msg = 'Monotonize method must be either "pava" or "scam"'
  )

  # ---- helper functions -----
  shift_right <- \(n,x){ if(n == 1) x else dplyr::lag(x, n, default = NA)}
  # function to simulate data for monotonize process
  generate_data <- \(dat, mod, no_sim=100) {
    link_inv <- family(mod)$linkinv
    n <- nrow(dat) - length(coef(mod))
    p <- (1 - ci)/2

    pred <- predict(mod, dat, se.fit = TRUE)  %>%
      as_tibble()  %>%
      select(fit, se.fit) %>%
      mutate(
        ymin = link_inv(fit + qt(p, n) * se.fit),
        ymax = link_inv(fit + qt(1 - p, n) * se.fit),
        y = link_inv(fit)
      )  %>%
      select(-se.fit,-fit)

    dat %>%
      bind_cols(pred) %>%
      pivot_longer(c(ymin, ymax, y),
                   names_to = "ys",
                   values_to = "prevalence")
  }
  # function to monotonize data using serosv pava function
  monotonize_data <- \(dat, grp){
    dat <- dat %>% arrange(mean_time)
    if(monotonize_method == "scam"){
      out <- tryCatch(
        {
          mod <- scam(prevalence ~ s(mean_time, bs = "mpi"), data=dat,family = betar)
          dat$prevalence <- predict(mod, list(mean_time = dat$mean_time), type = "response")
          dat
        },
        error = \(e) {
          e
        }
      )

      if (!("error" %in% class(out))) return(out)
    }

    dat %>%
      mutate(
        prevalence = serosv::pava(prevalence)$pai2
      )
  }
  # initialize model obj
  model <- list()

  # --- preprocess data ------
  check_input <- check_input(data)
  age_range <- range(data$age)
  age_grid <- seq(age_range[1], age_range[2], length.out = le)

  model$datatype <- check_input$type
  data <- check_input[c("age", "pos", "tot")] %>% as.data.frame() %>% bind_cols(data[c(time_col, grouping_col)])

  # ---- gam model for age-stratified prevalence for each group -----
  gam_mods <- data %>%
    group_by(.data[[grouping_col]]) %>% nest() %>%
    mutate(
      mod = map(data, \(dat){
        # handle potential error when dataset is small
        k <- if(length(unique(dat$age)) < 10) length(unique(dat$age)) - 1 else -1

        mod <- if(model$datatype == "linelisting") gam(pos ~ s(age, k=k), data = dat, family = binomial) else
          gam(cbind(pos, tot - pos) ~ s(age, k=k), data = dat, family = binomial)
        mod
      }),
      mean_time = map_dbl(data, \(dat){mean(dat[[time_col]])}) %>% as.Date()
    ) %>%
    ungroup()

  # ----- branching based on age_correct ---
  # if age_correct is TRUE: enforce monotonic increase in prevalence overtime within age group
  # otherwise: enforce monotonic increase in prevalence within cohort
  if(age_correct == FALSE){
    # simulate data + monotonize data using scam
    scam_out <- gam_mods %>%
      select(!!sym(grouping_col), mod, mean_time) %>%
      mutate(
        # simulate data to fit scam model
        sim_data = map(mod, \(mod){
          data.frame(age = age_grid) %>% generate_data(mod)
        })
      ) %>%
      select(-mod) %>% unnest(sim_data) %>%
      group_by(age, ys) %>%
      group_modify(monotonize_data) %>% ungroup()

    # modify monotonized data
    scam_out <- scam_out %>%
      pivot_wider(names_from = ys, values_from = prevalence) %>%
      group_by(!!sym(grouping_col), mean_time)  %>%
      nest()
  }else{
    dpy <- 365

    # simulate data to monotonize
    # return a data.frame of collection_time, age (at current collection time), cohort (age at first collection time)
    scam_data <- gam_mods %>%
      mutate(
        age = map(mean_time, \(.) {
          age_grid
        }),
        shift_no = (mean_time - min(mean_time)) / (dpy * mean(diff(age_grid))),
        cohort = map(shift_no, \(n) {
          shift_right(round(n), age_grid)
        }),
        sim_data = pmap(list(mod, age, cohort, mean_time),
                        \(mod, age, cohort, mean_time) {
                          data.frame(age = age, cohort = cohort) %>%
                            generate_data(mod)
                        })
      ) %>%
      select(!!sym(grouping_col), mean_time, sim_data) %>%
      unnest(sim_data)

    # ----- use scam model to monotonize cohort-stratifed prevalence over time----
    scam_out <- scam_data %>%
      filter(
        cohort < max(age) - diff(range(mean_time)) / dpy,
        !is.na(cohort)
      ) %>%
      group_by(cohort, ys) %>%
      group_modify(monotonize_data) %>%
      ungroup()


    # mapping to covert cohort to age
    cohort_age_mapping <- scam_data %>%
      select(!!sym(grouping_col), age, cohort) %>%
      unique()

    # map cohort from monotized data to age (at collection time)
    scam_out <- scam_out %>%
      left_join(
        cohort_age_mapping,
        by = join_by(
          !!sym(grouping_col) == !!sym(grouping_col), cohort == cohort, age == age
        )
      ) %>%
      pivot_wider(names_from = ys, values_from = prevalence) %>%
      group_by(!!sym(grouping_col), mean_time)  %>%
      nest()
  }

  # ------ Fit the monotonized data ------
  out <- scam_out %>%
    mutate(
      monotonized_mod = map(data, \(dat){
        # handle potential error when dataset is small
        k <- if(length(unique(dat$age)) < 10) length(unique(dat$age)) - 1 else -1
        # also when y range is small
        k <- if(diff(range(dat$y)) < 0.01) 3 else k

        gam(y ~ s(age, k=k), family = betar, data = dat)
      }),
      # also have model for smooth ci
      monotonized_ci_mod = map(data, \(dat){
        # handle potential error when dataset is small
        k <- if(length(unique(dat$age)) < 10) length(unique(dat$age)) - 1 else -1
        # also when y range is small
        k <- if(min(diff(range(dat$ymin)), diff(range(dat$ymax))) < 0.01) 3 else k

        list(
          "ymin" = gam(ymin ~ s(age, k=k), family = betar, data = dat),
          "ymax" = gam(ymax ~ s(age, k=k), family = betar, data = dat)
        )
      })
    ) %>%
    ungroup() %>%
    select(-data) %>%
    right_join(gam_mods,
               by = join_by(!!sym(grouping_col) == !!sym(grouping_col), mean_time == mean_time)) %>%
    select(-mean_time)

  # reformat output
  out <- out %>%
    # rename to follow the convention of other functions
    rename(
      df = data,
      info = mod,
      monotonized_info = monotonized_mod
    ) %>%
    # finally predict seroprevalence and foi for the input data
    mutate(
      sp = map2(df, info, \(dat, mod){
        predict(mod, list(age = dat$age), type="response")
      }),
      foi = map2(df, sp, \(dat, sp){
        est_foi(dat$age, sp)
      }),
      monotonized_sp = map2(df, monotonized_info, \(dat, mod){
        predict(mod, list(age = dat$age), type="response")
      }),
      monotonized_foi = map2(df, monotonized_sp, \(dat, sp){
        est_foi(dat$age, sp)
      })
    )

  model$out <- out
  model$grouping_col <- grouping_col
  model$age_correct <- age_correct


  class(model) <- "age_time_model"

  model
}

