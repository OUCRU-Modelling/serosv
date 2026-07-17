# compute the i * a^(i-1) matrix
X <- function(t, degree) {
  # X_matrix <- matrix(rep(1, length(t)), ncol = 1)
  # if (degree > 1) {
  #   for (i in 2:degree) {
  #     X_matrix <- cbind(X_matrix, i * t^(i-1))
  #   }
  # }
  -sapply(1:degree, function(i) i * t^max(i-1, 0))
}

#' Polynomial models
#'
#' @description Fit age-stratified seroprevalence data to serocatalytic models formulated as polynomials.
#'
#' @details
#' The seroprevalence is assumed to follow the general format
#' \deqn{
#' \pi(a)  = 1 - e^{-\Sigma_{i=1}^k \beta_i a^i}
#' }
#' Which implies the force of infection to be \eqn{\lambda(a) = \Sigma_{i=1}^k \beta_i i a^{i-1}}
#'
#' Where:
#'
#' - \eqn{\pi} is the seroprevalence at age \eqn{a}
#'
#' - \eqn{a} is the variable age
#'
#' - \eqn{k} is the degree of the polynomial
#'
#' The seroprevalence \eqn{\pi(a)} is fitted using a GLM with log link with
#' the linear predictor \eqn{\eta(a) = \Sigma_{i=1}^k \beta_i a^{i}}
#'
#' Muench (1934) model is equivalent to a degree 1 (\eqn{k=1}) linear predictor
#'
#' Griffith model is equivalent to a degree 2 (\eqn{k=2}) linear predictor
#'
#' Grenfell & Anderson (1985) suggested a higher order polynomials (\eqn{k \geq 3})
#'
#' Refer to section 6.1.1. of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' Grenfell, B. T., and R. M. Anderson. 1985. “The Estimation of
#' Age-Related Rates of Infection from Case Notifications and Serological Data.”
#' The Journal of Hygiene 95 (2): 419–36. \doi{https://doi.org/10.1017/s0022172400062859}.
#'
#' Muench, Hugo. 1934. “Derivation of Rates from Summation Data by the Catalytic Curve.”
#' Journal of the American Statistical Association 29 (185):
#' 25–38. \doi{https://doi.org/10.1080/01621459.1934.10502684}.
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR `age`, `status` (for linelisting data)
#' @param k  degree of the polynomial. (k=1 for Muench model, k=2 for Griffith model, k=3 for Grenfell model)
#' @param link link function (default link="log")
#' @param age_col name of the `age` column (default age_col="age")
#' @param pos_col name of the `pos` column (default pos_col="pos")
#' @param tot_col name of the `tot` column (default tot_col="tot")
#' @param status_col name of the `status` column (default status_col="status")
#' @param ... additional arguments to be passed to `glm()` function that fits the model
#'
#' @examples
#' data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ]
#' aggregated <- transform_data(data, stratum_col = "age", status_col="seropositive")
#'
#' # fit with aggregated data
#' model <- polynomial_model(aggregated, k = 1)
#' # fit with linelisting data
#' model <- polynomial_model(data,
#'     status_col = "seropositive",
#'     k = 1)
#' plot(model)
#'
#' @return a list of class polynomial_model with 5 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "glm" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'   \item{foi_mod}{function to compute FoI given a vector of age and estimated parameters}
#'
#' @export
polynomial_model <- function(data, k, link = "log",
                             age_col="age",pos_col="pos", tot_col="tot", status_col="status",
                             ...){
  model <- list()
  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  model$datatype <- data$type

  age <- data$age
  pos <- data$pos
  neg <- data$tot - pos

  df <- data.frame(cbind(age, pos,neg))

  # helper function to generate the polynomial given a k value
  # to be used for parameter selection if multiple values for k are given
  generate_polynomial <- function(k, df, link="log"){
    Age <- function(k){
      if(k>1){
        formula<- paste0("I","(",paste("age", 2:k,sep = "^"),")",collapse = "+")
        paste0("cbind(neg,pos)"," ~","-1+age+",formula)
      } else {
        paste0("cbind(neg,pos)"," ~","-1+age")
      }
    }


    tryCatch(
      mod <- glm(Age(k), family=binomial(link=link),df, ...),
      warning = function(w){
        warning(sprintf("glm warning for degree k=%d: %s", k, conditionMessage(w)), call. = FALSE)
        suppressWarnings(glm(Age(k), family = binomial(link = link), data = df))
      },
      error = function(e){
        warning(sprintf("glm failed for degree k=%d: %s", k, conditionMessage(e)), call. = FALSE)
      }
    )
  }

  # If a vector of values for k is provided -> select best value
  if(length(k) > 1){
    out <- nested_mod_selection(
      list("k" = k),
      model_fn = \(k, df){
        generate_polynomial(k, df, link=link)
      },
      dat = df
    )

    k <- out$best_par$k
    model$info <- out$mod
  }else{
    model$info <- generate_polynomial(k, df, link=link)
  }

  X <- X(age, k)
  model$sp <- 1 - model$info$fitted.values
  model$foi <- X%*%model$info$coefficients
  model$df <- data.frame(age=age, pos=pos, tot= pos + neg)
  # function to generate FoI given age and coefs
  model$foi_mod <- function(age, coefs){
    age_mat <- X(age, k)
    age_mat %*% coefs
  }
  model$k <- k
  class(model) <- "polynomial_model"
  model
}

# TODO: check if this can be generalized to other functions as well (e.g. fractional polynomial)
# function to return the best parameter of nested glm models using LRT
# par_range - list of parameters and its possible values
# model_fn - function to fit and return a model, must takes 2 arguments: par, df
#' @import tidyr
#' @importFrom purrr pmap keep
#' @importFrom stats anova
nested_mod_selection <- function(par_range, model_fn, dat, method="LRT"){
  # work around to resolve no visible binding note NOTE during check()
  `Pr(>Chi)` <- Deviance <- idx <- NULL

  # conditions for model to be usable
  is_usable_model <- function(mod){
    if (is.null(mod)) return(FALSE)
    if (!mod$converged) return(FALSE)
    dev <- tryCatch(deviance(mod), error = function(e) NA)
    if (is.na(dev) || is.infinite(dev)) return(FALSE)
    TRUE
  }

  # generate all combinations of parameters values
  par_combs <- tidyr::crossing(!!!par_range)

  # fit the model using specified parameter values
  mods_out <- par_combs %>%
    purrr::pmap(\(...){
      model_fn(..., df=dat)
    })

  mods_out <- purrr::keep(mods_out, is_usable_model)
  # perform LRT
  lrt_out <- do.call(
    anova,
    c(mods_out, list(test="LRT"))
  )

  # get the best model
  best_idx <- lrt_out %>%
    as.data.frame() %>%
    mutate(
      idx = 1:n()
    ) %>%
    filter(
      `Pr(>Chi)` < 0.05
    ) %>%
    arrange(
      Deviance
    ) %>%
    pull(idx)

  # handle scenario when the reference model (i.e., the first model) is in fact the best option
  # i.e., when the other parameter combinations do not result in statistically significant improvement
  if(length(best_idx)>1){
    best_idx <- best_idx[1]
  }else{
    best_idx <- 1
  }

  list(
    best_par = par_combs[best_idx, ] %>% as.list(),
    mod = mods_out[best_idx][[1]]
  )
}

