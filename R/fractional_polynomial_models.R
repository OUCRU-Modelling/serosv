is_monotone <- function(model) {
  (sum(diff(predict(model))<0)==0)
}

formulate <- function(p) {
  equation <- "cbind(pos,tot-pos)~"
  prev_term <- ""

  for (i in 1:length(p)) {
    if (i > 1 && p[i] == p[i-1]) {
      cur_term <- paste0("I(", prev_term, "*log(age))")
    } else if (p[i] == 0) {
      cur_term <- "log(age)"
    } else {
      cur_term <- paste0("I(age^", p[i], ")")
    }
    equation <- paste0(equation, "+", cur_term)
    prev_term <- cur_term
  }
  equation
}

#' Returns the powers of the fractional polynomial model which has the lowest deviance score.
#'
#' Return the best powers for a given degree
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR
#' `age`, `status` (for linelisting data)
#' @param p a powers sequence to be tested.
#' @param mc indicates if the returned model should be monotonic.
#' @param degree the degree of the model (i.e. number of power terms). Recommended to be <= 2.
#' @param link the link function. Defaulted to "logit".
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#'
#' @return list of 3 elements:
#'   \item{p}{The best power for fp model.}
#'   \item{deviance}{Deviance of the best fitted model.}
#'   \item{model}{The best model fitted}
#'
#' @examples
#' df <- hav_be_1993_1994
#' best_p <- find_best_fp_powers(
#' df,
#' p=seq(-2,3,0.1), mc=FALSE, degree=2, link="cloglog"
#' )
#' best_p
#'
#' @importFrom stats glm binomial as.formula
#'
#' @export
find_best_fp_powers <- function(data,
                                age_col="age",pos_col="pos", tot_col="tot", status_col="status",
                                p, mc, degree, link="logit"){
  data <- check_input(data, stratum_col=age_col, pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  glm_best <- NULL
  d_best <- NULL
  p_best <- NULL
  #----
  min_p <- 1
  max_p <- length(p)
  state <- rep(min_p, degree) # state=pointers (or indices) to current power values
  i <- degree
  #----

  # helper to get current powers from the state
  get_cur_p <- function(cur_state) {
    cur_p <- c()
    for (i in 1:degree) {
      cur_p <- c(cur_p, p[cur_state[i]])
    }
    cur_p
  }

  repeat {
    if (
      (i < degree && state[i] == max_p)
      || (i == degree && state[i] == max_p+1)
    ) {
      # stop when done looping through all the terms
      if (i-1 == 0) break

      if (state[i-1] < max_p) {
        # move on to the power of the current term
        state[i-1] <- state[i-1]+1
        for (j in i:degree) state[j] <- state[i-1]
        i <- degree
      } else {
        # move on to the next term
        i <- i-1
        next
      }
    }
    #------ iteration implementation -------
    p_cur <- get_cur_p(state)

    glm_cur <- glm(
      as.formula(formulate(p_cur)),
      family=binomial(link=link)
    )
    if (glm_cur$converged == TRUE) {
      # d_cur <- deviance(glm_cur)
      d_cur <- glm_cur$deviance
      if (is.null(glm_best) || d_cur < d_best) {
        if ((mc && is_monotone(glm_cur)) | !mc) {
          glm_best <- glm_cur
          d_best <- d_cur
          p_best <- p_cur
        }
      }
    }
    #---------------------------------------

    # stop after reaching the highest power for all terms
    if (sum(state != max_p) == 0) break

    state[i] <- state[i]+1 # increase power of the current term
  }
  return(list(p=p_best, deviance=d_best, model=glm_best))
}

#' A fractional polynomial model.
#'
#' @description Fractional polynomial model is a generalization of polynomial models
#' where the power of the terms can be fractions, allowing more flexibility and better
#' fit for data where asymptotic behavior is expected.
#'
#' @details
#' Instead of a polynomial, the linear predictor is now defined as
#' \deqn{
#'  \eta_m(a, \beta, p_1, p_2, ...,p_m) = \Sigma^m_{i=0} \beta_i H_i(a)
#' }
#' Where \eqn{m} is an integer, \eqn{p_1 \le p_2 \le... \le p_m} is a sequence of powers,
#' and \eqn{H_i(a)} is a transformation given by
#'
#' \deqn{
#' H_i = \begin{cases}
#' a^{p_i} & \text{ if } p_i \neq p_{i-1},
#' \\ H_{i-1}(a) \times log(a)  & \text{ if } p_i = p_{i-1},
#' \end{cases}
#' }
#'
#' Refers to section 6.2. of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR `age`, `status` (for linelisting data)
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param p the powers of the predictor.
#' @param link the link function for model. Defaulted to "logit".
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#' @importFrom stats predict as.formula
#'
#' @return a list of class fp_model with 5 items
#'   \item{datatype}{type of data used for fitting model (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{a fitted glm model}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso
#' [stats::glm()] for more information on glm object
#'
#' [polynomial_models()]
#'
#' @examples
#' df <- hav_be_1993_1994
#' model <- fp_model(
#'   df,
#'   p=c(1.5, 1.6), link="cloglog")
#' plot(model)
#'
#' @export
fp_model <- function(data,p,link="logit",
                     age_col="age",pos_col="pos", tot_col="tot", status_col="status") {
  model <- list()

  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  model$datatype <- data$type

  model$info <- glm(
    as.formula(formulate(p)),
    family=binomial(link=link)
  )
  model$sp  <- model$info$fitted.values
  model$foi <- est_foi(
    t=age,
    sp=model$info$fitted.values
  )
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "fp_model"
  model
}
