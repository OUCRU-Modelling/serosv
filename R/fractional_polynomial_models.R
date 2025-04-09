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

#' Returns the powers of the GLM fitted model which has the lowest deviance score.
#'
#' Refers to section 6.2.
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param p a powers sequence.
#' @param mc indicates if the returned model should be monotonic.
#' @param degree the degree of the model. Recommended to be <= 2.
#' @param link the link function. Defaulted to "logit".
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
find_best_fp_powers <- function(data, p, mc, degree, link="logit"){
  data <- check_input(data)
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  glm_best <- NULL
  d_best <- NULL
  p_best <- NULL
  #----
  min_p <- 1
  max_p <- length(p)
  state <- rep(min_p, degree)
  i <- degree
  #----

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
      if (i-1 == 0) break
      if (state[i-1] < max_p) {
        state[i-1] <- state[i-1]+1
        for (j in i:degree) state[j] <- state[i-1]
        i <- degree
      } else {
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
    if (sum(state != max_p) == 0) break
    state[i] <- state[i]+1
  }
  return(list(p=p_best, deviance=d_best, model=glm_best))
}

#' A fractional polynomial model.
#'
#' Refers to section 6.2.
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param p the powers of the predictor.
#' @param link the link function for model. Defaulted to "logit".
#'
#' @importFrom stats predict as.formula
#'
#' @return a list of class fp_model with 5 items
#'   \item{datatype}{type of data used for fitting model (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{a fitted glm model}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [stats::glm()] for more information on glm object
#'
#' @examples
#' df <- hav_be_1993_1994
#' model <- fp_model(
#'   df,
#'   p=c(1.5, 1.6), link="cloglog")
#' plot(model)
#'
#' @export
fp_model <- function(data,p,  link="logit") {
  model <- list()

  data <- check_input(data)
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
