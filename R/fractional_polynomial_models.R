is_monotone <- function(model) {
  (sum(diff(predict(model))<0)==0)
}

formulate <- function(powers) {
  equation <- "cbind(pos,tot-pos)~"
  prev_term <- ""

  for (i in 1:length(powers)) {
    if (i > 1 && powers[i] == powers[i-1]) {
      cur_term <- paste0("I(", prev_term, "*log(age))")
    } else if (powers[i] == 0) {
      cur_term <- "log(age)"
    } else {
      cur_term <- paste0("I(age^", powers[i], ")")
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
#' @param df the data to be fitted. MUST have `pos`, `tot`, and `age` columns.
#'
#' @param p a powers sequence.
#'
#' @param mc indicates if the returned model should be monotonic.
#'
#' @param degree the degree of the model. Recommended to be <= 2.
#'
#' @param link the link function. Defaulted to "logit".
#'
#' @examples
#' hbe_best <- hav_be_1993_1994 %>%
#'   find_best_fp_powers(p=seq(-2,3,0.1), mc=F, degree=2, link="cloglog")
#' hbe_best
#'
#' @export
find_best_fp_powers <- function(df, p, mc, degree, link="logit"){
  pos <- df$pos
  tot <- df$tot
  age <- df$age
  #----
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
      d_cur <- deviance(glm_cur)
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
  return(list(power=p_best, deviance=d_best, model=glm_best))
}

# assuming age is sorted!!
estimate_foi <- function(age, sp)
{
  dsp <- diff(sp)/diff(age)
  foi <- approx(
    (age[-1]+age[-length(age)])/2,
    dsp,
    age[c(-1,-length(age))]
  )$y/(1-p[c(-1,-length(age))])

  foi
}

#' A fractional polynomial model.
#'
#' Refers to section 6.2.
#'
#' @param p the powers of the predictor.
#'
#' @param link the link function for model. Defaulted to "logit".
#'
#' @examples
#' model <- fp_model(c(1.5, 1.6), link="cloglog")
#' model
#'
#' @export
fp_model <- function(p, link="logit") {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        as.formula(formulate(p)),
        family=binomial(link=link)
      )
      X <- generate_X(age)
      model$sp  <- model$info$fitted.values
      model$foi <- estimate_foi(age, model$info$fitted.values)
      model$df  <- df
      model
    })
  }
  model
}
