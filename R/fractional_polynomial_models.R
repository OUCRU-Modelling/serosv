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
#' @param degree the maximum degree (i.e. number of power terms) to search for the best model. Recommended to be <= 2.
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
#' @import dplyr tidyr
find_best_fp_powers <- function(data,
                                p, mc, degree, link="logit"){
  age <- data$age
  pos <- data$pos
  tot <- data$tot

  best_mod <- NULL # best model
  best_p <- NULL # best powers (p vector) for the given degree

  # Starting from the lowest degree
  # Get the best combinations of powers p
  # Try increment m and only accept when there is a statistically significant improvement (determined by LRT)
  for (curr_deg in 1:degree){
    # generate combinations of powers
    p_combis <- expand.grid(rep(list(p), curr_deg))
    # filter the rows that are increasing in order
    p_combis <- p_combis[apply(p_combis, 1, \(r){all(diff(r)>=0)}), ,drop=FALSE]
    # arrange by increasing order of power in first->second degree and so on
    p_combis <- p_combis[do.call(order,p_combis),,drop=FALSE]

    # generate formulas for all the combinations of powers
    # get model with lowest deviance
    best_dev <- 1e8 # best deviance for current degree
    curr_deg_mod <- NULL # best model for current degree
    curr_deg_p <- NULL # best powers for current degree

    for (row_id in 1:nrow(p_combis)){
      curr_p <- as.numeric(p_combis[row_id,])
      curr_mod <- glm(
        as.formula(formulate(curr_p)),
        family=binomial(link=link)
      )

      if(curr_mod$converged==TRUE){
        curr_dev <- curr_mod$deviance
        if (is.null(curr_deg_mod) || curr_dev < best_dev) {
          # make sure to only accept monotone model if specified
          if ((mc && is_monotone(curr_mod)) | !mc) {
            best_dev <- curr_dev
            curr_deg_mod <- curr_mod
            curr_deg_p <- curr_p
          }
        }
      }
    }

    # check if the best model with current degree is better than the last degree
    if(!is.null(curr_deg_mod)){
      if(!is.null(best_mod)){
        # perform LRT
        lrt_out <- anova(best_mod, curr_deg_mod, test="LRT")
        lrt_out <- as.data.frame(lrt_out)

        # Royston&Altman (1994) suggests significance level of 0.1
        if(lrt_out$`Pr(>Chi)`[2] < 0.1){
          best_mod <- curr_deg_mod
          best_p <- curr_deg_p
        }
      }else{
        best_mod <- curr_deg_mod
        best_p <- curr_deg_p
      }
    }
  }

  if(is.null(best_mod)) stop("Cannot find a converged model with the given degree and power")

  list(p=best_p, model=best_mod)
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
#' @param p is either:
#'   (1) a numeric vector specifying the powers to apply to the predictors, or
#'   (2) a named list with two elements, \code{"p_range"} and \code{"degree"}. \code{"p_range"} \
#'      is a sequence of powers and \code{"degree"} is the maximum degree to search over
#'   If (1), the supplied values are used directly as the powers for the predictors.
#'   If (2), the package searches for the best degree and power combinations
#' @param degree the degree of the model (i.e. number of power terms). Recommended to be <= 2.
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
fp_model <- function(data,p,monotonic=FALSE,link="logit",
                     age_col="age",pos_col="pos", tot_col="tot", status_col="status") {
  model <- list()

  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  # handle powers input here
  if(is.numeric(p)){
    model$info <- glm(
      as.formula(formulate(p)),
      family=binomial(link=link)
    )
  }else if(is.list(p) && all(c("p_range", "degree") %in% names(p))){
    out <- find_best_fp_powers(
      data = data.frame(age=age, pos=pos, tot=tot),
      p = p$p_range, degree = p$degree, mc = monotonic, link = link
    )
    model$info <- out$model
  }else{
    stop("Invalid value for `p`: either a numeric vector or a named list with
         2 elements `p_range` and `degree`")
  }


  model$sp  <- model$info$fitted.values
  model$foi <- est_foi(
    t=age,
    sp=model$info$fitted.values
  )
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "fp_model"
  model
}
