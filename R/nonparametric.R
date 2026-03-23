#' A local polynomial model.
#'
#' @description Fit the age-specific seroprevalence to a local polynomial model,
#' where the linear predictor is approximated locally at one particular age.
#'
#' @details
#' Consider a linear predictor \eqn{\eta(a)} approximated locally at one particular value \eqn{a_0}.
#'
#' For a general degree \eqn{p}, the linear predictor for a neighbor of \eqn{a_0}, labeled \eqn{a_i} is the Taylor approximation
#' \deqn{
#' \eta(a_i) = \eta(a_0) + \eta^{(1)}(a_0)(a_i - a_0) +
#' \frac{\eta^{(2)}(a_0)}{2}(a_i - a_0)^2 + ... + \frac{\eta^{(p)}(a_0)}{p!}(a_i - a_0)^p
#' }
#'
#' Where the estimator for the \eqn{k}-th derivative of \eqn{\eta(a_0)}, for \eqn{k = 0,1,…,p}
#' (degree of local polynomial) is as followed:
#' \deqn{
#'  \hat{\eta}^{(k)}(a_0) = k!\hat{\beta}_k(a_0)
#' }
#'
#' The estimator for the prevalence at age \eqn{a_0} is then given by
#' \deqn{
#' \hat{\pi}(a_0) = g^{-1}\{ \hat{\beta}_0(a_0) \}
#' }
#' Where \eqn{g} is the link function
#'
#' The estimator for the force of infection at age \eqn{a_0} by assuming \eqn{p \ge 1} is as followed
#' \deqn{
#' \hat{\lambda}(a_0) = \hat{\beta}_1(a_0) \delta \{ \hat{\beta}_0 (a_0) \}
#' }
#' Where \eqn{\delta \{ \hat{\beta}_0(a_0) \} = \frac{dg^{-1} \{ \hat{\beta}_0(a_0) \} } {d\hat{\beta}_0(a_0)}}
#'
#' Refer to section 7.1 and 7.2. of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR `age`, `status` (for linelisting data)
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#' @param nn Nearest neighbor component of the smoothing parameter.
#' Default value is 0.7, unless either h is provided, in which case the default is 0.
#' @param h The constant component of the smoothing parameter. Default: 0.
#' @param deg Degree of polynomial to use. Default: 2.
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#' @examples
#' df <- mumps_uk_1986_1987
#' model <- lp_model(
#'   df,
#'   nn=0.7, kern="tcub"
#'   )
#' plot(model)
#'
#' @importFrom locfit locfit lp crit crit<-
#' @importFrom graphics par
#' @importFrom stats fitted
#'
#' @return a list of class lp_model with 6 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{pi}{fitted locfit object for pi}
#'   \item{eta}{fitted locfit object for eta}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [locfit::locfit()] for more information on the fitted locfit object
#'
#' @export
lp_model <- function(data, kern="tcub", nn=0, h=0, deg=2,
                     age_col="age",pos_col="pos", tot_col="tot", status_col="status") {
  if (all(nn==0) & all(h==0))  {
    # default nn from lp()
    nn <- 0.7
    h <- 0
  }
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  y <- pos/tot

  if(length(nn) > 1 & length(h) > 1) stop("Tuning both `h` and `nn` may lead to parameter idenfifiability issues, please fix one of the parameters instead")

  if(length(nn) > 1 || length(h)>1){
    best_param <- best_lp_params(
      data=data,
      nn = nn,
      h = h,
      family="binomial",
      kern=kern
    )

    nn <- best_param$nn
    h <- best_param$h
  }

  # print(paste0("nn: ", nn))
  # print(paste0("h: ", h))

  model$info  <- locfit(y~lp(age, deg=deg, nn=nn, h=h), family="binomial", kern=kern)
  model$nn <- nn
  model$h <- h
  model$deg <- deg
  model$kern <- kern
  model$eta <- locfit(y~lp(age, deg=deg, nn=nn, h=h), family="binomial", kern=kern, deriv=1)
  model$sp  <- fitted(model$info)
  model$foi <- fitted(model$eta)*fitted(model$info) # λ(a)=η′(a)π(a)
  model$df  <- list(age=age, pos=pos, tot=tot)

  class(model) <- "lp_model"
  model
}

# function to return the best parameter of local polynomial model using GCV
# nn - range of values for nearest neighbor
# h - range of values for constant bandwidth
# if both nn and h are given, select either best nn or h, whichever gives the lowest GCV
#' @import tidyr dplyr locfit
best_lp_params <- function(data, nn=0, h=0, kern="tcub",deg=2, family="binomial"){
  # helper function to get df and GCV
  summary.gcvplot <- function(object, ...){
    z <- cbind(object$df, object$values)
    dimnames(z) <- list(NULL, c("df", object$cri))
    z
  }

  # helper function which return parameter value which gives the lowest gcv
  get_best_gcv <- function(nn_vals=0, h_vals=0){
    is_nn <- length(nn_vals) > 1 # check if we are tuning for nn or h
    par_vals <- if (is_nn) nn_vals else h_vals
    # generate parameters matrix
    alpha <-
      if (is_nn)
        cbind(nn_vals, rep(h_vals, length(nn_vals)))
      else
        cbind(rep(nn_vals, length(h_vals)), h_vals)

    # compute gcv
    gcv_out <- gcvplot(
      pos / tot ~ age,
      deg = deg,
      kern = kern,
      family = family,
      alpha = alpha,
      data=data
    )

    gcv_out <- cbind(par_vals, gcv_out$values)
    best_idx <- which.min(gcv_out[,2])

    c(
      gcv_out[best_idx, 1],
      gcv_out[best_idx, 2]
    )
  }

  res <- list()

  if(length(nn) > 1 & length(h) > 1) stop("Tuning both `h` and `nn` may lead to parameter idenfifiability issues, please fix one of the parameters instead")
  if(length(nn) > 1) res[c("h", "nn", "nn_gcv")] <- c(h, get_best_gcv(nn_vals=nn, h_vals=h))
  if(length(h) > 1) res[c("nn", "h", "h_gcv")] <- c(nn, get_best_gcv(nn_vals=nn, h_vals=h))

  res
}


