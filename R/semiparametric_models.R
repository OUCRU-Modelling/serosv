#' Penalized Spline model
#'
#' @description Fit age-specific seroprevalence to a semi-parametric model where
#' predictor is modeled with penalized splines. The penalized splines can be estimated by either
#' (1) penalized likelihood framework or (2) mixed model framework
#'
#' @details
#' In the semi-parametric model, the predictor is formulated as a penalized spline
#' with truncated power basis functions of degree \eqn{p}
#' and fixed knots \eqn{\kappa_1,\cdots, \kappa_k} as followed
#'
#' \deqn{
#' \eta(a_i) = \beta_0 + \beta_1a_i + \cdots + \beta_p a_i^p + \Sigma_{k=1}^ku_k(a_i - \kappa_k)^p_+
#' }{}
#'
#' Where:
#' \deqn{
#' (a_i - \kappa_k)^p_+ = \begin{cases}
#' 0, & a_i \le \kappa_k \\
#' (a_i - \kappa_k)^p, & a_i > \kappa_k
#' \end{cases}
#' }{}
#'
#' FOI can then be derived by
#'
#' \deqn{\hat{\lambda}(a_i) = [\hat{\beta_1} , 2\hat{\beta_2}a_i, \cdots,
#' p \hat{\beta} a_i ^{p-1} + \Sigma^k_{k=1} p \hat{u}_k(a_i - \kappa_k)^{p-1}_+] \delta(\hat{\eta}(a_i))
#' }{}
#'
#' Where \eqn{\delta(.)} is determined by the link function used in the model
#'
#' In matrix annotation, the mean structure model for \eqn{\eta(a_i)}{} becomes
#' \deqn{\eta = \textbf{X}\beta + \textbf{Zu}}{}
#'
#' Where \eqn{\eta = [\eta(a_i) \cdots \eta(a_N) ]^T}{}, \eqn{\beta = [\beta_0 \beta_1 \cdots \beta_p]^T}{},
#' and \eqn{\textbf{u} = [u_1 u_2 \cdots u_k]^T}{} are the regression with corresponding design matrices
#'
#' \deqn{
#' \textbf{X} = \begin{bmatrix}
#' 1 & a_1 & a_1^2 & \cdots & a_1^p \\
#' 1 & a_2 & a_2^2 & \cdots & a_2^p \\
#' \vdots & \vdots & \vdots & \dots & \vdots \\
#' 1 & a_N & a_N^2 & \cdots & a_N^p
#' \end{bmatrix}, \textbf{Z} = \begin{bmatrix}
#' (a_1 - \kappa_1 )_+^p & (a_1 - \kappa_2 )_+^p & \dots & (a_1 - \kappa_k)_+^p \\
#' (a_2 - \kappa_1 )_+^p & (a_2 - \kappa_2 )_+^p & \dots & (a_2 - \kappa_k)_+^p \\
#' \vdots & \vdots & \dots & \vdots \\
#' (a_N - \kappa_1 )_+^p & (a_N - \kappa_2 )_+^p & \dots & (a_N - \kappa_k)_+^p
#' \end{bmatrix}
#' }{}
#'
#' Under \bold{penalized likelihood framework}, the model is fitted by maximizing
#' the following likelihood
#'
#' \deqn{
#' \phi^{-1}[y^T(\textbf{X}\beta + \textbf{Zu} ) -  \textbf{1}^Tc(\textbf{X}\beta + \textbf{Zu} )] - \frac{1}{2}\lambda^2
#' \begin{bmatrix} \beta \\ \textbf{u} \end{bmatrix}^T D\begin{bmatrix} \beta \\ \textbf{u} \end{bmatrix}
#' }{}
#'
#' Where:
#'  \itemize{
#'     \item \eqn{X\beta + Zu} is the predictor
#'     \item \eqn{D} is a known semi-definite penalty matrix
#'     \item \eqn{y} is the response vector
#'     \item \eqn{\mathbf{1}} the unit vector, \eqn{c(.)} is determined by the link function used
#'     \item \eqn{\lambda} is the smoothing parameter (larger values -> smoother curves)
#'     \item \eqn{\phi} is the overdispersion parameter and equals 1 if there is no overdispersion
#'   }
#'
#' Under the \bold{mixed model} framework,
#' the model instead treats the coefficients \eqn{\textbf{u}}{} in the likelihood formulation
#' as random effects with \eqn{\textbf{u} \sim N(\textbf{0}, \boldsymbol{\sigma}^2_u \textbf{I})}{}
#'
#' Refer to section 8.1 and 8.2 of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#' @param data the input data frame, must either have columns for `age`, `pos`, `tot` (for aggregated data) OR
#' columns for `age`, `status` (for linelisting data)
#' @param s smoothing basis to use
#' @param sp smoothing parameter
#' @param link link function to use
#' @param framework which approach to fit the model ("pl" for penalized likelihood framework, "glmm" for generalized linear mixed model framework)
#' @param age_col name of the `age` column (default age_col="age").
#' @param pos_col name of the `pos` column (default pos_col="pos").
#' @param tot_col name of the `tot` column (default tot_col="tot").
#' @param status_col name of the `status` column (default status_col="status").
#'
#' @importFrom mgcv gam gamm
#' @importFrom stats binomial
#'
#' @return a list of class penalized_spline_model with 6 attributes
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{framework}{either pl or glmm}
#'   \item{info}{fitted "gam" model when framework is pl or "gamm" model when framework is glmm}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'
#' @seealso [mgcv::gam()], [mgcv::gamm()] for more information the fitted gam and gamm model
#'
#' @export
#'
#' @examples
#' data <- parvob19_be_2001_2003
#' data$status <- data$seropositive
#' model <- penalized_spline_model(data, framework="glmm")
#' model$info$gam
#' plot(model)
penalized_spline_model <- function(data,
                                   age_col="age",pos_col="pos", tot_col="tot", status_col="status",
                                   s = "bs", link = "logit", framework = "pl", sp = NULL){
  model <- list()

  data <- check_input(data, stratum_col=age_col,pos_col=pos_col, tot_col=tot_col, status_col=status_col)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  # s <- mgcv:::s
  neg <- tot - pos

  if (framework == "pl"){
    model$info <- if(data$type == "aggregated"){
        mgcv::gam(cbind(pos, neg) ~ s(age, bs = s, sp=sp), family = binomial(link = link))
      }else{
        mgcv::gam(pos ~ s(age, bs = s, sp=sp), family = binomial(link = link))
      }

    model$sp <- model$info$fitted.values
  }else if(framework == "glmm"){
    model$info <- if(data$type == "aggregated"){
        mgcv::gamm(cbind(pos, neg) ~ s(age, bs = s, sp=sp), family = binomial(link = link))
      }else{
        mgcv::gamm(pos ~ s(age, bs = s, sp=sp), family = binomial(link = link))
      }

    model$sp <- model$info$gam$fitted.values
  }else{
    stop(paste0('Invalid value for framework. Expected "pl" or "glmm", got ', framework))
  }

  # aggregate data after fitting for plotting
  model$df <- data.frame(age=age, pos = pos, tot = tot)
  model$foi <- est_foi(age, model$sp)
  model$framework <- framework

  class(model) <- "penalized_spline_model"
  model
}
