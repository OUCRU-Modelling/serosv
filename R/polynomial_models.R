X <- function(t, degree) {
  X_matrix <- matrix(rep(1, length(t)), ncol = 1)
  if (degree > 1) {
    for (i in 2:degree) {
      X_matrix <- cbind(X_matrix, i * t^(i-1))
    }
  }
  -X_matrix
}

#' Polynomial models
#'
#' @description Fit age-stratified seroprevalence data to serocatalytic models formulated as polynomials.
#' Supported models includes Muench (1934), Griffith, Grenfell & Anderson (1985)
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
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param k  degree of the model.
#' @param type name of method (Muench, Giffith, Grenfell).
#' @param link link function.
#'
#' @examples
#' data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ]
#' data$status <- data$seropositive
#' aggregated <- transform_data(data$age, data$seropositive, stratum_col = "age")
#'
#' # fit with aggregated data
#' model <- polynomial_model(aggregated, type = "Muench")
#' # fit with linelisting data
#' model <- polynomial_model(data, type = "Muench")
#' plot(model)
#'
#' @return a list of class polynomial_model with 5 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "glm" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#'
#' @export
polynomial_model <- function(data, k,type, link = "log"){
  model <- list()
  data <- check_input(data)
  model$datatype <- data$type

  Age <- data$age
  Pos <- data$pos
  Neg <- data$tot - Pos

  df <- data.frame(cbind(Age, Pos,Neg))
  if(missing(k)){
    k <- switch(type,
                "Muench" = 1 ,
                "Griffith" = 2,
                "Grenfell" = 3)}
  age <- function(k){
    if(k>1){
      formula<- paste0("I","(",paste("Age", 2:k,sep = "^"),")",collapse = "+")
      paste0("cbind(Neg,Pos)"," ~","-1+Age+",formula)
    } else {
      paste0("cbind(Neg,Pos)"," ~","-1+Age")
    }
  }
  model$info <- glm(age(k), family=binomial(link=link),df)
  X <- X(Age, k)
  model$sp <- 1 - model$info$fitted.values
  model$foi <- X%*%model$info$coefficients
  model$df <- list(age=Age, pos=Pos, tot= Pos + Neg)
  class(model) <- "polynomial_model"
  model
}

#' The Farrington (1990) model.
#'
#' @description Fit age-stratified seroprevalence data using the Farrington (1990) model, which assumes the
#' force of infection increases linearly with age and subsequently decreases exponentially.
#'
#' @details
#' The force of infection is defined as followed
#'
#' \deqn{
#' \lambda(a) = (\alpha a - \gamma)e^{-\beta a} + \gamma
#' }
#' Where \eqn{\gamma} is called the long term residual for FOI,
#' as \eqn{a \rightarrow \infty} , \eqn{\lambda (a) \rightarrow \gamma}
#'
#' The seroprevalence can thus be estimated using the non-linear model
#' \deqn{
#'  \pi(a) = 1 - exp\{ \frac{\alpha}{\beta}ae^{-\beta a} +
#'  \frac{1}{\beta}(\frac{\alpha}{\beta} -
#'  \gamma)(e^{-\beta a} - 1) -\gamma a \}
#' }
#'
#' Refer to section 6.1.2. of the the book by Hens et al. (2012) for further details.
#'
#' @references
#' Hens, Niel, Ziv Shkedy, Marc Aerts, Christel Faes, Pierre Van Damme,
#' and Philippe Beutels. 2012. Modeling Infectious Disease Parameters Based on
#' Serological and Social Contact Data: A Modern Statistical Perspective.
#' tatistics for Biology and Health. Springer New York.
#' \doi{https://doi.org/10.1007/978-1-4614-4072-7}.
#'
#'
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param start Named list of vectors or single vector.
#' Initial values for optimizer.
#' @param fixed Named list of vectors or single vector.
#' Parameter values to keep fixed during optimization.
#'
#' @return a list of class farrington_model with 5 items
#'   \item{datatype}{type of datatype used for model fitting (aggregated or linelisting)}
#'   \item{df}{the dataframe used for fitting the model}
#'   \item{info}{fitted "mle" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [stats4::mle()] for more information on the fitted mle object
#'
#' @examples
#' df <- rubella_uk_1986_1987
#' model <- farrington_model(
#'   df,
#'   start=list(alpha=0.07,beta=0.1,gamma=0.03)
#'   )
#' plot(model)
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(data, start, fixed=list())
{
  model <- list()

  # check input whether it is line-listing or aggregated data
  data <- check_input(data)
  age <- data$age
  pos <- data$pos
  tot <- data$tot
  model$datatype <- data$type

  farrington <- function(alpha,beta,gamma) {
    p=1-exp((alpha/beta)*age*exp(-beta*age)
            +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)-gamma*age)
    ll=pos*log(p)+(tot-pos)*log(1-p)
    return(-sum(ll))
  }

  model$info <- mle(farrington, fixed=fixed, start=start)
  alpha <- model$info@coef[1]
  beta  <- model$info@coef[2]
  gamma <- model$info@coef[3]
  model$sp <- 1-exp(
    (alpha/beta)*age*exp(-beta*age)
    +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
    -gamma*age)
  model$foi <- (alpha*age-gamma)*exp(-beta*age)+gamma
  model$df <- list(age=age, pos=pos, tot=tot)

  class(model) <- "farrington_model"
  model
}

