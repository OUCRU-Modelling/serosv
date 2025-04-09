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
#' Refers to section 6.1.1
#' @param data the input data frame, must either have `age`, `pos`, `tot` columns (for aggregated data) OR `age`, `status` for (linelisting data)
#' @param k  degree of the model.
#' @param type name of method (Muench, Giffith, Grenfell).
#' @param link link function.
#'
#' @examples
#' data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ]
#' data$status <- data$seropositive
#' aggregated <- transform_data(data$age, data$seropositive, heterogeneity_col = "age")
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
#' Refers to section 6.1.2.
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
#'   \item{info}{fitted "glm" object}
#'   \item{sp}{seroprevalence}
#'   \item{foi}{force of infection}
#' @seealso [stats::glm()] for more information on the fitted glm object
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

