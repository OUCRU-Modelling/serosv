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
#' @description
#' Either `status` (treated as line-listing dataset) or  `pos` & `tot` (treated as aggregated dataset) must be provided
#'
#' @param age the age vector.
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
#' @param k  degree of the model.
#' @param type name of method (Muench, Giffith, Grenfell).
#' @param link link function.
#'
#' @examples
#' data <- parvob19_fi_1997_1998[order(parvob19_fi_1997_1998$age), ]
#' aggregated <- transform_data(data$age, data$seropositive)
#'
#' # fit with aggregated data
#' model <- polynomial_model(aggregated$t, pos = aggregated$pos, tot = aggregated$tot, type = "Muench")
#' # fit with linelisting data
#' model <- polynomial_model(data$age, status = data$seropositive, type = "Muench")
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
polynomial_model <- function(age,k,type, pos=NULL,tot=NULL,status=NULL, link = "log"){
  stopifnot("Values for either `pos & tot` or `status` must be provided" = !is.null(pos) & !is.null(tot) | !is.null(status) )

  model <- list()
  Pos <- NULL
  Neg <- NULL

  # check input whether it is line-listing or aggregated data
  if (!is.null(pos) & !is.null(tot)){
    Pos <- as.numeric(pos)
    Neg <- as.numeric(tot) - Pos
    model$datatype <- "aggregated"
  }else{
    Pos <- as.numeric(status)
    Neg <- 1 - Pos
    model$datatype <- "linelisting"
  }

  Age <- as.numeric(age)

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
#' @param age the age vector.
#' @param pos the positive count vector (optional if status is provided).
#' @param tot the total count vector (optional if status is provided).
#' @param status the serostatus vector (optional if pos & tot are provided).
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
#'   df$age, pos = df$pos, tot = df$tot,
#'   start=list(alpha=0.07,beta=0.1,gamma=0.03)
#'   )
#' plot(model)
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(age, start,pos=NULL,tot=NULL,status=NULL, fixed=list())
{
  stopifnot("Values for either `pos & tot` or `status` must be provided" = !is.null(pos) & !is.null(tot) | !is.null(status) )
  model <- list()
  age <- as.numeric(age)

  # check input whether it is line-listing or aggregated data
  if (!is.null(pos) & !is.null(tot)){
    pos <- as.numeric(pos)
    tot <- as.numeric(tot)
    model$datatype <- "aggregated"
  }else{
    pos <- as.numeric(status)
    tot <- rep(1, length(pos))
    model$datatype <- "linelisting"
  }

  # print(head(pos))
  # print(head(tot))

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

