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
#'Refers to section 6.1.1
#' @param age the age vector
#' @param positive the positive vector
#' @param negative the negative vector
#' @param k  degree of the model
#' @param type name of method
#'
#' @examples
#' a <- hav_bg_1964
#' a$neg <- a$tot -a$pos
#' model <- polynomial_model(a$age,a$pos,a$neg, type = "Muench")
#'
polynomial_model <- function(age,positive,negative,k,type){
  Age <- as.numeric(age)
  Pos <- as.numeric(positive)
  Neg <- as.numeric(negative)
  df <- data.frame(cbind(Age, Pos,Neg))
  model <- list()
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
  model$info <- glm(age(k), family=binomial(link="log"),df)
  X <- X(Age, k)
  model$sp <- 1 - model$info$fitted.values
  model$foi <- X%*%model$info$coefficients
  model$df <- list(Age=Age, Pos=Pos, Tot= Pos + Neg)
  class(model) <- "polynomial_model"
  model
}

#' The Farrington (1990) model.
#'
#' Refers to section 6.1.2.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.
#' @param start Named list of vectors or single vector.
#' Initial values for optimizer.
#' @param fixed Named list of vectors or single vector.
#' Parameter values to keep fixed during optimization.
#'
#' @examples
#' df <- rubella_uk_1986_1987
#' model <- farrington_model(
#'   df$age, df$pos, df$tot,
#'   start=list(alpha=0.07,beta=0.1,gamma=0.03)
#'   )
#' plot(model)
#'
#' @importFrom stats4 mle
#'
#' @export
farrington_model <- function(age, pos, tot, start, fixed=list())
{
  farrington <- function(alpha,beta,gamma) {
    p=1-exp((alpha/beta)*age*exp(-beta*age)
            +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)-gamma*age)
    ll=pos*log(p)+(tot-pos)*log(1-p)
    return(-sum(ll))
  }

  model <- list()

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

