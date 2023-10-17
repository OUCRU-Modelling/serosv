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

