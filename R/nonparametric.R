
#' An isotonic regression estimator using Pavitt's algorithm.
#'
#' This function takes a vector of data as input and performs a Pavitt’s algorithm on it.
#' The Pavitt’s algorithm is used to smooth out a time series data by replacing any
#' decreasing values with the average of the two adjacent values.
#'
#' The function does the following:
#' 1.	Create two copies of the input data vector.
#' 2.	Iterate through the vector from the first element to the second to last element.
#' 3.	If the current element in the second copy is greater than the next element,
#'     replace the current and next elements in the second copy with their average.
#' 4.	Starting from the current element, iterate backwards through the first copy
#'     until the first element.
#' 5.	For each element in the first copy that is greater than the current element
#'     in the second copy, replace it with the average of all the elements between it
#'     and the current element in the second copy.
#' 6.	Return the smoothed data vector.
pavitt <- function(datai)
{
  pai1 <- pai2 <- datai
  N <- length(pai1)
  for(i in 1:(N - 1)) {
    if(pai2[i] > pai2[i + 1]) {
      pool <- (pai1[i] + pai1[i + 1])/2
      pai2[i:(i + 1)] <- pool
      k <- i + 1
      for(j in (k - 1):1) {
        if(pai2[j] > pai2[k]) {
          pool.2 <- sum(pai1[j:k])/length(pai1[j:k])
          pai2[j:k] <- pool.2
        }
      }
    }
  }
  pai2
}

library(locfit)
lp_model <- function(kern, alpha) {
  model <- list()
  model$parameters <- list(kern, alpha)
  model$fit <- function(df) {
    with(c(df, model$parameters), {
      sp <- pos/tot
      lpfit <- locfit(sp~age, deg=2, family="binomial", alpha=alpha, kern=kern)
      lpfitd1 <- locfit(sp~age, deriv=1, deg=2, family="binomial", kern=kern)
      model$sp <- fitted(lpfit)
      model$foi <- fitted(lpfitd1) * fitted(lpfit)
      model$df <- df
      model
    })
  }
  model
}

# alpha=0.7, deg=2, kernel=“tcub”
# params = c(kern="tcub", alpha=c(0, 15))
model <- lp_model(kern="tcub", alpha=0.7)
mumps_uk_1986_1987 %>%
  model$fit() %>%
  plot_p_foi_wrt_age()

