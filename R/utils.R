
#' Estimate force of infection
#'
#' @param t - time (in this case age) vector
#' @param sp - seroprevalence vector
#'
#' @importFrom stats approx
#'
#' @return computed foi vector
#' @export
est_foi <- function(t, sp)
{
  # handle duplicated t
  sp<-(sp[order(t)])[duplicated(sort(t))==F]
  t<-sort(unique(t))

  dsp <- diff(sp)/diff(t)
  foi <- approx(
    (t[-1]+t[-length(t)])/2,
    dsp,
    t[c(-1,-length(t))]
  )$y/(1-sp[c(-1,-length(t))])

  foi
}

#' Monotonize seroprevalence
#'
#' @param pos the positive count vector.
#' @param tot the total count vector.
#'
#' @importFrom stats approx
#'
#' @return computed list of 2 items pai1 for original values and pai2 for monotonized value
#' @export
pava<- function(pos=pos,tot=rep(1,length(pos)))
{
  gi<- pos/tot
  pai1 <- pai2 <- gi
  N <- length(pai1)
  ni<-tot
  for(i in 1:(N - 1)) {
    if(pai2[i] > pai2[i + 1]) {
      pool <- (ni[i]*pai1[i] + ni[i+1]*pai1[i + 1])/(ni[i]+ni[i+1])
      pai2[i:(i + 1)] <- pool
      k <- i + 1
      for(j in (k - 1):1) {
        if(pai2[j] > pai2[k]) {
          pool.2 <- sum(ni[j:k]*pai1[j:k])/(sum(ni[j:k]))
          pai2[j:k] <- pool.2
        }
      }
    }
  }
  return(list(pai1=pai1,pai2=pai2))
}

#' Generate a dataframe with `t`, `pos` and `tot` columns from
#' `t` and `seropositive` vectors.
#'
#' @param t the time vector.
#' @param spos the seropositive vector.
#'
#' @examples
#' df <- hcv_be_2006
#' hcv_df <- transform_data(df$dur, df$seropositive)
#' hcv_df
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr summarize
#' @import magrittr
#'
#' @return dataframe in aggregated format
#' @export
transform_data <- function(t, spos) {
  df <- data.frame(t, spos)
  df_agg <- df %>%
    group_by(t) %>%
    summarize(
      pos = sum(spos),
      tot = n()
    )
  df_agg
}

