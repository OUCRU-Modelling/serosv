
#' Estimate force of infection
#'
#' @param t time (in this case age) vector
#' @param sp seroprevalence vector
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

#' Aggregate data
#'
#' Generate a dataframe with `t`, `pos` and `tot` columns from
#' `t` and `seropositive` vectors.
#'
#' @param data a data frame with columns for age and serostatus
#' @param status_col name of the column for serostatus
#' @param stratum_col name of the column to stratify by (default to "age")
#'
#' @examples
#' df <- hcv_be_2006
#' hcv_df <- transform_data(df, stratum_col="dur", status_col="seropositive")
#' hcv_df
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr summarize
#' @importFrom magrittr %>%
#'
#' @return dataframe in aggregated format
#' @export
transform_data <- function(data, stratum_col="age", status_col="status") {
  df <- age <- status <-  NULL

  if( all(c(stratum_col, status_col) %in% names(data)) ) {
    df <- data.frame(
      age = data[[stratum_col]],
      status = data[[status_col]]
    )
  }else{
    stop(paste0(
      "Data must have `",
      stratum_col,
      "`, `",status_col ,"` columns"
    ))
  }


  df_agg <- df %>%
    group_by(age) %>%
    summarize(
      pos = sum(status),
      tot = n()
    )

  df_agg
}

# utility function to check input data
# return:
# - type of data (either linelisting or aggregated)
# - preprocessed pos and tot columns
#' @importFrom assertthat assert_that
check_input <- function(data, pos_col="pos",tot_col="tot",status_col="status",
                        stratum_col = "age"){
  assert_that(
    is.data.frame(data),
    msg = "Input must be a data.frame or tibble"
  )

  age <- NULL
  pos <- NULL
  tot <- NULL
  type <- NULL


  if( all(c(stratum_col, pos_col, tot_col) %in% colnames(data)) ){
    age <- as.numeric(data[[stratum_col]])
    pos <- as.numeric(data[[pos_col]])
    tot <- as.numeric(data[[tot_col]])
    type <- "aggregated"
  }else if( all(c(stratum_col, status_col) %in% colnames(data)) ){
    age <- as.numeric(data[[stratum_col]])
    pos <- as.numeric(data[[status_col]])
    tot <- rep(1, length(data[[status_col]]))
    type <- "linelisting"
  }else{
    stop(paste0(
      "Data must have `",
      stratum_col,
      "`, `", pos_col, "`, `", tot_col,"` columns for aggregated data OR `",
      stratum_col,
      "`, `",status_col ,"` columns for linelisting data"
    ))
  }

  list(
    age = age, pos = pos, tot = tot, type = type
  )
}



