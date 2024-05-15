# rename_col <- function(df, old, new)
# {
#   names(df)[names(df) == old] <- new
#   df
# }

est_foi <- function(t, sp)
{
  dsp <- diff(sp)/diff(t)
  foi <- approx(
    (t[-1]+t[-length(t)])/2,
    dsp,
    t[c(-1,-length(t))]
  )$y/(1-sp[c(-1,-length(t))])

  foi
}

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
#' hcv_df <- transform(df$dur, df$seropositive)
#' hcv_df
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr n
#' @importFrom dplyr summarize
#' @import magrittr
#'
#' @export
transform <- function(t, spos) {
  df <- data.frame(t, spos)
  df_agg <- df %>%
    group_by(t) %>%
    summarize(
      pos = sum(spos),
      tot = n()
    )
  df_agg
}

# detransform <- function(hcv_df) {
#   hcv_dur <- c()
#   hcv_sp <- c()
#
#   for (i in 1:nrow(hcv_df)) {
#     row <- hcv_df[i, ]
#     hcv_dur <- c(hcv_dur, rep(row$age, row$tot))
#     if (row$pos > 0) {
#       hcv_sp <- c(hcv_sp, rep(1, row$pos))
#     }
#     if (row$tot - row$pos > 0) {
#       hcv_sp <- c(hcv_sp, rep(0, row$tot - row$pos))
#     }
#   }
#   data.frame(
#     dur = hcv_dur,
#     seropositive = hcv_sp
#   )
# }

#' Plot the prevalence and force of infection against age.
#'
#' Refers to section 6.1.
#'
#' @param model a model returned from this package's function.
#'
#' @examples
#' model <- polynomial_model(degree=3)
#' hav_be_1993_1994 %>%
#'   model$fit() %>%
#'   plot_p_foi_wrt_age()
# plot_p_foi_wrt_age <- function(model)
# {
#   CEX_SCALER <- 4 # arbitrary number for better visual
#
#   with(model$df, {
#     par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
#     plot(
#       age,
#       pos/tot,
#       cex=CEX_SCALER*tot/max(tot),
#       xlab="age", ylab="seroprevalence",
#       xlim=c(0, max(age)), ylim=c(0,1)
#     )
#     lines(age, model$sp, lwd=2)
#     lines(age, model$foi, lwd=2, lty=2)
#     axis(side=4, at=round(seq(0.0, max(model$foi), length.out=3), 2))
#     mtext(side=4, "force of infection", las=3, line=2)
#   })
# }
