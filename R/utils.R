rename_col <- function(df, old, new)
{
  names(df)[names(df) == old] <- new
  df
}

estimate_foi <- function(age, sp)
{
  dsp <- diff(sp)/diff(age)
  foi <- approx(
    (age[-1]+age[-length(age)])/2,
    dsp,
    age[c(-1,-length(age))]
  )$y/(1-sp[c(-1,-length(age))])

  foi
}

X <- function(age, degree) {
  X_matrix <- matrix(rep(1, length(age)), ncol = 1)
  if (degree > 1) {
    for (i in 2:degree) {
      X_matrix <- cbind(X_matrix, i * age^(i-1))
    }
  }
  -X_matrix
}


#' Transform the hcv_be_2006 dataframe.
#'
#' @param df the hcv_be_2006 dataframe.
#'
#' @examples
#' hcv_df <- transform_hcv(hcv_be_2006)
#' hcv_df
#'
#' @importFrom dplyr group_by
#'
#' @export
transform_hcv <- function(t, spos) {
  # df$t <- ceiling(df$t)
  df <- data.frame(t, spos)
  df_agg <- df %>%
    group_by(t) %>%
    summarize(
      pos = sum(spos),
      tot = n()
    )
  names(df_agg)[names(df_agg) == "t"] <- "age"
  df_agg
}

detransform <- function(hcv_df) {
  hcv_dur <- c()
  hcv_sp <- c()

  for (i in 1:nrow(hcv_df)) {
    row <- hcv_df[i, ]
    hcv_dur <- c(hcv_dur, rep(row$age, row$tot))
    if (row$pos > 0) {
      hcv_sp <- c(hcv_sp, rep(1, row$pos))
    }
    if (row$tot - row$pos > 0) {
      hcv_sp <- c(hcv_sp, rep(0, row$tot - row$pos))
    }
  }
  data.frame(
    dur = hcv_dur,
    seropositive = hcv_sp
  )
}

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
plot_p_foi_wrt_age <- function(model)
{
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(model$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, model$sp, lwd=2)
    lines(age, model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(model$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}
