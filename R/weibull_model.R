
plot_p_foi_wrt_dur <- function(model)
{
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(c(model$df_dur, model$df_age), {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="exposure", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(dur, model$sp, lwd=2)
    lines(dur, model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(model$foi), length.out=10), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
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

weibull_model <- function()
{
  model <- list()
  model$parameters <- c()
  model$fit <- function(df) {
    df_de <- detransform(df)
    with(df_de, {
      model$info <- glm(
        seropositive~log(dur),
        family=binomial(link="cloglog")
      )
      b0 <- coef(model$info)[1]
      b1 <- coef(model$info)[2]
      model$foi <- exp(b0)*b1*exp(log(dur))^(b1-1)
      model$sp <- 1-exp(-exp(b0)*dur^b1)
      model$df_dur <- df_de
      model$df_age <- df
      model
    })
  }
  model
}

library(dplyr)
transform_hcv <- function(df) {
  # df$dur <- ceiling(df$dur)
  df_agg <- df %>%
    group_by(dur) %>%
    summarize(
      pos = sum(seropositive),
      tot = n()
    )
  names(df_agg)[names(df_agg) == "dur"] <- "age"
  df_agg
}
hcv_df <- transform_hcv(hcv_be_2006)
hcv_df

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

model_wei <- weibull_model()

hcv_df %>%
  model_wei$fit() %>%
  plot_p_foi_wrt_dur()

fitted_model <- model_wei$fit(hcv_df)
fitted_model$info

