#' Plot the prevalence and force of infection against age.
#'
#' Refers to section 6.1.2.
#'
#' @param model a model returned from this package's function.
#'
#' @examples
#' model <- weibull_model()
#' hcv_df %>%
#'   model$fit() %>%
#'   plot_p_foi_wrt_dur()
#'
#' @export
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

#' The Weibull model.
#'
#' Refers to section 6.1.2.
#'
#' @examples
#' hbe_best <- hav_be_1993_1994 %>%
#'   find_best_fp_powers(p=seq(-2,3,0.1), mc=F, degree=2, link="cloglog")
#' hbe_best
#'
#' @export
weibull_model_t <- function(t, spos)
{
  model <- list()

  model$info <- glm(
    spos~log(t),
    family=binomial(link="cloglog")
    )
  b0 <- coef(model$info)[1]
  b1 <- coef(model$info)[2]
  model$foi <- exp(b0)*b1*exp(log(t))^(b1-1)
  model$sp <- 1-exp(-exp(b0)*t^b1)
  model$df <- data.frame(t=t, spos=spos)

  class(model) <- "weibull_model"
  model
}

plot.weibull_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  df_ <- transform_hcv(x$df)

  with(c(x$df, df_), {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="exposure", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(t, model$sp, lwd=2)
    lines(t, model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(model$foi), length.out=10), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

df <- hcv_be_2006[order(hcv_be_2006$dur), ]
model <- weibull_model_t(
  t=df$dur,
  spos=df$seropositive)
plot(model)

