#' The Weibull model.
#'
#' Refers to section 6.1.2.
#'
#' @examples
#' df <- hcv_be_2006[order(hcv_be_2006$dur), ]
#' model <- weibull_model(
#'   t=df$dur,
#'   spos=df$seropositive
#'   )
#' plot(model)
#'
#' @export
weibull_model <- function(t, spos)
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

  df_ <- transform(x$df$t, x$df$spos)
  names(df_)[names(df_) == "t"] <- "exposure"

  with(c(x$df, df_), {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      exposure,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="exposure", ylab="seroprevalence",
      xlim=c(0, max(exposure)), ylim=c(0,1)
    )
    lines(t, model$sp, lwd=2)
    lines(t, model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(model$foi), length.out=10), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}
