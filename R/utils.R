rename_col <- function(df, old, new)
{
  names(df)[names(df) == old] <- new
  df
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
