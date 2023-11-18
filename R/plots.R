#### Polynomial model ####
#' plot() overloading for polynomial model
#'
#' @param x the polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.polynomial_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#### Non-linear ####

#### Farrington model ####
#' plot() overloading for Farrington model
#'
#' @param x the Farrington model object.
#' @param ... arbitrary params.
#'
#' @export
plot.farrington_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#### Weibull model ####
#' plot() overloading for Weibull model
#'
#' @param x the Weibull model object.
#' @param ... arbitrary params.
#'
#' @export
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

#### Fractional polynomial model ####
#' plot() overloading for fractional polynomial model
#'
#' @param x the fractional polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.fp_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age[c(-1,-length(age))], x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#### Non-parametric ####

#### Local polynomial model ####
#' plot() overloading for local polynomial model
#'
#' @param x the local polynomial model object.
#' @param ... arbitrary params.
#'
#' @export
plot.lp_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      age,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(age, x$sp, lwd=2)
    lines(age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

#### GCV values ####
#' Plotting GCV values with respect to different nn-s and h-s parameters.
#'
#' Refers to section 7.2.
#'
#' @param age the age vector.
#' @param pos the pos vector.
#' @param tot the tot vector.#'
#' @param nn_seq Nearest neighbor sequence.
#' @param h_seq Smoothing parameter sequence.
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' df <- mumps_uk_1986_1987
#' plot_gcv(
#'   df$age, df$pos, df$tot,
#'   nn_seq = seq(0.2, 0.8, by=0.1),
#'   h_seq = seq(5, 25, by=1)
#' )
#'
#' @import locfit
#' @import graphics
#'
#' @export
plot_gcv <- function(age, pos, tot, nn_seq, h_seq, kern="tcub", deg=2) {
  par(mfrow=c(1,2),lwd=2,las=1,cex.axis=1,cex.lab=1.1,mgp=c(2,0.5, 0),mar=c(3.1,3.1,3.1,3))

  y <- pos/tot
  res = cbind(nn_seq, summary(gcvplot(y~age, family="binomial", alpha=nn_seq)))
  plot(res[,1],res[,3],type="n",xlab="nn (% Neighbors)",ylab=" ")
  lines(res[,1],res[,3])
  mtext(side=2, "GCV", las=3, line=2.4, cex=0.9)

  h_seq_ <- cbind(rep(0, length(h_seq)), h_seq)
  res=cbind(h_seq_[,2],summary(gcvplot(y~age,family="binomial",alpha=h_seq_)))
  plot(res[,1],res[,3],type="n",xlab="h (Bandwidth)",ylab=" ")
  lines(res[,1],res[,3])
  mtext(side=2, "GCV", las=3, line=3, cex=0.9)
}
