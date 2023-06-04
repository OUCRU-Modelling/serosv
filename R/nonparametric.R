#' A local polynomial model.
#'
#' Refers to section 7.1. and 7.2.
#'
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#' @param nn Nearest neighbor component of the smoothing parameter.
#' Default value is 0.7, unless either h is provided, in which case the default is 0.
#' @param h The constant component of the smoothing parameter. Default: 0.
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' df <- mumps_uk_1986_1987
#' model <- lp_model(
#'   df$age, df$pos, df$tot,
#'   nn=0.7, kern="tcub"
#'   )
#' plot(model)
#'
#' @importFrom locfit locfit
#'
#' @export
# library(locfit)
lp_model <- function(age, pos, tot, kern="tcub", nn=0, h=0, deg=2) {
  if (missing(nn) & missing(h)) {
    nn <- 0.7 # default nn from lp()
  }

  model <- list()

  y <- pos/tot
  estimator <- lp(age, deg=deg, nn=nn, h=h)
  model$pi  <- locfit(y~estimator, family="binomial", kern=kern)
  model$eta <- locfit(y~estimator, family="binomial", kern=kern, deriv=1)
  model$sp  <- fitted(model$pi)
  model$foi <- fitted(model$eta)*fitted(model$pi) # λ(a)=η′(a)π(a)
  model$df  <- list(age=age, pos=pos, tot=tot)

  class(model) <- "lp_model"
  model
}

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
#' @importFrom locfit gcvplot
#'
#' @export
plot_gcv <- function(age, pos, tot, nn_seq, h_seq, kern="tcub", deg=2) {
  par(mfrow=c(1,2),lwd=2,las=1,cex.axis=1,cex.lab=1.1,mgp=c(2,0.5, 0),mar=c(3.1,3.1,3.1,3))

  with(df, {
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
  })
}

