#' A local polynomial model.
#'
#' Refers to section 7.1. and 7.2.
#'
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#'
#' @param nn Nearest neighbor component of the smoothing parameter.
#' Default value is 0.7, unless either h is provided, in which case the default is 0.
#'
#' @param h The constant component of the smoothing parameter. Default: 0.
#'
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' model <- lp_model()
#' mumps_uk_1986_1987 %>%
#'   model$fit() %>%
#'   plot_p_foi_wrt_age()
#'
#' @importFrom locfit locfit
#'
#' @export
# library(locfit)
lp_model <- function(kern="tcub", nn=0, h=0, deg=2) {
  if (is_missing(nn) & is_missing(h)) {
    nn <- 0.7 # default nn from lp()
  }
  model <- list()
  model$parameters <- list(kern, nn, h, deg)
  model$fit <- function(df) {
    with(c(df, model$parameters), {
      y <- pos/tot
      estimator <- lp(age, deg=deg, nn=nn, h=h)
      lpfit   <- locfit(y~estimator, family="binomial", kern=kern)
      lpfitd1 <- locfit(y~estimator, family="binomial", kern=kern, deriv=1)
      model$sp <- fitted(lpfit)
      model$foi <- fitted(lpfitd1) * fitted(lpfit)
      model$df <- df
      model
    })
  }
  model
}
windows(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(3, 0.5, 0))

model <- lp_model(h=14)
mumps_uk_1986_1987 %>%
  model$fit() %>%
  plot_p_foi_wrt_age()

#' Plotting GCV values with respect to different nn-s and h-s parameters.
#'
#' Refers to section 7.2.
#'
#' @param df The dataframe. MUST have `age`, `pos`, `tot` columns.
#'
#' @param nn_seq Nearest neighbor sequence.
#'
#' @param h_seq Smoothing parameter sequence.
#'
#' @param kern Weight function, default = "tcub".
#' Other choices are "rect", "trwt", "tria", "epan", "bisq" and "gauss".
#' Choices may be restricted when derivatives are required;
#' e.g. for confidence bands and some bandwidth selectors.
#'
#' @param deg Degree of polynomial to use. Default: 2.
#'
#' @examples
#' plot_gcv(
#'   mumps_uk_1986_1987,
#'   nn_seq = seq(0.2, 0.8, by=0.1),
#'   h_seq = seq(5, 25, by=1)
#' )
#'
#' @importFrom locfit gcvplot
#'
#' @export
plot_gcv <- function(df, nn_seq, h_seq, kern="tcub", deg=2) {
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
