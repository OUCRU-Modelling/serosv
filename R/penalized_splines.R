
# Penalized Likelihood Framework
plot.plf <- function(model, round.by=NULL, x.title=NULL, y.title=NULL, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual

  tdf <- transform(model$df$t, model$df$spos, round.by)
  names(tdf)[names(tdf) == "t"] <- "tt"

  with(c(model$df, tdf), {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      tt,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab=ifelse(is.null(x.title), "exposure", x.title),
      ylab=ifelse(is.null(y.title), "seroprevalence", y.title),
      xlim=c(0, max(tt)),
      ylim=c(min(model$foi, 0), 1),
      ...
    )
    lines(t, model$sp, lwd=2)
    lines(t[c(-1, -length(t))], model$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(
      min(model$foi, 0),
      ifelse(is.infinite(max(model$foi)), 1, max(model$foi)),
      length.out=5), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

### library Semipar: Ruppert et al. (2003)
library(SemiPar)
ps <- function(x, y, k=3, deg=2) {
  basis_fn <- function(a, b) { (a-b)^deg*(a>b) }
  expit <- function(a) { exp(a)/(1+exp(a)) }

  model <- list()

  n <- length(x)
  # Workaround as `spm()` get values from environment rather than passed in
  assign("y_", y, envir = globalenv())
  assign("x_", x, envir = globalenv())
  assign("k_", k, envir = globalenv())
  assign("d_", deg, envir = globalenv())
  # ----
  model$info <- spm(
    y_~f(x_, degree=d_, knots=default.knots(x_, k_), basis="trunc.poly"),
    family="binomial",
  )
  # clean up
  # for (i in c("y_", "x_", "k_", "d_")) {
  #   rm(i, envir = globalenv())
  # }
  # ----
  u_k <- outer(x, default.knots(x, k), basis_fn)*t(
    matrix(
      rep(model$info$fit$coefficients$random$dummy.group.vec.Handan, n),
      ncol=n
    )
  )
  model$semi_param_part <- apply(u_k, 1, sum)
  model$param_part <- (cbind(1, sapply(1:n, function(i) x^i))[, 1:(1+deg)]) %*%
    model$info$fit$coefficients$fixed
  model$sp <- expit(model$param_part + model$semi_param_part)
  model$foi <- est_foi(x, model$sp)
  model$df <- list(t=x, spos=y)
  class(model) <- "plf"
  model
}
#
# df <- vzv_be_2001_2003
# filter <- (df$age>0.5)&(df$age<40)&(!is.na(df$age))&!is.na(df$parvores)
# sdf <- df[filter,]
# sdf <- sdf[order(sdf$age), ]
# tps_fit <- tps(
#   x=sdf$age,
#   y=sdf$parvores,
#   k=3
# )
# plot(tps_fit, round.by=0)
# tps_fit$info$info$pen$knots

#----------------------------------------------------------------------
ss <- function(x, y, d=5, link="logit") {
  model <- list()
  model$info <- gam::gam(y~s(x, df=d),family=binomial(link=link))
  model$sp   <- model$info$fitted.values
  model$foi  <- est_foi(x, model$sp)
  model$df   <- list(t=x, spos=y)
  class(model) <- "plf"
  model
}

df <- parvob19_be_2001_2003
m <- gam::gam(df$seropositive~s(df$age, df=4),family=binomial(link="logit"))
f <- est.foi(df$age, m$fitted.values)
f$y

ss_fit <- ss(
  x=df$age,
  y=df$seropositive,
)
plot(ss_fit, round.by=0)

#----------------------------------------------------------------------
bss <- function(x, y, lambda=20, k=20, deg=3, m=2, link="logit") {
  model <- list()
  model$info <- pspline.fit(
    x.predicted=x, x.var=x, response=y, lambda=lambda, ps.intervals=k,
    degree=deg, order=m, link="logit", family="binomial"
  )
  model$sp   <- model$info$summary.predicted[,3]
  model$foi  <- est_foi(x, model$sp)
  model$df   <- list(t=x, spos=y)

  class(model) <- "plf"
  model
}

# source("R/pspline.R")
#
# bss_fit <- bss(
#   x=sdf$age,
#   y=sdf$parvores,
# )
# plot(bss_fit, round.by = 0)
#
#

#----------------------------------------------------------------------
# detach(package:gam)
library(mgcv)
crs <- function(x, y, link="logit") {
  model <- list()
  model$info <- mgcv::gam(y~s(x, bs="cr"),family=binomial(link=link))
  model$sp   <- predict(model$info, type="response")
  model$foi  <- est_foi(x, model$sp)
  model$df   <- list(t=x, spos=y)

  class(model) <- "plf"
  model
}


# df <- vzv_be_2001_2003
# df <- read.table("/Users/thanhlongb/Projects/oucru/the-book/RCodeBook/Chapter4/VZV-B19-BE.dat",header=T)
# filter <- (df$age>0.5)&(df$age<76)&(!is.na(df$age))&!is.na(df$parvores)
# sdf <- df[filter, ]
# sdf <- sdf[order(sdf$age), ]
# crs_fit <- crs(
#   x=sdf$age,
#   y=sdf$parvores,
# )
# plot(crs_fit, round.by=0)

#----------------------------------------------------------------------

tps <- function(x, y, link="logit") {
  model <- list()
  model$info <- mgcv::gam(y~s(x, bs="tp"),family=binomial(link=link))
  model$sp   <- predict(model$info, type="response")
  model$foi  <- est_foi(x, model$sp)
  model$df   <- list(t=x, spos=y)
  class(model) <- "plf"
  model
}

df <- parvob19_be_2001_2003
filter <- (df$age>0.5)&(df$age<76)&(!is.na(df$age))&!is.na(df$seropositive)
sdf <- df[filter,]
sdf <- sdf[order(sdf$age), ]
tps_fit <- tps(
  x=sdf$age,
  y=sdf$parvores,
)
# plot_cr(tps_fit, round.by=0)

