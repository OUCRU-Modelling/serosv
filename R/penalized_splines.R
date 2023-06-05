
# Penalized Likelihood Framework

parvovirus<-read.table("/Users/thanhlongb/Projects/oucru/the-book/RCodeBook/Chapter4/VZV-B19-BE.dat",header=T)
subset <- (parvovirus$age>0.5)&(parvovirus$age<40)&(!is.na(parvovirus$age))&!is.na(parvovirus$parvores)
parvovirus<-parvovirus[subset,]
y<-parvovirus$parvores[order(parvovirus$age)]
a<-parvovirus$age[order(parvovirus$age)]


knr <- 3
deg <- 2

# basis function
basis_fn <- function(x, y)
{
  (x-y)^deg*(x>y)
}

# not sure what this is
expit <- function(x)
{
  return(exp(x)/(1+exp(x)))
}

# estimator <- f(a,
#   degree=deg,
#   basis="trunc.poly",
#   knots=default.knots(a, knr)
# )

### library Semipar: Ruppert et al. (2003)
library(SemiPar)

vzv_be_2001_2003 <- read.table(
  "/Users/thanhlongb/Projects/oucru/the-book/RCodeBook/Chapter4/VZV-B19-BE.dat",
  header = TRUE)
vzv_be_2001_2003 <- vzv_be_2001_2003 %>%
  mutate(
    age = age, seropositive = VZVres,
    gender = sex, VZVmIUml = VZVmUIml, .keep = "none"
  ) %>%
  filter(!is.na(seropositive)) %>%
  relocate(gender, .after = seropositive)
# use_data(vzv_be_2001_2003, overwrite = TRUE)

vzv_be_2001_2003 <- vzv_be_2001_2003[order(vzv_be_2001_2003$age), ]

vzv_be_2001_2003

sp_model <- function(k, deg) {
  model <- list()
  model$parameters <- list(k=k, deg=deg)
  model$fit <- function(df) {
    y = df$seropositive
    with(c(df, model$parameters), {
      n <- length(age)
      # Workaround as `spm()` get values from environment rather than passed in
      assign("y_", seropositive, envir = globalenv())
      assign("a_", age, envir = globalenv())
      assign("k_", k, envir = globalenv())
      assign("d_", deg, envir = globalenv())
      # ----
      model$info <- spm(
        y_~f(a_, degree=d_, basis="trunc.poly", knots=default.knots(a_, k_)),
        family="binomial"
      )
      # clean up
      for (i in c("y_", "a_", "k_", "d_")) rm(i, envir = globalenv())
      # ----
      u_k <- outer(age, default.knots(age, k), basis_fn)*t(
        matrix(
          rep(model$info$fit$coefficients$random$dummy.group.vec.Handan, n),
          ncol=n
        )
      )
      model$semi_param_part <- apply(u_k, 1, sum)
      model$param_part <- (cbind(1, age, age^2, age^3)[, 1:(1+deg)]) %*% model$info$fit$coefficients$fixed
      model$sp <- expit(model$param_part + model$semi_param_part)
      model$df <- df
      model
    })
  }
  model
}

model <- sp_model(k=100, deg=2)
model_fitted <- model$fit(vzv_be_2001_2003)
model_fitted$foi <- foi.num(vzv_be_2001_2003$age, model_fitted$sp)

with(model_fitted, {
  par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
  plot(df$age, sp, ylim=c(0,1), xlab="age", type="l", lwd=2, ylab="sp")
  lines(foi$grid, foi$foi, lwd=2, lty=2)
  axis(side=4, at=round(seq(0.0, max(foi$foi), length.out=3), 2))
  mtext(side=4, "force of infection", las=3, line=2)
})

best_smooth_spline <- function(t, spos, d_seq, link) {
  best <- cbind()
  for (d in d_seq) {
    fit <- gam(spos~s(t, df=d),family=binomial(link=link))
    best <- cbind(best, list(
      BIC=unname(fit$deviance+log(length(fit$y))*(fit$nl.df+2)),
      fit=fit,
      d=d
    ))
  }

  best[, which.min(best[1,])]
}

library(gam)
best_ss <- best_smooth_spline(
  t=vzv_be_2001_2003$age,
  spos=vzv_be_2001_2003$seropositive,
  d_seq=seq(1, 10, by=0.5),
  link="cloglog"
  )
best_ss

detach(package:gam)
library(mgcv)
best_crs <- function(t, y, d_seq, link="logit") {
  best <- cbind()

  for (d in d_seq) {
    fit <- gam(y~s(t, bs="cr"),family=binomial(link=link))
    best <- cbind(best, list(
      BIC=unname(fit$deviance+log(length(fit$y))*sum(fit$edf)),
      fit=fit,
      d=d
    ))
  }

  # best[, which.min(best[1,])]
  best
}

df <- vzv_be_2001_2003
model <- gam(df$seropositive~s(df$age, bs="tp"),family=binomial(link="logit"))
model$edf

crs <- best_crs(
  t=vzv_be_2001_2003$age,
  y=vzv_be_2001_2003$seropositive,
  d_seq=seq(1, 10, by=0.5),
  link="cloglog"
)
crs$BIC

# library(gam)
best_ss <- best_smooth_spline(
  t=vzv_be_2001_2003$age,
  spos=vzv_be_2001_2003$seropositive,
  d_seq=seq(1, 10, by=0.5),
  link="cloglog"
)
best_ss


# library(splines)
# install.packages("pspline")
library(pspline)
library(face)
best_b_spline <- function(t, spos, lambda_seq, deg=3, m=2, k=20, link="logit") {
  best <- cbind()

  for (lambda in lambda_seq) {
    fit <- pspline.fit(
      response=spos, lambda=80,
      x.var=t, x.predicted=t, ps.intervals=k, degree=deg,
      order=m, link=link, family="binomial"
    )
    best <- cbind(best, list(
      BIC=unname(fit$dev+log((dim(fit$summary.predicted)[1]))*(fit$eff.df)),
      fit=fit,
      lambda=lambda
    ))
  }

  best[, which.min(best[1,])]
}

# df <- parvob19_be_2001_2003
# best_bs <- best_b_spline(
#   t=df$age,
#   spos=df$seropositive,
#   lambda_seq=seq(50, 200, by=10)
# )

