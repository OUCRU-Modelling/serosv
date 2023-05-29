
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
use_data(vzv_be_2001_2003, overwrite = TRUE)

head(vzv_be_2001_2003)

sp_model <- function(k, deg) {
  model <- list()
  model$parameters <- list(k=k, deg=deg)
  model$fit <- function(df) {
    with(c(df, model$parameters), {
      n <- length(age)
      # params naming must match for `spm()` to work
      y <- seropositive
      a <- age
      knr <- k
      # -----====
      model$info <- spm(
        y~f(a, degree=deg, basis="trunc.poly", knots=default.knots(a, knr)),
        family="binomial"
      )
      u_k <- outer(age, default.knots(age, k), basis_fn)*t(
        matrix(
          rep(model$info$fit$coefficients$random$dummy.group.vec.Handan, n),
          ncol=n
        )
      )
      model$semi_param_part <- apply(u_k, 1, sum)
      model$param_part <- (cbind(1, age, age^2, age^3)[, 1:(1+deg)]) %*% model$info$fit$coefficients$fixed
      model$sp <- expit(model$semi_param_part + model$semi_param_part)
      model$df <- df
      model
    })
  }
  model
}

model <- sp_model(k=3, deg=2)
model_fitted <- model$fit(vzv_be_2001_2003)


plot(a, sp, ylim=c(0,1), xlab="age", type="l", lwd=2, ylab=expression(pi(age)))
lines(a, expit(semi_param_part), lty=3, lwd=2)
lines(a, expit(param_part), lty=2, lwd=2)

foi.num<-function(x,p)
{
  grid<-sort(unique(x))
  pgrid<-(p[order(x)])[duplicated(sort(x))==F]
  dp<-diff(pgrid)/diff(grid)
  foi<-approx((grid[-1]+grid[-length(grid)])/2,dp,grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])
  return(list(grid=grid[c(-1,-length(grid))],foi=foi))
}

model_fitted$foi <- foi.num(vzv_be_2001_2003$age, model_fitted$sp)

length(expit(model_fitted$semi_param_part + model_fitted$param_part))
length(model_fitted$foi)

with(model_fitted, {
  par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
  plot(df$age, sp, ylim=c(0,1), xlab="age", type="l", lwd=2, ylab="sp")
  lines(foi$grid, foi$foi, lwd=2, lty=2)
  axis(side=4, at=round(seq(0.0, max(model$foi), length.out=3), 2))
  mtext(side=4, "force of infection", las=3, line=2)
})
