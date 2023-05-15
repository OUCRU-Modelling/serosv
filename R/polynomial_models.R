
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

polynomial_model <- function(degree) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        as.formula(linear_predictor_str()),
        family=binomial(link="log")
      )
      X <- generate_X(age)
      model$sp <- 1 - model$info$fitted.values
      model$foi <- X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  linear_predictor_str <- function() {
    formula <- "cbind(tot-pos, pos)~-1"
    for (i in 1:degree) {
      formula <- paste0(formula, "+I(age^", i, ")")
    }
    formula
  }
  generate_X <- function(age) {
    X <- matrix(rep(1, length(age)), ncol = 1)
    if (degree > 1) {
      for (i in 2:degree) {
        X <- cbind(X, i * age^(i-1))
      }
    }
    -X
  }
  model
}

model_poly <- polynomial_model(degree=3)

hav_be_1993_1994 %>%
  model_poly$fit() %>%
  plot_p_foi_wrt_age()

hcv_trans <- transform_hcv(hcv_be_2006)
wei <- model_wei$fit(hcv_trans)
summary(wei$foi)

summary(hcv_be_2006)
summary(wei$df_dur)

muench_model <- function() {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        cbind(tot-pos, pos)~-1+age,
        family=binomial(link="log")
      )
      X <- -matrix(rep(1,length(age)))
      model$sp <- 1 - model$info$fitted.values
      model$foi <- 5*X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  model
}

griffiths_model <- function(df) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        cbind(tot-pos,pos)~-1+age+I(age^2),
        family=binomial(link="log")
      )
      X <- -cbind(rep(1,length(age)),2*age)
      model$sp <- 1 - model$info$fitted.values
      model$foi <- 5*X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  model
}

grenfell_anderson_model <- function(df) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        cbind(tot-pos,pos)~-1+age+I(age^2)+I(age^3),
        family=binomial(link="log")
      )
      X <- -cbind(rep(1,length(age)),2*age,3*age^2)
      model$sp <- 1 - model$info$fitted.values
      model$foi <- 5*X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  model
}

library(stats4)
farrington_model <- function(parameters)
{
  model <- list()
  model$parameters <- parameters
  model$fit <- function(df) {
    with(c(df, parameters), {
      farrington <- function(alpha,beta,gamma) {
        p=1-exp((alpha/beta)*age*exp(-beta*age)
                +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)-gamma*age)
        ll=pos*log(p)+(tot-pos)*log(1-p)
        return(-sum(ll))
      }

      model$info <- mle(farrington, start=as.list(parameters))
      alpha <- coef(model$info)[1]
      beta  <- coef(model$info)[2]
      gamma <- coef(model$info)[3]
      model$sp <- 1-exp(
        (alpha/beta)*age*exp(-beta*age)
        +(1/beta)*((alpha/beta)-gamma)*(exp(-beta*age)-1)
        -gamma*age)
      model$foi <- (alpha*age-gamma)*exp(-beta*age)+gamma
      model$df <- df
      model
    })
  }
  model
}

parms <- c(alpha=1, beta=2, gamma=3)
parms_list <- as.list(parms)
parms_list$alpha

parms_list_2 <- list(alpha=1, beta=2, gamma=3)
parms_list_2


head(rubella)
a <- rubella %>%
  rename_col("age", "Age") %>%
  rename_col("pos", "Pos")
head(a)


model_far <- farrington_model(c(alpha=0.07,beta=0.1,gamma=0.03))
model_mue <- muench_model()
model_gre <- grenfell_anderson_model()
model_gri <- griffiths_model()

hav %>%
  model_wei$fit() %>%
  plot_p_foi_wrt_age()





#-----------------------------------


windows(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(2, 0.5, 0),mar=c(4.1,4.1,4.1,3))

names(rubella)

Age <- rubella$age
Pos <- rubella$pos
Tot <- rubella$tot

plot(Age,Pos/Tot,cex=0.018*Tot,xlab="age",xlim=c(0,45),ylim=c(0,1),ylab="seroprevalence")
alpha<-coef(model5)[1]
beta<-coef(model5)[2]
gamma<-coef(model5)[3]

p<-1-exp((alpha/beta)*Age*exp(-beta*Age)+(1/beta)*((alpha/beta)-gamma)*(exp(-beta*Age)-1)-gamma*Age)
lines(Age,p,lwd=2)
foi<-(alpha*Age-gamma)*exp(-beta*Age)+gamma
lines(Age,2*foi,lwd=2)

alpha<-coef(model5)[1]
beta<-coef(model5)[2]
gamma<-0

p<-1-exp((alpha/beta)*Age*exp(-beta*Age)+(1/beta)*((alpha/beta)-gamma)*(exp(-beta*Age)-1)-gamma*Age)
lines(Age,p,lwd=2,lty=2)
foi<-(alpha*Age-gamma)*exp(-beta*Age)+gamma
lines(Age,2*foi,lwd=2,lty=2)

axis(side=4,at=c(0.0,0.2,0.4),labels=c(0.00,0.1,0.2))
mtext(side=4,"force of infection", las=3,line=1.5)








model5=mle(farrington,start=list(alpha=0.07,beta=0.1,gamma=0.03))
summary(model5)
AIC(model5)

hav<-read.table("/Users/thanhlongb/Projects/oucru/the-book/RCodeBook/Chapter4/HAV-BUL.dat",header=T)
names(hav)
Age <- hav$Age
Pos <- hav$Pos
Tot <- hav$Tot

attach(hav)

windows(record=TRUE, width=5, height=5)
par(las=1,cex.axis=1.1,cex.lab=1.1,lwd=3,mgp=c(2, 0.5, 0),mar=c(4.1,4.1,4.1,3))

names(hav)[names(hav) == "Age"] <- "age"
names(hav)[names(hav) == "Pos"] <- "pos"
names(hav)[names(hav) == "Tot"] <- "tot"


head(hav)
model1 <- muench_model(hav)
model1$plot()
model3 <- griffiths_model(hav)
model4 <- grenfell_anderson_model(hav)

plot(Age,Pos/Tot,cex=0.1*Tot,xlab="age",xlim=c(0,86),ylim=c(0,1),ylab="seroprevalence")
lines(Age, model1$prevalence, lwd=2)
lines(Age, model1$foi, lwd=2)

lines(Age,1-model1$fitted.values,lwd=2)
lines(Age,1-model3$fitted.values,lwd=2,lty=2)
lines(Age,1-model4$fitted.values,lwd=2,lty=3)

X<--matrix(rep(1,length(Age)))
lines(Age,5*X%*%model1$coefficients,lwd=2)
X<--cbind(rep(1,length(Age)),2*Age)
lines(Age,5*X%*%model3$coefficients,lwd=2,lty=2)

X<--cbind(rep(1,length(Age)),2*Age,3*Age^2)
lines(Age,5*X%*%model4$coefficients,lwd=2,lty=3)
axis(side=4,at=c(0.0,0.2,0.4),labels=c(0.00,0.04,0.08))
mtext(side=4,"force of infection", las=3,line=2)

