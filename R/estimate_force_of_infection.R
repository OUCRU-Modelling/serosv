
fit_GLM <- function(age, pos, neg)
{
  lAge <- log(age)
  glm(cbind(pos, neg)~1,offset=lAge,family=binomial(link=cloglog))
}

hav<-read.table("/Users/thanhlongb/Projects/data-engineer/serosv/R/RCodeBook/Chapter3/HAV-BUL.dat",header=T)
age <- hav$Age
pos <- hav$Pos
tot <- hav$Tot
neg <- tot - pos
model <- fit_GLM(age, pos, neg)
summary(model)
coef(model)["(Intercept)"] # <- est. force of infection

### FIGURE 3.14
par(mfrow=c(1,1),cex.axis=1.2,cex.lab=1.2,lwd=3,las=1,mgp=c(3, 0.5, 0))
plot(age,pos/tot,cex=0.1*tot,xlab="age",ylab="seroprevalence",xlim=c(0,86),ylim=c(0,1))

beta <- 0.0505
pi <- 1-exp(-beta*age) # π
lines(age, pi)
