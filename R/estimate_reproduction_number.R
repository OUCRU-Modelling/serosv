install.packages("car")
library(car)

source("R/estimate_force_of_infection.R")

obj = coef(model)
names(obj)<-"Intercept"
deltaMethod(obj, "1/exp((Intercept))", vcov.=vcov(model))
deltaMethod(obj, "75*exp((Intercept))", vcov.=vcov(model))
