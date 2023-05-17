# same function but now with powergrid as parameter
# mc=monotonicity constraint, assumes x is ordered
pow<-seq(-2,3,0.1)

is_monotone <- function(model) {
  preds <- predict(model, type = "response")
  all(diff(preds) >= 0) || all(diff(preds) <= 0)
  # (sum(diff(predict(model))<0)==0)
}

search.fracpoly.twoR<-function(
  pos, tot, age, pow, mc
){
  glm_best <- NULL
  d_best <- NULL
  mistake <- NULL
  p_best <- NULL

  for (i in 1:(length(pow))){
    for (j in i:(length(pow))){
      if (pow[i]==0) {
        term1 <- log(age)
      } else {
        term1 <- age^pow[i]
      }
      if (pow[j]==pow[i]) {
        term2 <- term1*log(age)
      } else if (pow[j]==0) {
        term2 <- log(age)
      } else {
        term2 <- age^pow[j]
      }
      glm_ij <- glm(
        cbind(pos,tot-pos)~term1+term2,
        family=binomial(link="logit")
      )

      if (glm_ij$converged == FALSE) {
        mistake<-rbind(mistake, c(1,pow[i],pow[j]))
        next
      }

      # if () {
      #   glm_best <- glm_ij
      #   D <- deviance(glm_best)
      #   p_best <- c(pow[i],pow[j])
      #   next
      # }

      d_ij <- deviance(glm_ij)
      if (is.null(glm_best) || d_ij < d_best) {
        if ((mc && is_monotone(glm_ij)) | !mc) {
          glm_best <- glm_ij
          d_best <- d_ij
          p_best <- c(pow[i],pow[j])
        }
      }
    }
  }
  return(list(power=p_best, deviance=d_best, model=glm_best, mistake=mistake))
}

pow<-seq(-2,3,0.1)
# pow
hb93 <- hav_be_1993_1994
found <- search.fracpoly.twoR(
  pos=hb93$pos,
  tot=hb93$tot,
  age=hb93$age,
  pow=pow,
  mc=T
  )
found
is_monotone(found)

library(stats)
glm_fit <- glm(
  cbind(hb93$pos,hb93$tot-hb93$pos) ~ I(hb93$age^1)+I(hb93$age^1.6),
  family=binomial(link="logit")
)
glm_fit
deviance(glm_fit)
glm_fit$coefficients
is_monotone(glm_fit)

pred_glm_fit <- predict(glm_fit)
diff_pred_glm_fit <- diff(pred_glm_fit)
diff_pred_glm_fit < 0


foi.num<-function(x,p)
{
  grid<-sort(unique(x))
  pgrid<-(p[order(x)])[duplicated(sort(x))==F]
  dp<-diff(pgrid)/diff(grid)
  foi<-approx((grid[-1]+grid[-length(grid)])/2,dp,grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])
  return(list(grid=grid[c(-1,-length(grid))],foi=foi))
}

