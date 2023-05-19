is_monotone <- function(model) {
  preds <- predict(model, type = "response")
  all(diff(preds) >= 0) || all(diff(preds) <= 0)
  # alternative: (sum(diff(predict(model))<0)==0)
}

formulate <- function(powers) {
  equation <- "cbind(pos,tot-pos)~"
  prev_term <- ""

  for (i in 1:length(powers)) {
    if (i > 1 && powers[i] == powers[i-1]) {
      cur_term <- paste0("I(", prev_term, "*log(age))")
    } else if (powers[i] == 0) {
      cur_term <- "log(age)"
    } else {
      cur_term <- paste0("I(age^", powers[i], ")")
    }
    equation <- paste0(equation, "+", cur_term)
    prev_term <- cur_term
  }
  equation
}

find_best_frac <- function(df, p, mc, n){
  pos <- df$pos
  tot <- df$tot
  age <- df$age
  #----
  glm_best <- NULL
  d_best <- NULL
  p_best <- NULL
  #----
  min_p <- 1
  max_p <- length(p)
  state <- rep(min_p, n)
  i <- n
  #----

  get_cur_p <- function(cur_state) {
    cur_p <- c()
    for (i in 1:n) {
      cur_p <- c(cur_p, p[cur_state[i]])
    }
    cur_p
  }

  repeat {
    if (
      (i < n && state[i] == max_p)
      || (i == n && state[i] == max_p+1)
    ) {
      if (i-1 == 0) break
      if (state[i-1] < max_p) {
        state[i-1] <- state[i-1]+1
        for (j in i:n) state[j] <- state[i-1]
        i <- n
      } else {
        i <- i-1
        next
      }
    }
    #------ iteration implementation -------
    p_cur <- get_cur_p(state)

    glm_cur <- glm(
      as.formula(formulate(p_cur)),
      family=binomial(link="logit")
    )
    if (glm_cur$converged == TRUE) {
      d_cur <- deviance(glm_cur)
      if (is.null(glm_best) || d_cur < d_best) {
        if ((mc && is_monotone(glm_cur)) | !mc) {
          glm_best <- glm_cur
          d_best <- d_cur
          p_best <- p_cur
        }
      }
    }
    #---------------------------------------
    if (sum(state != max_p) == 0) break
    state[i] <- state[i]+1
  }
  return(list(power=p_best, deviance=d_best, model=glm_best))
}

# assuming age is sorted!!
estimate_foi <- function(age,p)
{
  dp<-diff(p)/diff(age)
  foi <- approx(
    (age[-1]+age[-length(age)])/2,
    dp,
    age[c(-1,-length(age))]
  )$y/(1-p[c(-1,-length(age))])
  return(list(age=age[c(-1,-length(age))],foi=foi))
}

fp_model <- function(p) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        as.formula(formulate(p)),
        family=binomial(link="logit")
      )
      X <- generate_X(age)
      model$sp  <- model$info$fitted.values
      model$foi <- est_foi(age, model$info$fitted.values)$foi
      model$df  <- df
      model
    })
  }
  model
}

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

hbe_best <- hav_be_1993_1994 %>%
  find_best_frac(p=seq(-2,3,0.1), mc=T, n=2)
hbe_best
model_hbe <- fp_model(hbe_best$power)
model_hbe_fitted <- model_hbe$fit(hav_be_1993_1994)
model_hbe_est <- estimate_foi(hav_be_1993_1994$age, model_hbe_fitted$info$fitted.values)

plot_p_foi_wrt_age(model_hbe_fitted)
lines(model_hbe_est$age, model_hbe_est$foi, lwd=2, lty=2)
axis(side=4, at=round(seq(0.0, max(model_hbe_est$foi), length.out=3), 2))
mtext(side=4, "force of infection", las=3, line=2)



