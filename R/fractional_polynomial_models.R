# same function but now with powergrid as parameter
# mc=monotonicity constraint, assumes x is ordered
pow<-seq(-2,3,0.1)

is_monotone <- function(model) {
  preds <- predict(model, type = "response")
  all(diff(preds) >= 0) || all(diff(preds) <= 0)
  # (sum(diff(predict(model))<0)==0)
}

# keep
formulate(c(-2,-2, 0, 0, 0))

formulate <- function(powers) {
  equation <- "cbind(pos,tot-pos)~-1"
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

# p_state <- rep(list(pow), 3)
# p_state[1]

seq(5, 1, -1)

best_model <- hb93 %>%
  find_best_frac(p=seq(-2,3,0.05), mc=T, n=2)
# model <- list()

foi.num<-function(x,p)
{
  grid<-sort(unique(x))
  pgrid<-(p[order(x)])[duplicated(sort(x))==F]
  dp<-diff(pgrid)/diff(grid)
  foi<-approx((grid[-1]+grid[-length(grid)])/2,dp,grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])
  return(list(grid=grid[c(-1,-length(grid))],foi=foi))
}

x <- hb93$age
p <- fitted$info$fitted.values
grid <- sort(unique(x))
pgrid<-(p[order(x)])[duplicated(sort(x))==F]
dp<-diff(pgrid)/diff(grid)
foi <- approx(
  (grid[-1]+grid[-length(grid)])/2,
  dp,
  grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])

grid[c(-1, -length(grid))]

help("approx")

hb93$age


model_fp <- fp_model(best_model$power)
fitted <- model_fp$fit(hb93)
model <- glm(cbind(hb93$pos,hb93$tot-hb93$pos) ~ I(hb93$age^1.25) + I(hb93$age^1.3), family=binomial(link="logit"))
foi <- foi.num(hb93$age, model$fitted.values)

fp_model <- function(p) {
  model <- list()
  model$fit <- function(df) {
    with(df, {
      model$info <- glm(
        as.formula(formulate(p)),
        family=binomial(link="logit")
      )
      X <- generate_X(age)
      print(X)
      model$sp <- 1 - model$info$fitted.values
      model$foi <- X%*%model$info$coefficients
      model$df <- df
      model
    })
  }
  generate_X <- function(age) {
    X <- matrix(nrow=length(age), ncol=length(p))
    if (length(p) > 1) {
      for (i in 1:length(p)) {
        X[,i] <- i * age^(p[i]-1)
      }
    }
    -X
  }
  model
}

m <- matrix(nrow = 10, ncol=5)
m[,2] <- 88
m
