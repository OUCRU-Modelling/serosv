#' A bivariate Dale model
#'
#' Refers to section 12.2
#'
#' @param age the age vector
#' @param y a dataframe of predictors (NN, NP, PN, PP), does not have to be in order
#' @param constant_or whether to set constant OR restriction or not
#' @param monotonized whether to monotonize seroprevalence or not
#' @param sparopt smoothing parameter for VGAM model
#'
#' @importFrom VGAM s vgam predictors binom2.or fitted deviance
#' @importFrom locfit expit
#'
#' @examples
#' data <- rubella_mumps_uk
#' model <- bivariate_dale_model(age = data$age, y = data[, c("NN", "NP", "PN", "PP")], monotonized=T)
#' plot(model, y1 = "Rubella", y2 = "Mumps", plot_type = "sp")
#'
#' @export

bivariate_dale_model <- function (age, y, constant_or = F, monotonized=T, sparopt=0){
  if(constant_or){
    constr <-  list(
      "(Intercept)" = diag(3), # intercept for all components
      "s(age, spar = sparopt)" = diag(3)[,1:2] # smoothing spline for fist 2 components
    )
  }else{
    constr <- NULL
  }

  model <- list()
  model$info <-  vgam(cbind(y$NN, y$NP, y$PN, y$PP) ~ s(age, spar=sparopt),
                      family=binom2.or(zero = NULL),
                      constraints = constr)

  model$df <- list(
    age = age,
    y = y
  )

  # --- return prevalence and conditional for marginal & conditional
  model$sp <- list(
    "y1m" = rowSums(fitted(model$info)[,c("10", "11")]), # pi1+ i.e. P(Y1 = positive)
    "y2m"= rowSums(fitted(model$info)[,c("01", "11")]), # pi+1 i.e. P(Y2 = positive)
    # P(Y1 = positive | Y2 = positive) = P(Y1 = positive | Y2 = positive)/ P(Y2 = positive)
    "y1condy2pos" = fitted(model$info)[,"11"]/rowSums(fitted(model$info)[,c("01", "11")]),
    # P(Y1 = positive | Y2 = negative) = P(Y1 = positive | Y2 = negative)/ P(Y2 = negative)
    "y1condy2neg" = fitted(model$info)[,"10"]/rowSums(fitted(model$info)[,c("10", "00")]),
    # P(Y2 = positive | Y1 = positive) = P(Y2 = positive & Y1 = positive) / P(Y1 = positive)
    "y2condy1pos" = fitted(model$info)[,"11"]/rowSums(fitted(model$info)[,c("10", "11")]),
    # P(Y2 = positive | Y1 = negative) = P(Y2 = positive & Y1 = negative) / P(Y1 = negative)
    "y2condy1neg" = fitted(model$info)[,"01"]/rowSums(fitted(model$info)[,c("00", "01")]),
    "joint" = fitted(model$info)[,"11"]
  )

  # ---- Monotonize sp
  if (monotonized){
    for (key in names(model$sp)){
      model$sp[[key]] <- pava(model$sp[[key]])$pai2
    }
  }

  # --- Compute foi
  for (key in names(model$sp)){
    model$foi[[key]] <- est_foi(age, model$sp[[key]])
  }

  # ---- Compute CI here
  suppressWarnings(
    {
      result <- bdm_bootstrapping(age, y, runs = 1000, constraint = constant_or)
    }
  )
  model$ci <- list(
    # each val is a data frame with 3 rows: val, lower bound, upper bound
    "y1m" = rbind(model$sp$y1m, expit(apply(result$y1m, 2, ci))),
    "y2m"= rbind(model$sp$y2m, expit(apply(result$y2m, 2, ci))),
    "or" = rbind(exp(predictors(model$info)[,"loglink(oratio)"]), exp(apply(result$or, 2, ci))),
    "pi11" = rbind(fitted(model$info)[,"11"], apply(result$pi11, 2, ci))
  )

  class(model) <- "bivariate_dale_model"

  model
}


#' Find best smoothing parameter for Bivariate Dale model
#'
#' @param age the age vector
#' @param y a dataframe of predictors (NN, NP, PN, PP), does not have to be in order
#' @param spar_seq sequence of smoothing parameters to be evaluated
#'
#' @importFrom VGAM s vgam deviance df.residual AIC binom2.or fitted
#'
#' @return data frame of optimal parameter for different criteria
#'
#' @export
find_best_bdm_param <- function(age, y, spar_seq = seq(0,0.1,0.001)){
  y<-cbind(y$NN, y$NP, y$PN, y$PP)
  a<-as.numeric(age)

  BICf<-function(fit){
    return(
      VGAM::deviance(fit)+log( dim(y)[1])*(3*dim(y)[1]-df.residual(fit))
    )
  }
  # TODO: compute UBRE for vgam
  # UBRE<-function(fit){
  #   return(
  #     VGAM::deviance(fit)+log(dim(y)[1])*(3*dim(y)[1]-df.residual(fit))
  #     )
  #   }
  out<-matrix(NA,ncol=3,nrow = length(spar_seq))


  for (i in 1:length(spar_seq)){
    skip <- FALSE
    suppressWarnings(
      {
        tryCatch({
          fit <- vgam(y~s(a,spar=spar_seq[i]),family=binom2.or(zero = NULL))
        }, error = function(e){
          skip <- TRUE
        })
      }
    )

    if(skip){
      # skip sparopt that returns error
      out[i,] <- c(spar_seq[i],0, 0)
      next
    }
    out[i,] <- c(spar_seq[i],
                 BICf(fit),
                 AIC(fit) # use built in funcition for AIC
                 # UBRE(fit)
    )
  }

  df.BIC <- out[which.min(out[,2]),1]
  df.AIC <- out[which.min(out[,3]),1]
  # df.UBRE <- out[which.min(out[,4]),1]
  val.BIC <- out[which.min(out[,2]),2]
  val.AIC <- out[which.min(out[,3]),3]
  # val.UBRE <- out[which.min(out[,4]),4]
  opt_df <- data.frame(
    "Methods" = c("BIC","AIC"),
    "Value"   = c(val.BIC,val.AIC),
    "sparopt" = c(df.BIC,df.AIC)
  )
  opt_df
}

#' Generate the (NN, NP, PN, PP) matrix from line-listing data
#'
#' @param data - dataframe being transformed
#' @param y1 - serostatus for y1
#' @param y2 - serostatus for y2
#' @param age - individual's age
#' @param discrete_age - round age so it is discrete
#'
#' @import dplyr magrittr
#' @return transformed matrix
#'
#' @examples
#' data <- vzv_parvo_be
#' data <- generate_quad_matrix(data, parvo_res, vzv_res, age)
#' @export
generate_quad_matrix <- function (data, y1, y2, age, discrete_age = F){
  # handle continuous age
  if(discrete_age) {data$age <- round(data$age)}

  data %>%
    dplyr::filter(!is.na(age), !is.na({{y1}}), !is.na({{y2}}) ) %>%
    group_by({{age}}) %>%
    mutate(
      age = as.numeric(age)
    ) %>%
    summarize(
      NN = sum({{y1}} == 0 & {{y2}} == 0),
      NP = sum({{y1}} == 0 & {{y2}} == 1),
      PN = sum({{y1}} == 1 & {{y2}} == 0),
      PP = sum({{y1}} == 1 & {{y2}} == 1),
      total = n()
    )
}

# --- Compute confidence interval
ci<-function(x){
  quantile(x,prob=c(0.025,0.975),na.rm=T)
}

# --- Perform bootstrapping to calculate CI
bdm_bootstrapping <- function(x, y, runs=1000, sparopt = 0, constraint=T){
  predictors <- c("y1m", "y2m", "or") # pi1 is pi1+ while pi2 is pi+1 in book
  quads <- c("pi00", "pi01", "pi10", "pi11")
  conditional_prob <- c("y1condy2", "y2condy1")

  output <- list()
  # --- Instantiate output
  #
  for (item in append(append(predictors, quads), conditional_prob)){
    output[[item]] <- matrix(nrow = runs, ncol = length(x))
  }

  for (i in 1:runs){
    # print(c("bootstrap ", i))
    sample_indices <- sample(c(1:length(x)), length(x), replace = T)
    sample_x <- x[sample_indices]
    sample_y <- y[sample_indices,]

    if(constraint){
      constr <- list(
        "(Intercept)" = diag(3), # intercept for all components
        "s(data_x, spar = sparopt)" = diag(3)[,1:2] # smoothing spline for fist 2 components
      )
    }else{
      constr <- NULL
    }

    bdm_fit <- vgam(
      cbind(sample_y$NN, sample_y$NP, sample_y$PN, sample_y$PP) ~ s(sample_x, spar=sparopt),
      family=binom2.or(zero = NULL),
      constraints = constr
    )
    # ---- calculate predictors
    for (pred in 1:length(predictors)){
      output[[ predictors[pred] ]][i, ] <- approx(x=sample_x,
                                                  y=predictors(bdm_fit)[,pred], xout=x)$y
    }
    # ---- calculate prob matrix
    for (quad in 1:length(quads)){
      output[[ quads[quad] ]][i,] <- approx(x=sample_x,y=fitted(bdm_fit)[,quad],xout=x)$y
    }
    # ---- calculate cond matrix
    output$y1condy2[i,] <- approx(x = sample_x,
                                  y=fitted(bdm_fit)[,"11"]/rowSums(fitted(bdm_fit)[,c("01","11")]), xout=x)$y
    output$y2condy1[i,] <- approx(x = sample_x,
                                  y=fitted(bdm_fit)[,"11"]/rowSums(fitted(bdm_fit)[,c("10","11")]), xout=x)$y
  }
  output
}
