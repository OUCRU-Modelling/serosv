library(locfit)

pavit<- function(pos=pos,tot=rep(1,length(pos)))
{
  gi<- pos/tot
  pai1 <- pai2 <- gi
  N <- length(pai1)
  ni<-tot
  for(i in 1:(N - 1)) {
    if(pai2[i] > pai2[i + 1]) {
      pool <- (ni[i]*pai1[i] + ni[i+1]*pai1[i + 1])/(ni[i]+ni[i+1])
      pai2[i:(i + 1)] <- pool
      k <- i + 1
      for(j in (k - 1):1) {
        if(pai2[j] > pai2[k]) {
          pool.2 <- sum(ni[j:k]*pai1[j:k])/(sum(ni[j:k]))
          pai2[j:k] <- pool.2
        }
      }
    }
  }
  return(list(pai1=pai1,pai2=pai2))
}

#' Returns the optimal smoothing parameter of GCV  
#'
#' Refers to section 9.3.3.
#' @param age the age vector
#' @param pos the positive vector
#' @param tot the total vector
#' @param alphagrid a alpha sequence
#' @param family family of gcvplot. Default to "binomal"
#' @param deg degree of gcvplot
#'
#' @examples
#' df <- hav_bg_1964
#' find_the_min_alpha(df$age,df$pos,df$tot,
#'            seq(0.2,2, by=0.05),family = "binomial", deg =2)

find_the_min_alpha <- function(age,pos,tot,alphagrid, family = "binomial",deg){
  neg <- tot - pos 
  a<-c(rep(age,pos),rep(age,neg))
  y<-c(rep(rep(1,length(age)),pos),rep(rep(0,length(age)),neg))
  gcvp<-gcvplot(y~a,family="binomial",alpha= alphagrid,deg = deg)
  alpha<-alphagrid[which.min(gcvp$values)]
  alpha
}



#' A Smooth then constrain model
#'
#' @param age the age vector
#' @param pos the positive vector
#' @param tot the total vector
#' @param alpha alpha of local fit
#' @param family family of local fit
#'
#' @examples
#' df <- hav_bg_1964
#' model <- stc(df$age,df$pos,df$tot,
#'             alpha = 0.35,family = "binomial")
stc <- function(age,pos,tot,alpha,family = "binomial"){
  neg <- tot - pos
  grid <- sort(age)
  model <- list()
  a<-c(rep(age,pos),rep(age,neg))
  y<-c(rep(rep(1,length(age)),pos),rep(rep(0,length(age)),neg))
  y<-y[order(a)]
  a<-a[order(a)]
  lpfit1 <- locfit(y~a,family= family ,alpha=alpha)
  lpfitd1 <- locfit(y~a,deriv=1,family= family,alpha=alpha)
  lpfoi1 <- fitted(lpfitd1)*fitted(lpfit1)
  model$sp <- pavit(pos=fitted(lpfit1))$pai2
  lpfoi2 <- apply(cbind(0,fitted(lpfitd1)),1,max)*model$sp
  model$foi <- apply(cbind(0,lpfoi2),1,max)
  model$df <- list(age=age, pos=pos, tot = tot, grid = grid ,a=a )
  class(model) <- "smooth_then_constrain_model"
  model
}

#' plot() overloading smooth then constrain model
#'
#' @param x the smooth then constrain model object
#' @param ... arbitrary params
#'
#'
#' @export
plot.smooth_then_constrain_model <- function(x, ...) {
  CEX_SCALER <- 4 # arbitrary number for better visual
  with(x$df, {
    par(las=1,cex.axis=1,cex.lab=1,lwd=2,mgp=c(2, 0.5, 0),mar=c(4,4,4,3))
    plot(
      grid,
      pos/tot,
      cex=CEX_SCALER*tot/max(tot),
      xlab="age", ylab="seroprevalence",
      xlim=c(0, max(age)), ylim=c(0,1)
    )
    lines(x$df$a,x$sp,lty=1)
    lines(x$df$a,x$foi,lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}

