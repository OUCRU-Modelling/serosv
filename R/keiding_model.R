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

foi.num<-function(x,p)
{
  grid<-sort(unique(x))
  pgrid<-(p[order(x)])[duplicated(sort(x))==F]
  dp<-diff(pgrid)/diff(grid)
  foi<-approx((grid[-1]+grid[-length(grid)])/2,dp,grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])
  return(list(grid=grid[c(-1,-length(grid))],foi=foi))
}

#' Keiding model
#'
#' @param age the age vector
#' @param pos the positive vector
#' @param tot the total vector
#' @param kernel kernel-based estimate
#' @param bw bandwidth
#'
#'
#' @examples
#' df <- hav_bg_1964
#' model <- keiding_model(df$age,df$pos,df$tot, kernel = "normal", bw = 30)
#' plot(model)

keiding_model <- function(age, pos, tot,kernel ="normal", bw){
  grid <- sort(age)
  model <- list()
  xx <- pavit(pos=pos,tot=tot)
  foi.k1<-foi.num(grid,xx$pai2)$foi
  foi.k1[is.na(foi.k1)]<-0
  foi.k1[foi.k1>10]<-0
  age.k1<-foi.num(grid,xx$pai2)$grid
  fit.k1<- ksmooth(age.k1,foi.k1,kernel=kernel,bandwidth=bw,n.points=length(age.k1))
  model$age <- fit.k1$x
  model$foi <- fit.k1$y
  model$sp <- 1-exp(-cumsum(c(age.k1[1],diff(age.k1))*model$foi))
  model$df <- list(age=age, pos=pos, tot=tot,grid = grid)
  class(model) <- "keiding_model"
  model
  }


#' plot() overloading Keiding model
#'
#' @param x the keiding model object.
#' @param ... arbitrary params
#'
#' 
#' @export
plot.keiding_model <- function(x, ...) {
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
    lines(x$age, x$sp, lwd=2)
    lines(x$age, x$foi, lwd=2, lty=2)
    axis(side=4, at=round(seq(0.0, max(x$foi), length.out=3), 2))
    mtext(side=4, "force of infection", las=3, line=2)
  })
}





