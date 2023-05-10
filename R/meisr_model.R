#' MSEIR model
#'
#' Refers to section 3.4.
#'
#' @param a age sequence
#' @param gamma time in maternal class.
#' @param lambda time in susceptible class.
#' @param sigma time in latent class.
#' @param nu time in infected class.
#'
#' @examples
#' model <- mseir_model(
#'   a=seq(from=1,to=20,length=500), # age range from 0 -> 20 yo
#'   gamma=1/0.5, # 6 months in the maternal antibodies
#'   lambda=0.2,  # 5 years in the susceptible class
#'   sigma=26.07, # 14 days in the latent class
#'   nu=36.5      # 10 days in the infected class
#' )
#' model
#'
#' @export
mseir_model <- function(a, gamma, lambda, sigma, nu)
{
  ma  <- exp(-gamma*a)
  sa  <- (gamma/(gamma-lambda))*(exp(-lambda*a)-exp(-gamma*a))
  ea  <- ((lambda*gamma)/(gamma-lambda))*
    (
      ((exp(-sigma*a)-exp(-lambda*a))/(lambda-sigma))
      -((exp(-sigma*a)-exp(-gamma*a))/(gamma-sigma))
    )
  ia  <- (sigma*lambda*gamma)*
    (
      ((exp(-nu*a)-exp(-sigma*a))/((lambda-sigma)*(gamma-sigma)*(sigma-nu)))
      +((exp(-nu*a)-exp(-lambda*a))/((lambda-gamma)*(lambda-sigma)*(lambda-nu)))
      +((exp(-nu*a)-exp(-gamma*a))/((gamma-lambda)*(gamma-sigma)*(gamma-nu)))
    )
  data.frame(
    a =  c(0, a),
    ma = c(1, ma),
    sa = c(0, sa),
    ea = c(0, ea),
    ia = c(0, ia),
    ra = c(0, ma - sa - ea - ia)
  )
}

# par(mfrow=c(3,2),lwd=2,las=1,cex.axis=1.1,cex.lab=1.1,mgp=c(3, 0.5, 0))
# a <- seq(from=1,to=20,length=500) # age range from 0 -> 20 yo
# model <- mseir_model(
#   a,
#   gamma=1/0.5, # 6 months in the maternal antibodies
#   lambda=0.2,  # 5 years in the susceptible class
#   sigma=26.07, # 14 days in the latent class
#   nu=36.5      # 10 days in the infected class
# )
# model
#
# plot(m$a,m$ma,type="l",xlab="Age",ylab="M(a)",pch=0.5)
# title("a:Proportion of host with maternal antibodies",adj=0,cex=0.35)
# plot(m$a,m$sa,type="l",xlab="Age",ylab="S(a)",pch=0.5)
# title("b:Proportion of susceptibles",adj=0,cex=0.35)
# plot(m$a,m$ea,type="l",xlab="Age",ylab="E(a)",pch=0.5)
# title("c:Proportion of host in the latent class",adj=0,cex=0.35)
# plot(m$a,m$ia,type="l",xlab="Age",ylab="I(a)",pch=0.5)
# title("d:Proportion of infected",adj=0,cex=0.35)
# plot(m$a,m$ra,type="l",xlab="Age",ylab="R(a)",pch=0.5)
# title("e:Proportion host in the immune class",adj=0,cex=0.35)
# plot(m$a,(m$ma+m$ia+m$ra),type="l",xlab="Age",ylab="seroprevalence",pch=0.5)
# title("f:Proportion of sero-positive",adj=0,cex=0.35)
