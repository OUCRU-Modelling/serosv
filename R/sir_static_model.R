# install.packages("deSolve")
require(deSolve)

# states:
#		s: % susceptible
#		i: % infected
#		r: % recovered
# parameters:
#   lambda: natural death rate
#		nu: recovery rate
SIR_static <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)),
     {
       ds <- -lambda*s
       di <- lambda*s - nu*i
       dr <- nu*i
       list(c(ds, di, dr))
     }
  )
}

get_SIR_static <- function(ages, state, parameters)
{
  as.data.frame(ode(y=state,times=ages,func=SIR_static,parms=parameters))
}

plot_SIR_static <- function(df)
{
  par(mfrow=c(1,3),cex.axis=1.2,cex.lab=1.2,lwd=3,las=1,mgp=c(3, 0.5, 0))

  plot(
    df$time,
    df$s,
    xlab="age", ylab="% susceptible", type="l", lwd=2
  )
  plot(
    df$time,
    df$i,
    xlab="age", ylab="% infected", type="l", lwd=2
  )
  plot(
    df$time,
    df$r,
    xlab="age", ylab="% recovered (seroprevalence)", type="l", lwd=2
  )
}

# state <- c(s=0.99,i=0.01,r=0)
# parameters <- c(
#   lambda = 0.05,
#   nu=1/(14/365) # 2 weeks to recover
# )
# ages<-seq(0, 90, by=0.01)
#
# df = get_SIR_static(ages, state, parameters)
# plot_SIR_static(df)

