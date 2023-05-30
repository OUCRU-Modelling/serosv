foi.num<-function(x,p) {
  grid<-sort(unique(x))
  pgrid<-(p[order(x)])[duplicated(sort(x))==F]
  dp<-diff(pgrid)/diff(grid)
  foi<-approx((grid[-1]+grid[-length(grid)])/2,dp,grid[c(-1,-length(grid))])$y/(1-pgrid[c(-1,-length(grid))])
  return(list(grid=grid[c(-1,-length(grid))],foi=foi))
}
