"pspline.fit"<-
  function(response, x.var, ps.intervals = 8, wts = NULL, degree = 3, order = 3,
           link = "default", family = "gaussian", m.binomial = NULL, r.gamma =
             NULL, lambda = 0, x.predicted = NULL, ridge.adj = 0.0001)
  {
    # Function pspline.fit: univariate smoother using P-splines.
    # Input: x.var= explanatory variable on abcissae.
    # Input: response= response variable.
    # Input: family=gaussian, binomial, poisson, Gamma distribution.
    # Input: wts= vector of weights; default is vector of ones.
    # Input: m.binomial=vector of binomial trials. Default is 1 vector.
    # Input: r.gamma=vector of gamma shape parameters. Default is 1 vector.
    # Input: link= link function (identity, log, sqrt, logit, probit, cloglog, loglog, recipical).
    # Input: ps.intervals= number of intervals for B-splines. Default=8.
    # Input: degree= degree of B-splines. Default=3.
    # Input: order= order of difference penalty. Default=3.
    # Input: lambda= smoothness regulalizing parameter ( >= 0). Default=0.
    # Input: x.predicted=a list of x variables for prediction and twice stderr limits.
    # Result: a scatterplot of (response, x.var) with smoothed fit and se bands.
    # Output: A list: including AIC= deviance + 2*trace(Hat), dispers.parm, etc.
    #
    # Reference: Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and
    #            penalties (with comments and rejoinder). Statistical Science, 11(2): 89-121.
    #
    #
    # Support functions: pspline.checker(), pspline.fitter(), pspline.predictor()
    #
    #
    # (c) 1995 Paul Eilers & Brian Marx
    #
    y <- response
    x <- x.var
    if(missing(wts)) {
      wts <- rep(1, length(y))
    }
    parms <- pspline.checker(family, link, degree, order, ps.intervals,
                             lambda, ridge.adj, wts)
    family <- parms$family
    link <- parms$link
    q <- parms$degree
    d <- parms$order
    ridge.adj <- parms$ridge.adj
    lambda <- parms$lambda
    ndx <- parms$ps.intervals
    wts <- parms$wts
    if(missing(m.binomial)) {
      m.binomial <- rep(1, length(y))
    }
    if(missing(r.gamma)) {
      r.gamma <- rep(1, length(y))
    }
    n <- length(y)
    xl <- min(x)
    xr <- max(x)
    xmax <- xr + 0.01 * (xr - xl)
    xmin <- xl - 0.01 * (xr - xl)
    dx <- (xmax - xmin)/ndx
    knots <- seq(xmin - q * dx, xmax + q * dx, by = dx)
    b <- spline.des(knots, x, q + 1, 0 * x)$design
    n.col <- ncol(b)
    if(d < 0) {
      d <- min(3, (n.col - 1))
      warning(paste("penalty order cannot be negative: have used", d)
      )
    }
    if((d - n.col + 1) > 0) {
      d <- n.col - 1
      warning(paste("penalty order was too large: have used", d))
    }
    if(ridge.adj > 0) {
      nix.ridge <- rep(0, n.col)
      p.ridge <- sqrt(ridge.adj) * diag(rep(1, n.col))
    }
    p <- diag(n.col)
    if(d != 0) {
      for(j in 1:d) {
        p <- diff(p)
      }
    }
    p <- sqrt(lambda) * p
    nix <- rep(0, n.col - d)
    b <- as.matrix(b)
    ps.fit <- pspline.fitter(family, link, n.col, m.binomial, r.gamma, y, b,
                             p, p.ridge, nix, nix.ridge, ridge.adj, wts)
    mu <- ps.fit$mu
    coef <- ps.fit$coef
    w <- ps.fit$w
    e <- 1e-009
    h <- hat(ps.fit$f$qr, intercept = F)[1:n]
    trace <- sum(h) - 1
    if(family == "binomial") {
      dev <- 2 * sum((y + e) * log((y + e)/(mu + e)) + (m.binomial -
                                                          y + e) * log((m.binomial - y + e)/(m.binomial - mu + e)
                                                          ))
      dispersion.parm <- 1
    }
    if(family == "poisson") {
      dev <- 2 * sum(y * log(y + e) - y - y * log(mu) + mu)
      dispersion.parm <- 1
    }
    if(family == "Gamma") {
      dev <- -2 * sum(r.gamma * (log((y + e)/mu) - ((y - mu)/mu)))
      ave.dev <- dev/n
      dispersion.parm <- (ave.dev * (6 + ave.dev))/(6 + 2 * ave.dev)
    }
    if(family == "gaussian") {
      dev <- sum(ps.fit$f$residuals^2)
      dispersion.parm <- dev/(n - trace)
    }
    aic <- dev + 2 * trace
    x.seq <- seq(xl, xr, length = 50)
    b.seq <- spline.des(knots, x.seq, q + 1, 0 * x.seq)$design
    w.aug <- c(w, (nix + 1))
    yhat <- b.seq %*% as.vector(ps.fit$coef)
    half.meat <- sqrt(c(w)) * b
    meat <- t(half.meat) %*% half.meat
    if(ridge.adj > 0) {
      bread <- solve(meat + t(p) %*% p + t(p.ridge) %*% p.ridge)
    }
    if(ridge.adj == 0) {
      bread <- solve(meat + t(p) %*% p)
    }
    half.sw <- half.meat %*% bread
    var.beta <- t(half.sw) %*% half.sw
    var.yhat <- b.seq %*% var.beta %*% t(b.seq)
    stdev.yhat <- as.vector(sqrt(diag(var.yhat)))
    stdev.yhat <- sqrt(dispersion.parm) * stdev.yhat
    pivot <- 2 * stdev.yhat
    upper <- yhat + pivot
    lower <- yhat - pivot
    summary.yhat <- cbind(lower, yhat, upper)
    if(link == "logit") {
      summary.yhat <- 1/(1 + exp( - summary.yhat))
    }
    if(link == "probit") {
      summary.yhat <- apply(summary.yhat, c(1, 2), pnorm)
    }
    if(link == "cloglog") {
      summary.yhat <- (1 - exp( - exp(summary.yhat)))
    }
    if(link == "loglog") {
      summary.yhat <- exp( - exp( - summary.yhat))
    }
    if(link == "sqrt") {
      summary.yhat <- summary.yhat^2
    }
    if(link == "log") {
      summary.yhat <- exp(summary.yhat)
    }
    if(link == "recipical") {
      summary.yhat <- 1/(summary.yhat)
    }
    # if(family == "binomial" && mean(m.binomial) != 1) {
    #   matplot(x.seq, summary.yhat, type = "l", lty = c(2, 1, 2), xlab
    #           = "regressor", ylab = "estimated mean", main =
    #             "P-spline fit with twice std error bands")
    # }
    # if(mean(m.binomial) == 1) {
    #   matplot(x.seq, summary.yhat, type = "l", lty = c(2, 1, 2), xlab
    #           = "regressor", ylab = "estimated mean", ylim = c(min(
    #             min(y), min(summary.yhat[, 1])), max(max(y), max(
    #               summary.yhat[, 3]))), main =
    #             "P-spline fit with twice std error bands")
    #   matpoints(x, y, type = "p", pch = "O")
    # }
    ps.predict <- NULL
    if(!missing(x.predicted)) {
      ps.predict <- pspline.predictor(x.predicted, knots, link, coef,
                                      q, var.beta, dispersion.parm)
    }
    llist <- list()
    llist$family <- family
    llist$link <- link
    llist$ps.intervals <- ndx
    llist$order <- d
    llist$degree <- q
    llist$lambda <- lambda
    llist$aic <- aic
    llist$deviance <- dev
    llist$eff.df <- trace
    llist$df.resid <- n - trace
    llist$dispersion.param <- dispersion.parm
    llist$summary.predicted <- ps.predict$summary.pred
    llist$coef <- coef
    llist
  }
"pspline.fitter"<-
  function(family, link, n.col, m.binomial, r.gamma, y, b, p, p.ridge, nix,
           nix.ridge, ridge.adj, wts, ...)
  {
    coef.est <- rep(1, ncol(b))
    if(family == "binomial") {
      mu <- (y + 0.5 * m.binomial)/2
    }
    if(family == "Gamma" || family == "poisson") {
      mu <- (y + 3)
    }
    if(family == "gaussian") {
      mu <- rep(mean(y), length(y))
    }
    it <- 0
    repeat {
      if(it == 0) {
        if(link == "identity") {
          eta <- mu
        }
        if(link == "log") {
          eta <- log(mu)
        }
        if(link == "sqrt") {
          eta <- sqrt(mu)
        }
        if(link == "logit") {
          eta <- log(mu/(m.binomial - mu))
        }
        if(link == "recipical") {
          eta <- 1/mu
        }
        if(link == "probit") {
          eta <- qnorm(mu/m.binomial)
        }
        if(link == "cloglog") {
          eta <- log( - log(1 - mu/m.binomial))
        }
        if(link == "loglog") {
          eta <-  - log( - log(mu/m.binomial))
        }
      }
      it <- it + 1
      if(it > 25)
        break
      if(link == "identity") {
        mu <- eta
        h.prime <- 1
      }
      if(link == "log") {
        mu <- exp(eta)
        h.prime <- mu
      }
      if(link == "sqrt") {
        mu <- eta^2
        h.prime <- 2 * eta
      }
      if(link == "logit") {
        mu <- m.binomial/(1 + exp( - eta))
        h.prime <- mu * (1 - mu/m.binomial)
      }
      if(link == "recipical") {
        mu <- 1/eta
        h.prime <-  - (mu^2)
      }
      if(link == "probit") {
        mu <- m.binomial * pnorm(eta)
        h.prime <- m.binomial * dnorm(eta)
      }
      if(link == "cloglog") {
        mu <- m.binomial * (1 - exp( - exp(eta)))
        h.prime <- (m.binomial) * exp(eta) * exp( - exp(eta))
      }
      if(link == "loglog") {
        mu <- m.binomial * exp( - exp( - eta))
        h.prime <- m.binomial * exp( - eta) * exp( - exp( - eta
        ))
      }
      if(family == "gaussian") {
        w <- rep(1, length(y))
      }
      if(family == "poisson") {
        w <- h.prime^2/mu
      }
      if(family == "binomial") {
        w <- h.prime^2/(mu * (1 - mu/m.binomial))
      }
      if(family == "Gamma") {
        w <- (r.gamma * h.prime^2)/mu^2
      }
      u <- (y - mu)/h.prime + eta
      if(ridge.adj > 0) {
        f <- lsfit(rbind(b, p, p.ridge), c(u, nix, nix.ridge),
                   wt = c(wts, nix + 1, nix.ridge + 1) * c(w, (nix +
                                                                 1), (nix.ridge + 1)), intercept = F)
      }
      if(ridge.adj == 0) {
        f <- lsfit(rbind(b, p), c(u, nix), wt = c(wts, nix + 1) *
                     c(w, (nix + 1)), intercept = F)
      }
      coef.old <- coef.est
      coef.est <- as.vector(f$coef)
      d.coef <- max(abs((coef.est - coef.old)/coef.old))
      if(d.coef < 1e-008)
        break
      # print(c(it, d.coef))
      eta <- b %*% coef.est
    }
    if(it > 24) {
      warning(paste("parameter estimates did NOT converge in 25 iterations"
      ))
    }
    llist <- list(coef = coef.est, mu = mu, f = f, w = w * wts)
    return(llist)
  }
"pspline.checker"<-
  function(family, link, degree, order, ps.intervals, lambda, ridge.adj, wts)
  {
    if(link == "default" && family == "gaussian") {
      link <- "identity"
    }
    if(link == "default" && family == "poisson") {
      link <- "log"
    }
    if(link == "default" && family == "binomial") {
      link <- "logit"
    }
    if(link == "default" && family == "Gamma") {
      link <- "log"
    }
    if(family != "binomial" && family != "gaussian" && family != "poisson" &&
       family != "Gamma") {
      warning(paste("Improper FAMILY option. Choose: gaussian, poisson, binomial or Gamma"
      ))
    }
    if((family == "binomial") && (link != "logit" && link != "probit" &&
                                  link != "cloglog" && link != "loglog")) {
      warning(paste("Improper LINK option with family=binomial. Choose: logit, probit, loglog, cloglog"
      ))
    }
    if((family == "Gamma") && (link != "log" && link != "recipical" && link !=
                               "identity")) {
      warning(paste("Improper LINK option with family=Gamma. Choose: recipical, log, identity"
      ))
    }
    if((family == "poisson") && (link != "log" && link != "sqrt" && link !=
                                 "identity")) {
      warning(paste("Improper LINK option with family=poisson. Choose: log, sqrt, identity"
      ))
    }
    if((family == "gaussian") && (link != "identity")) {
      warning(paste("Improper LINK option with family=gaussian. Choose: identity"
      ))
    }
    if(degree < 0) {
      degree <- 1
      warning(paste("degree must be non-neg integer: have used 1"))
    }
    if(order < 0) {
      order <- 0
      warning(paste("order must be non-neg integer: have used 0"))
    }
    if(ps.intervals < 2) {
      ps.intervals <- 2
      warning(paste("ps.intervals must be positive integer, > 1: have used 2"
      ))
    }
    if(lambda < 0) {
      lambda <- 0
      warning(paste("lambda cannot be negative: have used 0"))
    }
    if(ridge.adj < 0) {
      ridge.adj <- 0
      warning(paste("ridge.adj cannot be negative: have used 0"))
    }
    if(min(wts) < 0) {
      warning(paste("At least one weight entry is negative"))
    }
    llist <- list(family = family, link = link, degree = degree, order =
                    order, ps.intervals = ps.intervals, lambda = lambda, ridge.adj
                  = ridge.adj, wts = wts)
    return(llist)
  }
"pspline.predictor"<-
  function(x.predicted, knots, link, coef, q, var.beta, dispersion.parm, ...)
  {
    b.pred <- spline.des(knots, x.predicted, q + 1, 0 * x.predicted)$design
    eta.pred <- b.pred %*% as.vector(coef)
    b.pred <- as.matrix(b.pred)
    if(length(x.predicted) > 1) {
      var.pred <- (b.pred) %*% var.beta %*% t(b.pred)
    }
    if(length(x.predicted) == 1) {
      var.pred <- t(b.pred) %*% var.beta %*% (b.pred)
    }
    stdev.pred <- as.vector(sqrt(diag(var.pred)))
    stdev.pred <- sqrt(dispersion.parm) * stdev.pred
    pivot <- as.vector(2 * stdev.pred)
    upper <- eta.pred + pivot
    lower <- eta.pred - pivot
    summary.pred <- cbind(lower, eta.pred, upper)
    if(link == "logit") {
      summary.pred <- 1/(1 + exp( - summary.pred))
    }
    if(link == "probit") {
      summary.pred <- apply(summary.pred, c(1, 2), pnorm)
    }
    if(link == "cloglog") {
      summary.pred <- (1 - exp( - exp(summary.pred)))
    }
    if(link == "loglog") {
      summary.pred <- exp( - exp( - summary.pred))
    }
    if(link == "sqrt") {
      summary.pred <- summary.pred^2
    }
    if(link == "log") {
      summary.pred <- exp(summary.pred)
    }
    if(link == "recipical") {
      summary.pred <- summary.predd <- 1/(summary.pred)
      summary.pred <- summary.predd[, 3:1]
    }
    summary.pred <- as.matrix(summary.pred)
    dimnames(summary.pred) <- list(NULL, c("-2std_Lower", "Predicted",
                                           "+2std_Upper"))
    llist <- list(summary.pred = summary.pred)
    return(llist)
  }
