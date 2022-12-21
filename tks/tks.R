# TODO:
# * add lag phase to growth and/or kill functions e.g. Kx*(1-exp(-ax*t)) (x = g or k) 
# * uncouple TeX from /user/kzth541/bib/
# * add AIC goodness of fit: http://en.wikipedia.org/wiki/Akaike_information_criterion
# * add BIC goodness of fit: http://en.wikipedia.org/wiki/Bayesian_information_criterion
# * override initial guess for N0, kg from .opt file using values from fit.Kg()
# * better parallel execution on cluster
# * error propagation to D_over_Kg via .cov matrix
# * alternatives to optim()? see http://cran.r-project.org/web/views/Optimization.html
#
# BUGS:  * null device 1

# Population mass balance:  N' = Kg*N*(1-N/Nmax) - r(C)*N
# C: antimicrobial agent concentration
# r(C): microbial kill-rate coefficient (non-decreasing function of C)
#
# http://en.wikipedia.org/wiki/Bernoulli_differential_equation
# Bernoulli differential equation: y' + P(x)*y = Q(x)*y^n
# Population mass balance:         N' = Kg*N*(1-N/Nmax) - r(C)*N (has exactly the same form as Bernoulli)
# N' + (rkill-Kg)*N = -(Kg/Nmax)*N^2 -> P = r(C) - Kg; Q = -Kg/Nmax; n = 2
#
# C: antimicrobial agent concentration
# r(C): microbial kill-rate coefficient (non-decreasing function of C)
#
# solving Bernoulli:
# substitution: w = 1/(y^(n-1)) gives
# w'/(1-n) + P*w = Q i.e. solvable first order linear equation
# M(x) = exp((1-n)*integral(P*dx)); [ M=exp(Kg*t) for C=0 ]
# w = (1/M) * (integral((1-n)*Q*M*dx) + const); since n=2 and w=1/y
# y = M / (integral((1-n)*Q*Mdx) + const)
# y = M / ((Kg/Nmax)*integral(Mdx) + const)

library(XML)     # for xml parsing

tam.model <- function(name='PA.Levofloxacin') {
   #
   # returns parameters for various Tam's models
   #
   #cat('tam.model: name =', name, '\n')
   if (name == 'VT2005') {
      src <- list(doi='http://dx.doi.org/10.1093/jac/dki086', journal='JAC 2005, 55, 699')
      C    <- c(0,    0.25, 1,    4,    16,    64)
      Kg   <- c(0.31, 0.61, 0.31, 2.31, 2.85, 2.58)
      Kk   <- c(4.08, 9.60, 4.08, 3.12, 3.41, 3.9)
      C50k <- c(3.53, 0.44, 3.53, 0.74, 2.61, 2.90)
      H    <- c(1.48, 3.82, 1.48, 1.67, 7.74, 0.79)
      Nmax <- c(2.39, 9.99, 2.29, 9.99, 9.99, 9.99)*1e9
      N0   <- c(1.49, 2.50, 1.49, 3.03, 2.62, 3.14)*1e8 # inoculum
      beta <- c(6.64, 3.79, 6.64, 4.48, 4.85, 8.70)
      tau  <- c(11,   21,   11,   6.9,   5.1,  1.4)/1000
      df <- data.frame(C, Kg, Kk, C50k, H, Nmax, N0, beta, tau)
      parms <- c(N0=1e8,   Nmax=10^9.4,
                 Kg=0.5,   Kk=5,       # h^(-1)
                 C50k=0.7,             # mg/L
                 H=3.7,    beta=42, tau=0.01)
      v <- list(agent='Meropenem', bacteria='P.aeruginosa ATCC 27853', agent='Meropenem') # df=df)
   } else if (name == 'VT2007') {
      # parameters for: "Levofloxacin x P. aeruginosa ATCC 27853" from the paper
      # "Modeling of microbial population responses to time-periodic concentrations of antimicrobial agents."
      # Nikolaou, M., A. N. Schilling, G. Vo, K. T. Chang, and V. H. Tam Ann. Biomed Eng. 2007, 35:1458–1470. 
      # http://dx.doi.org/10.1007/s10439-007-9306-x
      parms <- c(N0=1e8,   Nmax=10^9.4,
                 Kg=0.22,  Kk=13.5,     Kb=13.5,   KA=0.7, # h^(-1)
                 C50=0.86, C50A=0.45,   C50b=39.4,         # mg/L
                 H=1.7,    HA=10,       Hb=24)
      src <- list(doi='http://dx.doi.org/10.1007/s10439-007-9306-x', journal='Annals of Biomedical Engineering 2007, 35, 1458')
      v <- list(agent='Levofloxacin', 'bacteria'='P.aeruginosa ATCC 27853', MIC=2)
      if (0) {
        cat('name =', name, '\n')
        cat('parms =', parms, '\n')
        cat('names(parms) =', names(parms), '\n')
      }
   } else if (name == 'SIMPLE') {
      parms <- c(N0=1e8,   Nmax=10^9.4,
                 Kg=0.5,   Kk=5,       # h^(-1)
                 C50k=0.7,             # mg/L
                 H=3.7)
      v <- list()
      src <- NULL     
   }
   v$name <- name
   v$src <- src
   v$parms <- parms
   v$idx <- 1:length(parms)
   v$dt <- 2^(-6)
   v$solver <- 'bernoulli'
   return (v)
}

rkill.SIMPLE <- function(C, t, model=NULL, debug=0) {
   if (length(C) != length(t)) {
     cat('rkill.SIMPLE: ERROR: length(C) != length(t); length(C), length(t) =', length(C), length(t), '\n')
     stop('length(C) != length(t)')
   }
   p <- as.list(model$parms)
   if (is.null(C)) stop('C is NULL')
   rC <- p$Kk*C^p$H/(C^p$H + p$C50k^p$H + .tiny)
   v <- list(rC=rC, Kg=p$Kg, Nmax=p$Nmax, N0=p$N0)
   return(v)
}

rkill.VT2005 <- function(C, t, model=NULL, debug=0) {
   if (length(C) != length(t)) {
     cat('rkill.VT2005: ERROR: length(C) != length(t); length(C), length(t) =', length(C), length(t), '\n')
     stop('length(C) != length(t)')
   }
   p <- as.list(model$parms)
   if (is.null(C)) stop('C is NULL')
   alpha <- 1 + p$beta*(1 - exp(-C*t*p$tau))
   rC <- p$Kk*C^p$H/(C^p$H + (alpha*p$C50k)^p$H + .tiny)
   v <- list(rC=rC, alpha=alpha, Kg=p$Kg, Nmax=p$Nmax, N0=p$N0)
   #v <- list(rC=rC)
   return(v)
}

rkill.VT2007 <- function(C, t, model=NULL, debug=0) {
   if (is.null(C)) stop('C is NULL')
   if (length(C) != length(t)) {
     cat('rkill.VT2007: ERROR: length(C) != length(t); length(C), length(t) =', length(C), length(t), '\n')
     stop('length(C) != length(t)')
   }
   p <- as.list(model$parms)
   b <- p$Kb*C^p$Hb/(C^p$Hb + p$C50b^p$Hb + .tiny)
   R <- p$Kk*C^p$H/(C^p$H + p$C50^p$H + .tiny) - b
   A <- p$KA*C^p$HA/(C^p$HA + p$C50A^p$HA + .tiny)
   rC <- R*exp(-A*t) + b
   v <- list(rC=rC, b=b, R=R, A=A, Kg=p$Kg, Nmax=p$Nmax, N0=p$N0)
   #v <- list(rC=rC)
   return(v)
}

fit.Kg <- function(t, N) {
  # solution for ODE: dN/dt = Kg*(1 - N/Nmax)*N
  # for C=0: N(t) = Nx*N0*exp(Kg*t)/(N0*exp(Kg*t) + Nx-N0); dN/dt = Nx*N0*(Nx-N0)*Kg*exp(Kg*t)/(N0*exp(Kg*t) - Nx-N0)^2
  # t <- seq(from=0, to=24, by=0.5); Nt <- Nx*N0*exp(Kg*t)/(N0*exp(Kg*t) + Nx-N0); plot(t, Nt, type='l')
  # Kg <- 0.593292; N0 <- 10^6.872358; Nx <- 10^8.607037
  #
  if (length(t) != length(N)) {
    cat('fit.Kg: ERROR: length(t) != length(N); length(t), length(N) =', length(t), length(N), '\n')
    stop('length(t) != length(N)')
  }
  m <- !is.na(N) # guard against NA values
  N <- N[m]; t <- t[m]
  N0 <- N[1]; Nmax <- N[length(N)] # <<< initial estimate
  Kg <- log(N[1]/N[2])/(t[1]-t[2]) # << initial estimate
  if ((Nmax-N0) < 1) {
    cat('fit.Kg: WARNING: (Nmax-N0) < 1; Nmax, N0 =', Nmax, N0, '\n')
    #stop('Nmax-N0 < 1')
  }
  par <- c(N0=N0, Nmax=Nmax, Kg=Kg)
  fn <- function(par, t, N) {
    N0 <- par['N0']; Nmax <- par['Nmax']
    v <- 0
    for (i in 1:length(t)) {
      e <- exp(par['Kg']*t[i])
      d <- Nmax*N0*e/(N0*e+Nmax-N0) - N[i]
      v <- v + d^2
    }
    return (log10(v))
  }
  opt <- optim(par, fn, gr=NULL, t, N, method='BFGS')
  if (0) {
    Nx <- opt$par['Nmax']; N0 <- opt$par['N0']; Kg <- opt$par['Kg']
    Nt <- Nx*N0*exp(Kg*t)/(N0*exp(Kg*t) + Nx-N0)
    plot(t, log10(N), type='p')
    lines(t, log10(Nt))
  }
  if (opt$convergence != 0) {
    cat('fit.Kg: ERROR: opt$convergence != 0\n')
    stop('opt$convergence != 0')
  }
  return (opt$par)
}
if (0) {
  N <- c(6.72, 7.97, 8.74, 8.82, 3.32)
  t <- c(0.0, 2.0, 4.0, 6.0, 8.0, 24.0)
  v <- fit.Kg(t, N)
}


get.tiny <- function() {
   for (i in 1024:2048) {
     x <- 2^(-i)
     if (x==0) break
     tiny <- x; it <- i
   }
   #return (list(i=it, tiny=tiny))
   return (tiny)
}

get.huge <- function() {
   for (i in 512:1024) {
     x <- 2^i
     if (x==Inf) break
     huge <- x; it <- i
   }
   #return (list(i=it, huge=huge))
   return (huge)
}

get.precision <- function() {
   for (i in 32:64) {
     p <- 2^(-i)
     x <- 1 + p
     if ((x-1) == 0) { break
     } else {
       .precision <- p; it <- i
     }
   }
   #return (list(i=it, precision=precision))
   return (.precision)
}

# .tiny [1] 4.940656e-324
# .HUGE [1] 8.988466e+307
# .Machine$double.eps [1] 2.220446e-16

dN.dt <- function(t, N, model) {
    if (model$debug > 1) cat ('dN.dt.vt: N=', N, 't=', t, '\n')
    p <- model$parms
    v <- rkill(model$C, t, model)
    if (model$debug > 1) cat("dN.dt: p['Kg'], N, p['Nmax'] =", p['Kg'], N, p['Nmax'], '\n')
    dN <- p['Kg']*(1-N/p['Nmax'])*N - v$rC * N
    if (length(dN) != 1) {
      cat('dN.dt: ERROR: length(N) =', length(N), '\n')
      stop('length(dN) != 1')
    }
    names(dN) <- 'dN/dt'
    if (model$debug) {
      global.dN.dt.count <<- global.dN.dt.count + 1
      global.dN.dt[global.dN.dt.count,] <<- t(model$parms)
    }
    return (list(dN=dN))
}

P.calc <- function(t, model, debug=0) {
   v <- rkill(model$C, t, model, debug)
   v$P <- v$rC - v$Kg
   return(v)
}

M.calc <- function(P, dt=1, n=2) {
   m <- is.finite(P)
   if (length(P) != length(P[m])) {
     cat('M.calc: ERROR: P =', P, '\n')
   }
   iP <- i.diffinv(P, dt)
   v <- exp((1-n)*iP)
   return (v)
}

i.diffinv <- function(f, dt, xi=0) { # indefinite integral...
   m <- is.finite(f)
   if (length(f) != length(f[m])) {
     cat('i.diffinv: ERROR: f =', f, '\n')
   }
   s <- diffinv(f, xi=xi)
   n <- length(f)
   a <- (s[1:n] + s[2:(n+1)])*dt/2
   a <- a - a[1]
   return (a)
}

C.calc <- function(times, C.max=1, injection.times=0, t.half=Inf, debug=0) {
   if (debug) {
      cat('C.calc: C.max, t.half =', C.max, t.half, '\n')
      cat('C.calc: injection.times =', injection.times, '\n')
      cat('C.calc: times =', time, '\n')
   }
   v <- numeric(length(times))
   i <- 0
   for (t in times) {
      i <- i + 1
      s <- 0
      for (it in injection.times) {
         dt <- t - it
         if (dt >= 0) s <- s + exp(-dt/t.half)
      }
      v[i] <- s
   }
   return (v*C.max)
}

solve.bernoulli <- function(model, times=NULL, debug=0) {
   #if (is.null(dt)) dt <- (max(times)-min(times))/(length(times)-1)
   if (is.null(times)) times <- model$times
   #v <- P.calc(model$times, model, debug)
   if (0) {
     cat('solve.bernoulli: length(times), times =', length(times), times[1:min(16, length(times))], '\n')
     cat('solve.bernoulli: length(model$times), model$times =', length(model$times), times[1:min(16, length(model$times))], '\n')
   }
   v <- P.calc(times, model, debug)
   m <- is.finite(v$P)
   if (length(v$P) != length(v$P[m])) {
     cat('solve.bernoulli: ERROR: v$P =', v$P, '\n')
     cat('solve.bernoulli: ERROR: model$parms =', model$parms, '\n')
   }
   M <- M.calc(v$P, model$dt)
   m <- is.finite(M)
   if (length(M) != length(M[m])) {
     cat('solve.bernoulli: ERROR: M =', M, '\n')
     cat('solve.bernoulli: ERROR: model$parms =', model$parms, '\n')
   }
   iM <- i.diffinv(M, model$dt)
   Nmax <- model$parms['Nmax']; N0 <- model$parms['N0']
   N <- Nmax * M / (v$Kg*iM + Nmax/N0) 
   #m <- match(model$exp.data$times, model$times)
   #m <- match(model$exp.data$times, times)
   #df <- data.frame(times=model$times[m], yp=N[m])
   df <- data.frame(times=model$times, yp=N)
   #msg <- sprintf('solve.bernoulli: nrow(df), ncol(df) = %d %d', nrow(df), ncol(df), '\n')
   #stop(msg)
   return (df)
}
# df <- solve.bernoulli(model, times, C)

set.device <- function(fn, device, width=4, height=4, units='in', pointsize=12) {
   fn <- sprintf('%s.%s', fn, device)
   if (device=='png') {
     png(filename=fn, width=width, height=height, units=units, pointsize=pointsize, res=360)
   } else if (device=='eps') {
     postscript(fn, width=width, height=height, pointsize=pointsize, horizontal=F, onefile = F, paper = 'special')
   }
}

# p <- as.numeric(parms[1,]); names(p) <- names(parms)
# sumsq(log(p), model, rmse.only=FALSE)
sumsq <- function(par, model, times=NULL, rmse.only=TRUE, debug=0, ...) {
  debug <- 0
  #cat('sumsq: class(par), par =', class(par), par, '\n')
  #cat('sumsq: =', class(par), par, '\n')
  #if (is.null(model$idx)) model$idx <- 1:length(model$parms)
  if (is.null(model$idx)) stop('is.null(model$idx)')
  #n <- names(model$parms)
  model$parms[model$idx] <- exp(par)
  #print(model$parms)
  #stop()
  if (outbound.parms(model)) {
    v <- sum(par^2)
    if (model$debug) cat('sumsq: after outbound.parms: v =', v, '\n')
    if (rmse.only) return (v)
  }
  if (is.null(times)) times <- model$exp.data$times
  nc <- length(model$exp.data$C)
  #yp.all <- model$exp.data$log10CFU
  yp.all <- matrix(nrow=nrow(model$exp.data$log10CFU), ncol=length(model$times))
  if (debug) {
    cat('sumsq: length(times), length(model$times) =', length(times), length(model$times), '\n')
    cat('sumsq: dim(yp.all) =', dim(yp.all), '\n')
    cat('sumsq: length(model$times) =', length(model$times), '\n')
    cat('sumsq: length(times) =', length(times), '\n')
    flush(stdout())
    #stop()
  }
  #dt <- (max(model$times)-min(model$times))/(length(times)-1)
  rmse <- np <- 0
  ym <- match(model$exp.data$times, model$times)
  for (i in 1:nc) {
    yo <- model$exp.data$log10CFU[i,]
    #if (model$debug) cat('sumsq: i, C =', i, model$C, '\n')
    if (model$solver=='ode') {
       model$C <- model$exp.data$C[i]
       #o <- ode(y=model$parms['N0'], times=model$exp.data$times, func=dN.dt, parms=model, verbose=F)
       o <- ode(y=model$parms['N0'], times=model$times, func=dN.dt, parms=model, verbose=F)
       df <- as.data.frame(o)
       if (ncol(df) > 2) { y <- apply(df[,-1],1,sum)
       } else { y <- df[,2]}
    } else if (model$solver=='bernoulli') {
       model$C <- C.calc(model$times, C.max=model$exp.data$C[i])
       df <- solve.bernoulli(model, times=model$times, debug=0)
    }
    #m <- match(model$exp.data$times, model$times)
    #y <- df[m,2]
    y <- df[,2]
    #msg <- sprintf('sumsq: length(y) = %d', length(y))
    #stop(msg)
    if (debug) {
       cat('sumsq: length(y) =', length(y), '\n')
       cat('sumsq: yo =', yo, '\n')
       cat('sumsq: ym =', ym, '\n')
       cat('sumsq: dim(yp.all) =', dim(yp.all), '\n')
       #stop()
    }
    y[y<1] <- 1
    yp <- log10(y)
    yp.all[i,] <- yp
    np <- np + length(yo)
    #rmse <- rmse + sum((yo-yp)^2, na.rm=T)
    rmse <- rmse + sum((yo-yp[ym])^2, na.rm=T)
  }
  rmse <- sqrt(rmse/np)
  if (debug) {
    cat('sumsq: RMSE=', rmse, 'model$parms=', model$parms, '\n')
    cat('sumsq: length(times), length(model$times) = ', length(times), length(model$times), '\n')
    cat('sumsq: model$times = ', model$times[1:min(16, length(model$times))], '\n')
    cat('sumsq: times = ', times[1:min(16, length(times))], '\n')
    cat('sumsq: ym =', ym, '\n')
    stop()
  }
  if (rmse.only) {
    v <- rmse
  } else {
    v <- list(rmse=rmse, yo=model$exp.data$log10CFU, yp=yp.all)
  }
  return (v)
}
# v.paper <- sumsq(log(par.paper), model, rmse.only=F)

plot.model <- function(model, vp=NULL, fn=NULL, device=NULL, main=NULL, cex=1.2, exp.lines=TRUE) {
  t <- model$exp.data$times
  C <- model$exp.data$C
  y <- model$exp.data$log10CFU
  xlim <- range(t); xlim[2] <- 1.4*xlim[2]
  ylim <- c(0, max(y, na.rm=T))
  #ylim <- range(y, na.rm=T)
  nr <- nrow(y)
  cols <- rainbow(nr)
  if (is.null(main)) main <- sprintf('%s:%s / %s', model$bacteria['id'], model$agent['name'], model$name)
  if (!is.null(fn) & !is.null(device)) set.device(fn, device)
  lg <- character(nr)
  plot(t, t, xlim=xlim, ylim=ylim, main=main, type='n', xlab='time [h]', ylab='log10(CFU/mL)')
  for (i in 1:nr) {
    if (exp.lines) lines(t, y[i,], col=cols[i], lwd=2)
    lg[i] <- sprintf('%sx', as.character(C[i]))
    points(t, y[i,], col=cols[i])
  }
  xL <- 0.75*xlim[2]; yL <- 0.95*ylim[2]; dy <- (ylim[2]-ylim[1])/20
  if (!is.null(vp)) {
    #for (i in 1:nr) lines(t, vp$yp[i,], col=cols[i], lwd=2)
    for (i in 1:nr) lines(model$times, vp$yp[i,], col=cols[i], lwd=2)
    xt <- xL; yt <- yL + dy
    text(xt, yt, sprintf('RMSE = %4.2f', vp$rmse), adj=c(0,0), cex=cex)
  }
  legend(xL, yL, legend=lg, col=cols, lwd=rep(2,nr), bg='grey', cex=cex, adj=c(0.2,0.5))
  if (!is.null(fn) & !is.null(device)) dev.off(dev.cur())
}

shake.parms <- function(model, par=NULL, C.max=128, t.max=24, sd=1) {
    #cat('shake.parms: sd =', sd, '\n')
    if (is.null(par)) {
      #cat('shake.parms: class(model$parms)=', class(model$parms), '\n')
      par <- log(model$parms)[model$idx]
    }
    pn <- names(par)
    par <- rnorm(length(par), mean=par, sd=sd)
    names(par) <- pn
    m <- par > model$upper; par[m] <- model$upper[m]
    m <- par < model$lower; par[m] <- model$lower[m]
    return (par)
}

shake.data <- function(x, sd, preserve.zeros=TRUE) {
  # x <- model$exp.data$log10CFU
  n <- length(x)
  m <- x == 0
  x <- x + rnorm(n, sd=sd)
  if (preserve.zeros) x[m] <- 0
  return (x)
}

# v <- par.limits(model)
par.limits <- function(model, C.max=128, t.max=24, factor=2) {
   # limits for log(model$parms)
   # model$parms
   #cat('par.limits: class(model$parms)=', class(model$parms), '\n')
   if (is.null(model$upper)) {
     upper <- rep(+5, length(model$parms))
     lower <- rep(-8, length(model$parms))
   } else {
     upper <- model$upper
     lower <- model$lower
   }
   names(lower) <- names(model$parms)
   names(upper) <- names(model$parms)
   lower['N0'] <- log(10^3);    upper['N0'] <- log(10^8)
   lower['Nmax'] <- log(10^5);  upper['Nmax'] <- log(10^11)
   #
   # set limits for N0, Nmax and Kg from experimental data
   #
   #t <- model$exp.data$times; N <- 10^(model$exp.data$log10CFU[1,])
   #p0 <- fit.Kg(t, N)
   #lower[names(p0)] <- log(p0/factor)
   #upper[names(p0)] <- log(p0*factor)
   L <- lower[model$idx]
   U <- upper[model$idx]
   if (0) {
   cat('par.limits: model$idx =', model$idx, '\n')
   cat('par.limits: L =', L, '\n')
   cat('par.limits: U =', U, '\n')
   cat('par.limits: names(U) =', names(U), '\n')
   #stop()
   }
   return (list(lower=L, upper=U))
}

outbound.parms <- function(model) {
  #cat('outbound.parms: class(model$parms)=', class(model$parms), '\n')
  par <- log(model$parms)[model$idx]
  w.lower <- which(par < model$lower)
  w.upper <- which(par > model$upper)
  err.lower <- err.upper <- 0
  err.msg <- c('lower', 'upper')
  if (length(w.lower) > 0) err.lower <- 1
  if (length(w.upper) > 0) err.upper <- 1
  err <- err.lower + err.upper
  if (err) {
    if (model$debug) {
      if (err.lower) cat('outbound.parms: ERROR: lower bound violated\n')
      if (err.upper) cat('outbound.parms: ERROR: upper bound violated\n')
      cat('outbound.parms: names(par) =', names(par), '\n')
      cat('outbound.parms: par =', par, '\n')
      cat('outbound.parms: model$lower =', model$lower, '\n')
      cat('outbound.parms: model$upper =', model$upper, '\n')
      cat('outbound.parms: w.upper, names(par)[w.upper], par[w.upper] =', w.upper, names(par)[w.upper], par[w.upper], '\n')
      cat('outbound.parms: w.lower, names(par)[w.lower], par[w.lower] =', w.lower, names(par)[w.lower], par[w.lower], '\n')
      #stop('outbound.parms: bounds violated')
    }
  }
  return (err)
}
# outbound.parms(model, lower, upper)

par.sens <- function(model, opts, sd=NULL, nrep=NULL, minimizer='optim', method='BFGS', maxit=1e4, trace=0, device=NULL, verbose=TRUE) {
  #
  time.usage <- data.frame(t(c(0,0,0))); names(time.usage) <- c('user', 'system', 'elapsed')
  #cat('par.sens: class(model$parms)=', class(model$parms), '\n')
  # TBD: need to add code to respect preoptimized N0, Nmax and Kg!!!
  if (is.null(opts$optimized.parms)) {
    par0 <- log(model$parms)[model$idx]
    origin <- 'DEFAULT'
  } else {
    cat('opts$optimized.parms =', opts$optimized.parms, '\n')
    x <- read.table(opts$optimized.parms, header=T)
    i <- sample((1:nrow(x)), 1)
    p0 <- x[i, names(model$parms)]
    origin <- sprintf('%d:%s.%s', i, as.character(x[i,c('Bug')]), as.character(x[i,c('Drug')]))
    p0 <- as.numeric(p0)
    names(p0) <- names(model$parms)
    par0 <- log(p0)[model$idx]
    cat('i, nrow(x) =', i, nrow(x), '\n')
  }
  cat('par.sens: starting parameter values: origin, names(par0), exp(par0) =', origin, names(par0), exp(par0), '\n')
  #stop()
  par <- par0 + rnorm(length(model$idx), sd=sd$opt)
  if (0) {
    cat('par.sens: sd$exp =', sd$exp, '\n')
    x <- var(par-par0)
    cat('par.sens: var(par-par0)=', x, '\n')
    stop()
  }
  df.0 <- df.1 <- as.data.frame(t(c(0, par)))
  names(df.0) <- names(df.1) <- c('RMSE', names(par))
  if (0) {
    cat('par.sens: dim(df.0)=', dim(df.0), '\n')
    cat('par.sens: names(par)=', names(par), '\n')
    cat('par.sens: model$idx=', model$idx, '\n')
    cat('par.sens: model$lower=\n')
    print(model$lower)
    cat('par.sens: model$upper=\n')
    print(model$upper)
    stop()
  }
  optim.control <- list(trace=trace, maxit=maxit)
  nlminb.control <- list(trace=trace, eval.max=maxit, iter.max=maxit)
  opt.list <- par.lst <- model.lst <- list()
  data <- model$exp.data$log10CFU
  fail.count <- 0
  ij <- 0
  #cat('as.numeric(sd)=', as.numeric(sd), '\n')
  for (i in 1:nrep) {
    t0 <- proc.time()
    ij <- ij + 1
    model$exp.data$log10CFU <- shake.data(data, sd$exp)
    if (nrow(df.1) > 0) {
    #if (0) {
      o <- order(df.1$RMSE)
      #m <- match('RMSE', names(df.1))
      #p <- df.1[o[1],-m]
      p <- df.1[o[1], 2:ncol(df.1)]
      pn <- names(p)
      p <- as.numeric(p)
      names(p) <- pn
      par.lst[[ij]] <- par <- shake.parms(model, par=p, sd=sd$opt)
    } else {
      par.lst[[ij]] <- par <- shake.parms(model, par=NULL, sd=sd$opt)
    }
    cat('par.sens: ij=', ij, ': par =', par, '\n')
    ssq0 <- sumsq(par, model) # main='BEFORE optim')
    df.0[ij,] <- c(log(ssq0), par)
    # other options: nlm, nlminb
    if (minimizer=='optim') {
      v <- optim(par, sumsq, gr=NULL, model, method=method, control=optim.control)
      par <- v$par; value <- v$value; code <- v$convergence
    } else if (minimizer=='nlm') {
      v <- nlm(sumsq, par, model, stepmax=1, iterlim=maxit)
      par <- v$estimate; value <- v$minimum; code <- v$code
    } else if (minimizer=='nlminb') {
      v <- nlminb(par, sumsq, gradient = NULL, hessian = NULL, model, control=nlminb.control,
                  lower=model$lower, upper=model$upper)
      par <- v$par; value <- v$objective; code <- v$convergence
    } else {
      stop('unknown minimizer')
    }
    if (code == 0) { # TBD: change code from nlm and nlminb 
      opt.list[[ij]] <- v
      df.1[ij,] <- x <- c(log(value), par)
      if (!prod(is.finite(x))) {
        cat('WARNING: par.sens: optimization succeded, but: code, x =', code, x, '\n')
      }
    } else {
      fail.count <- fail.count + 1
      cat('WARNING: optimization failed: fail.count, code =', fail.count, code, '\n')
    }
    t1 <- proc.time()
    time.usage[ij,] <- (t1-t0)[1:3]
    ssq1 <- sumsq(par, model) #  main='AFTER optim')
    threshold <- (2^5)*.precision
    if (abs(value-ssq1) > threshold) {
       cat('par.sens: WARNING: abs(value-ssq1) >ssq1,  value-ssq1, threshold, .precision =', ssq1, value-ssq1, threshold, .precision, '\n')
       #stop('abs(value-ssq1) > 2*.precision')
    }
    if (verbose) cat('par.sens: before and after optim: i, ssq0, ssq1 =', i, ssq0, ssq1, '\n')
  }
  o <- order(df.1$RMSE)
  return (list(df.0=exp(df.0[o,]), df.1=exp(df.1[o,]), time.usage=time.usage, opt=opt.list, par=par.lst, minimizer=minimizer))
  #return (list(df.0=df.0, df.1=df.1, time.usage=time.usage, opt=opt, par=par))
} # end-of-par.sens

get.matrix <- function(str, ncol) {
  if (!is.null(str)) {
    y <- strsplit(str, split=' ')[[1]]
    m <- y != ""
    z <- t(matrix(data=as.numeric(y[m]), ncol=ncol))
  } else {
    z <- NULL
  }
  return (z)
}

to.numbers <- function(s) {
  # converts string to numbers
  v <- as.numeric(strsplit(s, split=' ')$values)
  v <- v[is.finite(v)]
  return (v)
}

#> names(model)
# [1] "name"     "src"      "agent"    "bacteria" "MIC"      "parms"    "exp.data" "debug"    "idx"      "solver"   "dt"       "times"
#get.xml.data <- function(input, dt=0.5, solver='bernoulli', name='VT2005', debug=0) {
get.xml.data <- function(input, opts) {
  xml <- xmlTreeParse(input)
  x <- xmlToList(xml)
  times <- to.numbers(x$time['values'])
  C <- to.numbers(x$concentration['values'])
  q <- all.equal(unique(C), C)
  err <- FALSE
  if (!is.logical(q)) {
    cat('get.xml.data: ERROR: C duplicates', q, '\n')
    err <- TRUE
  }  
  q <- all.equal(unique(times), times)
  if (!is.logical(q)) {
    cat('get.xml.data: ERROR: times duplicates', q, '\n')
    err <- TRUE
  }
  #if (err) stop('get.xml.data: ERROR: duplicates in C or times')
  ncol <- length(times)
  nrow <- length(C)
  time.units <- x$time['units']
  C.units <- x$concentration['units']
  MIC <- as.numeric(x$MIC['values'])
  MIC.units <- x$MIC['units']
  bacteria.id <- x$bacteria['id']
  if (bacteria.id == '') bacteria.id <- x$bacteria['strain']
  agent <- x$agent['name']
  options(warn=0) # to avoid: Error in as.vector(data): (converted from warning) NAs introduced by coercion
  Mean <- get.matrix(x$Mean['values'], ncol=nrow)
  #stop('aaa')
  CFU <- get.matrix(x$CFU['values'], ncol=nrow)
  StDev <- get.matrix(x$StDev['values'], ncol=nrow)
  options(warn=2)
  if (is.null(Mean)) {
    Mean <- log10(CFU)
    m <- !is.finite(Mean)
    Mean[m] <- 0
  }
  exp.data <- list(times=times, C=C, CFU=CFU, log10CFU=Mean, StDev=StDev)
  t.begin <- min(exp.data$times); t.end <- max(exp.data$times)
  n <- (t.end-t.begin)/opts$dt
  tm <- tam.model(opts$model.name)
  model <- list(publication=x$publication, agent=x$agent, bacteria=x$bacteria, exp.data=exp.data, dt=opts$dt, solver=opts$solver, name=opts$model.name,
                parms=tm$parms, fixed=opts$fixed, t.begin=t.begin, t.end=t.end, times=t.begin + (t.end-t.begin)*(0:n)/n, debug=opts$debug)

  m <- match(opts$fixed, names(model$parms))
  model$idx <- setdiff(1:length(model$parms), m)
  if (0) {
    cat('get.xml.data: opts$fixed =', opts$fixed, '\n')
    cat('get.xml.data: model$idx =', model$idx, '\n')
    cat('get.xml.data: n, model$times =', length(model$imes), model$times[1:min(16, length(model$times))], '\n')
    stop()
  }
  model$lower <- opts$lower
  model$upper <- opts$upper
  if (!is.null(model$lower) & !is.null(model$upper)) {
    x <- 0.5*(model$lower+model$upper)
    if (prod(names(x) == names(model$parms))) {
      model$parms <- x
    } else {
      cat('ERROR: get.xml.data: names(model$parms)=', names(model$parms), '\n')
      cat('ERROR: get.xml.data: names(model$lower)=', names(model$lower), '\n')
      cat('ERROR: get.xml.data: names(model$upper)=', names(model$upper), '\n')
      stop('parms names mismatch')
    }
  }
  v <- par.limits(model, C.max=max(model$exp.data$C), t.max=max(model$exp.data$times))
  model$lower <- v$lower; model$upper <- v$upper
  if (0) {
     cat('par.sens: after par.limits: model$lower =', model$lower, '\n')
     cat('par.sens: after par.limits: model$upper =', model$upper, '\n')
     stop()
  }

  t <- model$exp.data$times; N <- 10^(model$exp.data$log10CFU[1,])
  p0 <- fit.Kg(t, N)
  model$parms[names(p0)] <- p0
  model$lower[names(p0)] <- log(p0/opts$limit.factor)
  model$upper[names(p0)] <- log(p0*opts$limit.factor)

  if (0) {
    cat('get.xml.data: names(p0) =', names(p0), '\n')
    cat('get.xml.data: p0 =', p0, '\n')
    cat('get.xml.data: names(model$parms) =', names(model$parms), '\n')
    cat('get.xml.data: model$parms =', model$parms, '\n')
    #stop()
    cat('get.xml.data: log10(model$parms) =', log10(model$parms), '\n')
    cat('get.xml.data: log10(exp(model$upper)) =', log10(exp(model$upper)), '\n')
    cat('get.xml.data: log10(exp(model$lower)) =', log10(exp(model$lower)), '\n')
    stop()
  }

  return (model)
}

#if (0)
#  model <- get.xml.data(input, opts)
#  stop('***')
#}

make.scatterplot <- function(yo, yp, fn=NULL, device=NULL, adj=0) {
  if (length(yo) != length(yp)) {
    cat('ERROR: make.scatterplot: length(yo) != length(yp); length(yo), length(yp)=', length(yo), length(yp), '\n')
    stop('length(yo) != length(yp)')
  }
  xlim <- ylim <- range(yo, yp, na.rm=T)
  #m <- match(model$exp.data$times, model$times)
  #mse <- sum((yo-yp[m])^2, na.rm=T)/length(yo)
  mse <- sum((yo-yp)^2, na.rm=T)/length(yo)
  var.o <- var(as.numeric(yo), na.rm=T)
  q2 <- 1 - mse/var.o
  if (!is.null(fn)) set.device(fn, device) 
  #cat('make.scatterplot: m =', m, '\n')
  #cat('make.scatterplot: length(yo), length(yp[m]) =', length(yo), length(yp[m]), '\n')
  plot(yo, yp, xlab='observed log10(CFU/mL)', ylab='fitted log10(CFU/mL)', xlim=xlim, ylim=ylim)
  abline(0, 1)
  dx <- dy <- (xlim[2]-xlim[1])/10
  xt <- xlim[1] + 0*dx; yt <- ylim[2] - dy
  text(xt, yt, 'q2 = 1 - MSE/Var(Y.obs)', adj=adj)
  text(xt, yt-1*dy, sprintf('Var(Y.obs) = %4.2f', var.o), adj=adj)
  text(xt, yt-2*dy, sprintf('MSE = %4.2f', mse), adj=adj)
  text(xt, yt-3*dy, sprintf('q2 = %4.2f', q2), adj=adj)
  if (!is.null(fn)) dev.off(dev.cur())
}

is.valid <- function(m) {
  #
  #
  #
  i <- 0; err <- list()
  if (is.null(m)) {i <- i+1; err[[i]] <- 'm is NULL'}
  if (is.null(m$dt)) {i <- i+1; err[[i]] <- 'm$dt is NULL'}
  if (is.null(m$solver)) {i <- i+1; err[[i]] <- 'm$solver is NULL'}
  if (length(err) > 0) {
    for (x in err) cat(x, '\n')
    stop('invalid model')
  }
  return (err)
}
# err <- is.valid(model)

cluster.parms <- function(x, cut.RMSE=NULL, method='kmeans', nc=6, plt=NULL, fn=NULL, device=NULL, ...) {
  #
  #
  #
  if (!is.null(cut.RMSE)) { # only rows with x$RMSE < cut.RMSE will be used
    m <- x$RMSE < cut.RMSE
    x <- x[m,]
  }
  if (nrow(x) < 3) {
    cat('cluster.parms: WARNING: to few rows for clustering, returning x\n')
    #return (NULL)
    return (x)
  }
  nc <- min(nc, nrow(x)-1)
  if (method == 'pam') {
    require(cluster)
    v <- pam(x, nc)
    cl <- v$medoids
  } else if (method == 'kmeans') {
    v <- kmeans(x, nc)
    cl <- v$centers
  } else {
    cat('cluster.parms: ERROR: unknown method =', method, '\n')
    stop('unknown method')
  }
  if (is.numeric(plt)) {
    if (!is.null(fn) & !is.null(device)) set.device(fn, device)
    plot(x[,plt])
    points(cl[,plt], col='red', pch=16)
    if (!is.null(fn) & !is.null(device)) dev.off(dev.cur())
  }
  cl <- as.data.frame(cl)
  o <- order(cl$RMSE)
  return (cl[o,])
}
# cluster.parms(x, plt=4:5)

find.mostlikely <- function(p) {
  #
  # finds most likely values of p
  #
  v <- as.numeric(p[1,])
  names(v) <- names(p)
  for (i in 1:ncol(p)) {
    if (nrow(p) > 1) {
      d <- density(p[,i])
      max.d <- max(d$y)
      m <- match(max.d, d$y)
      v[i] <- d$x[m]
    } else {
      v[i] <- p[1,i]
    }
  }
  return (v)
}

add.row <- function(x) {
  # if there is only one row make it two
  if (nrow(x) < 2) {
    if (nrow(x) < 1) {
      cat('ERROR: zero rows in data\n')
      stop('ERROR: zero rows in data')
    } else {
      x <- rbind(x,x)
    }
  }
  return (x)
}

stats <- function(x, cut.RMSE=NULL, delta.RMSE=NULL, debug=0) {
  #if (is.null(fn)) fn <- sprintf('%s.txt', model$name)
  #x <- read.table(fn, header=TRUE)
  if (class(sd) != 'function') stop('class(sd) != "function"')
  if (is.null(cut.RMSE)) cut.RMSE <- min(x$RMSE, na.rm=T) + delta.RMSE
  props <- list(Mean=mean, Median=median, StDev=sd, Min=min, Max=max) # Likely=find.mostlikely)
  nc <- ncol(x); nr <- length(props)+1
  if (debug) {
    cat('stats: class(x), dim(x) =', class(x), dim(x), '\n')
    cat('stats: nr, nc=', nr, nc, '\n')
    #cat('stats: DEBUG\n')
    #stop()
  }
  st <- as.data.frame(matrix(rep(0,nc*nr), ncol=nc))
  names(st) <- names(x)
  row.names(st) <- c(names(props), 'Most.Likely')
  if (!is.null(cut.RMSE)) {
    o <- order(x$RMSE)
    if (!prod(o == 1:nrow(x))) {
      cat('stats: WARNING: unsorted x?\n')
      stop('unsorted x?')
      x <- x[o,]
    }
    m <- x$RMSE < cut.RMSE
    x <- x[m,]
  }
  if (nrow(x) == 0) { # x <- x.orig # needed when all x$RMSE values are the same
    cat('stats: ERROR: nrow(x) == 0\n')
    stop('nrow(x) == 0')
  }
  for (i in 1:length(props)) st[i,] <- apply(x, 2, props[[i]])
  st[nrow(st),] <- exp(find.mostlikely(log(x))) # exp...log... to avoid negative values!
  return (list(stats=st, x=x))
}
# st <- stats(x)

generate.report <- function(model, opts, data=NULL, parms.file=NULL, template=NULL, prefix='prefix', debug=1) {
  cat('generate.report: START\n')
  err <- is.valid(model)
  if (length(err) > 0) {
    cat('generate.report: ERROR: model is not valid:', '\n')
    for (x in err) cat(x, '\n')
    stop()
  }
  if (is.null(data)) {
    df <- read.table(parms.file, header=T)
    cat('optimized parameters read from:', output, '\n')
  } else {
    df <- data
    cat('optimized parameters passed from opt$x\n')
  }
  n0 <- nrow(df)
  .data <<- df <- na.omit(df) # use <<- to pass .data to Sweave through globalEnv
  if (n0 != nrow(df)) cat('WARNING: generate.report: NA values in data: n0, nrow(df) =', n0, nrow(df), '\n')
  st <- stats(df, cut.RMSE=opts$cut.RMSE, delta.RMSE=opts$delta.RMSE)
  df <- st$x
  df.mean <- mean(df)
  df.cov <- cov(df)
  df.median <- apply(df, 2, median)
  x <- rbind(Mean=df.mean, Median=df.median, df.cov) # output Mean and Covariance
  write.table(signif(x, digits=4), file=sprintf('%s.cov', prefix), , quote=F, sep='\t', row.names=T)
  #n1 <- nrow(df)
  if (is.null(template)) {
    cat('generate.report: ERROR: null template:', '\n')
    stop()
  }
  stamp <- time.stamp()
  #target <- sprintf('%s-%s', model$name, stamp)
  target <- prefix
  cmd <- sprintf('cp %s %s', template, target)
  system(cmd)
  #if (is.null(opt)) {
  #  cat('generate.report: opt is not null\n')
  #  stop()
  #}
  Sweave(target)
  m <- create.makefile(target)
  cat('generate.report: target =', target, '\n')
  cat('generate.report: m =', m, '\n')
  #stop()
  system(sprintf('make -f %s pdf >& %s.log ', m, m))
  #system(sprintf('make -f %s pdf', m))
}

# m <- create.makefile(target)
create.makefile <- function(target) {
  x <- list()
  x[[1]] <- sprintf('TARGET=%s\n', target)
  x[[2]] <- 'export BSTINPUTS = /user/kzth541/bib/styles/achemso\n'
  x[[3]] <- 'export BIBINPUTS = /user/kzth541/bib\n'
  #x[[3]] <- 'export BIBINPUTS = $(HOME)/bib\n'
  x[[4]] <- 'export TEXINPUTS = ./:$(BSTINPUTS):\n'
  x[[5]] <- 'dvi:\n'
  x[[6]] <- '\tlatex $(TARGET).tex; bibtex $(TARGET); latex $(TARGET).tex; latex $(TARGET).tex\n'
  x[[7]] <- 'ps:\tdvi\n'
  x[[8]] <- '\tdvips -f < $(TARGET).dvi > $(TARGET).ps\n'
  x[[9]] <- 'pdf:\tps\n'
  x[[10]] <- '\tps2pdf $(TARGET).ps $(TARGET).pdf\n'
  x[[11]] <- '\trm -f $(TARGET) $(TARGET).{aux,bbl,blg,dvi,ps,tex} $(TARGET)*.eps Makefile.$(TARGET)\n'
  fn <- sprintf('Makefile.%s', target)
  cat(x[[1]], file=fn)
  for (i in 2:length(x)) cat(x[[i]], file=fn, append=T)
  return (fn)
}

time.stamp <- function() {
  x <- Sys.time()
  date <- substr(x, 1, 10)
  time <- substr(x, 12, 19)
  substr(time, 3, 3) <- '.'
  substr(time, 6, 6) <- '.'
  stamp <- sprintf('%s.%s', date, time)
  #return (stamp)
  return (date)
}

add.fixed <- function(model, x) {
  #
  # recombine optimized parms with fixed from model$parms
  # sort order in the same as model$parms
  #
  nr <- nrow(x)
  m <- setdiff(1:length(model$parms), model$idx)
  df <- t(as.data.frame(rep(as.data.frame(model$parms[m]), nr)))
  df <- as.data.frame(df)
  names(df) <- names(model$parms[m])
  row.names(df) <- 1:nr
  df <- cbind(x, df)
  m <- match(names(model$parms), names(df))
  df <- cbind(RMSE=x$RMSE, df[,m])
  return (df)
}

print.options <- function(opts) {
  for (x in names(opts)) {
    cat('option: ', x, '=', as.character(opts[[x]]), '\n')
  }
}

# MAIN

args <- commandArgs()
#cat('class(args), length(args) =', class(args), length(args), '\n')
#cat('args =', args, '\n')
if (length(args) < 2) stop() # this should allow to source() all code above

models <- list(Tam2007=list(name='PA.Levofloxacin', url='http://dx.doi.org/10.1007/s10439-007-9306-x'),
               Tam2007=list(name='PA.Meropenem', url='http://dx.doi.org/10.1007/s10439-007-9306-x'))

model.names <- c("PA.Levofloxacin", 'VT2007', 'PA.Meropenem', 'VT2005', 'SIMPLE')

opts <- list()

opts$user <- NULL
opts$firstName <- NULL
opts$lastName <- NULL
opts$model.name <- model.name <- "VT2005"
#opts$parms <- parms <- 'NONE' # starting values for parameter optimization
#opts$fixed <- fixed <- 'NONE' # parameters to fix
opts$parms <- parms <- NULL # starting values for parameter optimization
opts$fixed <- fixed <- NULL # parameters to fix
opts$lower <- lower <- NULL # parameters lower bounds
opts$upper <- upper <- NULL # parameters upper bounds
opts$method <- method <- 'BFGS'
opts$minimizer <- minimizer <- 'optim'
opts$maxit <- maxit <- 1e3
opts$prefix <- prefix <- 'prefix'
opts$solver <- solver <- 'bernoulli' # 'ode'
opts$dt <- dt <- 2^(-6)
opts$StDev <- StDev <- list(exp=0.5, opt=0.01)
opts$nrep <- nrep <- 2; # list(exp=1, opt=1)
opts$verbose <- verbose <- FALSE
opts$debug <- debug <- FALSE
opts$device <- device <- 'png'
opts$report.flag <- report.flag <- 0
opts$report.template <- report.template <- '$HOME/apps/tks/R/tks.Rnw'
opts$cut.RMSE <- cut.RMSE <- NULL
opts$delta.RMSE <- delta.RMSE <- 0.2
opts$seed <- seed <- as.integer(runif(1)*1e6)
opts$limit.factor <- 2

usage <- function() {
  cat('Usage: Rscript tks.R -input infile -output outfile [options] 
                     -help      : prints this usage
                     -input     : xml file []
                     -debug     : turn on debugging output
                     -output    : TAB separated data frame with optimized parameters []
                     -solver    : ODE solver [bernoulli|ode]
                     -model     : analytical model name [SIMPLE|VT2005|VT2007]
                     -parms     : starting values for the parameters [e.g. "c(N0=7461000,Nmax=4.05e+08, Kg=0.5932,Kk=13.97,C50k=1.192,H=4.207,beta=37.46,tau=0.02189)"]
                     -fixed     : parameters to fix [e.g. "c(\'N0\', \'Nmax\', \'Kg\')"]
                     -lower     : parameters lower bound [e.g. "c(N0=1e6, Nmax=1e8, Kg=0.1)"]
                     -upper     : parameters upper bound [e.g. "c(N0=1e7, Nmax=1e9, Kg=1.0)"]
                     -prefix    : prefix used to tag generated output files [prefix]
                     -optimized : file with historical optimized parameters [$TKS_HOME/data/models/MODEL-NAME.opt]
                     -dt        : time step used for "bernoulli" solver [2^(-6)]
                     -maxit     : maximum number of optimization steps [1e3]
                     -nrep      : number of restarts perturbing experimental data and params with gaussian noise [1]
                     -sde       : StDev to use for gaussian noise sample to perturb experimental data [0.5]
                     -sdo       : StDev to use for gaussian noise sample to perturb initial values for optimized parameters [0.01]
                     -device    : plot file format [png|eps]
                     -report    : (0|1|2) 0: - no report, 1: report only, 2: opt and report uses PREFIX.txt file)
                     -cutRMSE   : factor to use par values with RMSE < cutRMSE [use all values]
                     -template  : template for Sweave [$HOME/apps/tks/$version/tks.Rnw]
                     -seed      : for random number generator
                     -user      : user name (default: system $USER)
                     -firstName : first name
                     -lastName  : last name
                     -limit     : limiting multiplicative factor used for lower/upper parameters limit on the log scale [2]
                     -verbose   : turn on verbose output 
     \n')
  stop('STOP')
}

for (a in args) {
  if (a == '-help')      usage()
  if (a == '-input')     opts$input <- input <- as.character(args[which(args==a)+1])
  if (a == '-output')    opts$output <- output <- as.character(args[which(args==a)+1])
  if (a == '-solver')    opts$solver <- ode.solver <- as.character(args[which(args==a)+1])
  if (a == '-model')     opts$model.name <- model.name <- as.character(args[which(args==a)+1])
  if (a == '-parms')     parms <- as.character(args[which(args==a)+1])
  if (a == '-lower')     lower <- as.character(args[which(args==a)+1])
  if (a == '-upper')     upper <- as.character(args[which(args==a)+1])
  if (a == '-fixed')     fixed <- as.character(args[which(args==a)+1])
  if (a == '-prefix')    opts$prefix <- prefix <- as.character(args[which(args==a)+1])
  if (a == '-device')    opts$device <- device <- as.character(args[which(args==a)+1])
  if (a == '-template')  opts$report.template <- report.template <- as.character(args[which(args==a)+1])
  if (a == '-optimized') opts$optimized.parms <- optimized.parms <- as.character(args[which(args==a)+1])
  if (a == '-dt')        opts$dt <- dt <- as.numeric(args[which(args==a)+1])
  if (a == '-maxit')     opts$maxit <- maxit <- as.numeric(args[which(args==a)+1])
  if (a == '-nrep')      opts$nrep <- nrep <- as.numeric(args[which(args==a)+1])
  if (a == '-sde')       StDev$exp <- as.numeric(args[which(args==a)+1])
  if (a == '-sdo')       StDev$opt <- as.numeric(args[which(args==a)+1])
  if (a == '-report')    opts$report.flag <- report.flag <- as.numeric(args[which(args==a)+1])
  if (a == '-cutRMSE')   opts$cut.RMSE <- cut.RMSE <- as.numeric(args[which(args==a)+1])
  if (a == '-deltaRMSE') opts$delta.RMSE <- delta.RMSE <- as.numeric(args[which(args==a)+1])
  if (a == '-seed')      opts$seed <- seed <- as.numeric(args[which(args==a)+1])
  if (a == '-user')      opts$user <- as.character(args[which(args==a)+1])
  if (a == '-firstName') opts$firstName <- as.character(args[which(args==a)+1])
  if (a == '-lastName')  opts$lastName <- as.character(args[which(args==a)+1])
  if (a == '-limit')     opts$limit.factor <- as.numeric(args[which(args==a)+1])
  if (a == '-debug')     opts$debug <- debug <- TRUE
  if (a == '-verbose')   opts$verbose <- verbose <- TRUE
}

if (is.null(opts$user)) opts$user <- paste(opts$firstName, opts$lastName)
if (length(opts$user) == 0) opts$user <- NULL

options(warn=2)
cat('length(opts$user), opts$user =', length(opts$user), opts$user, '\n')
#stop('AAA')
cat('parms =', parms, '\n')
cat('fixed =', fixed, '\n')

parse.parms <- function(x) {
  if (!is.null(x)) {
    e <- parse(text=x)
    v <- eval(e)
  } else {
    v <- NULL
  }
  return (v)
}

opts$parms <- parse.parms(parms)
opts$lower <- parse.parms(lower)
opts$upper <- parse.parms(upper)
#cat ('1: opts$optimized ="', opts$optimized, '"\n')
#cat ('1: opts$optimized.parms ="', opts$optimized.parms, '"\n')
if (opts$optimized.parms == 'NONE') opts$optimized.parms <- NULL
#cat ('2: opts$optimized.parms ="', opts$optimized.parms, '"\n')
#stop()

#print.options(opts)
#stop('lll')

cat('before print: parms =', parms, '\n')
#print(parms)
#stop()
opts$StDev <- StDev
set.seed(seed)

#cat('StDev$exp =', StDev$exp, '\n')
#stop()

if (model.name == 'PA.Meropenem') model.name <- 'VT2005'
if (model.name == 'PA.Levofloxacin') model.name <- 'VT2007'

i <- 0; errs <- list()
if (!is.element(model.name, model.names)) {i <- i + 1; errs[[i]] <- sprintf('ERROR: invalid model.name= %s\n', model.name)}

n <- length(errs)
if (n > 0) {
  for (i in 1:n) cat(errs[[i]])
  stop()
}

if (opts$solver == 'ode') library(deSolve)

if (1) {
  if (1) print.options(opts)
  #stop()
  t0 <- proc.time()
  .tiny <- get.tiny()
  .HUGE <- get.huge()
  .precision <- get.precision()
  model <- get.xml.data(input, opts)
  
  # TBD : for some strange reasons this gives 
  # Error in diffinv.vector(x, lag, differences, xi) :
  # NA/NaN/Inf in foreign function call (arg 1)
  # Calls: par.sens ... i.diffinv -> diffinv -> diffinv.default -> diffinv.vector -> .C
  # setup rkill appropriate for the model, e.g. rkill <- rkill.VT2005
  e <- parse(text=sprintf('rkill <- rkill.%s', model$name)); eval(e)

  #t.begin <- min(model$exp.data$times); t.end <- max(model$exp.data$times)
  #n <- (t.end-t.begin)/dt
  #model$times <- t.begin + (t.end-t.begin)*(0:n)/n
  # initial guess for parameters
  if (class(opts$parms) == 'numeric') {
    if (!prod( names(model$parms) == names(opts$parms))) {
      cat('ERROR: names(model$parms) != names(opts$parms)\n')
      cat('ERROR: names(model$parms) =', names(model$parms), '\n')
      cat('ERROR: names(opts$parms) =', names(opts$parms), '\n')
      stop('wrong initial parmas')
    }
    if (length(match(opts$fixed, names(opts$parms))) < length(opts$fixed)) {
        cat('ERROR: some parameters selected to fix are not in parameter list\n')
        cat('ERROR: names(opts$parms) =', names(opts$parms), '\n')
        cat('ERROR: opts$fixed =', opts$ix, '\n')
        stop('some parameters selected to fix are not in parameter list')
    }
    model$parms <- opts$parms
    model$fixed <- opts$fixed
  } else {
    model$parms <- model$parms*0 + 1
    model$parms[1] <- 10^mean(model$exp.data$log10CFU[,1], na.rm=TRUE)
    model$parms[2] <- 10^max(model$exp.data$log10CFU, na.rm=TRUE)
  }
  #verbose <- 1
  cat ('report.flag =', report.flag, '\n')
  #cat ('model$name =', model$name, '\n')
  if (report.flag==1) {
    generate.report(model, opts, data=NULL, parms.file=opts$output, template=report.template, prefix=prefix)
    stop('report only')
  }
  opt <- par.sens(model, opts, sd=StDev, nrep=nrep, maxit=maxit, minimizer=minimizer, method=method, device=device, verbose=verbose)
  df <- add.fixed(model, opt$df.1)
  opt$x <- x <- df
  write.table(signif(x, digits=4), file=output, quote=F, sep='\t', row.names=F)
  #cat('END.\n')
  if (report.flag==2) {
    generate.report(model, opts, data=df, parms.file=NULL, template=report.template, prefix=prefix)
  }
  t1 <- proc.time()
  cat('time: user system elapsed =', (t1-t0)[1:3], '\n')
}
