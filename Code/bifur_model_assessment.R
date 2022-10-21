########################################################################################
## Preamble ##
########################################################################################
bifur.mod.scripts <- list.files(path ="Code/bifur_models", pattern="*.R",full.names = T) 
purrr::walk(bifur.mod.scripts, source) # source silently

load("Data/wrangled_genus_plank_data.Rdata")

########################################################################################
## Perform Assessments ##
########################################################################################
kin_density <- cbind(seq_along(as.numeric(kin_mth_dat$Date[kin_mth_dat$Date >1980])),rowSums(kin_mth_dat[kin_mth_dat$Date >1980,2:50]))

kinOU <- fit_bif_mod(kin_density,model = "OU",method = "L-BFGS-B",lower = 0)
kinSN <- fit_bif_mod(kin_density,model = "LSN",method = "L-BFGS-B",lower = 0)
kinTC <- fit_bif_mod(kin_density,model = "LTC",method = "L-BFGS-B",lower = 0)
kinPF <- fit_bif_mod(kin_density,model = "LPF",method = "L-BFGS-B",lower = 0)

observedSN <- -2*(kinOU$loglik - kinSN$loglik)
observedTC<- -2*(kinOU$loglik - kinTC$loglik)
observedPF <- -2*(kinOU$loglik - kinPF$loglik)
knitr::kable(cbind(observedSN,observedTC,observedPF))

ui <- matrix(c(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)),ncol=4,byrow = T)
ci <- c(0,0,0,0)
ui %*% p - ci

tt <- compare(A,B)
reps <- pbmcapply::pbmclapply(1:50, function(i) compare(A, B),mc.cores=8)
reps <- pbmcapply::pbmclapply(1:50, function(i) compare(A, C),mc.cores=8)

lr <- lik_ratios(reps)

ggplot(lr) + 
  geom_density(aes(value, fill=simulation), alpha=0.4) + 
  geom_vline(aes(xintercept=observedc)) + theme_bw()

compare <- function(A, B){
  done <- 0
  while(!done){
    Asim <- simulate(A)
    AfitA <- update(A, X = Asim)
    BfitA <- update(B, X = Asim)
    done <- AfitA$convergence && BfitA$convergence
  }
  done <- 0
  while(!done){
    Bsim <- simulate(B)
    BfitB <- update(B, X = Bsim)
    AfitB <- update(A, X = Bsim)
    done <- AfitB$convergence && BfitB$convergence
  }
  list(AA = AfitA, BB = BfitB, AB = AfitB, BA = BfitA)
}





f_closure <- function(X, setmodel){
  function(p){
    n <- length(X[,1])
    out <- -sum(dc.gauss(setmodel, x = X[2:n,2], x0 = X[1:(n-1),2], to=X[1:(n-1),1],
                         t1=X[2:n,1], p, log=T))
    if(abs(out) == Inf | is.nan(out))
      out <- 1e19
    out
  }
}

dc.gauss  <- function(setmodel, x, x0, to, t1, pars, log = FALSE){
  P <- setmodel(x0, to, t1, pars)
  dnorm(x, mean=P$Ex, sd=sqrt(P$Vx), log=log)
}

constOU <- function(Xo, to, t1, pars){
  Dt <- t1 - to
  Ex <- pars["theta"]*(1 - exp(-pars["Ro"] * Dt)) + Xo *
    exp(-pars["Ro"] * Dt) 
  Vx <- 0.5 * pars["sigma"] ^ 2 * 
    (1 - exp(-2 * pars["Ro"] * Dt)) / pars["Ro"]
  if(pars['Ro'] < 0 ) Vx <- rep(Inf, length(Xo)) 
  if(pars['sigma'] < 0 ) Vx <- rep(Inf, length(Xo)) 
  return(list(Ex = Ex, Vx = Vx))
}



LSN <- function(Xo, to, t1, pars){
  
  R <- function(t, pars){pars[1] + pars[2]*t }
  
  check <- any(R(t1,pars) < 0)
  if(is.na(check) | check | pars[3] < 0){
    Ex <- Xo
    Vx <- rep(Inf,length(Xo))
  } else {
    
    
    moments <- function(t,y,p){ 
      sqrtR <- sqrt(R(t,pars)) 
      yd1 <- 2*sqrtR*(sqrtR+pars[3] - y[1]) 
      yd2 <- -2*sqrtR*y[2] + p[4]^2*(sqrtR+pars[3])
      list(c(yd1=yd1, yd2=yd2))
    }
    jacfn <- function(t,y,p){
      sqrtR <- sqrt(R(t,pars)) 
      c(
        -2*sqrtR, 0,
        0, -2*sqrtR
      )}
    ## The apply calls needed to work with vector inputs as Xo (whole timeseries)
    times <- matrix(c(to, t1), nrow=length(to))
    out <- lapply(1:length(Xo), function(i){
      deSolve::lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, 
                     parms=pars, jacfunc=jacfn) 
    })
    Ex <- sapply(1:length(Xo), function(i) out[[i]][2,2]) # times are in rows, cols are time, par1, par2
    Vx <- sapply(1:length(Xo), function(i) out[[i]][2,3])
  }
  # Handle badly defined parameters by creating very low probability returns,
  # needed particularly for the L-BFGS-B bounded method, which can't handle Infs 
  # and does rather poorly... Worth investigating a better way to handle this.  
  # Note errors can apply to some of the timeseries, possibly by estimates of m 
  # that force system into bifurcation on later parameters (for which terms  
  # become undefined)
  
  if (pars[4] <= 0){
    warning(paste("sigma=",pars[4]))
    Vx <- rep(Inf, length(Xo))
  }
  if (any(Vx < 0, na.rm=TRUE)){
    warning(paste("discard negative Vx,  "))
    Vx[Vx<0] <- Inf
  }
  return(list(Ex=Ex, Vx=Vx))
}


stability_model <- function(X, model=c("LSN", "OU"), p = NULL, ..., 
                            store_data=TRUE){
  model <- match.arg(model)
  # reformat time series objects into proper data frames
  if(is(X, "ts"))
    X <- data.frame(as.numeric(time(X)), X@.Data)
  # if time values are not provided
  else if(is.null(dim(X)))
    X <- data.frame(1:length(X), X)
  
  
  if(!is.null(p)){ ## got everything? then rock & roll
    f1 <- switch(model, 
                 LSN = f_closure(X, LSN),
                 OU = f_closure(X, constOU))
    o <- optim(p, f1, ...)
    
  } else if(is.null(p)){ ## oh, need p? try:
    p <- c(Ro=1/max(time(X[,1])), theta=mean(X[,2]), sigma=sd(X[,2]))
    f2 <- f_closure(X, constOU)
    o <- optim(p, f2, ...)
    
    # if model is "OU", we're done.  otherwise:
    if(model=="LSN"){
      f3 <- f_closure(X, LSN) # switch to the LSN model
      p_est <- o$par  # & use the OU estimated pars as starting guess
      # but rescale them to the new definitions:
      Ro <- as.numeric(p_est[1]^2)
      theta <- as.numeric(p_est[2]+p_est[1])
      sigma <- as.numeric(abs(p_est[3]/sqrt(2*p_est[1]+ p_est[2])))
      p <- c(Ro=Ro, m=0, theta=theta, sigma=sigma)
      
      ## now fit the LSN model
      o <- optim(p, f3, ...)
    }
  }
  names(X) <- c("time", "value")
  ## Collect the results and we're done
  if(!store_data) # remove the data object to save space?
    X <- NULL
  out <- list(X=X, pars=o$par, model=model, loglik = -o$value,
              convergence=(o$convergence==0) )
  class(out) <- c("gauss", "list")
  out
}

rc.gauss <- function(setmodel, n=1, x0, to, t1, pars){
  P <- setmodel(x0, to, t1, pars)
  rnorm(n, mean=P$Ex, sd=sqrt(P$Vx))
}

update.gauss <- function(object, ...){
  stability_model(..., model = object$model, p = object$pars)
}

simulate.gauss <- function(object, nsim = 1, seed = NULL, ...){
  if(object$model == "LSN"){
    setmodel <- LSN
  }else if(object$model == "LTC"){ 
    setmodel <- LTC
  }else if(object$model == "LPF"){ 
    setmodel <- LPF
  }else if(object$model == "OU"){ 
    setmodel <- constOU
  }
  time <- object$X[,1]
  N <- length(time)
  value <- sapply(1:nsim, function(j){
    X <- numeric(N - 1)
    X[1] <- object$X[1,2]
    for(i in 1:(N - 1) ){
      X[i + 1] <- rc.gauss(setmodel, 1, x0 = X[i], to = time[i], 
                           t1 = time[i + 1], object$pars)
    }
    X
  })
  data.frame(time, value)
}

lik_ratios <- function(reps){
  dat <- sapply(reps, function(rep){
    null <- -2*(rep$AA$loglik - rep$BA$loglik)
    test <- -2*(rep$AB$loglik - rep$BB$loglik)
    c(null=null, test=test)
  })
  dat <- reshape2::melt(dat)
  names(dat) <- c("simulation", "rep", "value") 
  dat
}
