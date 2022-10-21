fit_bif_mod <- function(X, model=c("LSN","LTC","LPF","OU"), p = NULL, ..., 
                            store_data=TRUE){
  model <- match.arg(model,choices = c("LSN","LTC","LPF","OU"))
  # reformat time series objects into proper data frames
  if(is(X, "ts")){
    X <- data.frame(as.numeric(time(X)), X@.Data)
  }
  # if time values are not provided
  else if(is.null(dim(X))){
    X <- data.frame(1:length(X), X)
  }
  
  if(!is.null(p)){ ## got everything? then rock & roll
    f1 <- switch(model, 
                 LSN = f_closure(X, LSN),
                 LTC = f_closure(X, LTC),
                 LPF = f_closure(X, LPF),
                 OU = f_closure(X, constOU))
    o <- optim(p, f1, ...)
    
  } else if(is.null(p)){ ## oh, need p? try:
    p <- c(Ro=1/max(time(X[,1])), theta=mean(X[,2]), sigma=sd(X[,2]))
    f2 <- f_closure(X, constOU)
    o <- optim(p, f2, ...)
    
    # if model is "OU", we're done.  otherwise:
    if(model != "OU"){
      
    if(model=="LSN"){
      f3 <- f_closure(X, LSN) # switch to the LSN model
    }else if(model=="LTC"){
        f3 <- f_closure(X, LTC) # switch to the LTC model
      }else if(model=="LPF"){
        f3 <- f_closure(X, LPF) # switch to the LPF model
      }
      p_est <- o$par  # & use the OU estimated pars as starting guess
      # but rescale them to the new definitions:
      Ro <- as.numeric(p_est[1]^2)
      theta <- as.numeric(p_est[2]+p_est[1])
      sigma <- as.numeric(abs(p_est[3]/sqrt(2*p_est[1]+ p_est[2])))
      p <- c(Ro=Ro, m=0, theta=theta, sigma=sigma)
      
      ## now fit the bifurcating model
      o <- optim(p, f3, method = "L-BFGS-B",lower = c(NA,0,NA,0))
      o2 <- optim(p, f3,method = "L-BFGS-B",lower = c(NA,0,NA,0))
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


# an internal helper function
f_closure <- function(X, setmodel){
  function(p){
    n <- length(X[,1])
    out <- -sum(dc.gauss(setmodel, X[2:n,2], X[1:(n-1),2], to=X[1:(n-1),1],
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

rc.gauss <- function(setmodel, n=1, x0, to, t1, pars){
  P <- setmodel(x0, to, t1, pars)
  rnorm(n, mean=P$Ex, sd=sqrt(P$Vx))
}

update.gauss <- function(object, ...){
  fit_bif_mod(..., model = object$model, p = object$pars)
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

