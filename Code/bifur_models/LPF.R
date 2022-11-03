###########  Linearized Pitchfork (LTC) Model ##################
#' Linearized supercritical Pitchfork model
#'
#' Estimate the mean and variance of the linearized supercritical 
#' pitchfork process (in which the control constant 
#' changes over time)
#' @param Xo initial condition
#' @param to initial time
#' @param t1 final time
#' @param pars numeric of parameters named Ro, theta, and sigma
#' @return a list with values Ex, the expected x value, and 
#' Vx, the expected variance
#' @keywords internal
superLPF <- function(Xo, to, t1, pars){
  
  R <- function(t, p){p[1] + p[2]*t }
  
  check <- any(R(t1,pars) < 0)
  if(is.na(check) | check | pars[3] < 0){
    Ex <- Xo
    Vx <- rep(Inf,length(Xo))
  } else {
    moments <- function(t,y,p){ 
      r_t <- R(t,p)
      sqrtR <- sqrt(r_t) 
      yd1 <- 2*r_t*(y[1]-p[3]-sqrtR)
      yd2 <- -2*r_t*r_t*(y[1] - p[3]-sqrtR)*y[2] + p[4]^2*(p[3] + sqrtR)
      list(c(yd1=yd1, yd2=yd2))
    }
    jacfn <- function(t,y,p){
      r_t <- R(t,p)
      c(
        2*r_t, 0,
        0,2*r_t
      )}
    ## The apply calls needed to work with vector inputs as Xo (whole timeseries)
    times <- matrix(c(to, t1), nrow=length(to))
    out <- lapply(1:length(Xo), function(i){
      deSolve::lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, 
            parms=pars,jacfunc = jacfn)
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


LPF <- function(Xo, to, t1, pars){
  
  R <- function(t, p){p[1] + p[2]*t }
  
  check <- any(R(t1,pars) < 0)
  if(is.na(check) | check | pars[3] < 0){
    Ex <- Xo
    Vx <- rep(Inf,length(Xo))
  } else {
    moments <- function(t,y,p){ 
      r_t <- R(t,p)
      sqrtR <- sqrt(r_t) 
      yd1 <- 2*r_t*(y[1]-p[3]-sqrtR)
      yd2 <- -2*r_t*r_t*(y[1] - p[3]-sqrtR)*y[2] + p[4]^2*(p[3] + sqrtR)
      list(c(yd1=yd1, yd2=yd2))
    }
    jacfn <- function(t,y,p){
      r_t <- R(t,p)
      c(
        2*r_t, 0,
        0,2*r_t
      )}
    ## The apply calls needed to work with vector inputs as Xo (whole timeseries)
    times <- matrix(c(to, t1), nrow=length(to))
    out <- lapply(1:length(Xo), function(i){
      tryCatch( {R.utils::withTimeout(deSolve::lsoda(y=c(xhat=Xo[i], sigma2=0), times=times[i,], func=moments, 
                                                     parms=pars,jacfunc = jacfn)
                                      ,timeout = 0.1,onTimeout = "error")},
                TimeoutException = function(ex) matrix(rep(NA,6),nrow = 2))
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