#' Significance Testing of Rolling Window Early Warning Signals
#'
#' A function for identifying whether a warning has been generated from rolling early warning signal data using permutation tests.
#'
#' @param data A dataframe where the first column is an equally spaced time vector and the second column is the time series to be assessed.
#' @param metrics String vector of early warning signal metrics to be assessed. Options include: "ar1", "cv", "SD", "acf", "rr", "dr", "skew", "kurt", and "trait".
#' @param winsize Numeric value. If method = "rolling", defines the window size of the rolling window as a percentage of the time series length.
#' @param perm_method String dictating the pseudo-randomisation technique to be used. Options include: "arima" (sampled from predictions of an ARIMA model), "red.noise" (red noise process using data mean, variance and autocorrelation coef) or "replacement" (sampled from observed data with replacement).
#' @param iter. Numeric value. The number of permutations.
#' @returns A dataframe modifying \code{ews.data} with whether a threshold has been crossed `$threshold.crossed` and if that matches the expectations of \code{outcome} (\code{$prediction}).
#'

perm_rollEWS<- function(data, metrics, winsize = 50, perm.meth = "arima", iter = 500, variate ="uni",...){
  
  red.noise.ts <- function(ts,lag,length){
    x  <- rep(NA, length)
    ts <- rnorm(length,mean = mean(ts,na.rm=T),sd=sd(ts,na.rm=T))
    x[1] <- ts[1]
    for(i in seq(from = 2, to = length, by =1)){
      x[i] <- lag*x[i-1] + ((1-lag^2)^0.5)*ts[i]
    }
    return(x) #create rednoise process from a time series' mean, variance and autocorrelation coef
  }
  
  data <- as.data.frame(data) 
  
  if(variate == "uni"){
    
  
  if(perm.meth == "arima"){
  fit <- forecast::auto.arima(stats::as.ts(data[,2]),allowdrift =F,seasonal = T)
  perm.df <- data

  for(perm in seq_len(iter)){ # for each permutation randomly sample from confidence interval at each time point
    new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
    perm.df[,new_col_name] <- arima.sim(n = length(data[,1]), as.list(coef(fit)), sd = sqrt(fit$sigma2))
  }
  
  }
  
  if(perm.meth == "red.noise"){
    
  d1.ar1 <- tryCatch({stats::arima(stats::as.ts(data[,2]), order = c(1, 0, 0),optim.control = list(maxit = 1000),method="ML")$coef[1]
  },error = function(err){return(0)})
  
  perm.df <- data
  
  for (perm in seq_len(iter)) {
    new_col_name <- paste0("perm_", perm)
    perm.df[,new_col_name] <- red.noise.ts(data[,2],lag=d1.ar1,length=length(data[,1])) #permuted autocorrelated surrogates
  }  
  }
  
  if(perm.meth == "replacement"){
    ## Create permutation order
    perm.df <- data
    
    for (perm in seq_len(iter)){ # for each permutation randomly sample from ts with replacement at each time point
      new_col_name <- paste0("perm_", perm) # new column and colname for each permutation
      perm.df[,new_col_name] <- sample(x = data[,2],replace = F,size = length(data[,1]))
      while( any(rle(perm.df[,new_col_name])$lengths >= floor(length(perm.df[,new_col_name])*(winsize/100)))){ #if ts contains many zeroes or runs to prevent svd, prevent the sampling with replacement 
      #while(mean( perm.df[,new_col_name] == min( perm.df[,new_col_name])) > 0.47 ){
          
        perm.df[,new_col_name] <- sample(x = data[,2],replace = F,size = length(data[,1]))
      }
      }
  
  }
  
  true.ews <- EWSmethods::uniEWS(perm.df[,c(1,2)],metrics =metrics,method = "rolling",ggplotIt = F,
                                               winsize = winsize)
  
  # perm.ews <- pbapply::pblapply(3:dim(perm.df)[2],FUN = function(i){
  #   ews.roll <- EWSmethods::uniEWS(perm.df[,c(1,i)],metrics =metrics,method = "rolling",ggplotIt = F,
  #                             winsize = winsize)
  #   return(ews.roll$cor)
  # },...)
  
  perm.ews <- pbmcapply::pbmclapply(3:dim(perm.df)[2],FUN = function(i){
    ews.roll <- tryCatch({EWSmethods::uniEWS(perm.df[,c(1,i)],metrics =metrics,method = "rolling",ggplotIt = F,
                                   winsize = winsize)},error = function(err){
                                     return(list(cor = -Inf))
                                   })
    return(ews.roll$cor)
  },...)
  
  }else if(variate == "multi"){

    if(perm.meth == "arima"){
      
      perm.ls <- lapply(1:iter,FUN = function(perm){
      
      sapply(2:dim(data)[2], function(i){
        fit <- forecast::auto.arima(stats::as.ts(data[,i]),allowdrift =F,seasonal = T)
        arima.sim(n = length(data[,1]), as.list(coef(fit)), sd = sqrt(fit$sigma2))
      })
      })
    }
    
    if(perm.meth == "red.noise"){
      
      perm.ls <- lapply(1:iter,FUN = function(perm){
        
        sapply(2:dim(data)[2], function(i){
          d1.ar1 <- tryCatch({stats::arima(stats::as.ts(data[,i]), order = c(1, 0, 0),optim.control = list(maxit = 1000),method="ML")$coef[1]
          },error = function(err){return(0)})
          red.noise.ts(ts = data[,i],lag=d1.ar1,length=length(data[,1]))
        })
      })
      
    }
    
    if(perm.meth == "replacement"){

      perm.ls <- lapply(1:iter,FUN = function(perm){
        sapply(2:dim(data)[2], function(i){
          s_out <- sample(x = data[,i],replace = F,size = length(data[,1]))
          while(mean(s_out == min(s_out)) > 0.47 | any(rle(s_out)$lengths >= floor(length(s_out)*0.5))){ #if ts contains many zeroes or runs to prevent svd, prevent the sampling with replacement 
            s_out <- sample(x = data[,i],replace = F,size = length(data[,1]))
          }
          return(s_out)
        })
      })
      
    }
    
    true.ews <- EWSmethods::multiEWS(data,metrics =metrics,method = "rolling",ggplotIt = F,
                                   winsize = winsize)
    
    # perm.ews <-  pbapply::pblapply(1:length(perm.ls),FUN = function(i){
    #   ews.roll <- EWSmethods::multiEWS(cbind(data[,1],perm.ls[[i]]),metrics =metrics,method = "rolling",ggplotIt = F,
    #                                  winsize = winsize)
    #   return(ews.roll$cor)
    # },...)
    
    # perm.ews <-  pbmcapply::pbmclapply(1:length(perm.ls),FUN = function(i){
    #   ews.roll <- EWSmethods::multiEWS(cbind(data[,1],perm.ls[[i]]),metrics =metrics,method = "rolling",ggplotIt = F,
    #                                    winsize = winsize)
    #   return(ews.roll$cor)
    # },...)
    
    perm.ews <-  pbmcapply::pbmclapply(1:length(perm.ls),FUN = function(i){
      #print(i)
      ews.roll <- tryCatch({EWSmethods::multiEWS(cbind(data[,1],perm.ls[[i]]),metrics =metrics,method = "rolling",ggplotIt = F,
                                       winsize = winsize)},error = function(err){
                                         return(list(cor = -Inf))
                                       })
      return(ews.roll$cor)
    },...)
    
  }
  perm.ews <- do.call("rbind",perm.ews)
  
  out.ews <- sapply(seq_len(length(perm.ews)), FUN = function(x){
    ifelse(!all(is.na(perm.ews[,x])),ecdf(perm.ews[,x])(true.ews$cor[x]),NA)
  })
  
  true.ews$cor["perm_pvalue",] <- out.ews
  
  return(true.ews)
  
}
