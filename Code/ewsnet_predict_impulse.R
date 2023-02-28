
ewsnet_predict_impulse <- function(x, scaling = TRUE, ensemble = 25, envname){
  
  if(!envname %in% (reticulate::conda_list()$name)){
    warning("Call 'ewsnet_init()' before attempting to use ewsnet_predict(), or check your spelling of envname")
  }else{
    
    if(!is.vector(x)){
      stop('Time series is not a vector')
    }
    
    if(!is.numeric(x) ){
      stop('Time series is not numeric')
    }
    
    wd <- getwd() #get working directory so it can be reset when Python alters the directory
    EWSNet <- NULL # global variable to be populated by Python
    
    ensemble <- match.arg(paste(ensemble), choices = paste(1:25))
    

    if(isTRUE(scaling)){
      scaling_string <- paste("Scaled")
    }else if(isFALSE(scaling)){
      scaling_string <- paste("Unscaled")
    }
    
    directory_string = paste(c("directory_string = '", wd,"'"),collapse = "")
    
    reticulate::py_run_string(directory_string)
    reticulate::py_run_string("import os")
    reticulate::py_run_string("os.chdir(directory_string)")
  
    reticulate::source_python("python/src/inference/ewsNET_generic_testing.py")

    ewsnet_obj <- EWSNet(ensemble = as.integer(25), weight_dir = paste(c(directory_string,"python/weights/Pretrained",scaling_string),collapse = "/"), prefix = "", suffix = ".h5")
    
    if(isTRUE(scaling)){
      pred <- ewsnet_obj$predict(data_scaling(x),ensemble_subset = as.integer(ensemble))
    }else if(isFALSE(scaling)){
      pred <- ewsnet_obj$predict(x,ensemble_subset = as.integer(ensemble))
    }
    
    out <- data.frame("pred" = pred[[1]],
                      "no_trans_prob" = pred[[2]]$`No Transition`,
                      "smooth_trans_prob" = pred[[2]]$`Smooth Transition`,
                      "critical_trans_prob" = pred[[2]]$`Critical Transition`)
    
    setwd(wd) # reset working directory
   
    return(out)
  }
  
}

data_scaling <- function(x){
  
  x_min <- min(x,na.rm = T)
  s <- (x-x_min)/(max(x,na.rm = T)-x_min)
  
  return(s+1)
  
}