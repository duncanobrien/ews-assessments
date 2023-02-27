require(dplyr)
require(EWSmethods)
require(foreach)
load("Data/wrangled_genus_plank_data.Rdata")

lapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,
            wind_yr_dat,wash_yr_dat,leve_yr_dat,UZ_yr_dat,mon_yr_dat,
            kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,
            wind_mth_dat,wash_mth_dat,leve_mth_dat,UZ_mth_dat,mon_mth_dat),
       function(x){
         
         if(any(stringr::str_length(x[,1])>4)){
           x <- x |>
             select(-c(pca1,pca2))|>
             select(where(function(x){all(rle(x)$lengths[rle(x)$values == 0] < 12)}))
         }else if(all(stringr::str_length(x[,1])==4)){
           x <- x |>
             select(-c(pca1,pca2))|>
             select(where(function(x){all(rle(x)$lengths[rle(x)$values == 0] < 2)}))
         }
         return(x)
       }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat")) |>
  list2env(env = .GlobalEnv) #replace global env objects with wrangled

rm(list = ls()[!ls() %in% c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
                            "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat")])


################################################################################################################
## Multivariate EWS Assessment Expanding ##
################################################################################################################

exp_multi_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                          wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                          kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                          wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                     FUN = function(x){
                                       
                                       sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                         mean(x[,i] == min(x[,i])) <= 0.47
                                       })]
                                       
                                       out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                         if(j != "none"){
                                           sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                           sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                             sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                           })
                                         }
                                         
                                         out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                           
                                           if(k != "none"){
                                             
                                             if(any(stringr::str_length(sub_dat[,1])>4)){
                                               sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                               sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                               sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                               
                                               if(k == "decompose"){
                                                 sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                               }
                                             }
                                           }
                                           
                                           out.det <- EWSmethods::multiEWS(sub_dat,
                                                                           metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                           method = "expanding",ggplotIt = F,
                                                                           burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                           
                                           out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                             metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                         "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                             method = "expanding",ggplotIt = F,
                                                                             burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                           
                                           return(data.frame(rbind(out.det$raw,out.undet$raw)))
                                           
                                         } ) |>
                                           `names<-`(c("none","average","decompose","stl")) |>
                                           data.table::rbindlist(idcol="deseason_meth")
                                         
                                         return(out_inner)
                                       }) |>
                                         `names<-`(c("none","linear","gaussian","loess")) |>
                                         data.table::rbindlist(idcol="detrend_meth")
                                       
                                       return(out)
                                     }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "expanding",
    troph_level = "phytoplankton") |>
  select(-name)


exp_multi_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                        wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                        kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                        wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                   FUN = function(x){
                                     
                                     sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                       mean(x[,i] == min(x[,i])) <= 0.47
                                     })]
                                     
                                     out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                       if(j != "none"){
                                         sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                         sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                           sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                         })
                                       }
                                       
                                       out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                         
                                         if(k != "none"){
                                           
                                           if(any(stringr::str_length(sub_dat[,1])>4)){
                                             sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                             sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                             sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                             
                                             if(k == "decompose"){
                                               sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                             }
                                           }
                                         }
                                         
                                         out.det <- EWSmethods::multiEWS(sub_dat,
                                                                         metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                         method = "expanding",ggplotIt = F,
                                                                         burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                         
                                         out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                           metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                       "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                           method = "expanding",ggplotIt = F,
                                                                           burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                         
                                         return(data.frame(rbind(out.det$raw,out.undet$raw)))
                                         
                                       } ) |>
                                         `names<-`(c("none","average","decompose","stl")) |>
                                         data.table::rbindlist(idcol="deseason_meth")
                                       
                                       return(out_inner)
                                     }) |>
                                       `names<-`(c("none","linear","gaussian","loess")) |>
                                       data.table::rbindlist(idcol="detrend_meth")
                                     
                                     return(out)
                                     
                                   }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "expanding",
    troph_level = "zooplankton") |>
  select(-name)

exp_multicomp <- rbind(exp_multi_phyto,exp_multi_zoo)

################################################################################################################
## Multivariate EWS Assessment Rolling ##
################################################################################################################

roll_multi_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                           wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                           kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                           wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                      FUN = function(x){
                                        
                                        sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                          mean(x[,i] == min(x[,i])) <= 0.47
                                        })]
                                        
                                        out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                          if(j != "none"){
                                            sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                            sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                              sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                            })
                                          }
                                          
                                          out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                            
                                            if(k != "none"){
                                              
                                              if(any(stringr::str_length(sub_dat[,1])>4)){
                                                sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                                sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                                sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                                
                                                if(k == "decompose"){
                                                  sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                                }
                                              }
                                            }
                                            out.det <- EWSmethods::multiEWS(sub_dat,
                                                                            metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                            method = "rolling",ggplotIt = F,
                                                                            winsize = 50)
                                            
                                            out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                              metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                          "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                              method = "rolling",ggplotIt = F,
                                                                              winsize = 50)
                                            
                                            return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                            
                                          } ) |>
                                            `names<-`(c("none","average","decompose","stl")) |>
                                            data.table::rbindlist(idcol="deseason_meth")
                                          
                                          return(out_inner)
                                        }) |>
                                          `names<-`(c("none","linear","gaussian","loess")) |>
                                          data.table::rbindlist(idcol="detrend_meth")
                                        
                                        return(out)
                                      }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling",
    troph_level = "phytoplankton") |>
  select(-name)


roll_multi_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                         wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                         kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                         wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                    FUN = function(x){
                                      
                                      sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                        mean(x[,i] == min(x[,i])) <= 0.47
                                      })]
                                      
                                      out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                        if(j != "none"){
                                          sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                          sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                            sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                          })
                                        }
                                        
                                        out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                          
                                          if(k != "none"){
                                            
                                            if(any(stringr::str_length(sub_dat[,1])>4)){
                                              sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                              sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                              sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                              
                                              if(k == "decompose"){
                                                sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                              }
                                            }
                                          }
                                          
                                          out.det <- EWSmethods::multiEWS(sub_dat,
                                                                          metrics = c("meanAR", "maxAR", "meanSD", "maxSD"),
                                                                          method = "rolling",ggplotIt = F,
                                                                          winsize = 50)
                                          
                                          out.undet <- EWSmethods::multiEWS(sub_dat,
                                                                            metrics = c("eigenMAF", "mafAR", "mafSD",
                                                                                        "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                                            method = "rolling",ggplotIt = F,
                                                                            winsize = 50)
                                          
                                          return(data.frame(cbind(out.det$cor,out.undet$cor)))
                                          
                                        } ) |>
                                          `names<-`(c("none","average","decompose","stl")) |>
                                          data.table::rbindlist(idcol="deseason_meth")
                                        
                                        return(out_inner)
                                      }) |>
                                        `names<-`(c("none","linear","gaussian","loess")) |>
                                        data.table::rbindlist(idcol="detrend_meth")
                                      
                                      return(out)
                                      
                                    }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling",
    troph_level = "zooplankton") |>
  select(-name)

roll_multicomp <- rbind(roll_multi_phyto,roll_multi_zoo)

################################################################################################################
## Permuted Multivariate EWS Assessment Rolling ##
################################################################################################################
phyto_ls <- list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                 wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                 kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                 wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

perm_roll_multi_phyto <- lapply(seq_along(phyto_ls),
                                FUN = function(i){
                                  
                                  x <- phyto_ls[[i]]
                                  
                                  sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                    mean(x[,i] == min(x[,i])) <= 0.47
                                  })]
                                  
                                  out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                    
                                    if(j != "none"){
                                      sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                      sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                        sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                      })
                                    }
                                    
                                    out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                      
                                      if(k != "none"){
                                        
                                        if(any(stringr::str_length(sub_dat[,1])>4)){
                                          sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                          sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                          sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                          
                                          if(k == "decompose"){
                                            sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                          }
                                        }
                                      }
                                      print(paste(names(phyto_ls)[i],j,k,collapse = "_"))
                                      out.det <- perm_rollEWS(data = sub_dat,
                                                              metrics = c("meanAR", "maxAR", "meanSD", "maxSD",
                                                                          "eigenMAF", "mafAR", "mafSD",
                                                                          "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                              perm.meth = "replacement",iter=500,
                                                              winsize = 50,variate = "multi",mc.cores=5)
                                      
                                      # out.undet <- perm_rollEWS(sub_dat,metrics = c("eigenMAF", "mafAR", "mafSD",
                                      #                                               "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                      #                           perm.meth = "replacement",iter=500,
                                      #                           winsize = 50,variate = "multi"
                                      
                                      return(data.frame(out.det$cor))
                                      
                                    } ) |>
                                      `names<-`(c("none","average","decompose","stl")) |>
                                      data.table::rbindlist(idcol="deseason_meth")
                                    
                                    return(out_inner)
                                  }) |>
                                    `names<-`(c("none","linear","gaussian","loess")) |>
                                    data.table::rbindlist(idcol="detrend_meth")
                                  
                                  return(out)
                                }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling",
    troph_level = "phytoplankton") |>
  select(-name)

zoo_ls <- list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
               wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
               kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
               wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

perm_roll_multi_zoo <- lapply(seq_along(zoo_ls),
                              FUN = function(i){
                                
                                x <- zoo_ls[[i]]
                                
                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                  mean(x[,i] == min(x[,i])) <= 0.47
                                })]
                                
                                out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                  if(j != "none"){
                                    sub_dat <- EWSmethods::detrend_ts(sub_dat,method = j,span=0.5)   
                                    sub_dat[,-1] <- sapply(colnames(sub_dat)[-1],FUN = function(i){
                                      sub_dat[,i] <- sub_dat[,i] + abs(min(sub_dat[,i])) #set 0 as minimum
                                    })
                                  }
                                  
                                  out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                    
                                    if(k != "none"){
                                      
                                      if(any(stringr::str_length(sub_dat[,1])>4)){
                                        sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                        sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = k,order = "ymd")
                                        sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
                                        
                                        if(k == "decompose"){
                                          sub_dat <- na.omit(sub_dat) #remove NAs at head and tail of df
                                        }
                                      }
                                    }
                                    
                                    print(paste(names(zoo_ls)[i],j,k,collapse = "_"))
                                    out.det <- perm_rollEWS(data = sub_dat,
                                                            metrics = c("meanAR", "maxAR", "meanSD", "maxSD",
                                                                        "eigenMAF", "mafAR", "mafSD",
                                                                        "pcaAR", "pcaSD", "eigenCOV", "maxCOV", "mutINFO"),
                                                            perm.meth = "replacement",iter=500,
                                                            winsize = 50,variate = "multi",mc.cores=5)
                                    
                                    
                                    return(data.frame(out.det$cor))
                                    
                                  } ) |>
                                    `names<-`(c("none","average","decompose","stl")) |>
                                    data.table::rbindlist(idcol="deseason_meth")
                                  
                                  return(out_inner)
                                }) |>
                                  `names<-`(c("none","linear","gaussian","loess")) |>
                                  data.table::rbindlist(idcol="detrend_meth")
                                
                                return(out)
                                
                              }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling",
    troph_level = "zooplankton") |>
  select(-name)

load("Results/perm_ews_raw_data.RData")
perm_roll_multicomp <- rbind(perm_roll_multi_phyto,perm_roll_multi_zoo)

################################################################################################################
## Univariate EWS Assessment Expanding ##
################################################################################################################

exp_uni_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                        wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                        kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                        wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                   FUN = function(x){
                                     
                                     out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                       
                                       if(j != "none"){
                                         x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                         x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                           x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                         })
                                       }
                                       
                                       out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                         
                                         if(k != "none"){
                                           if(any(stringr::str_length(x[,1])>4)){
                                             x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                             x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                             x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                             
                                             if(k == "decompose"){
                                               x <- na.omit(x) #remove NAs at head and tail of df
                                             }
                                           }
                                         }
                                         
                                         time_vec <- colnames(x)[1]
                                         
                                         foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                           
                                           out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                     burn_in = ceiling(dim(x)[1]*0.50))
                                           
                                           data.frame(cbind(data_source = paste(i),
                                                            out))
                                         }
                                       }) |>
                                         `names<-`(c("none","average","decompose","stl")) |>
                                         data.table::rbindlist(idcol="deseason_meth")
                                       
                                       return(out_inner)
                                     }) |>
                                       `names<-`(c("none","linear","gaussian","loess")) |>
                                       data.table::rbindlist(idcol="detrend_meth")
                                     return(out)
                                     
                                   }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "expanding",
    troph_level = "phytoplankton") |>
  select(-name)

exp_uni_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                      wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                      kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                      wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                 FUN = function(x){
                                   
                                   
                                   out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                     
                                     if(j != "none"){
                                       x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                       x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                         x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                       })
                                     }
                                     
                                     out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                       
                                       if(k != "none"){
                                         if(any(stringr::str_length(x[,1])>4)){
                                           x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                           x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                           x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                           
                                           if(k == "decompose"){
                                             x <- na.omit(x) #remove NAs at head and tail of df
                                           }
                                         }
                                       }
                                       
                                       time_vec <- colnames(x)[1]
                                       
                                       foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                         
                                         out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                   burn_in = ceiling(dim(x)[1]*0.50))
                                         
                                         data.frame(cbind(data_source = paste(i),
                                                          out))
                                       }
                                     }) |>
                                       `names<-`(c("none","average","decompose","stl")) |>
                                       data.table::rbindlist(idcol="deseason_meth")
                                     
                                     return(out_inner)
                                   }) |>
                                     `names<-`(c("none","linear","gaussian","loess")) |>
                                     data.table::rbindlist(idcol="detrend_meth")
                                   return(out)
                                   
                                 }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = "Kinneret",
         res =  case_when(
           grepl("yr", name) ~ "Yearly",
           grepl("mth", name) ~ "Monthly"),
         method = "univariate EWS",
         computation = "expanding",
         troph_level = "zooplankton") |>
  select(-name)

exp_unicomp <- rbind(exp_uni_phyto,exp_uni_zoo)

################################################################################################################
## Univariate EWS Assessment Rolling ##
################################################################################################################

roll_uni_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                         wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                         kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                         wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                    FUN = function(x){
                                      
                                      out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                        
                                        if(j != "none"){
                                          x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                          x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                            x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                          })
                                        }
                                        
                                        out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                          
                                          if(k != "none"){
                                            
                                            if(any(stringr::str_length(x[,1])>4)){
                                              x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                              x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                              x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                              
                                              if(k == "decompose"){
                                                x <- na.omit(x) #remove NAs at head and tail of df
                                              }
                                            }
                                          }
                                          
                                          time_vec <- colnames(x)[1]
                                          
                                          foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                            
                                            out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                      winsize = 50)
                                            
                                            data.frame(cbind(data_source = paste(i),
                                                             out$cor))
                                          }
                                        }) |>
                                          `names<-`(c("none","average","decompose","stl")) |>
                                          data.table::rbindlist(idcol="deseason_meth")
                                        
                                        return(out_inner)
                                      }) |>
                                        `names<-`(c("none","linear","gaussian","loess")) |>
                                        data.table::rbindlist(idcol="detrend_meth")
                                      return(out)
                                      
                                    }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling",
    troph_level = "phytoplankton") |>
  select(-name)


roll_uni_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                       wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                       kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                       wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                  FUN = function(x){
                                    
                                    out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                      
                                      if(j != "none"){
                                        x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                        x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                          x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                        })
                                      }
                                      
                                      out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                        
                                        if(k != "none"){
                                          
                                          if(any(stringr::str_length(x[,1])>4)){
                                            x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                            x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                            x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                            
                                            if(k == "decompose"){
                                              x <- na.omit(x) #remove NAs at head and tail of df
                                            }
                                          }
                                        }
                                        
                                        time_vec <- colnames(x)[1]
                                        
                                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                          
                                          out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                    winsize = 50)
                                          
                                          data.frame(cbind(data_source = paste(i),
                                                           out$cor))
                                        }
                                      }) |>
                                        `names<-`(c("none","average","decompose","stl")) |>
                                        data.table::rbindlist(idcol="deseason_meth")
                                      
                                      return(out_inner)
                                    }) |>
                                      `names<-`(c("none","linear","gaussian","loess")) |>
                                      data.table::rbindlist(idcol="detrend_meth")
                                    return(out)
                                    
                                  }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling",
    troph_level = "zooplankton") |>
  select(-name)


roll_unicomp <- rbind(roll_uni_phyto,roll_uni_zoo)

################################################################################################################
## Permuted Univariate EWS Assessment Rolling ##
################################################################################################################
phyto_ls <- list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                 wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                 kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                 wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

perm_roll_uni_phyto <- lapply(seq_along(phyto_ls),
                              FUN = function(l){
                                
                                x <- phyto_ls[[l]]
                                
                                out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                  
                                  if(j != "none"){
                                    x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                    x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                      x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                    })
                                  }
                                  
                                  out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                    
                                    if(k != "none"){
                                      
                                      if(any(stringr::str_length(x[,1])>4)){
                                        x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                        x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                        x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                        
                                        if(k == "decompose"){
                                          x <- na.omit(x) #remove NAs at head and tail of df
                                        }
                                      }
                                    }
                                    
                                    time_vec <- colnames(x)[1]
                                    print(paste(names(phyto_ls)[l],j,k,collapse = "_"))
                                    
                                    foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                      
                                      out <-  perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                           winsize = 50, perm.meth = "replacement", variate = "uni",
                                                           iter = 500,mc.cores=5)   
                                      
                                      data.frame(cbind(data_source = paste(i),
                                                       out$cor))
                                    }
                                  }) |>
                                    `names<-`(c("none","average","decompose","stl")) |>
                                    data.table::rbindlist(idcol="deseason_meth")
                                  
                                  return(out_inner)
                                }) |>
                                  `names<-`(c("none","linear","gaussian","loess")) |>
                                  data.table::rbindlist(idcol="detrend_meth")
                                return(out)
                                
                              }) |>
  `names<-`(c("wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling",
    troph_level = "phytoplankton") |>
  select(-name)

zoo_ls <- list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
               wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
               kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
               wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))

perm_roll_uni_zoo <- lapply(seq_along(zoo_ls),
                            FUN = function(l){
                              
                              x <- zoo_ls[[l]]
                              
                              out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                
                                if(j != "none"){
                                  x <- EWSmethods::detrend_ts(x,method = j,span = 0.5) 
                                  x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                    x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                  })
                                }
                                
                                out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                  
                                  if(k != "none"){
                                    
                                    if(any(stringr::str_length(x[,1])>4)){
                                      x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                      x <- EWSmethods::deseason_ts(x,increment = "month", method = k,order = "ymd")
                                      x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                      
                                      if(k == "decompose"){
                                        x <- na.omit(x) #remove NAs at head and tail of df
                                      }
                                    }
                                  }
                                  
                                  time_vec <- colnames(x)[1]
                                  
                                  print(paste(names(zoo_ls)[l],j,k,collapse = "_"))
                                  
                                  foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                    
                                    out <-  perm_rollEWS(data = x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),
                                                         winsize = 50, perm.meth = "replacement", variate = "uni",
                                                         iter = 500,mc.cores=5)   
                                    
                                    data.frame(cbind(data_source = paste(i),
                                                     out$cor))
                                  }
                                }) |>
                                  `names<-`(c("none","average","decompose","stl")) |>
                                  data.table::rbindlist(idcol="deseason_meth")
                                
                                return(out_inner)
                              }) |>
                                `names<-`(c("none","linear","gaussian","loess")) |>
                                data.table::rbindlist(idcol="detrend_meth")
                              return(out)
                              
                            }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling",
    troph_level = "zooplankton") |>
  select(-name)


perm_roll_unicomp <- rbind(perm_roll_uni_phyto,perm_roll_uni_zoo)

################################################################################################################
## EWSNet ##
################################################################################################################
ewsnet_init("EWSNET_env", auto=T)

ewsnet_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                       wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                       kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                       wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                  FUN = function(x){
                                    
                                    out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                      
                                      if(j != "none"){
                                        x <- detrend_ts(x,method = j,span = 0.5) 
                                        x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                          x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                        })
                                      }
                                      
                                      out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                        
                                        if(k != "none"){
                                          if(any(stringr::str_length(x[,1])>4)){
                                            x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                            x <- deseason_ts(x,increment = "month", method = k,order = "ymd")
                                            x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                            
                                            if(k == "decompose"){
                                              x <- na.omit(x) #remove NAs at head and tail of df
                                            }
                                          }
                                        }
                                        
                                        time_vec <- colnames(x)[1]
                                        
                                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                          
                                          ews.tmp.scaled <- EWSmethods::ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                          ews.tmp.unscaled <- EWSmethods::ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                          
                                          data.frame(data_source = paste(i),
                                                     model_ensemble = "25",
                                                     predictionScaled = ews.tmp.scaled$pred,
                                                     no_transScaled = ews.tmp.scaled$no_trans_prob,
                                                     smth_transScaled =  ews.tmp.scaled$smooth_trans_prob,
                                                     crt_transScaled = ews.tmp.scaled$critical_trans_prob,
                                                     predictionUnscaled = ews.tmp.unscaled$pred,
                                                     no_transUnscaled = ews.tmp.unscaled$no_trans_prob,
                                                     smth_transUnscaled =  ews.tmp.unscaled$smooth_trans_prob,
                                                     crt_transUnscaled = ews.tmp.unscaled$critical_trans_prob)
                                          
                                        } 
                                      }) |>
                                        `names<-`(c("none","average","decompose","stl")) |>
                                        data.table::rbindlist(idcol="deseason_meth")
                                      
                                      return(out_inner)
                                    }) |>
                                      `names<-`(c("none","linear","gaussian","loess")) |>
                                      data.table::rbindlist(idcol="detrend_meth")
                                    return(out)
                                    
                                  }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "ML",
    troph_level = "phytoplankton") |>
  select(-name)


ewsnet_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                     wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                     kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                     wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                FUN = function(x){
                                  
                                  out <- lapply(c("none","linear","gaussian","loess"),function(j){
                                    
                                    if(j != "none"){
                                      x <- detrend_ts(x,method = j,span = 0.5) 
                                      x[,-1] <- sapply(colnames(x)[-1],FUN = function(i){
                                        x[,i] <- x[,i] + abs(min(x[,i])) #set 0 as minimum
                                      })
                                    }
                                    
                                    out_inner <- lapply(c("none","average","decompose","stl"), function(k){
                                      
                                      if(k != "none"){
                                        if(any(stringr::str_length(x[,1])>4)){
                                          x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                          x <- deseason_ts(x,increment = "month", method = k,order = "ymd")
                                          x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                          
                                          if(k == "decompose"){
                                            x <- na.omit(x) #remove NAs at head and tail of df
                                          }
                                        }
                                      }
                                      
                                      time_vec <- colnames(x)[1]
                                      
                                      foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                        
                                        ews.tmp.scaled <- EWSmethods::ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                        ews.tmp.unscaled <- EWSmethods::ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                                        
                                        data.frame(data_source = paste(i),
                                                   model_ensemble = "25",
                                                   predictionScaled = ews.tmp.scaled$pred,
                                                   no_transScaled = ews.tmp.scaled$no_trans_prob,
                                                   smth_transScaled =  ews.tmp.scaled$smooth_trans_prob,
                                                   crt_transScaled = ews.tmp.scaled$critical_trans_prob,
                                                   predictionUnscaled = ews.tmp.unscaled$pred,
                                                   no_transUnscaled = ews.tmp.unscaled$no_trans_prob,
                                                   smth_transUnscaled =  ews.tmp.unscaled$smooth_trans_prob,
                                                   crt_transUnscaled = ews.tmp.unscaled$critical_trans_prob)
                                        
                                      } 
                                    }) |>
                                      `names<-`(c("none","average","decompose","stl")) |>
                                      data.table::rbindlist(idcol="deseason_meth")
                                    
                                    return(out_inner)
                                  }) |>
                                    `names<-`(c("none","linear","gaussian","loess")) |>
                                    data.table::rbindlist(idcol="detrend_meth")
                                  return(out)
                                  
                                }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat","wash_yr_dat","leve_yr_dat","UZ_yr_dat","mon_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat","wash_mth_dat","leve_mth_dat","UZ_mth_dat","mon_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere",
    grepl("UZ", name) ~ "Upper Zurich",
    grepl("mon", name) ~ "Monona",
    grepl("leve", name) ~ "Loch Leven",
    grepl("wash", name) ~ "Washington"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "ML",
    troph_level = "zooplankton") |>
  select(-name)

ewsnet_comp <- rbind(ewsnet_phyto,ewsnet_zoo) |>
  mutate(model = "Original")

################################################################################################################
## Combine ##
################################################################################################################

#ewsnet_comp <- read.csv("/Users/ul20791/Downloads/ewsnet_comp.csv")[,-1]

save(ewsnet_comp,exp_multicomp,roll_unicomp,roll_multicomp,perm_roll_unicomp,perm_roll_multicomp, file = "Results/ews_raw_data_a.RData")
save(exp_uni_phyto,file = "Results/ews_raw_data_b.RData")
save(exp_uni_zoo,file = "Results/ews_raw_data_c.RData")
