require(dplyr)
require(EWSmethods)
require(foreach)
load("Data/wrangled_genus_plank_data.Rdata")
source("Code/extract_ews_pred_fn.R")

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
  list2env(env = .GlobalEnv)

################################################################################################################
## Multivariate EWS Assessment Expanding ##
################################################################################################################

exp_multi_detrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                                   wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                                   kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                                   wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                              FUN = function(x){
                                                
                                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                  mean(x[,i] == min(x[,i])) <= 0.47
                                                })]
                                                
                                                sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                                
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


exp_multi_detrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                 wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                 kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                 wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                            FUN = function(x){
                                              
                                              sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                mean(x[,i] == min(x[,i])) <= 0.47
                                              })]
                                              
                                              sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                              
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

exp_multi_deseason_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                                  wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                                  kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                                  wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                             FUN = function(x){
                                               
                                               sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                 mean(x[,i] == min(x[,i])) <= 0.47
                                               })]
                                               
                                               sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                               
                                               if(any(stringr::str_length(sub_dat[,1])>4)){
                                                 sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                                 sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = "average",order = "ymd")
                                                 sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
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


exp_multi_deseason_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                           FUN = function(x){
                                             
                                             sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                               mean(x[,i] == min(x[,i])) <= 0.47
                                             })]
                                             
                                             sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                             
                                             if(any(stringr::str_length(sub_dat[,1])>4)){
                                               sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                               sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = "average",order = "ymd")
                                               sub_dat[,1] <- as.numeric(zoo::as.yearmon(sub_dat[,1])) #convert from date back to numeric for rbind
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

exp_multi_nodetrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                                   wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                                   kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                                   wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                              FUN = function(x){
                                                
                                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                  mean(x[,i] == min(x[,i])) <= 0.47
                                                })]
                                                
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


exp_multi_nodetrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                 wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                 kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                 wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                            FUN = function(x){
                                              
                                              sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                mean(x[,i] == min(x[,i])) <= 0.47
                                              })]
                                              
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

exp_multi_detrendews <- rbind(exp_multi_detrend_phyto,exp_multi_detrend_zoo)
exp_multi_deseasonews <- rbind(exp_multi_deseason_phyto,exp_multi_deseason_zoo)
exp_multi_nodetrendews <- rbind(exp_multi_nodetrend_phyto,exp_multi_nodetrend_zoo)


lake_outcome_troph <- subset(exp_multi_detrendews, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)|>
  dplyr::mutate(label = ifelse(outcome == "trans",1,0))


expmultiews_detrend <- extract_ews_pred(ews.data = exp_multi_detrendews,
                                         sensitivity = 2,
                                         outcome = lake_outcome_troph,
                                         method = "expanding") |>
  mutate(processing = "detrend")

expmultiews_deseason <- extract_ews_pred(ews.data = exp_multi_deseasonews,
                                          sensitivity = 2,
                                          outcome = lake_outcome_troph,
                                          method = "expanding") |>
  mutate(processing = "deseason")

expmultiews_nodetrend <- extract_ews_pred(ews.data = exp_multi_nodetrendews,
                                           sensitivity = 2,
                                           outcome = lake_outcome_troph,
                                           method = "expanding")|>
  mutate(processing = "no detrend")


exp_multiews <- rbind(expmultiews_nodetrend,expmultiews_detrend,expmultiews_deseason)|>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") 

exp_multi_roc <- lapply(unique(exp_multiews$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 5)) |> `names<-`(c("fpr","tpr","auc","res","processing"))
  
  for(j in unique(exp_multiews$processing)){
    
    for(i in unique(exp_multiews$res)){
      sub_dat <- subset(exp_multiews,metric.code == x & res == i & processing == j) |>
        dplyr::select(lake,metric.code,res,processing,prop_run,label) |>
        na.omit()
      
      pred <- ROCR::prediction(sub_dat$prop_run, sub_dat$label)
      auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
      per <- ROCR::performance(pred,"tpr","fpr")
      out.dat <- rbind(out.dat,
                       data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i,"processing" = j))
      
    }
  }
  return(out.dat)
}) |>
  `names<-`(unique(exp_multiews$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "multivariate",
                computation = "expanding")


ggplot(exp_multi_roc,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_abline(slope=1) + 
  geom_line() + facet_grid(res~processing) + theme_bw()

ggplot(exp_multi_roc,
       aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + facet_wrap(~res) + theme_bw()

################################################################################################################
## Multivariate EWS Assessment Rolling ##
################################################################################################################

roll_multi_detrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                             wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                             kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                             wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                        FUN = function(x){
                                          
                                          sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                            mean(x[,i] == min(x[,i])) <= 0.47
                                          })]
                                          
                                          sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           

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


roll_multi_detrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                 wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                 kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                 wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                              FUN = function(x){
                                                
                                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                  mean(x[,i] == min(x[,i])) <= 0.47
                                                })]
                                                
                                                sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                                
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

roll_multi_deseason_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                                   wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                                   kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                                   wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                              FUN = function(x){
                                                
                                                sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                  mean(x[,i] == min(x[,i])) <= 0.47
                                                })]
                                                
                                                sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                                
                                                if(any(stringr::str_length(sub_dat[,1])>4)){
                                                  sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                                  sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = "average",order = "ymd")
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


roll_multi_deseason_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                 wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                 kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                 wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                            FUN = function(x){
                                              
                                              sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                mean(x[,i] == min(x[,i])) <= 0.47
                                              })]
                                              
                                              sub_dat <- EWSmethods::detrend_ts(sub_dat,method = "gaussian",span = 0.75)                                           
                                              
                                              if(any(stringr::str_length(sub_dat[,1])>4)){
                                                sub_dat[,1] <- zoo::as.Date(zoo::as.yearmon(sub_dat[,1]))
                                                sub_dat <- EWSmethods::deseason_ts(sub_dat,increment = "month", method = "average",order = "ymd")
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

roll_multi_nodetrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                                    wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                                    kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                                    wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                               FUN = function(x){
                                                 
                                                 sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                   mean(x[,i] == min(x[,i])) <= 0.47
                                                 })]
                                                 
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


roll_multi_nodetrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                                  wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                                  kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                                  wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                             FUN = function(x){
                                               
                                               sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                                 mean(x[,i] == min(x[,i])) <= 0.47
                                               })]
                                               
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

roll_multi_detrendews <- rbind(roll_multi_detrend_phyto,roll_multi_detrend_zoo)
roll_multi_deseasonews <- rbind(roll_multi_deseason_phyto,roll_multi_deseason_zoo)
roll_multi_nodetrendews <- rbind(roll_multi_nodetrend_phyto,roll_multi_nodetrend_zoo)


lake_outcome_troph <- subset(roll_multi_detrendews, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)|>
  dplyr::mutate(label = ifelse(outcome == "trans",1,0))


rollmultiews_detrend <- extract_ews_pred(ews.data = roll_multi_detrendews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling") |>
  mutate(processing = "detrend") 

rollmultiews_deseason <- extract_ews_pred(ews.data = roll_multi_deseasonews,
                                     sensitivity = 0.7,
                                     outcome = lake_outcome_troph,
                                     method = "rolling") |>
  mutate(processing = "deseason")

rollmultiews_nodetrend <- extract_ews_pred(ews.data = roll_multi_nodetrendews,
                                      sensitivity = 0.7,
                                      outcome = lake_outcome_troph,
                                      method = "rolling")|>
  mutate(processing = "no detrend")


roll_multiews <- rbind(rollmultiews_nodetrend,rollmultiews_detrend,rollmultiews_deseason)|>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") 

roll_multi_roc <- lapply(unique(roll_multiews$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 5)) |> `names<-`(c("fpr","tpr","auc","res","processing"))
  
  for(j in unique(roll_multiews$processing)){
    
    for(i in unique(roll_multiews$res)){
      sub_dat <- subset(roll_multiews,metric.code == x & res == i & processing == j) |>
        dplyr::select(-data_source) |>
        na.omit()
      
      pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
      auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
      per <- ROCR::performance(pred,"tpr","fpr")
      out.dat <- rbind(out.dat,
                       data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i,"processing" = j))
      
    }
  }
  return(out.dat)
}) |>
  `names<-`(unique(roll_multiews$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "multivariate",
                computation = "rolling")

ggplot(roll_multi_roc,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_abline(slope=1) + 
  geom_line() + facet_grid(res~processing) + theme_bw()

ggplot(roll_multi_roc,
       aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + facet_wrap(~res) + theme_bw()

################################################################################################################
## Univariate EWS Assessment Rolling ##
################################################################################################################

roll_detrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                             wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                             kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                             wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                        FUN = function(x){
                                          
                                          x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                          
                                          time_vec <- colnames(x)[1]
                                          
                                          foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                            
                                            out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                      winsize = 50)
                                            
                                            data.frame(cbind(data_source = paste(i),
                                                             out$cor))
                                            
                                          } 
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

roll_detrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                             wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                             kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                             wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                        FUN = function(x){
                                          
                                          x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                          
                                          time_vec <- colnames(x)[1]
                                          
                                          foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                            
                                            out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                      winsize = 50)
                                            
                                            data.frame(cbind(data_source = paste(i),
                                                             out$cor))
                                            
                                          } 
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


roll_deseason_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                              wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                              kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                              wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                  FUN = function(x){
                                    
                                    x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                    
                                    if(any(stringr::str_length(x[,1])>4)){
                                      x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                      x <- EWSmethods::deseason_ts(x,increment = "month", method = "average",order = "ymd")
                                    }
                                    
                                    time_vec <- colnames(x)[1]
                                    
                                    foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                      
                                      out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                winsize = 50)
                                      
                                      data.frame(cbind(data_source = paste(i),
                                                       out$cor))
                                      
                                    } 
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

roll_deseason_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                            wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                            kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                            wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                         FUN = function(x){
                                           
                                           x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                           
                                           if(any(stringr::str_length(x[,1])>4)){
                                             x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                             x <- EWSmethods::deseason_ts(x,increment = "month", method = "average",order = "ymd")
                                           }
                                           
                                           time_vec <- colnames(x)[1]
                                           
                                           foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                             
                                             out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                       winsize = 50)
                                             
                                             data.frame(cbind(data_source = paste(i),
                                                              out$cor))
                                             
                                           } 
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

roll_nodetrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                               wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                               kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                               wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                  FUN = function(x){
                                    
                                    time_vec <- colnames(x)[1]
                                    
                                    foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                      
                                      out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                winsize = 50)
                                      
                                      data.frame(cbind(data_source = paste(i),
                                                       out$cor))
                                      
                                    } 
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

roll_nodetrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                             wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                             kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                             wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                          FUN = function(x){
                                            
                                            time_vec <- colnames(x)[1]
                                            
                                            foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                              
                                              out <- EWSmethods::uniEWS(data =x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                                        winsize = 50)
                                              
                                              data.frame(cbind(data_source = paste(i),
                                                               out$cor))
                                              
                                            } 
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

roll_uni_detrendews <- rbind(roll_detrend_phyto,roll_detrend_zoo)
roll_uni_deseasonews <- rbind(roll_deseason_phyto,roll_deseason_zoo)
roll_uni_nodetrendews <- rbind(roll_nodetrend_phyto,roll_nodetrend_zoo)


lake_outcome_troph <- subset(roll_uni_detrendews, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)|>
  dplyr::mutate(label = ifelse(outcome == "trans",1,0))


rollews_detrend <- extract_ews_pred(ews.data = roll_uni_detrendews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling") |>
  mutate(processing = "detrend")

rollews_deseason <- extract_ews_pred(ews.data = roll_uni_deseasonews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling") |>
  mutate(processing = "deseason")

rollews_nodetrend <- extract_ews_pred(ews.data = roll_uni_nodetrendews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling")|>
  mutate(processing = "no detrend")


rollews <- rbind(rollews_nodetrend,rollews_detrend,rollews_deseason)|>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") 

roll_uni_roc <- lapply(unique(rollews$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 5)) |> `names<-`(c("fpr","tpr","auc","res","processing"))
  
  for(j in unique(rollews$processing)){
  
  for(i in unique(rollews$res)){
    sub_dat <- subset(rollews,metric.code == x & res == i & processing == j)
    
    pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i,"processing" = j))
    
  }
  }
  return(out.dat)
}) |>
  `names<-`(unique(rollews$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "univariate",
                computation = "rolling")

ggplot(roll_uni_roc,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_abline(slope=1)+
  geom_line() + facet_grid(res~processing) + theme_bw()

ggplot(roll_uni_roc,
       aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + facet_wrap(~res) + theme_bw()

################################################################################################################
## Univariate EWS Assessment Expanding ##
################################################################################################################

exp_detrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                             wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                             kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                             wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                        FUN = function(x){
                                          
                                          x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                          
                                          time_vec <- colnames(x)[1]
                                          
                                          foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                            
                                            out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                      burn_in = ceiling(dim(x)[1]*0.50))
                                            
                                            data.frame(cbind(data_source = paste(i),
                                                             out))
                                            
                                          } 
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

exp_detrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                           wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                           kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                           wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                      FUN = function(x){
                                        
                                        x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                        
                                        time_vec <- colnames(x)[1]
                                        
                                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                          
                                          out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                    burn_in = ceiling(dim(x)[1]*0.50))
                                          
                                          data.frame(cbind(data_source = paste(i),
                                                           out))
                                          
                                        } 
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
    troph_level = "zooplankton") |>
  select(-name)


exp_deseason_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                              wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                              kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                              wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                         FUN = function(x){
                                           
                                           x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                           
                                           if(any(stringr::str_length(x[,1])>4)){
                                             x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                             x <- EWSmethods::deseason_ts(x,increment = "month", method = "average",order = "ymd")
                                             x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                             }
                                           
                                           time_vec <- colnames(x)[1]
                                           
                                           foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                             
                                             out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                       burn_in = ceiling(dim(x)[1]*0.50))
                                             
                                             data.frame(cbind(data_source = paste(i),
                                                              out))
                                             
                                           } 
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

exp_deseason_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                            wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                            kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                            wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                       FUN = function(x){
                                         
                                         x <- EWSmethods::detrend_ts(x,method = "gaussian",span = 0.75)                                           
                                         
                                         if(any(stringr::str_length(x[,1])>4)){
                                           x[,1] <- zoo::as.Date(zoo::as.yearmon(x[,1]))
                                           x <- EWSmethods::deseason_ts(x,increment = "month", method = "average",order = "ymd")
                                           x[,1] <- as.numeric(zoo::as.yearmon(x[,1])) #convert from date back to numeric for rbind
                                            }
                                         
                                         time_vec <- colnames(x)[1]
                                         
                                         foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                           
                                           out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                     burn_in = ceiling(dim(x)[1]*0.50))
                                           
                                           data.frame(cbind(data_source = paste(i),
                                                            out))
                                           
                                         } 
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
    troph_level = "zooplankton") |>
  select(-name)

exp_nodetrend_phyto <- pbapply::pblapply(list(kin_yr_dat[,1:14],kas_yr_dat[,1:25],LZ_yr_dat[,1:45],mad_yr_dat[,1:32],
                                               wind_yr_dat[,1:18],wash_yr_dat[,1:7],leve_yr_dat[,1:4],UZ_yr_dat[,1:59],mon_yr_dat[,1:45],
                                               kin_mth_dat[,1:14],kas_mth_dat[,1:25],LZ_mth_dat[,1:36],mad_mth_dat[,1:17],
                                               wind_mth_dat[,1:12],wash_mth_dat[,1:7],leve_mth_dat[,1:4],UZ_mth_dat[,1:58],mon_mth_dat[,1:29]),
                                          FUN = function(x){
                                            
                                            time_vec <- colnames(x)[1]
                                            
                                            foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                              
                                              out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                        burn_in = ceiling(dim(x)[1]*0.50))
                                              
                                              data.frame(cbind(data_source = paste(i),
                                                               out))
                                              
                                            } 
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

exp_nodetrend_zoo <- pbapply::pblapply(list(kin_yr_dat[,c(1,15:32)],kas_yr_dat[,c(1,26:34)],LZ_yr_dat[,c(1,46:51)],mad_yr_dat[,c(1,33:38)],
                                             wind_yr_dat[,c(1,20:23)],wash_yr_dat[,c(1,8:15)],leve_yr_dat[,c(1,5:7)],UZ_yr_dat[,c(1,60:65)],mon_yr_dat[,c(1,46:53)],
                                             kin_mth_dat[,c(1,15:22)],kas_mth_dat[,c(1,26:34)],LZ_mth_dat[,c(1,37:41)],mad_mth_dat[,c(1,18:22)],
                                             wind_mth_dat[,c(1,14:17)],wash_mth_dat[,c(1,8:13)],leve_mth_dat[,c(1,5:6)],UZ_mth_dat[,c(1,59:63)],mon_mth_dat[,c(1,30:36)]),
                                        FUN = function(x){
                                          
                                          time_vec <- colnames(x)[1]
                                          
                                          foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                            
                                            out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                      burn_in = ceiling(dim(x)[1]*0.50))
                                            
                                            data.frame(cbind(data_source = paste(i),
                                                             out))
                                            
                                          } 
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
    troph_level = "zooplankton") |>
  select(-name)

exp_uni_detrendews <- rbind(exp_detrend_phyto,exp_detrend_zoo)
exp_uni_deseasonews <- rbind(exp_deseason_phyto,exp_deseason_zoo)
exp_uni_nodetrendews <- rbind(exp_nodetrend_phyto,exp_nodetrend_zoo)

lake_outcome_troph <- subset(exp_uni_detrendews, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)|>
  dplyr::mutate(label = ifelse(outcome == "trans",1,0))


expews_detrend <- extract_ews_pred(ews.data = exp_uni_detrendews,
                                    sensitivity = 2,
                                    outcome = lake_outcome_troph,
                                    method = "expanding") |>
  mutate(processing = "detrend")

expews_deseason <- extract_ews_pred(ews.data = exp_uni_deseasonews,
                                     sensitivity = 2,
                                     outcome = lake_outcome_troph,
                                     method = "expanding") |>
  mutate(processing = "deseason")

expews_nodetrend <- extract_ews_pred(ews.data = exp_uni_nodetrendews,
                                      sensitivity = 2,
                                      outcome = lake_outcome_troph,
                                      method = "expanding")|>
  mutate(processing = "no detrend")


expews <- rbind(expews_nodetrend,expews_detrend,expews_deseason)|>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference")

exp_uni_roc <- lapply(unique(expews$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 5)) |> `names<-`(c("fpr","tpr","auc","res","processing"))
  
  for(j in unique(expews$processing)){
    
    for(i in unique(expews$res)){
      sub_dat <- subset(expews,metric.code == x & res == i & processing == j)
      
      pred <- ROCR::prediction(sub_dat$prop_run, sub_dat$label)
      auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
      per <- ROCR::performance(pred,"tpr","fpr")
      out.dat <- rbind(out.dat,
                       data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i,"processing" = j))
      
    }
  }
  return(out.dat)
}) |>
  `names<-`(unique(expews$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "univariate",
                computation = "expanding")

ggplot(exp_uni_roc,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_abline(slope=1)+
  geom_line() + facet_grid(res~processing) + theme_bw()

ggplot(exp_uni_roc,
       aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + facet_wrap(~res) + theme_bw()


################################################################################################################
## Combine ##
################################################################################################################

tot_roc <- rbind(exp_uni_roc,exp_multi_roc,roll_uni_roc,roll_multi_roc)

ggpubr::ggarrange(
  ggplot(subset(tot_roc, variate == "univariate"),
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_abline(slope=1)+
  geom_line() + 
  ggh4x::facet_nested(res ~computation + processing,scales = "free") +
    ggtitle("a) univariate") +
  theme_bw(),
  ggplot(subset(tot_roc, variate == "multivariate"),
         aes(x=fpr,y=tpr,col=metric.code)) + 
    geom_abline(slope=1)+
    geom_line() + 
    ggh4x::facet_nested(res ~computation + processing,scales = "free") +
    ggtitle("b) multivariate") +
    theme_bw(),
  ncol=2
  )

ggplot(tot_roc,
       aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + facet_grid(computation~res) + theme_bw()


ggpubr::ggarrange(
  ggplot(subset(tot_roc, variate == "univariate"),
         aes(x=processing,y=auc)) + 
    geom_hline(yintercept=0.5,col="grey") + 
    geom_text(aes(label=metric.code)) + 
    facet_grid(computation~res) + 
    ggtitle("a) univariate") +
    theme_bw(),
  ggplot(subset(tot_roc, variate == "multivariate"),
    aes(x=processing,y=auc)) + 
  geom_hline(yintercept=0.5,col="grey") + 
  geom_text(aes(label=metric.code)) + 
  facet_grid(computation~res) + 
  ggtitle("b) multivariate") +
  theme_bw(),
  ncol=2
)
