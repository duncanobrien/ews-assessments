################################################################################################################
## Preamble ##
################################################################################################################
require(EWSmethods)
require(tidyverse)
require(foreach) 
require(pbapply)

EWSmethods::ewsnet_init("EWSNET_env", auto=T)
reticulate::py_config() #check for 'forced by RETICULATE_PYTHON
source("Code/perm_rollEWS_fn.R")
# save(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
#      kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat,
#     file =  "Data/wrangled_plank_data.Rdata" )

################################################################################################################
## Load Data ##
################################################################################################################

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kinneret/Data/kinneret_environmental_data.R")
#phyo, zoo + combined raw concentrations at week, month and year groupings

source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/zurich_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Zurich/Data/LZ_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Madison/Data/madison_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Windermere/Data/windermere_environmental_data.R")

source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/kasumigaura_plankton_data.R")
source("/Users/ul20791/Desktop/Academia/PhD/Data/Kasumigaura/Data/Kasumigaura_environmental_data.R")

################################################################################################################
## Wangle Data ##
#prewrangled load("Data/wrangled_plank_data.Rdata")
################################################################################################################

kin_yr_dat <- plank_env.data.yr |> #drop environmentals
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[25])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(`2-Anabaena`:`Anuraeopsis fissa`,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(`2-Anabaena`:`Anuraeopsis fissa`) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  dplyr::select(c(Date:`Anuraeopsis fissa`,pca1,pca2))|>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(Date)) |>
  as.data.frame()

kin_mth_dat <- plank_env.data.mth |> 
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[288])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(`2-Anabaena`:`Anuraeopsis fissa`,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(`2-Anabaena`:`Anuraeopsis fissa`) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(Date,`2-Anabaena`:`Anuraeopsis fissa`,pca1,pca2))|>
  mutate(Date = base::as.Date(zoo::as.Date(Date)))|>
  deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date)))

kas_yr_dat <- plank_env.kasyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[16])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Acanthoceras_zachariasii:Thermocyclops_taihokuensis,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Acanthoceras_zachariasii:Thermocyclops_taihokuensis) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(date:Thermocyclops_taihokuensis,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

kas_mth_dat <- plank_env.kasmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[185])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Acanthoceras_zachariasii:Thermocyclops_taihokuensis,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Acanthoceras_zachariasii:Thermocyclops_taihokuensis) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(date,Acanthoceras_zachariasii:Thermocyclops_taihokuensis,pca1,pca2))|>
  mutate(date = base::as.Date(zoo::as.Date(date)))|>
  deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date)))

LZ_yr_dat <- plank_env.LZyrdata |>
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Microcystis_sp:Polyphemus_pediculus,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Microcystis_sp:Polyphemus_pediculus) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(date,Microcystis_sp:Polyphemus_pediculus,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

LZ_mth_dat <- plank_env.LZmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[length(date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Microcystis_sp:Polyphemus_pediculus,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Microcystis_sp:Polyphemus_pediculus) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(date,Microcystis_sp:Polyphemus_pediculus,pca1,pca2))|>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date)))

mad_yr_dat <- plank_env.madyrdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[15]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2))|>
  select(c(date:TROPOCYCLOPS_PRASINUS_MEXICANUS,pca1,pca2))|>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(date)) |>
  as.data.frame()

mad_mth_dat <- plank_env.madmthdata |> 
  dplyr::filter(as.numeric(date) <= as.numeric(date[168]))|>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS) %>%
                                   dplyr::select_if(colSums(.) != 0))) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  dplyr::select(c(date,`Achnanthes minutissima`:TROPOCYCLOPS_PRASINUS_MEXICANUS,pca1,pca2))|>
  mutate(date = base::as.Date(zoo::as.Date(date))) |>
  deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date))) |>
  dplyr::select(-Aphanothece.clathrata) #randomly dropped species due to bug in svd

wind_yr_dat <- phyto_env.windyrdata |>
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:TotCyclopoids) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(Date:TotCyclopoids,pca1,pca2)) |>
  mutate(across(c(-Date), ~ .x +  abs(min(.x))),
         Date = as.numeric(Date)) |>
  as.data.frame()

wind_mth_dat <- phyto_env.windmthdata |>
  dplyr::filter(as.numeric(Date) <= as.numeric(Date[length(Date)*0.85])) |>
  dplyr::select_if(~ !is.numeric(.) || sum(.) != 0) |>
  mutate(across(Asterionella:TotCyclopoids,~log1p(.x))) |>
  nest(data = everything())|>
  dplyr::mutate(tmp = purrr::map(data, ~.x |> 
                                   dplyr::select(Asterionella:TotCyclopoids) %>%
                                   dplyr::select_if(colSums(.) != 0) )) |>
  dplyr::mutate(pca1 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,1]),
                pca2 = purrr::map(tmp, ~prcomp(.x,scale. = T)$x[,2])) |>
  dplyr::select(-tmp)|>
  unnest(c(data,pca1,pca2)) |>
  select(c(Date,Asterionella:TotCyclopoids,pca1,pca2))|>
  mutate(Date = base::as.Date(zoo::as.Date(Date))) |>
  deseason_ts(increment = "month",order = "ymd",method="stl")  |>
  mutate(across(c(-date), ~ .x +  abs(min(.x))),
         date = as.numeric(zoo::as.yearmon(date)))

################################################################################################################
## EWSNet Assessment ##
################################################################################################################
list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
     kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat)

lake_ewsnet <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
                           kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat),
                      FUN = function(x){
                        
                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                          
                          ews.tmp.scaled <- ewsnet_predict(c(x[,i]),scaling = T, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                          ews.tmp.unscaled <- ewsnet_predict(c(x[,i]),scaling = F, ensemble = 25,envname = "EWSNET_env") #only run on pre-regime shift
                          
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
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "EWSNet",
    computation = "ML") |>
  select(-name)

write.csv(lake_ewsnet,file = "Results/lake_results/lake_ewsnet.csv")

################################################################################################################
## Multivariate EWS Assessment ##
################################################################################################################

lake_roll_multi_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
                           kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat),
                      FUN = function(x){
                        
                        # if(any(colMeans(x == 0) <= 0.5)){
                        #   
                        #   x[,-1] <- x[,-1] + 1
                        #   
                        # }
                        sub_dat <- x[sapply(colnames(x),FUN = function(i){
                           mean(x[,i] == min(x[,i])) <= 0.49
                        })]
                        #sub_dat <- x[colMeans(x == 0) <= 0.49] #ensures not 50% of timeseries is the same value to break rolling window
                    
                    out <- EWSmethods::multiEWS(sub_dat,method = "rolling",ggplotIt = F,
                                          winsize = 50)
                    
                    return(data.frame(out$cor))
                        
                      }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling") |>
  select(-name)

write.csv(lake_roll_multi_ews,file = "Results/lake_results/lake_roll_multi_ews.csv")

################################################################################################################
## Multivariate EWS Assessment Expanding ##
################################################################################################################

lake_exp_multi_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
                                              kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat),
                                         FUN = function(x){
                                           
                                           # if(any(colMeans(x == 0) <= 0.5)){
                                           #   
                                           #   x[,-1] <- x[,-1] + 1
                                           #   
                                           # }
                                           sub_dat <- x[sapply(colnames(x),FUN = function(i){
                                             mean(x[,i] == min(x[,i])) <= 0.49
                                           })]
                                           #sub_dat <- x[colMeans(x == 0) <= 0.49] #ensures not 50% of timeseries is the same value to break rolling window
                                           
                                           out <- EWSmethods::multiEWS(sub_dat,method = "expanding",ggplotIt = T,
                                                                       burn_in  = ceiling(dim(sub_dat)[1]*0.50))
                                           
                                           return(out$raw)
                                           
                                         }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "multivariate EWS",
    computation = "rolling") |>
  select(-name)

write.csv(lake_exp_multi_ews,file = "Results/lake_results/lake_exp_multi_ews.csv")

out <- EWSmethods::multiEWS(sub_dat,method = "expanding",ggplotIt = T,
                            burn_in  = ceiling(dim(sub_dat)[1]*0.50))


################################################################################################################
## Univariate EWS Assessment ##
################################################################################################################
list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
     kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat)

lake_roll_uni_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
                                       kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat),
                      FUN = function(x){
                        
                        foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                          
                          time_vec <- colnames(x)[1]
                          out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "rolling",ggplotIt = F,
                                                      winsize = 50)
                          
                          data.frame(cbind(data_source = paste(i),
                                    out$cor))
                          
                        } 
                        
                      }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "rolling") |>
  select(-name)

write.csv(lake_roll_uni_ews,file = "Results/lake_results/lake_roll_uni_ews.csv")

################################################################################################################
## Univariate EWS Assessment Expanding ##
################################################################################################################
list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
     kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat)

lake_exp_uni_ews <- pbapply::pblapply(list(kin_yr_dat,kas_yr_dat,LZ_yr_dat,mad_yr_dat,wind_yr_dat,
                                       kin_mth_dat,kas_mth_dat,LZ_mth_dat,mad_mth_dat,wind_mth_dat),
                                  FUN = function(x){
                                    
                                    foreach(i=colnames(x[,-1]), .combine='rbind',.verbose = F) %do%{
                                      
                                      time_vec <- colnames(x)[1]
                                      out <- EWSmethods::uniEWS(x[,c(time_vec,i)],metrics = c("ar1","SD","skew"),method = "expanding",ggplotIt = F,
                                                                burn_in = ceiling(dim(x)[1]*0.25))
                                      
                                      data.frame(cbind(data_source = paste(i),
                                                       out))
                                      
                                    } 
                                    
                                  }) |>
  `names<-`(c("kin_yr_dat","kas_yr_dat","LZ_yr_dat","mad_yr_dat","wind_yr_dat",
              "kin_mth_dat","kas_mth_dat","LZ_mth_dat","mad_mth_dat","wind_mth_dat"))  |> 
  data.table::rbindlist(idcol="name") |> 
  mutate(lake = case_when(
    grepl("kin", name) ~ "Kinneret",
    grepl("kas", name) ~ "Kasumigaura",
    grepl("LZ", name) ~ "Lower Zurich",
    grepl("mad", name) ~ "Mendota",
    grepl("wind", name) ~ "Windermere"),
    res =  case_when(
      grepl("yr", name) ~ "Yearly",
      grepl("mth", name) ~ "Monthly"),
    method = "univariate EWS",
    computation = "expanding") |>
  select(-name)
  
  write.csv(lake_exp_uni_ews,file = "Results/lake_results/lake_exp_uni_ews.csv")
