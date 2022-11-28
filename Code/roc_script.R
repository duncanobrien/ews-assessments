##########################################################################################
# Preamble
##########################################################################################
require(ROCR)
require(Matrix)
require(tidyverse)
source("Code/extract_ews_pred_fn.R")

genus_lake_roll_uni_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_phyto.csv")[,-1],
                                 read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_zoo.csv")[,-1])

genus_lake_roll_multi_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_phyto.csv")[,-1],
                                   read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_zoo.csv")[,-1])

genus_lake_ewsnet <- rbind(read.csv(file = "Results/lake_results/genus_lake_ewsnet_phyto.csv")[,-1],
                           read.csv(file = "Results/lake_results/genus_lake_ewsnet_zoo.csv")[,-1])

genus_lake_exp_multi_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_exp_multi_ews_phyto.csv")[,-1], #drop write.csv() introduced X column
                                  read.csv(file = "Results/lake_results/genus_lake_exp_multi_ews_zoo.csv")[,-1])

genus_lake_exp_uni_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_exp_uni_ews_phyto.csv.gz")[,-1],
                                read.csv(file = "Results/lake_results/genus_lake_exp_uni_ews_zoo.csv.gz")[,-1])


lake_outcome_troph <- subset(genus_lake_roll_uni_ews, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T) |>
  dplyr::mutate(label = ifelse(outcome == "trans",1,0))

##########################################################################################
# Extract EWS success rates
##########################################################################################

rolluniews_diff_df <- extract_ews_pred(ews.data = genus_lake_roll_uni_ews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling") |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") |>
  dplyr::select(data_source,lake,res,troph_level,metric.code,metric.score,label) 

rollmultiews_diff_df <- extract_ews_pred(ews.data = genus_lake_roll_multi_ews,
                                    sensitivity = 0.7,
                                    outcome = lake_outcome_troph,
                                    method = "rolling") |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") |>
  dplyr::select(lake,res,troph_level,metric.code,metric.score,label) |>
  na.omit()

ewsnet_diff_df <- extract_ews_pred(ews.data =  genus_lake_ewsnet,
                                   sensitivity = 0.7,
                                   outcome = lake_outcome_troph,
                                   method = "ML") |>
  dplyr::filter(pos_outcome == "critical")|>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") |>
  dplyr::select(lake,res,troph_level,scaling,prob,label) 


expuniews_diff_df <-extract_ews_pred(ews.data = genus_lake_exp_uni_ews,
                                     sensitivity = 2,
                                     outcome = lake_outcome_troph,
                                     method = "expanding") |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") |>
  dplyr::select(lake,res,troph_level,metric.code,longest_run,prop_run,label) |>
  dplyr::rename("metric.score" = prop_run)


expmultiews_diff_df <- extract_ews_pred(ews.data = genus_lake_exp_multi_ews,
                                        sensitivity = 2,
                                        outcome = lake_outcome_troph,
                                        method = "expanding") |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_")) |>
  dplyr::left_join(lake_outcome_troph,by = "reference") |>
  dplyr::select(lake,res,troph_level,metric.code,longest_run,prop_run,label) |>
  dplyr::rename("metric.score" = prop_run)


##########################################################################################
# Fit ROCs
##########################################################################################

roll_multi_roc <- lapply(unique(rollmultiews_diff_df$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 4)) |> `names<-`(c("fpr","tpr","auc","res"))
  for(i in unique(rollmultiews_diff_df$res)){
    sub_dat <- subset(rollmultiews_diff_df,metric.code == x & res == i)
    
    pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i))
    
  }
  return(out.dat)
  }) |>
  `names<-`(unique(rollmultiews_diff_df$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code")|>
  dplyr::mutate(variate = "multivariate",
                computation = "rolling")

roll_uni_roc <- lapply(unique(rolluniews_diff_df$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 4)) |> `names<-`(c("fpr","tpr","auc","res"))
  for(i in unique(rollmultiews_diff_df$res)){
    sub_dat <- subset(rolluniews_diff_df,metric.code == x & res == i)
    
    pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i))
    
  }
  return(out.dat)
}) |>
  `names<-`(unique(rolluniews_diff_df$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "univariate",
                computation = "rolling")


ewsnet_roc <- lapply(unique(ewsnet_diff_df$scaling),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 4)) |> `names<-`(c("fpr","tpr","auc","res"))
  for(i in unique(ewsnet_diff_df$res)){
    sub_dat <- subset(ewsnet_diff_df,scaling == x & res == i)
    
    pred <- ROCR::prediction(sub_dat$prob, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i))
    
  }
  return(out.dat)
}) |>
  `names<-`(unique(ewsnet_diff_df$scaling)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "univariate",
                computation = "ML")

exp_multi_roc <- lapply(unique(expmultiews_diff_df$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 4)) |> `names<-`(c("fpr","tpr","auc","res"))
  for(i in unique(expmultiews_diff_df$res)){
    sub_dat <- subset(expmultiews_diff_df,metric.code == x & res == i)
    
    pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i))
    
  }
  return(out.dat)
}) |>
  `names<-`(unique(expmultiews_diff_df$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code")|>
  dplyr::mutate(variate = "multivariate",
                computation = "expanding")
  

exp_uni_roc <- lapply(unique(expuniews_diff_df$metric.code),function(x){
  
  out.dat <- data.frame(matrix(NA,nrow = 0,ncol = 4)) |> `names<-`(c("fpr","tpr","auc","res"))
  for(i in unique(expmultiews_diff_df$res)){
    sub_dat <- subset(expuniews_diff_df,metric.code == x & res == i)
    
    pred <- ROCR::prediction(sub_dat$metric.score, sub_dat$label)
    auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
    per <- ROCR::performance(pred,"tpr","fpr")
    out.dat <- rbind(out.dat,
                     data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],"auc" = auc,"res" = i))
    
  }
  return(out.dat)
}) |>
  `names<-`(unique(expuniews_diff_df$metric.code)) |>
  data.table::rbindlist(idcol = "metric.code") |>
  dplyr::mutate(variate = "univariate",
                computation = "expanding")

##########################################################################################
# Visualise
##########################################################################################

ggpubr::ggarrange(ggplot(rbind(roll_multi_roc,roll_uni_roc) ,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_line() + facet_grid(res~variate) + theme_bw() + ggtitle("a) rolling computation"),
  ggplot(rbind(exp_multi_roc,exp_uni_roc) ,
       aes(x=fpr,y=tpr,col=metric.code)) + 
  geom_line() + facet_grid(res~variate) + theme_bw() + ggtitle("b) expanding computation"),
  ncol=2,common.legend = T)



ggplot(rbind(roll_multi_roc,roll_uni_roc,ewsnet_roc,exp_multi_roc,exp_uni_roc),
       aes(x=computation,y=auc)) + 
  geom_text(aes(label=metric.code)) + facet_grid(variate~res) + theme_bw()

