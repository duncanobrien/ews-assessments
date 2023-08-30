##########################################################################################
# Preamble
##########################################################################################
require(ROCR)
require(caret)
require(Matrix)
require(tidyverse)
source("Code/extract_ews_pred_fn.R")

load(file = "ews_raw_data_a.RData")
load(file = "ews_raw_data_b.RData")
load(file = "ews_raw_data_c.RData")
exp_unicomp <- rbind(exp_uni_phyto,exp_uni_zoo) #merge separated expanding univariate EWSs as file too large for Github

load(file = "Results/ews_raw_data_ewsnet.RData")

##########################################################################################
# Extract EWS success
##########################################################################################

lake_outcome_troph <- subset(ewsnet_comp, select = c(lake,res,troph_level)) |>
  dplyr::mutate(outcome = case_when(
    grepl("Kinneret", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Kasumigaura", lake) & troph_level == "zooplankton" ~ "trans",
    grepl("Washington", lake) & troph_level == "phytoplankton" ~ "trans",
    grepl("Monona", lake) & troph_level == "zooplankton" ~ "trans",
    TRUE ~ "no.trans"
  )) |>
  dplyr::mutate(reference = paste(lake,res,troph_level,sep="_"))|>
  dplyr::select(reference,outcome)|>
  dplyr::distinct(reference,.keep_all = T)

suc_ewsnet_max_comp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(ewsnet_comp,detrend_meth == x & deseason_meth == j)|>
                       dplyr::select(-c(detrend_meth,deseason_meth)),
                     sensitivity = NULL,
                     outcome = lake_outcome_troph,
                     method = "ML") |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_exp_multicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    out.mth <- extract_ews_pred(ews.data = subset(exp_multicomp,detrend_meth == x & deseason_meth == j & res == "Monthly")|>
                                  select(-c(detrend_meth,deseason_meth)),
                                sensitivity = 2,
                                outcome = lake_outcome_troph,
                                method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    out.yr <- extract_ews_pred(ews.data = subset(exp_multicomp,detrend_meth == x & deseason_meth == j & res == "Yearly")|>
                                 select(-c(detrend_meth,deseason_meth)),
                               sensitivity = 1,
                               outcome = lake_outcome_troph,
                               method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    return(rbind(out.mth,out.yr))
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_exp_unicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    out.mth <- extract_ews_pred(ews.data = subset(exp_unicomp,detrend_meth == x & deseason_meth == j & res == "Monthly")|>
                                  select(-c(detrend_meth,deseason_meth)),
                                sensitivity = 2,
                                outcome = lake_outcome_troph,
                                method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    out.yr <- extract_ews_pred(ews.data = subset(exp_unicomp,detrend_meth == x & deseason_meth == j & res == "Yearly")|>
                                 select(-c(detrend_meth,deseason_meth)),
                               sensitivity = 1,
                               outcome = lake_outcome_troph,
                               method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
    
    return(rbind(out.mth,out.yr))
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_roll_perm_multicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(perm_roll_multicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.05,
                     outcome = lake_outcome_troph,
                     method = "rolling",
                     surrogate = TRUE) |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_roll_perm_unicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(perm_roll_unicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.05,
                     outcome = lake_outcome_troph,
                     method = "rolling",
                     surrogate = TRUE) |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

##########################################################################################
# Merge dataframes
##########################################################################################

all_ews_data <- as.data.table(suc_ewsnet_max_comp) %>%
  .[,.SD[1], by = c("data_source","troph_level","scaling","lake","res","detrend_meth","deseason_meth")] %>%
  .[,c("data_source","troph_level","scaling","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","prob")] %>%
  setnames(old = "scaling",new= "metric.code") %>%
  setnames(old = "prob",new= "metric_score") %>%
  rbind(as.data.table(suc_exp_multicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","prop_run")] %>% setnames(old = "prop_run",new= "metric_score")) %>%
  rbind(as.data.table(suc_exp_unicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","prop_run")] %>% setnames(old = "prop_run",new= "metric_score")) %>%
  rbind(as.data.table(suc_roll_perm_multicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","obs")] %>% setnames(old = "obs",new= "metric_score")) %>%
  rbind(as.data.table(suc_roll_perm_unicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth","obs")] %>% setnames(old = "obs",new= "metric_score")) %>%
  .[,prediction := ifelse(prediction %in% c("match","prior"),
                          1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML","multivariate_rolling","multivariate_expanding"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,variate := factor(variate, levels = c("univariate","multivariate"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]

computation_data <- all_ews_data %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()

##########################################################################################
# Fit ROCs
##########################################################################################
data_grid <- expand.grid("method_code" = unique(computation_data$method_code),
                         "res" = unique(computation_data$res),
                         "detrend_meth" = unique(computation_data$detrend_meth),
                         "deseason_meth" = unique(computation_data$deseason_meth))

ews_roc <- lapply(1:NROW(data_grid),function(x){
  
  dat <- subset(computation_data,
                 method_code == data_grid$method_code[x] &
                 detrend_meth == data_grid$detrend_meth[x] &
                 deseason_meth == data_grid$deseason_meth[x] &
                 res == data_grid$res[x]) |>
    dplyr::mutate(metric.code = as.character(metric.code)) |>
    tidyr::drop_na(metric_score)
  
  metrics <- c("all",unique(dat$metric.code))
  out.dat <- lapply(metrics,function(j){
    
    if(j != "all"){
      sub_dat <- subset(dat,metric.code == j) |>
        dplyr::mutate(prediction = case_when(prediction == 1 & outcome == "trans" ~ "trans",
                                      prediction != 1 & outcome == "trans" ~ "no.trans",
                                      prediction == 1 & outcome == "no.trans" ~ "no.trans",
                                      prediction != 1 & outcome == "no.trans" ~ "trans")) |>
        dplyr::mutate(dplyr::across(c(prediction,outcome), ~ifelse(.x == "trans",1,0))) |>         #convert from correct prediction flag to actual prediction label
        dplyr::mutate(dplyr::across(c(prediction,outcome),~ factor(.x,levels = c(0,1))))
      }else{
      sub_dat <- dat |>
        dplyr::mutate(metric.code = "all") |>
        dplyr::mutate(prediction = case_when(prediction == 1 & outcome == "trans" ~ "trans",
                                      prediction != 1 & outcome == "trans" ~ "no.trans",
                                      prediction == 1 & outcome == "no.trans" ~ "no.trans",
                                      prediction != 1 & outcome == "no.trans" ~ "trans"))|>
        dplyr::mutate(dplyr::across(c(prediction,outcome), ~ifelse(.x == "trans",1,0))) |>         #convert from correct prediction flag to actual prediction label
        dplyr::mutate(dplyr::across(c(prediction,outcome),~ factor(.x,levels = c(0,1))))
    }
      
      pred <- tryCatch(ROCR::prediction(sub_dat$metric_score, as.numeric(sub_dat$outcome)),
                       error = function(err){NA})
      # pred <- tryCatch(ROCR::prediction(as.numeric(sub_dat$prediction), as.numeric(sub_dat$outcome)),
      #                  error = function(err){NA})
      if(!inherits(pred,"prediction")){
        out <- data.frame("fpr" = NA,"tpr" = NA,
                          "auc" = NA, "pr_auc" = NA,
                          "f1_binary" = NA,
                          "f1_cont" = NA,
                          "bal_acc" = NA) |>
          cbind(sub_dat |> dplyr::select(method,computation,variate,method_code,detrend_meth,deseason_meth,res, metric.code) |>
                  dplyr::distinct()) 
      }else{
      auc <- ROCR::performance(pred,measure = "auc")@y.values[[1]]
      pr_auc <- ROCR::performance(pred,measure = "aucpr")@y.values[[1]]
      
      per <- ROCR::performance(pred,"tpr","fpr")
      f1 <- na.omit(tail(ROCR::performance(pred,"f")@y.values[[1]],n=1))

      c_mat <-caret::confusionMatrix(data = sub_dat$prediction,
                             reference = sub_dat$outcome)
      
      out <- data.frame("fpr" = per@x.values[[1]],"tpr" = per@y.values[[1]],
                        "auc" = auc, "pr_auc" = pr_auc,
                        "f1_binary" = c_mat$byClass[["F1"]],
                        "f1_cont" = f1,
                        "bal_acc" = c_mat$byClass[["Balanced Accuracy"]]) |>
        cbind(sub_dat |> dplyr::select(method,computation,variate,method_code,detrend_meth,deseason_meth,res, metric.code) |>
                dplyr::distinct()) 
      }
      return(out)
    
  }) |>
    data.table::rbindlist()
  return(out.dat)
}) |>
  data.table::rbindlist()

saveRDS(ews_roc,"Results/ROC/ews_roc.rds")

##########################################################################################
# F1-statistic + AUC equivalents of Figures 3 and 4
##########################################################################################

plot_data <- rbind(subset(ews_roc, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding"),
                   subset(ews_roc, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML")) |>
  rbind(subset(ews_roc, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding"),
        subset(ews_roc, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML")) |>
  dplyr::mutate(variate = ifelse(grepl("ML",method_code),"EWSNet",ifelse(grepl("multivariate",method_code),"Multivariate","Univariate")),
                res = factor(res,levels = c("Monthly","Yearly"))) |>
  dplyr::select(metric.code,method_code,f1_cont,f1_binary,auc,pr_auc,bal_acc,res,variate) |>
  dplyr::distinct()

auc_fig3 <- ggplot(subset(plot_data,metric.code == "all"),
                   aes(x=auc,y=method_code,col=variate))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(size = 4) +
  facet_grid(~res) +
  scale_color_manual(values=c("#529928","#5d3099","#bfbd3d"))+
  labs(x="Area Under Curve", y = "Early warning signal method") +
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

auc_fig4 <- ggplot(subset(plot_data,metric.code != "all") |>
         dplyr::mutate(reference = paste(variate,method_code,sep = "_")),
       aes(x=auc,y=metric.code))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(data = subset(plot_data |>
                             dplyr::mutate(method_code = as.character(method_code),
                                           method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                           reference = paste(variate,method_code,sep = "_")),
                           method_code == "EWSNet" & metric.code != "all"), 
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  geom_point(data= subset(plot_data |>
                            mutate(method_code = as.character(method_code),
                                   method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                   reference = paste(variate,method_code,sep = "_")
                            ), method_code != "EWSNet" & metric.code != "all"),
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  scale_color_manual(values=c("#529928","#160B24","#5d3099","#403F14","#bfbd3d"),
                     labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling","Univariate\nexpanding","Univariate\nrolling"),
                     guide = "none")+
  scale_fill_manual(values=c("#529928","#160B24","#5d3099","#403F14","#529928","#bfbd3d"),
                    labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling", "Univariate\nexpanding","Univariate\nrolling"))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  labs(x="Area Under Curve", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation")+
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("Results/ROC/AUC_fig.pdf",
       ggpubr::ggarrange(auc_fig3 + ggtitle("Figure 5"),
                         auc_fig4 + ggtitle("Figure 6"),
                         labels="AUTO",font.label = list(face="plain"),
                         common.legend = F,ncol=2,widths = c(1,1.5),heights = c(0.75,1)),
       width = 8,height = 5)

f1_fig3 <- ggplot(subset(plot_data,metric.code == "all"),
       aes(x=f1_binary,y=method_code,col=variate))  +
  geom_point(size=4) +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  facet_grid(~res) +
  scale_color_manual(values=c("#529928","#5d3099","#bfbd3d"))+
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  labs(x="F1 statistic", y = "Early warning signal method") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

f1_fig4 <- ggplot(subset(plot_data,metric.code != "all"),
                  aes(x=f1_binary,y=metric.code))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(data = subset(plot_data |>
                             dplyr::mutate(method_code = as.character(method_code),
                                           method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                           reference = paste(variate,method_code,sep = "_")),
                           method_code == "EWSNet" & metric.code != "all"), 
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  geom_point(data= subset(plot_data |>
                            mutate(method_code = as.character(method_code),
                                   method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                   reference = paste(variate,method_code,sep = "_")
                            ), method_code != "EWSNet" & metric.code != "all"),
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  scale_color_manual(values=c("#529928","#160B24","#5d3099","#403F14","#bfbd3d"),
                     labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling","Univariate\nexpanding","Univariate\nrolling"),
                     guide = "none")+
  scale_fill_manual(values=c("#529928","#160B24","#5d3099","#403F14","#529928","#bfbd3d"),
                    labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling", "Univariate\nexpanding","Univariate\nrolling"))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  labs(x="F1 statistic", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation")+
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave("Results/ROC/F1_fig.pdf",
       ggpubr::ggarrange(f1_fig3 + ggtitle("Figure 5"),
                         f1_fig4 + ggtitle("Figure 6"),
                         labels="AUTO",font.label = list(face="plain"),
                         common.legend = F,ncol=2,nrow=1,widths = c(1,1.5)),
       width = 8,height = 5)

bal_acc_fig3 <- ggplot(subset(plot_data,metric.code == "all"),
                      aes(x=bal_acc,y=method_code,col=variate))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(size = 4) +
  facet_grid(~res) +
  scale_color_manual(values=c("#529928","#5d3099","#bfbd3d"))+
  labs(x="Balanced accuracy", y = "Early warning signal method") +
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

bal_acc_fig4 <- ggplot(subset(plot_data,metric.code != "all") |>
                        dplyr::mutate(reference = paste(variate,method_code,sep = "_")),
                      aes(x=bal_acc,y=metric.code))  +
  geom_vline(xintercept = 0.5,linetype = "dashed", colour="grey50") + 
  geom_point(data = subset(plot_data |>
                             dplyr::mutate(method_code = as.character(method_code),
                                           method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                           reference = paste(variate,method_code,sep = "_")),
                           method_code == "EWSNet" & metric.code != "all"), 
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  geom_point(data= subset(plot_data |>
                            mutate(method_code = as.character(method_code),
                                   method_code = ifelse(grepl("ML",method_code),"EWSNet",method_code),
                                   reference = paste(variate,method_code,sep = "_")
                            ), method_code != "EWSNet" & metric.code != "all"),
             aes(fill = reference),size = 4, shape = 21,alpha = 0.8) +
  scale_color_manual(values=c("#529928","#160B24","#5d3099","#403F14","#bfbd3d"),
                     labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling","Univariate\nexpanding","Univariate\nrolling"),
                     guide = "none")+
  scale_fill_manual(values=c("#529928","#160B24","#5d3099","#403F14","#bfbd3d"),
                    labels = c("EWSNet","Multivariate\nexpanding","Multivariate\nrolling", "Univariate\nexpanding","Univariate\nrolling"))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  scale_x_continuous(breaks = c(0.2,0.5,0.8),limits = c(0,1)) +
  labs(x="Balanced accuracy", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation")+
  theme_bw() +
  theme(#legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())

ggsave("Results/ROC/bal_acc_fig.pdf",
       ggpubr::ggarrange(bal_acc_fig3 + ggtitle("Figure 5"),
                         bal_acc_fig4 + ggtitle("Figure 6"),
                         labels="AUTO",font.label = list(face="plain"),
                         common.legend = F,ncol=2,widths = c(1,1.5),heights = c(0.75,1)),
       width = 8,height = 5)
