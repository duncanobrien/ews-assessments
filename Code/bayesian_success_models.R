##########################################################################################
# Preamble
##########################################################################################

require(Matrix)
require(tidybayes)
require(ggplot2)
require(data.table)
require(tidyverse)
require(brms)
require(patchwork)

load(file = "Results/ews_raw_data_a.RData")
load(file = "Results/ews_raw_data_b.RData")
load(file = "Results/ews_raw_data_c.RData")
exp_unicomp <- rbind(exp_uni_phyto,exp_uni_zoo) #merge separated expanding univariate EWSs as file too large for Github

load(file = "Results/ews_raw_data_ewsnet.RData")
source("Code/extract_ews_pred_fn.R")

##########################################################################################
# Wrangle Data
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
                       select(-c(detrend_meth,deseason_meth)),
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
                     sensitivity = 0.5,
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
                     sensitivity = 0.5,
                     outcome = lake_outcome_troph,
                     method = "rolling",
                     surrogate = TRUE) |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()


#all_ews_data <- as.data.table(suc_ewsnet_comp) %>%
all_ews_data <- as.data.table(suc_ewsnet_max_comp) %>%
  .[,.SD[1], by = c("data_source","troph_level","scaling","lake","res","detrend_meth","deseason_meth")] %>%
  .[,c("data_source","troph_level","scaling","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")] %>%
  setnames(old = "scaling",new= "metric.code") %>%
  rbind(as.data.table(suc_exp_multicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  rbind(as.data.table(suc_exp_unicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  #rbind(as.data.table(suc_roll_multicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  #rbind(as.data.table(suc_roll_unicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  rbind(as.data.table(suc_roll_perm_multicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  rbind(as.data.table(suc_roll_perm_unicomp)[,.SD[1], by = c("data_source","troph_level","metric.code","lake","res","detrend_meth","deseason_meth")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction","detrend_meth","deseason_meth")]) %>%
  .[,prediction := ifelse(prediction %in% c("match","prior"),
                       1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML","multivariate_rolling","multivariate_expanding"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,variate := factor(variate, levels = c("univariate","multivariate"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]

ts_balance <- copy(all_ews_data) %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,data_ref := paste(data_source,lake,sep ="_")] %>%
  .[,.(data_ref = length(unique(data_ref))),by="outcome"] 

overall_table <- lapply(unique(paste(all_ews_data$detrend_meth,all_ews_data$deseason_meth,
                                     sep = "_")), function(x){
                                       
                                       lapply(unique(all_ews_data$res),function(k){        
                                         gtsummary::tbl_strata(subset(all_ews_data,
                                                                      detrend_meth == gsub("_.*","",x) & 
                                                                        deseason_meth ==  gsub(".*_","",x) &
                                                                        res == k) |>
                                                                 dplyr::select(method_code,prediction) |>
                                                                 dplyr::mutate(prediction = paste(prediction)),
                                                               method_code,.tbl_fun =tbl_summary) |>
                                           as.data.frame() %>%
                                           rbind(colnames(.),.) |>
                                           t() |>
                                           as.data.frame() |>
                                           dplyr::slice(-1) |>
                                           `colnames<-`(c("Number of Trials","drop","Number of failures","Number of successes")) |>
                                           dplyr::select(-drop) |>
                                           dplyr::mutate(`Method code` = levels(all_ews_data$method_code),
                                                         Resolution = k, .before = `Number of Trials`) |>
                                           dplyr::mutate(`Detrend method` = gsub("_.*","",x), 
                                                         `Deseason method` = gsub(".*_","",x),.before = `Number of Trials`) |>
                                           dplyr::mutate(`Number of Trials` = gsub("\\D+","",`Number of Trials`))
                                       }) |> data.table::rbindlist()
                                     }) |> data.table::rbindlist() |>
  dplyr::arrange(Resolution,`Method code`)

metric_table <- lapply(unique(paste(all_ews_data$detrend_meth,all_ews_data$deseason_meth,
                                    sep = "_")), function(x){
                                      
                                      lapply(paste(unique(all_ews_data$method_code),rep(unique(all_ews_data$res),5),sep = ":"),function(k){
                                        
                                        sub_dat <- subset(all_ews_data,
                                                          detrend_meth == gsub("_.*","",x) & 
                                                            deseason_meth ==  gsub(".*_","",x) &
                                                            method_code == gsub(":.*","",k) &
                                                            res == gsub(".*:","",k)) |>
                                          dplyr::select(metric.code,prediction) |>
                                          dplyr::mutate(prediction = paste(prediction))
                                        
                                        gtsummary::tbl_strata(sub_dat,
                                                              metric.code,.tbl_fun =tbl_summary) |>
                                          as.data.frame() %>%
                                          rbind(colnames(.),.) |>
                                          t() |>
                                          as.data.frame() |>
                                          dplyr::slice(-1) |>
                                          `colnames<-`(c("Number of Trials","drop","Number of failures","Number of successes")) |>
                                          dplyr::select(-drop) |>
                                          dplyr::mutate(Indicator = unique(sub_dat$metric.code), .before = `Number of Trials`) |>
                                          dplyr::mutate(`Method code` = gsub(":.*","",k),
                                                        Resolution =  gsub(".*:","",k), 
                                                        .before = `Number of Trials`) |>
                                          dplyr::mutate(`Detrend method` = gsub("_.*","",x), 
                                                        `Deseason method` = gsub(".*_","",x),.before = `Number of Trials`) |>
                                          dplyr::mutate(`Number of Trials` = gsub("\\D+","",`Number of Trials`))
                                      }) |> data.table::rbindlist()
                                    }
) |> data.table::rbindlist() |>
  dplyr::arrange(Resolution,Indicator)

computation_data <- all_ews_data %>%
  .[,reference := paste(lake,res,troph_level,sep="_")] %>% #trials in terms of assessed time series
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,offset := length(unique(data_source)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth") ] %>% #trials in terms of assessed time series
  .[,offset2 := length(unique(data_source))*length(unique(metric.code)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth")] %>% #trials in terms of assessed time series AND metrics 
  .[,.(total_success = sum(prediction),
       offset = unique(offset2),
       #ts_length = unique(ts_length),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_"),
       #weights = unique(weights),
       outcome = unique(outcome),
       #weights = unique(ifelse(outcome == "trans",(ts_balance$data_ref[1]),(ts_balance$data_ref[2]))),
       weights = unique(ifelse(outcome == "trans",(ts_balance$data_ref[2]),(ts_balance$data_ref[1])))
       ),
    by = c("lake", "res", "method_code","troph_level","detrend_meth","deseason_meth")] %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()

metric_data <- all_ews_data %>%
  .[,offset := length(unique(data_source)), by =c("metric.code","lake","res","troph_level","detrend_meth","deseason_meth","method_code") ] %>% #trials in terms of assessed time series
  .[,.(total_success = sum(prediction),
       offset = unique(offset),
       #ts_length = unique(ts_length),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_")),by = c("lake", "res", "metric.code","troph_level","detrend_meth","deseason_meth","method_code")] %>%
  .[,indicator := paste(metric.code,sub(".*_", "",method_code),sep = "_")] %>%
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()

its <- 10000
thn <- 0.0005*its
wrmup <- 0.1*its

##########################################################################################
# Compare Detrend/Deseason Methods
##########################################################################################

ews_mod_detrend_mth_uni_roll <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                     data = computation_data[computation_data$res == "Monthly" & grepl("univariate_rolling",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                       mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_rolling_none_none"))
                                     ,
                                     iter = 10000,
                                     thin =thn,
                                     warmup = 2000,
                                     prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                              prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                     family = binomial(link = "logit"), 
                                     chains = 4,
                                     backend = "cmdstanr",
                                     control = list(adapt_delta =0.99, max_treedepth = 20,step_size = 0.01),
                                     seed = 12345, cores = 4,sample_prior = TRUE)
#linear_stl
ews_mod_detrend_mth_uni_exp <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                          data = computation_data[computation_data$res == "Monthly" & grepl("univariate_expanding",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                           mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_expanding_none_none"))
                                         ,
                                          iter = 10000,
                                          thin =thn,
                                          warmup = 2000,
                                          prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                                   prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                          family = binomial(link = "logit"), 
                                          chains = 4,
                                         backend = "cmdstanr",
                                          control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                          seed = 12345, cores = 4,sample_prior = TRUE)
#linear-decompose
ews_mod_detrend_mth_multi_roll <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                            data = computation_data[computation_data$res == "Monthly" & grepl("multivariate_rolling",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                              mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "multivariate_rolling_none_none"))
                                            ,
                                            iter = 10000,
                                            thin =thn,
                                            warmup = 2000,
                                            prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                                     prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                            family = binomial(link = "logit"), 
                                            chains = 4,
                                            backend = "cmdstanr",
                                            control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                            seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian_none
ews_mod_detrend_mth_multi_exp <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                      data = computation_data[computation_data$res == "Monthly" & grepl("multivariate_expanding",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                        mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "multivariate_expanding_none_none"))
                                      ,
                                      iter = 10000,
                                      thin =thn,
                                      warmup = 2000,
                                      prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                               prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                      family = binomial(link = "logit"), 
                                      chains = 4,
                                      backend = "cmdstanr",
                                      control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                      seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian-average
ews_mod_detrend_mth_ml <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                      data = computation_data[computation_data$res == "Monthly" & grepl("ML",computation_data$method_code)  & computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington"),] |>
                                      mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_ML_none_none"))
                                    ,
                                      iter = 10000,
                                      thin =thn,
                                      warmup = 2000,
                                      prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                               prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                      family = binomial(link = "logit"), 
                                      chains = 4,
                                      backend = "cmdstanr",
                                      control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                      seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian-none
saveRDS(ews_mod_detrend_mth_uni_roll,file = "Results/supplementary_info/ews_mod_detrend_mth_uni_roll.rds")
saveRDS(ews_mod_detrend_mth_uni_exp,file = "Results/supplementary_info/ews_mod_detrend_mth_uni_exp.rds")
saveRDS(ews_mod_detrend_mth_multi_roll,file = "Results/supplementary_info/ews_mod_detrend_mth_multi_roll.rds")
saveRDS(ews_mod_detrend_mth_multi_exp,file = "Results/supplementary_info/ews_mod_detrend_mth_multi_exp.rds")
saveRDS(ews_mod_detrend_mth_ml,file = "Results/supplementary_info/ews_mod_detrend_mth_ml.rds")

dat_detrend_trials_alt <-  ews_mod_detrend_mth_uni_roll |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_combo_code", "", .variable),
         .variable = gsub("b_Intercept", "univariate_rolling_none_none", .variable) ) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  rbind(ews_mod_detrend_mth_uni_exp |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_combo_code", "", .variable),
                 .variable = gsub("b_Intercept", "univariate_expanding_none_none", .variable) ) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5))) |>
  rbind(ews_mod_detrend_mth_multi_roll |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_combo_code", "", .variable),
                 .variable = gsub("b_Intercept", "multivariate_rolling_none_none", .variable) ) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5))) |>
  rbind(ews_mod_detrend_mth_multi_exp |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_combo_code", "", .variable),
                 .variable = gsub("b_Intercept", "multivariate_expanding_none_none", .variable) ) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5))) |>
  rbind(ews_mod_detrend_mth_ml |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_combo_code", "", .variable),
                 .variable = gsub("b_Intercept", "univariate_ML_none_none", .variable) ) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5))) |>
mutate(res = "Monthly",
         reference = ifelse(grepl("univariate_rolling_none_none|multivariate_rolling_none_none",.variable),1,0)) |>
  separate(.variable,c("variate","computation","detrend_meth","deseason_meth"),sep = "_")

ggplot(subset(dat_detrend_trials_alt,.width == 0.95),
       aes(y = computation, x = .value)) +
  geom_rect(data=subset(dat_detrend_trials_alt, .width == 0.95 & reference == 1), 
            aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.5)+
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  ggh4x::facet_nested(res + detrend_meth~variate + deseason_meth,scales = "free_y") +
  
  
  geom_errorbar(aes(xmin=.lower, xmax=.upper,col=computation),
                width=0, alpha=0.9, size=1.3)+
  geom_point(size=2,aes(fill=computation),shape=21) +
  scale_fill_manual(values=c("#B4D5DE",
                             "#DEC5BF",
                             "#917C67"))+
  scale_color_manual(values=c("#B4D5DE",
                              "#DEC5BF",
                              "#917C67"))+
  labs(x="Posterior estimate", y = "Processing method") +
  coord_cartesian(xlim = c(-1.5,1.5))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank())


##########################################################################################
# Compare Computation Methods
##########################################################################################
#ews_mod_method_mth<- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 + (1||lake/outcome) ), 
ews_mod_method_mth <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 +(1|ID1|lake) + (1|ID2|outcome) + (1|ID3|lake:outcome)), 
#ews_mod_method_mth <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 + offset(weights) + outcome ), 
                                data = rbind(subset(computation_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
                                             subset(computation_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
                                             subset(computation_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling"),
                                             subset(computation_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding"),
                                             subset(computation_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML"))|>
                                      dplyr::mutate(outcome = factor(outcome,levels = c("trans","no.trans"))),
                                iter = its,
                                thin = thn,
                                warmup = wrmup,
                                prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                         prior(normal(1, 1),class = sd)),
                                #prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5)),
                                family = binomial(link = "logit"), 
                                chains = 4,
                                backend = "cmdstanr",
                                control = list(adapt_delta = .99, max_treedepth = 20),
                                seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_method_mth,file = "Results/ews_models/computation_models/ews_mod_method_mth.rds")

ews_mod_method_yr <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 +  (1|ID1|lake) + (1|ID2|outcome) + (1|ID3|lake:outcome)), 
#ews_mod_method_yr<- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 + offset(weights) + outcome), 
#ews_mod_method_yr<- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 +  (1|lake/outcome)), 
                              data =  rbind(subset(computation_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling"),
                                            subset(computation_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding"),
                                            subset(computation_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling"),
                                            subset(computation_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding"),
                                            subset(computation_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML")) |>
                                dplyr::mutate(outcome = factor(outcome,levels = c("trans","no.trans"))),
                              iter = its,
                              thin = thn,
                              warmup = wrmup,
                              prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                       prior(normal(1, 1),class = sd)),
                              #prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5)),
                              family = binomial(link = "logit"), 
                              chains = 4,
                              backend = "cmdstanr",
                              control = list(adapt_delta = .99, max_treedepth = 20),
                              seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_method_yr,file = "Results/ews_models/computation_models/ews_mod_method_yr.rds")

dat_method_trials <-  ews_mod_method_mth |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_method_code", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate",ifelse(grepl("ML",.variable),"EWSNet","univariate")),
         res = "Monthly") |>
  rbind(ews_mod_method_yr |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_method_code", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate",ifelse(grepl("ML",.variable),"EWSNet","univariate")),
                 res = "Yearly")) |>
    mutate(.variable = forcats::fct_reorder(.variable,.value,.fun = mean,.desc = FALSE)) #reorder variable in to increasing ability
    
dat_method_halfeye <- ews_mod_method_mth |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_method_code", "", .variable)) |> 
  mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate",ifelse(grepl("ML",.variable),"EWSNet","univariate")),
         res = "Monthly") |>
  rbind( ews_mod_method_yr |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_method_code", "", .variable)) |> 
           mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate",ifelse(grepl("ML",.variable),"EWSNet","univariate")),
                  res = "Yearly")) |>
  mutate(.variable = factor(.variable, levels = levels(dat_method_trials$.variable))) #reorder variable in to increasing average ability based on dat_method_trials


ggsave(ggplot(dat_method_trials,aes(y = .variable, x = .value)) +
         geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
         tidybayes::stat_slab(data= dat_method_halfeye, alpha=0.5, aes(fill=variate)) +
         # geom_errorbar(data = subset(dat_metric_trials,.width == 0.95), 
         #               aes(xmin=.lower, xmax=.upper,col=variate), 
         #                width=0, alpha=0.9, size=1.3)+
         # geom_point(size=4,aes(fill=variate),shape=21) +
         scale_fill_manual(values=c("#529928",
                                    "#5d3099",
                                    "#bfbd3d"))+
         scale_color_manual(values=c("#529928",
                                     "#5d3099",
                                     "#bfbd3d"))+
         tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                                       interval_size_range = c(0.8, 2)
                                       ,colour = "black" ) +
         labs(x="Probability of correct prediction", y = "Early warning signal method") +
         #xlim(-3,3)+
         scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
         scale_y_discrete()+
         coord_cartesian(xlim = c(-3,3))+
         facet_wrap(~res)+
         theme_bw()+
         theme(legend.position = "none",
               panel.grid.minor = element_blank(),
               #panel.grid.major = element_line(colour = "grey80"),
               #plot.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF'),
               #panel.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF')
               ),
       filename = "Figures/figure_3.pdf",width = 6,height = 5)

##########################################################################################
# Compare Indicators' True Positive Rate
##########################################################################################

ews_mod_ind_mth_true <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                     data =  rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "trans"),
                                                   subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "trans"),
                                                   subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling" & outcome == "trans"),
                                                   subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding" & outcome == "trans"),
                                                   subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                     iter = its,
                                     thin = thn,
                                     warmup = wrmup,
                                     prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                              prior(normal(1, 1),class = sd)),
                                     family = binomial(link = "logit"), 
                                     chains = 4,
                                     backend = "cmdstanr",
                                     control = list(adapt_delta = .975, max_treedepth = 20),
                                     seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_ind_mth_true,file = "Results/ews_models/indicator_models/ews_mod_ind_mth_true.rds")

ews_mod_ind_yr_true<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 +(1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                   data =rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling" & outcome == "trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding" & outcome == "trans"),
                                               subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   backend = "cmdstanr",
                                   control = list(adapt_delta = .975, max_treedepth = 20),
                                   seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_ind_yr_true,file = "Results/ews_models/indicator_models/ews_mod_ind_yr_true.rds")

dat_ind_true_trials <-  ews_mod_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_ind_yr_true |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = forcats::fct_reorder(.variable,.value,.fun = mean,.desc = FALSE))  #reorder variable in to increasing ability

dat_ind_true_halfeye <- ews_mod_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_ind_yr_true |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))

ind_true_plot3 <- ggplot(data = dat_ind_true_trials |>   
                           mutate(.variable = sub("_.*", "",.variable),
                                  method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                  reference = paste(variate,method_code,sep = "_")),
                         aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye |>
                                      mutate( method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")),
                                    method_code == "EWSNet"), alpha=0.5, aes(fill=reference),normalize = "panels") +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye |>
                                      mutate( method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")
                                      ), variate != "EWSNet"),
                       aes(fill=reference,group=method_code),alpha=0.65,position = position_dodge(width=0.65),normalize = "panels") +
  scale_fill_manual(values=c("#529928", "#160B24","#5d3099","#403F14","#bfbd3d"),
                    labels = c("EWSNet", "Multivariate\nexpanding","Multivariate\nrolling",  "Univariate\nexpanding","Univariate\nrolling"))+
  tidybayes::geom_pointinterval(data = ~subset(.,method_code %in% c("expanding","rolling")),aes(xmin = .lower, xmax = .upper,group = method_code),colour = "black",position = position_dodge(width=0.65),interval_size_range = c(0.4, 1.2)) +
  tidybayes::geom_pointinterval(data = ~subset(.,method_code == "EWSNet"),aes(xmin = .lower, xmax = .upper),colour = "black",interval_size_range = c(0.4, 1.2)) +
  labs(x="True positive prediction probability", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation") +
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4,4))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines')) 


##########################################################################################
# Compare Indicators' True Negative Rate
##########################################################################################

ews_mod_ind_mth_false <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                      data =  rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "no.trans"),
                                                    subset(metric_data, res  == "Monthly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "no.trans"),
                                                    subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling" & outcome == "no.trans"),
                                                    subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding" & outcome == "no.trans"),
                                                    subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                      iter = its,
                                      thin = thn,
                                      warmup = wrmup,
                                      prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                               prior(normal(1, 1),class = sd)),
                                      family = binomial(link = "logit"), 
                                      chains = 4,
                                      backend = "cmdstanr",
                                      control = list(adapt_delta = .99, max_treedepth = 20),
                                      seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_ind_mth_false,file = "Results/ews_models/indicator_models/ews_mod_ind_mth_false.rds")

ews_mod_ind_yr_false<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|ID1|lake) + (1|ID2|method_code) + (1|ID3|lake:method_code)), 
                                    data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "stl" & method_code == "univariate_rolling" & outcome == "no.trans"),
                                                 subset(metric_data, res  == "Yearly" & detrend_meth == "linear" & deseason_meth == "decompose" & method_code == "univariate_expanding" & outcome == "no.trans"),
                                                 subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "multivariate_rolling" & outcome == "no.trans"),
                                                 subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code == "multivariate_expanding" & outcome == "no.trans"),
                                                 subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                    iter = its,
                                    thin = thn,
                                    warmup = wrmup,
                                    prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                             prior(normal(1, 1),class = sd)),
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    backend = "cmdstanr",
                                    control = list(adapt_delta = .99, max_treedepth = 20,step_size = 0.01),
                                    seed = 12345, cores = 4,sample_prior = TRUE)
saveRDS(ews_mod_ind_yr_false,file = "Results/ews_models/indicator_models/ews_mod_ind_yr_false.rds")

dat_ind_false_trials <-  ews_mod_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_ind_yr_false |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable,levels = levels(dat_ind_true_trials$.variable)))  #reorder variable in to increasing ability based upon true postive results

dat_ind_false_halfeye <- ews_mod_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_ind_yr_false |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))

ind_false_plot3 <- ggplot(data = dat_ind_false_trials |>   
         mutate(.variable = sub("_.*", "",.variable),
                method_code = ifelse(method_code == "ML","EWSNet",method_code),
                reference = paste(variate,method_code,sep = "_")),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye |>
                                      mutate( .variable = sub("_.*", "",.variable),
                                              method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")),
                                    method_code == "EWSNet"), alpha=0.5, aes(fill=reference),normalize = "panels") +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye |>
                                      mutate( .variable = sub("_.*", "",.variable),
                                              method_code = ifelse(method_code == "ML","EWSNet",method_code),
                                              reference = paste(variate,method_code,sep = "_")
                                      ), variate != "EWSNet"),
                       aes(fill=reference,group=method_code),alpha=0.65,position = position_dodge(width=0.65),normalize = "panels") +

  scale_fill_manual(values=c("#529928", "#160B24","#5d3099","#403F14","#bfbd3d"),
                    labels = c("EWSNet", "Multivariate\nexpanding","Multivariate\nrolling",  "Univariate\nexpanding","Univariate\nrolling"))+
  tidybayes::geom_pointinterval(data = ~subset(.,method_code %in% c("expanding","rolling")),aes(xmin = .lower, xmax = .upper,group = method_code),colour = "black",position = position_dodge(width=0.65),interval_size_range = c(0.4, 1.2)) +
  tidybayes::geom_pointinterval(data = ~subset(.,method_code == "EWSNet"),aes(xmin = .lower, xmax = .upper),colour = "black",interval_size_range = c(0.4, 1.2)) +
  labs(x="True negative prediction probability", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation") +
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4,4))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines')) 

ggsave(ind_true_plot3 +
  ind_false_plot3 +
  patchwork::plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A'),
  filename = "Figures/figure_4.pdf",width = 8,height = 6)
