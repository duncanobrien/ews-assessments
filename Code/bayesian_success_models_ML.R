##########################################################################################
# Preamble
##########################################################################################
require(Matrix)
require(tidybayes)
require(ggplot2)
require(data.table)
require(tidyverse)
require(brms)

load(file = "Results/ews_raw_data_a.RData")
load(file = "Results/ews_raw_data_b.RData")
load(file = "Results/ews_raw_data_c.RData")
exp_unicomp <- rbind(exp_uni_phyto,exp_uni_zoo) #merge separated expanding univariate EWSs as file too large for Github

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


suc_ewsnet_comp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(ewsnet_comp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.7,
                     outcome = lake_outcome_troph,
                     method = "ML") |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

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
    
    extract_ews_pred(ews.data = subset(exp_multicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 2,
                     outcome = lake_outcome_troph,
                     method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()


suc_exp_unicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(exp_unicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 2,
                     outcome = lake_outcome_troph,
                     method = "expanding") |>
      mutate(detrend_meth = x,deseason_meth = j)
  }) |>
    data.table::rbindlist()
}) |>
  data.table::rbindlist()

suc_roll_multicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(roll_multicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.5,
                     outcome = lake_outcome_troph,
                     method = "rolling") |>
      mutate(detrend_meth = x,deseason_meth = j)
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


suc_roll_unicomp <- lapply(c("none","linear","gaussian","loess"),function(x){
  
  lapply(c("none","average","decompose","stl"),function(j){
    
    extract_ews_pred(ews.data = subset(roll_unicomp,detrend_meth == x & deseason_meth == j)|>
                       select(-c(detrend_meth,deseason_meth)),
                     sensitivity = 0.5,
                     outcome = lake_outcome_troph,
                     method = "rolling") |>
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
  .[,outcome := ifelse(prediction %in% c("match","prior"),
                       1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML","multivariate_rolling","multivariate_expanding"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,variate := factor(variate, levels = c("univariate","multivariate"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]


computation_data <- all_ews_data %>%
  .[,offset := length(unique(data_source)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth") ] %>% #trials in terms of assessed time series
  .[,offset2 := length(unique(data_source))*length(unique(metric.code)), by =c("method_code","lake","res","troph_level","detrend_meth","deseason_meth")] %>% #trials in terms of assessed time series AND metrics 
  .[,.(total_success = sum(outcome),
       offset = unique(offset2),
       #ts_length = unique(ts_length),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_")),by = c("lake", "res", "method_code","troph_level","detrend_meth","deseason_meth")] %>%
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  .[,detrend_meth := factor(detrend_meth, levels = c("none","linear","loess","gaussian"))] %>%
  .[,deseason_meth := factor(deseason_meth, levels = c("none","average","decompose","stl"))] %>%
  as.data.frame()

metric_data <- all_ews_data %>%
  .[,offset := length(unique(data_source)), by =c("metric.code","lake","res","troph_level","detrend_meth","deseason_meth","method_code") ] %>% #trials in terms of assessed time series
  .[,.(total_success = sum(outcome),
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
ews.mod.detrend.mth.uni.roll <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|lake/outcome)), 
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
                                     control = list(adapt_delta =0.99, max_treedepth = 20,stepsize = 0.01),
                                     seed = 12345, cores = 4,sample_prior = TRUE)
#linear_stl
ews_mod_detrend_mth_uni_exp <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|lake/outcome)), 
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
                                          control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                          seed = 12345, cores = 4,sample_prior = TRUE)
#none-decompose
ews_mod_detrend_mth_multi_exp <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|lake/outcome)), 
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
                                      control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                      seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian-average
ews_mod_detrend_mth_multi_roll <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|lake/outcome)), 
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
                                           control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                           seed = 12345, cores = 4,sample_prior = TRUE)
#gaussian_none
ews_mod_detrend_mth_ml <- brms::brm(brms::bf(total_success | trials(offset) ~ combo_code +  (1|lake/outcome)), 
                                      data = computation_data[computation_data$res == "Monthly" & grepl("ML",computation_data$method_code)  & !(computation_data$lake %in% c("Kinneret", "Kasumigaura", "Monona", "Washington")),] |>
                                      mutate(combo_code = relevel(factor(paste(method_code,detrend_meth,deseason_meth,sep = "_")), ref = "univariate_ML_none_none"))
                                    ,
                                      iter = 10000,
                                      thin =thn,
                                      warmup = 2000,
                                      prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                               prior(normal(1, 1),class = sd)), #1.2 sd as corresponds to 1% and 99% success rate
                                      family = binomial(link = "logit"), 
                                      chains = 4,
                                      control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                      seed = 12345, cores = 4,sample_prior = TRUE)
#none-none/none-decompose
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
# Compare Indicators (Gaussian detrend and Average deseason)
##########################################################################################

ews_mod_metric_mth_cond<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code/outcome) ), 
                                    data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML"),
                                                 subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML")),
                                    iter = its,
                                    thin = thn,
                                    warmup = wrmup,
                                    prior= c(prior(normal(0, 1.2), class = b),
                                             prior(normal(1, 1),class = sd)),
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    control = list(adapt_delta = .975, max_treedepth = 20),
                                    seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_metric_yr_cond<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code/outcome) ), 
                                   data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML"),
                                                subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   control = list(adapt_delta = .975, max_treedepth = 20),
                                   seed = 12345, cores = 4,sample_prior = TRUE)
#bayestestR::describe_posterior(ews.mod.metric.mth.cond, ci = 0.95, test="none")

dat_metric_trials <-  ews_mod_metric_mth_cond |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"ML","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_metric_yr_cond |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"ML","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable))

dat_metric_halfeye <- ews_mod_metric_mth_cond |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"ML","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_metric_yr_cond |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"ML","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable))

all_metric_plot <- ggplot(dat_metric_trials,aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= dat_metric_halfeye, alpha=0.5, aes(fill=variate)) +
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
  labs(x="Posterior estimate", y = "Early warning signal metric") +
  coord_cartesian(xlim = c(-3.75,3.75))+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  facet_wrap(~res)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank())


##########################################################################################
# Compare Computation Methods
##########################################################################################

ews_mod_method_mth <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 + (1|lake/outcome) ), 
                                data = rbind(subset(computation_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML"),
                                      subset(computation_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML")),
                                iter = its,
                                thin = thn,
                                warmup = wrmup,
                                prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                         prior(normal(1, 1),class = sd)),
                                family = binomial(link = "logit"), 
                                chains = 2,
                                control = list(adapt_delta = .99, max_treedepth = 20),
                                seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_method_yr<- brms::brm(brms::bf(total_success | trials(offset) ~ method_code - 1 + (1|lake/outcome) ), 
                              data =  rbind(subset(computation_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML"),
                                            subset(computation_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML")),
                              iter = its,
                              thin = thn,
                              warmup = wrmup,
                              prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                       prior(normal(1, 1),class = sd)),
                              family = binomial(link = "logit"), 
                              chains = 2,
                              control = list(adapt_delta = .99, max_treedepth = 20),
                              seed = 12345, cores = 4,sample_prior = TRUE)

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
       filename = "/Users/ul20791/Downloads/eg_bayes_fig3.pdf",width = 6,height = 5)

##########################################################################################
# Compare Metrics' True Positive Rate
##########################################################################################

ews_mod_metric_mth_true <- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code) ), 
                                     data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "trans"),
                                                  subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                     iter = its,
                                     thin = thn,
                                     warmup = wrmup,
                                     prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                              prior(normal(1, 1),class = sd)),
                                     family = binomial(link = "logit"), 
                                     chains = 4,
                                     control = list(adapt_delta = .99, max_treedepth = 20),
                                     seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_metric_yr_true<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code) ), 
                                   data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "trans"),
                                                subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   control = list(adapt_delta = .99, max_treedepth = 20),
                                   seed = 12345, cores = 4,sample_prior = TRUE)

dat_metric_true_trials <-  ews_mod_metric_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_metric_yr_true |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable)) |>
  mutate(.variable = forcats::fct_reorder(.variable,.value,.fun = mean,.desc = FALSE))  #reorder variable in to increasing ability

dat_metric_true_halfeye <- ews_mod_metric_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_metric_yr_true |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_metric_true_trials$.variable))) #reorder variable in to increasing average ability based on dat_metric_true_trials


metric_true_plot <- ggplot(dat_metric_true_trials |>   
                             mutate(.variable = sub("_.*", "",.variable)),aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= dat_metric_true_halfeye, alpha=0.5, aes(fill=variate)) +
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
                                interval_size_range = c(0.4, 1.2)
                                ,colour = "black" ) +
  labs(x="True positive prediction probability", y = "Early warning signal indicator") +
  #coord_cartesian(xlim = c(-4.5,4.5))+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  facet_wrap(~res)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())  

##########################################################################################
# Compare Metrics' False Negative Rate
##########################################################################################
ews_mod_metric_mth_false <- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code) ), 
                                     data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "no.trans"),
                                                  subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                     iter = its,
                                     thin = thn,
                                     warmup = wrmup,
                                     prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                              prior(normal(1, 1),class = sd)),
                                     family = binomial(link = "logit"), 
                                     chains = 4,
                                     control = list(adapt_delta = .99, max_treedepth = 20),
                                     seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_metric_yr_false<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code - 1 + (1|lake/method_code) ), 
                                   data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "no.trans"),
                                                subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                   seed = 12345, cores = 4,sample_prior = TRUE)


dat_metric_false_trials <-  ews_mod_metric_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_metric_yr_false |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable)) |>
  mutate(.variable = factor(.variable,levels = levels(dat_metric_true_trials$.variable)))  #reorder variable in to increasing ability based upon true postive results

dat_metric_false_halfeye <- ews_mod_metric_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_metric_yr_false |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_metric_true_trials$.variable))) #reorder variable in to increasing average ability based on dat_metric_true_trials



metric_false_plot <- ggplot(dat_metric_false_trials |>
                              mutate(.variable = sub("_.*", "",.variable))
                            ,aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= dat_metric_false_halfeye, alpha=0.5, aes(fill=variate)) +
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
                                interval_size_range = c(0.4, 1.2)
                                ,colour = "black") +
  labs(x="True negative prediction probability", y = "Early warning signal indicator") +
  #coord_cartesian(xlim = c(-4.5,4.5))+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  facet_wrap(~res)+
  #facet_grid(method_code~res,scales = "free_y")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank())      

require(patchwork)

ggsave(metric_true_plot +
  metric_false_plot +
  patchwork::plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'a') &
    theme(plot.tag = element_text(face = 'bold')
          # ,legend.background = element_rect(fill = '#EFEFEF'),
          # panel.grid.major = element_line(colour = "grey80"),
          # plot.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF'),
          # panel.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF')
          ) &
    labs(fill="EWS method"),
  filename = "/Users/ul20791/Downloads/eg_bayes_fig2.pdf",width = 10,height = 6)


##########################################################################################
# Compare Indicators' True Positive Rate
##########################################################################################

ews_mod_ind_mth_true <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|lake/method_code) ), 
                                     data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "trans"),
                                                  subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                     iter = its,
                                     thin = thn,
                                     warmup = wrmup,
                                     prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                              prior(normal(1, 1),class = sd)),
                                     family = binomial(link = "logit"), 
                                     chains = 4,
                                     control = list(adapt_delta = .975, max_treedepth = 20),
                                     seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_ind_yr_true<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|lake/method_code) ), 
                                   data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "trans"),
                                                subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "trans")),
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   control = list(adapt_delta = .975, max_treedepth = 20),
                                   seed = 12345, cores = 4,sample_prior = TRUE)

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




ind_true_plot <- ggplot(dat_ind_true_trials |>   
                             mutate(.variable = sub("_.*", "",.variable)),aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= dat_ind_true_halfeye, alpha=0.5, aes(fill=variate)) +
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
tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),interval_size_range = c(0.4, 1.2),colour = "black" ) +
  labs(x="True positive prediction probability", y = "Early warning signal indicator") +
  #coord_cartesian(xlim = c(-4.5,4.5))+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  #facet_wrap(~res)+
  facet_grid(method_code~res,scales = "free_y")+
  
  theme_bw()+
  theme(panel.grid.minor = element_blank())  



ind_true_plot2 <- ggplot(data = dat_ind_true_trials |>   
         mutate(.variable = sub("_.*", "",.variable),
                computation = sub("*._", "",method_code),
                computation = ifelse(computation == "ML","EWSNet",computation)) |>
         filter(.width == 0.95),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye, method_code == "ML"), alpha=0.5, aes(fill=variate)) +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye, method_code == "expanding" & variate == "univariate"), aes(fill=variate),alpha=0.5,) +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye, method_code == "expanding" & variate == "multivariate"), aes(fill=variate),alpha=0.5) +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye, method_code == "rolling" & variate == "univariate"), aes(fill=variate),alpha=0.5) +
  tidybayes::stat_slab(data= subset(dat_ind_true_halfeye, method_code == "rolling" & variate == "multivariate"), aes(fill=variate),alpha=0.5) +
  scale_fill_manual(values=c("#529928", "#5d3099","#bfbd3d"))+
  tidybayes::geom_pointinterval(data = ~subset(.,computation %in% c("expanding","EWSNet")),aes(xmin = .lower, xmax = .upper,alpha = computation),colour = "black") +
  tidybayes::geom_pointinterval(data = ~subset(.,computation == "rolling"),aes(xmin = .lower, xmax = .upper,alpha = computation),colour = "grey") +
  labs(x="True positive prediction probability", y = "Early warning signal indicator",
       colour = "EWS method",fill = "EWS method",alpha = "Computation") +
  scale_alpha_manual(breaks = c("rolling","expanding"), values=rep(0.75,3),labels = c("rolling","expanding/\nEWSNet"))+
  guides(alpha = guide_legend(override.aes = list(colour = c("grey","black")) ) )+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank()) 

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
# Compare Indicators' False Negative Rate
##########################################################################################
ews_mod_ind_mth_false <- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|lake/method_code) ), 
                                      data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "no.trans"),
                                                   subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                      iter = its,
                                      thin = thn,
                                      warmup = wrmup,
                                      prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                               prior(normal(1, 1),class = sd)),
                                      family = binomial(link = "logit"), 
                                      chains = 4,
                                      control = list(adapt_delta = .99, max_treedepth = 20),
                                      seed = 12345, cores = 4,sample_prior = TRUE)

ews_mod_ind_yr_false<- brms::brm(brms::bf(total_success | trials(offset) ~ indicator - 1 + (1|lake/method_code) ), 
                                    data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" & outcome == "no.trans"),
                                                 subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML" & outcome == "no.trans")),
                                    iter = its,
                                    thin = thn,
                                    warmup = wrmup,
                                    prior= c(prior(normal(0, 1.2), class = b,lb= -5.5,ub=5.5),
                                             prior(normal(1, 1),class = sd)),
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    control = list(adapt_delta = .99, max_treedepth = 20,stepsize = 0.01),
                                    seed = 12345, cores = 4,sample_prior = TRUE)


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


ind_false_plot <- ggplot(dat_ind_false_trials |>
                              mutate(.variable = sub("_.*", "",.variable))
                            ,aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= dat_ind_false_halfeye, alpha=0.5, aes(fill=variate)) +
  # geom_errorbar(data = subset(dat_metric_trials,.width == 0.95), 
  #               aes(xmin=.lower, xmax=.upper,col=variate), 
  #                width=0, alpha=0.9, size=1.3)+
  # geom_point(size=4,aes(fill=variate),shape=21) +
  scale_fill_manual(values=c("#529928",
                                      "#5d3099",
                                      "#bfbd3d"))+
scale_color_manual(values=c("#529928","#5d3099","#bfbd3d"))+
tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),interval_size_range = c(0.4, 1.2),colour = "black") +
  labs(x="True negative prediction probability", y = "Early warning signal indicator") +
  #coord_cartesian(xlim = c(-4.5,4.5))+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  #facet_wrap(~res)+
  facet_grid(method_code~res,scales = "free_y")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank())      

require(patchwork)

ggsave(ind_true_plot +
         ind_false_plot +
         patchwork::plot_layout(guides = 'collect') + 
         plot_annotation(tag_levels = 'a') &
         theme(plot.tag = element_text(face = 'bold')
               # ,legend.background = element_rect(fill = '#EFEFEF'),
               # panel.grid.major = element_line(colour = "grey80"),
               # plot.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF'),
               # panel.background = element_rect(fill = '#EFEFEF',colour = '#EFEFEF')
         ) &
         labs(fill="EWS method"),
       filename = "/Users/ul20791/Downloads/eg_bayes_fig2.pdf",width = 10,height = 6)




ind_false_plot2 <- ggplot(data = dat_ind_false_trials |>   
         mutate(.variable = sub("_.*", "",.variable),
                computation = sub("*._", "",method_code),
                computation = ifelse(computation == "ML","EWSNet",computation)) |>
         filter(.width == 0.95),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye, method_code == "ML"), alpha=0.5, aes(fill=variate)) +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye, method_code == "expanding" & variate == "univariate"), aes(fill=variate),alpha=0.5) +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye, method_code == "expanding" & variate == "multivariate"), aes(fill=variate),alpha=0.5) +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye, method_code == "rolling" & variate == "univariate"), aes(fill=variate),alpha=0.5) +
  tidybayes::stat_slab(data= subset(dat_ind_false_halfeye, method_code == "rolling" & variate == "multivariate"), aes(fill=variate),alpha=0.5) +
  scale_fill_manual(values=c("#529928", "#5d3099","#bfbd3d"))+
  tidybayes::geom_pointinterval(data = ~subset(.,computation %in% c("expanding","EWSNet")),aes(xmin = .lower, xmax = .upper,alpha = computation),colour = "black") +
  tidybayes::geom_pointinterval(data = ~subset(.,computation == "rolling"),aes(xmin = .lower, xmax = .upper,alpha = computation),colour = "grey") +
  labs(x="True positive prediction probability", y = "Early warning signal indicator",
       colour = "EWS method",fill = "EWS method",alpha = "Computation") +
  scale_alpha_manual(breaks = c("rolling","expanding"), values=rep(0.75,3),labels = c("rolling","expanding/\nEWSNet"))+
  guides(alpha = guide_legend(override.aes = list(colour = c("grey","black")) ) )+
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4.5,4.5))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank())      


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
  labs(x="True positive prediction probability", y = "Early warning signal indicator",
       fill = "EWS method",colour = "Computation") +
  scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
  coord_cartesian(xlim = c(-4,4))+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'lines')) 

ind_true_plot3 +
  ind_false_plot3 +
  patchwork::plot_layout(guides = 'collect') + 
  plot_annotation(tag_levels = 'A')

##########################################################################################
# Compare Indicators (Lakes)
##########################################################################################

ews_mod_metric_mth_lake<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code:lake - 1 + (1|method_code/outcome) ), 
                                    data = rbind(subset(metric_data, res  == "Monthly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" ),
                                                 subset(metric_data, res  == "Monthly" & detrend_meth == "none" & deseason_meth == "none" & method_code == "univariate_ML")),                                    
                                    iter = its,
                                    thin = thn,
                                    warmup = wrmup,
                                    prior= c(prior(normal(0, 1.2), class = b),
                                             prior(normal(1, 1),class = sd)),
                                    family = binomial(link = "logit"), 
                                    chains = 4,
                                    control = list(adapt_delta = .975, max_treedepth = 20),
                                    seed = 12345, cores = 4)

ews_mod_metric_yr_lake<- brms::brm(brms::bf(total_success | trials(offset) ~ metric.code:lake - 1 + (1|method_code/outcome) ), 
                                   data = rbind(subset(metric_data, res  == "Yearly" & detrend_meth == "gaussian" & deseason_meth == "average" & method_code != "univariate_ML" ),
                                                subset(metric_data, res  == "Yearly" & detrend_meth == "none" & deseason_meth == "decompose" & method_code == "univariate_ML")),                                       
                                   iter = its,
                                   thin = thn,
                                   warmup = wrmup,
                                   prior= c(prior(normal(0, 1.2), class = b),
                                            prior(normal(1, 1),class = sd)),
                                   family = binomial(link = "logit"), 
                                   chains = 4,
                                   control = list(adapt_delta = .99, max_treedepth = 20),
                                   seed = 12345, cores = 4)
#bayestestR::describe_posterior(ews.mod.metric.mth.cond, ci = 0.95, test="none")


brms::mcmc_plot(ews_mod_metric_mth_lake, 
                type = "areas",
                #type = "intervals",
                prob = 0.95)
inv_logit_scaled(fixef(ews_mod_metric_mth_lake)[,1])
inv_logit_scaled(fixef(ews_mod_metric_mth_lake)[,1])


dat_metric_lake_trials <-  ews_mod_metric_mth_lake |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_metric_yr_lake |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable)) |>
  separate(.variable,c(".variable","lake"),sep = ":") |>
  mutate(lake = gsub("lake","",lake))


dat_metric_lake_halfeye <- ews_mod_metric_mth_lake |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_metric_yr_lake |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_metric.code", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","multivariate")),
                  res = "Yearly"))|>
  separate(.variable,c(".variable","lake"),sep = ":") |>
  mutate(.variable = gsub("P"," + ",.variable))|>
  mutate(lake = gsub("lake","",lake))

ggplot(subset(dat_metric_lake_trials,.width == 0.95),
       aes(y = .variable, x = .value
       )) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  #tidybayes::stat_slab(data= dat_metric_lake_halfeye, alpha=0.5, aes(fill=variate)) +
  geom_errorbar(aes(xmin=.lower, xmax=.upper,col=variate),
                width=0, alpha=0.9, size=1.3)+
  geom_point(aes(fill=variate),size=2,shape=21) +
  scale_fill_manual(values=c("#529928",
                                      "#5d3099",
                                      "#bfbd3d"))+
                                        scale_color_manual(values=c("#529928",
                                                                             "#5d3099",
                                                                             "#bfbd3d"))+
                                                                               # tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),
  #                               interval_size_range = c(0.8, 2)
  #                               ,colour = "black" ) +
  labs(x="Probability of correct prediction", y = "Early warning signal metric") +
  coord_cartesian(xlim = c(-3.1,3.1))+
  facet_grid(lake~res)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

