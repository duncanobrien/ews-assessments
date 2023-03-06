########################################################################################################################
# Figure 4 #
########################################################################################################################
require(tidyverse)
require(patchwork)
require(tidybayes)

ews_mod_ind_mth_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_true.rds")
ews_mod_ind_yr_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_true.rds")
ews_mod_ind_mth_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_false.rds")
ews_mod_ind_yr_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_false.rds")

#as the back convertion from log odds to probabilities is a sigmoidal function, add scaling function
#to x axis to represent this. Causes skew of posterior probability distributions on plot
# labels <- seq(-4,4,2)
# breaks <-  seq(-4,4,2)

# labels <- c(-4,-2,-1.4,-0.4,0,0.4,1.4,2,4)
# breaks <- c(-4,-2,-1.4,-0.4,0,0.4,1.4,2,4)

labels <- c(-4,round(brms::logit_scaled(c(0.25,0.5,0.75)),1),4) #ensure breaks @ 0.25 probability intervals
breaks <- c(-4,round(brms::logit_scaled(c(0.25,0.5,0.75)),1),4)

#back transformation function from log odds to probabilities via inverse logit
inv_logit_perc <- scales::trans_new("inv_logit_perc",
                        transform = function(x){suppressWarnings(brms::inv_logit_scaled(x))},
                        inverse = function(x){suppressWarnings(brms::logit_scaled(x))})

############ 
#Extract posterior draws
############ 
dat_ind_true_trials <-  ews_mod_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_ind_yr_true |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = forcats::fct_reorder(.variable,.value,.fun = mean,.desc = FALSE))  #reorder variable in to increasing ability

dat_ind_true_halfeye <- ews_mod_ind_mth_true |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_ind_yr_true |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))

dat_ind_false_trials <-  ews_mod_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  tidybayes::median_qi(.width = c(.95, .8, .5)) |>
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
         res = "Monthly") |>
  rbind(ews_mod_ind_yr_false |>
          tidybayes::gather_draws(`b.*`,regex = T) |>
          mutate(.variable = gsub("b_indicator", "", .variable)) |> 
          tidybayes::median_qi(.width = c(.95, .8, .5)) |>
          mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
                 res = "Yearly")) |>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable,levels = levels(dat_ind_true_trials$.variable)))  #reorder variable in to increasing ability based upon true postive results

dat_ind_false_halfeye <- ews_mod_ind_mth_false |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_indicator", "", .variable)) |> 
  mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
         res = "Monthly") |>
  rbind( ews_mod_ind_yr_false |>
           tidybayes::gather_draws(`b.*`,regex = T) |>
           mutate(.variable = gsub("b_indicator", "", .variable)) |> 
           mutate(variate = ifelse(grepl("ar1|^SD|skew",.variable),"Univariate",ifelse(grepl("scaled|unscaled",.variable),"EWSNet","Multivariate")),
                  res = "Yearly"))|>
  mutate(.variable = gsub("P"," + ",.variable),
         method_code = gsub(".*_", "",.variable)) |>
  mutate(.variable = factor(.variable, levels = levels(dat_ind_true_trials$.variable))) |> #reorder variable in to increasing average ability based on dat_metric_true_trials
  mutate(.variable = sub("_.*", "",.variable))


############ 
#Create individual panels
############ 
ind_true_plot <- ggplot(data = dat_ind_true_trials |>   
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
  #scale_x_continuous(labels = function(i){round(brms::inv_logit_scaled(i),1)})+
  #coord_cartesian(xlim = c(-4,4))+
  scale_x_continuous(labels =  ifelse(is.na(labels),"",round(brms::inv_logit_scaled(labels),1)), breaks = breaks)+
  coord_trans(x = inv_logit_perc)+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.size = unit(1.5, 'lines')) 

ind_false_plot <- ggplot(data = dat_ind_false_trials |>   
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
  #scale_x_continuous(labels = function(i){round(brms::inv_logit_scaled(i),1)})+
  #coord_cartesian(xlim = c(-4,4))+
  scale_x_continuous(labels =  ifelse(is.na(labels),"",round(brms::inv_logit_scaled(labels),1)), breaks = breaks)+
  coord_trans(x = inv_logit_perc)+
  facet_grid(variate~res,scales = "free_y",space = "free")+
  theme_bw()+
  theme(axis.title.y =  element_blank(),
        axis.text.y =  element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.key.size = unit(1.5, 'lines'))

############ 
#Create figure
############ 
ggsave(ind_true_plot +
         ind_false_plot +
         patchwork::plot_layout(guides = 'collect') + 
         plot_annotation(tag_levels = 'A'),
       filename = "Figures/figure_4.pdf",width = 9,height = 6)
