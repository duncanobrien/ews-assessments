########################################################################################################################
# Figure 3 #
########################################################################################################################
require(tidyverse)
require(patchwork)
require(tidybayes)

ews_mod_method_mth <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_mth.rds")
ews_mod_method_yr <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_yr.rds")

############ 
#Extract posterior draws
############ 
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

############ 
#Create figure
############ 
#as the back convertion from log odds to probabilities is a sigmoidal function, add scaling function
#to x axis to represent this. Causes skew of posterior probability distributions on plot
labels <- c(-4,round(brms::logit_scaled(c(0.25,0.5,0.75)),1),4) #ensure breaks @ 0.25 probability intervals
breaks <- c(-4,round(brms::logit_scaled(c(0.25,0.5,0.75)),1),4)

#back transformation function from log odds to probabilities via inverse logit
inv_logit_perc <- scales::trans_new("inv_logit_perc",
                                    transform = function(x){suppressWarnings(brms::inv_logit_scaled(x))},
                                    inverse = function(x){suppressWarnings(brms::logit_scaled(x))})


ggsave(ggplot(dat_method_trials,aes(y = .variable, x = .value)) +
         geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
         tidybayes::stat_slab(data= dat_method_halfeye, alpha=0.5, aes(fill=variate)) +
         scale_fill_manual(values=c("#529928",
                                             "#5d3099",
                                             "#bfbd3d"))+
                                               scale_color_manual(values=c("#529928","#5d3099","#bfbd3d"))+
         tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),interval_size_range = c(0.8, 2),colour = "black" ) +
         labs(x="Probability of correct prediction", y = "Early warning signal method") +
         #scale_x_continuous(labels = function(i){round(inv_logit_scaled(i),1)})+
         #coord_cartesian(xlim = c(-3,3))+
         scale_x_continuous(labels = round(brms::inv_logit_scaled(labels),1), breaks = breaks)+
         coord_trans(x = inv_logit_perc,xlim = c(-4,4))+
         scale_y_discrete(labels = c("univariate rolling", "multivariate rolling", "EWSNet", "multivariate expanding","univariate expanding"))+
         facet_wrap(~res)+
         theme_bw()+
         theme(legend.position = "none",
               panel.grid.minor = element_blank()),
       filename = "Figures/figure_3.pdf",width = 5,height = 4)
