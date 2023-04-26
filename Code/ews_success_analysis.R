##########################################################################################
# Preamble
##########################################################################################

require(tidybayes)
require(ggplot2)
require(data.table)
require(Matrix)
require(tidyverse)
require(brms)

source("Code/extract_ews_pred_fn.R")
#load("Data/wrangled_genus_plank_data.Rdata")

genus_lake_exp_multi_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_exp_multi_ews_phyto.csv")[,-1], #drop write.csv() introduced X column
                                  read.csv(file = "Results/lake_results/genus_lake_exp_multi_ews_zoo.csv")[,-1])

genus_lake_exp_uni_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_exp_uni_ews_phyto.csv.gz")[,-1],
                                read.csv(file = "Results/lake_results/genus_lake_exp_uni_ews_zoo.csv.gz")[,-1])

genus_lake_roll_uni_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_phyto.csv")[,-1],
                                 read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_zoo.csv")[,-1])

genus_lake_roll_multi_ews <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_phyto.csv")[,-1],
                                   read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_zoo.csv")[,-1])

genus_lake_ewsnet <- rbind(read.csv(file = "Results/lake_results/genus_lake_ewsnet_phyto.csv")[,-1],
                           read.csv(file = "Results/lake_results/genus_lake_ewsnet_zoo.csv")[,-1])

genus_lake_roll_uni_ews_perm <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_perm_phyto.csv")[,-1],
                                      read.csv(file = "Results/lake_results/genus_lake_roll_uni_ews_perm_zoo.csv")[,-1])

genus_lake_roll_multi_ews_perm <- rbind(read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_perm_phyto.csv")[,-1],
                                        read.csv(file = "Results/lake_results/genus_lake_roll_multi_ews_perm_zoo.csv")[,-1])
##########################################################################################
# Extract EWS success
##########################################################################################

# lake_outcome_whole <- data.frame("reference" = unique(paste(genus_lake_roll_uni_ews_perm$lake,genus_lake_roll_uni_ews_perm$res,sep = "_"))) |>
#                       dplyr::mutate("outcome" = 
#                         case_when(
#                           grepl("Kinneret", reference) ~ "trans",
#                           grepl("Kasumigaura", reference)  ~ "trans",
#                           grepl("Washington", reference) ~ "trans",
#                           grepl("Monona", reference) ~ "trans",
#                           TRUE ~ "no.trans"
#                         )) 

lake_outcome_troph <- subset(rbind(genus_lake_exp_multi_ews), select = c(lake,res,troph_level)) |>
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


ewsnet_diff_df <- extract_ews_pred(ews.data =  genus_lake_ewsnet,
                                   sensitivity = 0.7,
                                   outcome = lake_outcome_troph,
                                   method = "ML")

rollews_diff_df <- extract_ews_pred(ews.data = genus_lake_roll_uni_ews,
                                          sensitivity = 0.7,
                                          outcome = lake_outcome_troph,
                                          method = "rolling")

rollews_diff_perm_df <- extract_ews_pred(ews.data = genus_lake_roll_uni_ews_perm,
                                         sensitivity = 0.95,
                                         outcome = lake_outcome_troph,
                                         method = "rolling",surrogate = T)

expews_diff_df <- extract_ews_pred(ews.data = genus_lake_exp_uni_ews,
                                   sensitivity = 2,
                                   outcome = lake_outcome_troph,
                                   method = "expanding")


rollmultiews_diff_df <- extract_ews_pred(ews.data = genus_lake_roll_multi_ews,
                                         sensitivity = 0.7,
                                         outcome = lake_outcome_troph,
                                         method = "rolling")

rollmultiews_diff_perm_df <- extract_ews_pred(ews.data = genus_lake_roll_multi_ews_perm,
                                              sensitivity = 0.95,
                                              outcome = lake_outcome_troph,
                                              method = "rolling",
                                              surrogate = T)

expmultiews_diff_df <- extract_ews_pred(ews.data = genus_lake_exp_multi_ews,
                                        sensitivity = 2,
                                        outcome = lake_outcome_troph,
                                        method = "expanding")

ggplot(expews_diff_df |>
         group_by(data_source,metric.code,lake,res) |>
         slice_head(n=1) |>
         mutate(miss.fac = ifelse(prediction %in% c("miss","unknown"), 
                                  "Undetected","Detected")), 
       aes(x=metric.code, y = as.numeric(interval.prior))) +
  geom_point(aes(col=data_source,
                 fill = miss.fac,
                 shape = miss.fac), size = 3, alpha = 0.5)+
  scale_shape_manual(values = c(16, NA),name = "Warning")+
  facet_grid(lake~res,scales = "fixed") +
  theme(legend.position = "none")


ggplot(expews_diff_df |>
         group_by(data_source,metric.code,lake,res) |>
         slice_head(n=1) |>
         group_by(lake,res,metric.code,prediction) |>
         summarise(count=n()), 
       aes(x=prediction, y = count,group = metric.code)) +
  geom_linerange(aes(ymin = 0, ymax = count),position = position_dodge(width = 0.75),size=0.25)+
  geom_point(aes(col=metric.code), size = 3,position = position_dodge(width = 0.75),alpha = 0.75)+
  scale_color_manual(values = c("#6886c4",
                                "#bfbd3d",
                                "#5d3099",
                                "#69c756",
                                "#e281fe",
                                "#6ca181",
                                "#76c3ef"))+ 
  facet_grid(lake~res) +
  theme_bw()

pal <- c("#6886c4",
         "#bfbd3d",
         "#5d3099",
         "#69c756",
         "#e281fe",
         "#6ca181",
         "#76c3ef",
         "#d06329",
         "#90676f",
         "#ce5c6e",
         "#5d4216",
         "black",
         "grey")

##########################################################################################
# Fit Models
##########################################################################################

all_ews_data <- as.data.table(ewsnet_diff_df) %>%
  .[,.SD[1], by = c("data_source", "scaling","lake","res")] %>%
  .[,c("data_source","troph_level","scaling","lake","res","method", "computation","prediction")] %>% #slice to first row of each group
  setnames(old = "scaling",new= "metric.code") %>%
  #.[metric.code == "scaled",] %>%
  #rbind(as.data.table(rollews_diff_df)[,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","metric.code","lake","res","method", "computation","prediction")]) %>%
  rbind(as.data.table(rollews_diff_perm_df)[,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction")]) %>%
  rbind(as.data.table(expews_diff_df)[,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction")]) %>%
  rbind(as.data.table(expmultiews_diff_df)[,data_source := NA][,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction")]) %>%
  #rbind(as.data.table(rollmultiews_diff_df)[,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","metric.code","lake","res","method", "computation","prediction")]) %>%
  rbind(as.data.table(rollmultiews_diff_perm_df)[,.SD[1], by = c("data_source", "metric.code","lake","res")][,c("data_source","troph_level","metric.code","lake","res","method", "computation","prediction")]) %>%
  .[,outcome := ifelse(prediction %in% c("match","prior"),
                       1,0)] %>%
  .[,variate := ifelse(grepl("multivariate",method),"multivariate","univariate")]%>%
  # .[,ts_length := ifelse(lake == "Kinneret" & res == "Yearly",dim(kin_yr_dat)[1],
  #                        ifelse(lake == "Kinneret" & res == "Monthly",dim(kin_mth_dat)[1],
  #                               ifelse(lake == "Kasumigaura" & res == "Yearly",dim(kas_yr_dat)[1],
  #                                      ifelse(lake == "Kasumigaura" & res == "Monthly",dim(kas_mth_dat)[1],
  #                                             ifelse(lake == "Windermere" & res == "Yearly",dim(wind_yr_dat)[1],
  #                                                    ifelse(lake == "Windermere" & res == "Monthly",dim(wind_mth_dat)[1],
  #                                                           ifelse(lake == "Mendota" & res == "Yearly",dim(mad_yr_dat)[1],
  #                                                                  ifelse(lake == "Mendota" & res == "Monthly",dim(mad_mth_dat)[1],
  #                                                                         ifelse(lake == "Lower Zurich" & res == "Yearly",dim(LZ_yr_dat)[1],dim(LZ_mth_dat)[1])))))))))] %>%
  # .[,ts_length := scale(as.numeric(ts_length))] %>%
  .[,method_code := paste(variate,computation,sep = "_")] %>%
  .[,method_code := factor(method_code, levels = c("univariate_rolling","univariate_expanding","univariate_ML","multivariate_rolling","multivariate_expanding"))] %>%
  .[,computation := factor(computation, levels = c("rolling","expanding","ML"))] %>%
  .[,variate := factor(variate, levels = c("univariate","multivariate"))] %>%
  .[,res := factor(res, levels = c("Yearly","Monthly"))]

computation_data <- all_ews_data %>%
  .[,offset := length(unique(data_source)), by =c("method_code","lake","res","troph_level") ] %>% #trials in terms of assessed time series
  .[,offset2 := length(unique(data_source))*length(unique(metric.code)), by =c("method_code","lake","res","troph_level")] %>% #trials in terms of assessed time series AND metrics 
  .[,.(total_success = sum(outcome),
       offset = unique(offset2),
       #ts_length = unique(ts_length),
       variate = unique(variate),
       reference = paste(lake,res,troph_level,sep="_")),by = c("lake", "res", "method_code","troph_level")] %>%
  merge(.,as.data.table(lake_outcome_troph),by="reference") %>%
  as.data.frame()

ews.mod.trials.mth <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code-1), 
                        data = computation_data[computation_data$res == "Monthly",],
                        iter = 10000,
                        thin = 0.0005*10000,
                        warmup = 0.1*10000,
                        prior= c(prior(normal(0, 1), class = b)),
                        family = binomial(link = "logit"), 
                        chains = 4,
                        control = list(adapt_delta = .975, max_treedepth = 20),
                        seed = 12345, cores = 4,sample_prior = TRUE)

ews.mod.trials.yr <- brms::brm(brms::bf(total_success | trials(offset) ~ method_code -1), 
                            data = computation_data[computation_data$res == "Yearly",],
                            iter = 10000,
                            thin = 0.0005*10000,
                            warmup = 0.1*10000,
                            prior= c(prior(normal(0, 1), class = b)),
                            family = binomial(link = "logit"), 
                            chains = 4,
                            control = list(adapt_delta = .975, max_treedepth = 20),
                            seed = 12345, cores = 4,sample_prior = TRUE)

brms::mcmc_plot(ews.mod.trials.mth, 
                type = "areas",
                #type = "intervals",
                prob = 0.95)

summary(ews.mod.trials.mth)
inv_logit_scaled(fixef(ews.mod.trials.mth)[,1])

dat_trials<-  ews.mod.trials.mth |>
  tidybayes::gather_draws(`b.*`,regex = T) |>
  mutate(.variable = gsub("b_method_code", "", .variable),
         .variable = gsub("b_Intercept", "univariate_rolling", .variable)) |> 
  #filter(.variable != "b_logoffset")
  filter(!grepl("*offset*",.variable)) |>
  mutate(res = "Monthly")|>
  rbind( ews.mod.trials.yr |>
               tidybayes::gather_draws(`b.*`,regex = T) |>
               mutate(.variable = gsub("b_method_code", "", .variable),
                      .variable = gsub("b_Intercept", "univariate_rolling", .variable)) |> 
               #filter(.variable != "b_logoffset")
               filter(!grepl("*offset*",.variable))|>
               mutate(res = "Yearly"))

bayesian.p1 <- ggplot(data = ews.mod.trials.mth |>
         tidybayes::gather_draws(`b.*`,regex = T) |>
         mutate(.variable = gsub("b_method_code", "", .variable),
                .variable = gsub("b_Intercept", "univariate_rolling", .variable)) |> 
         #filter(.variable != "b_logoffset")|>
         filter(!grepl("*offset*",.variable)) |>
         tidybayes::median_qi(.width = c(.95, .8, .5)) |>
         mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate","univariate"),
                res = "Monthly") |>
         rbind(ews.mod.trials.yr |>
                     tidybayes::gather_draws(`b.*`,regex = T) |>
                     mutate(.variable = gsub("b_method_code", "", .variable),
                            .variable = gsub("b_Intercept", "univariate_rolling", .variable)) |> 
                     #filter(.variable != "b_logoffset")|>
                     filter(!grepl("*offset*",.variable)) |>
                     tidybayes::median_qi(.width = c(.95, .8, .5)) |>
                     mutate(variate = ifelse(grepl("multivariate",.variable),"multivariate","univariate"),
                            res = "Yearly")),
       aes(y = .variable, x = .value)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
  tidybayes::stat_slab(data=dat_trials, alpha=0.5, aes(fill=.variable)) +
  scale_fill_manual(values=c("#5d3099",
                             "black",
                             "#5d3099",
                             "#bfbd3d",
                             "black"))+
  tidybayes::geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                                interval_size_range = c(0.8, 2)
                                ,colour = "black" ) +
  labs(x="Posterior estimate", y = "Early warning signal method") +
  facet_wrap(~res)+
  xlim(-1.6,1.6)+
  theme_bw()+
  theme(legend.position = "none")

ggsave(bayesian.p1, 
       filename = "Results/lake_results/model_estimates.pdf",
       width = 6,height = 4, dpi=144)


#####################################################################################
# Raw Assessment Counts
#####################################################################################
p1 <- ggplot(genus_lake_roll_uni_ews |>
         pivot_longer(c(ar1:skew),names_to = "metric.code",values_to = "tau") |>
           mutate(res = factor(res,levels = c("Yearly","Monthly"))),
       aes(x=lake, y=tau,col=metric.code)) + geom_boxplot() + 
  facet_wrap(~res) + theme_bw() +  ylab("Tau") +
  scale_x_discrete(labels = ~ stringr::str_wrap(.x, width = 10)) +
  guides(col=guide_legend(title="EWS metric")) +
  xlab("Lake system") 

p1 <- ggplot( subset(brms_data,method_code =="univariate_rolling") |>
                group_by(lake,res,metric.code) |>
                summarise(count = sum(outcome)/offset),
              aes(x=lake, y=count,col=metric.code)) + 
  geom_linerange(aes(ymin = 0, ymax = count),position = position_dodge(width = 0.75),size=0.25)+
  geom_point(position = position_dodge(width = 0.75)) + theme_bw() + 
  ylab("Proportion correct") + 
  xlab("Lake system") +
  scale_x_discrete(labels = ~ stringr::str_wrap(.x, width = 10)) +
  scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10),
                     values = scales::hue_pal()(3),
                     name = "EWS metric") +
  facet_wrap(~res)

p2 <- ggplot( subset(brms_data,method_code =="univariate_expanding") |>
          group_by(lake,res,metric.code) |>
         summarise(count = sum(outcome)/offset),
       aes(x=lake, y=count,col=metric.code)) + 
  geom_linerange(aes(ymin = 0, ymax = count),position = position_dodge(width = 0.75),size=0.25)+
  geom_point(position = position_dodge(width = 0.75)) + theme_bw() + 
  ylab("Proportion correct") + 
  xlab("Lake system") +
  scale_x_discrete(labels = ~ stringr::str_wrap(.x, width = 10)) +
  scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10),
                     values = scales::hue_pal()(7),
                     name = "EWS metric") +
  facet_wrap(~res)

p3 <- ggplot(ewsnet_diff_df |>
               mutate(res = factor(res,levels = c("Yearly","Monthly"))),
       aes(x=lake, y=prob,col=pos_outcome)) + geom_boxplot() + 
  scale_x_discrete(labels = ~ stringr::str_wrap(.x, width = 10)) +
  guides(col=guide_legend(title="Prediction")) +
facet_wrap(~res) + theme_bw() + ylab("Probability") +
  xlab("Lake system") 

require(patchwork)
ggsave((p1/ p2 / p3) + plot_layout(guides = 'collect'), 
       filename = "Results/lake_results/raw_ews_predictions.pdf",
       width = 9,height = 7, dpi=144)


