########################################################################################################################
# Model summaries #
########################################################################################################################

require(brms)
require(tidyverse)

############ 
#Section: Early warning signal overall ability
############ 
#load summary tables
tableS6 <- read.csv("Results/supplementary_info/supplementary_tables/table_S6.csv")
tableS7 <- read.csv("Results/supplementary_info/supplementary_tables/table_S7.csv")

multi_exp_mth <- brms::inv_logit_scaled(subset(tableS6,Parameter == "multivariate_expanding")[,c("Median","CI_low","CI_high")])
multi_exp_yr <- brms::inv_logit_scaled(subset(tableS7,Parameter == "multivariate_expanding")[,c("Median","CI_low","CI_high")])

uni_exp_mth <- brms::inv_logit_scaled(subset(tableS6,Parameter == "univariate_expanding")[,c("Median","CI_low","CI_high")])
uni_exp_yr <- brms::inv_logit_scaled(subset(tableS7,Parameter == "univariate_expanding")[,c("Median","CI_low","CI_high")])

uni_roll_yr <- brms::inv_logit_scaled(subset(tableS7,Parameter == "univariate_rolling")[,c("Median","CI_low","CI_high")])

############ 
#Section: Individual indicator ability
############ 
#load summary tables
tableS10 <- read.csv("Results/supplementary_info/supplementary_tables/table_S10.csv")
tableS11 <- read.csv("Results/supplementary_info/supplementary_tables/table_S11.csv")
tableS12 <- read.csv("Results/supplementary_info/supplementary_tables/table_S12.csv")
tableS13 <- read.csv("Results/supplementary_info/supplementary_tables/table_S13.csv")

ewsnet_mth_true_unsc <- brms::inv_logit_scaled(subset(tableS10,Parameter == "unscaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_yr_true_unsc <- brms::inv_logit_scaled(subset(tableS11,Parameter == "unscaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_mth_false_unsc <- brms::inv_logit_scaled(subset(tableS12,Parameter == "unscaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_yr_false_unsc <- brms::inv_logit_scaled(subset(tableS13,Parameter == "unscaled_ML")[,c("Median","CI_low","CI_high")])

ewsnet_mth_true_sc <- brms::inv_logit_scaled(subset(tableS10,Parameter == "scaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_yr_true_sc <- brms::inv_logit_scaled(subset(tableS11,Parameter == "scaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_mth_false_sc <- brms::inv_logit_scaled(subset(tableS12,Parameter == "scaled_ML")[,c("Median","CI_low","CI_high")])
ewsnet_yr_false_sc <- brms::inv_logit_scaled(subset(tableS13,Parameter == "scaled_ML")[,c("Median","CI_low","CI_high")])

ewsnet_mth_true_summary <- tableS10 |>
  dplyr::mutate(method = gsub(".*_","",Parameter),
                method = case_when(
                  method == "ML" ~ method,
                  grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                  TRUE ~ paste("univariate",method,sep="_")
                )) |>
  dplyr::mutate(across(c(Median,CI_low,CI_high),~ brms::inv_logit_scaled(.x))) |>
  dplyr::reframe(median = median(Median),sd = sd(Median),
          .by= "method")

ewsnet_yr_true_summary <- tableS11 |>
  dplyr::mutate(method = gsub(".*_","",Parameter),
                method = case_when(
                  method == "ML" ~ method,
                  grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                  TRUE ~ paste("univariate",method,sep="_")
                )) |>
  dplyr::mutate(across(c(Median,CI_low,CI_high),~ brms::inv_logit_scaled(.x))) |>
  dplyr::reframe(median = median(Median),sd = sd(Median),
          .by= "method")

ewsnet_mth_false_summary <- tableS12 |>
  dplyr::mutate(method = gsub(".*_","",Parameter),
                method = case_when(
                  method == "ML" ~ method,
                  grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                  TRUE ~ paste("univariate",method,sep="_")
                )) |>
  dplyr::mutate(across(c(Median,CI_low,CI_high),~ brms::inv_logit_scaled(.x))) |>
  dplyr::reframe(median = median(Median),sd = sd(Median),
                 .by= "method")

ewsnet_yr_false_summary <- tableS13 |>
  dplyr::mutate(method = gsub(".*_","",Parameter),
                method = case_when(
                  method == "ML" ~ method,
                  grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                  TRUE ~ paste("univariate",method,sep="_")
                )) |>
  dplyr::mutate(across(c(Median,CI_low,CI_high),~ brms::inv_logit_scaled(.x))) |>
  dplyr::reframe(median = median(Median),sd = sd(Median),
                 .by= "method")

all_indicator <- rbind(tableS10 |>
                         dplyr::mutate(method = gsub(".*_","",Parameter),
                                       method = case_when(
                                         method == "ML" ~ method,
                                         grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                                         TRUE ~ paste("univariate",method,sep="_")
                                       )) |>
                         dplyr::mutate(res = "Monthly", ts = "true"),
                       tableS9 |>
                         dplyr::mutate(method = gsub(".*_","",Parameter),
                                       method = case_when(
                                         method == "ML" ~ method,
                                         grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                                         TRUE ~ paste("univariate",method,sep="_")
                                       )) |>
                         dplyr::mutate(res = "Yearly", ts = "true"),
                       tableS10 |>
                          dplyr::mutate(method = gsub(".*_","",Parameter),
                                        method = case_when(
                                          method == "ML" ~ method,
                                          grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                                          TRUE ~ paste("univariate",method,sep="_")
                                        )) |>
                         dplyr::mutate(res = "Monthly", ts = "false"),
                       tableS11 |>
                         dplyr::mutate(method = gsub(".*_","",Parameter),
                                       method = case_when(
                                         method == "ML" ~ method,
                                         grepl("eigen|maf|pca|mean|max|mut",Parameter) ~ paste("multivariate",method,sep="_"),
                                         TRUE ~ paste("univariate",method,sep="_")
                                       )) |>
                         dplyr::mutate(res = "Yearly", ts = "false")) |>
  dplyr::mutate(Parameter = gsub("_.*","",Parameter)) |>
  dplyr::mutate(across(c(Median,CI_low,CI_high),~ brms::inv_logit_scaled(.x)))|>
  group_by(Parameter,method) |>
  filter(all(Median >= 0.48))
