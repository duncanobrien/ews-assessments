########################################################################################################################
# Supplementary Tables 1 - 5 #
########################################################################################################################

require(brms)
require(bayestestR)
require(tidyverse)

ews_mod_detrend_mth_uni_roll <- readRDS(file = "Results/supplementary_info/ews_mod_detrend_mth_uni_roll.rds")
ews_mod_detrend_mth_uni_exp <- readRDS(file = "Results/supplementary_info/ews_mod_detrend_mth_uni_exp.rds")
ews_mod_detrend_mth_multi_roll <- readRDS(file = "Results/supplementary_info/ews_mod_detrend_mth_multi_roll.rds")
ews_mod_detrend_mth_multi_exp <- readRDS(file = "Results/supplementary_info/ews_mod_detrend_mth_multi_exp.rds")
ews_mod_detrend_mth_ml <- readRDS(file = "Results/supplementary_info/ews_mod_detrend_mth_ml.rds")

############ 
#Extract model coefficients
############ 
#drop first row (as is the reference none-none detrending and deseasoning method) and third column (as is jsut stating the credible interval value - i.e. 95%)
tableS1 <- bayestestR::describe_posterior(ews_mod_detrend_mth_uni_roll, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_combo_codeunivariate_rolling_","",Parameter),
         Parameter = gsub("_","-",Parameter)) |>
  dplyr::mutate(Parameter = factor(Parameter,levels = c("linear-none","loess-none","gaussian-none","none-average","none-decompose","none-stl",
                                                 "linear-average", "loess-average","gaussian-average","linear-decompose","loess-decompose","gaussian-decompose",
                                                 "linear-stl","loess-stl","gaussian-stl"))) |>
  dplyr::arrange(Parameter) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS2 <- bayestestR::describe_posterior(ews_mod_detrend_mth_uni_exp, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_combo_codeunivariate_expanding_","",Parameter),
                Parameter = gsub("_","-",Parameter)) |>
  dplyr::mutate(Parameter = factor(Parameter,levels = c("linear-none","loess-none","gaussian-none","none-average","none-decompose","none-stl",
                                                        "linear-average", "loess-average","gaussian-average","linear-decompose","loess-decompose","gaussian-decompose",
                                                        "linear-stl","loess-stl","gaussian-stl"))) |>
  dplyr::arrange(Parameter) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS3 <- bayestestR::describe_posterior(ews_mod_detrend_mth_multi_roll, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_combo_codemultivariate_rolling_","",Parameter),
                Parameter = gsub("_","-",Parameter)) |>
  dplyr::mutate(Parameter = factor(Parameter,levels = c("linear-none","loess-none","gaussian-none","none-average","none-decompose","none-stl",
                                                        "linear-average", "loess-average","gaussian-average","linear-decompose","loess-decompose","gaussian-decompose",
                                                        "linear-stl","loess-stl","gaussian-stl"))) |>
  dplyr::arrange(Parameter) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS4 <- bayestestR::describe_posterior(ews_mod_detrend_mth_multi_exp, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_combo_codemultivariate_expanding_","",Parameter),
                Parameter = gsub("_","-",Parameter)) |>
  dplyr::mutate(Parameter = factor(Parameter,levels = c("linear-none","loess-none","gaussian-none","none-average","none-decompose","none-stl",
                                                        "linear-average", "loess-average","gaussian-average","linear-decompose","loess-decompose","gaussian-decompose",
                                                        "linear-stl","loess-stl","gaussian-stl"))) |>
  dplyr::arrange(Parameter) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS5 <- bayestestR::describe_posterior(ews_mod_detrend_mth_ml, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_combo_codeunivariate_ML_","",Parameter),
                Parameter = gsub("_","-",Parameter)) |>
  dplyr::mutate(Parameter = factor(Parameter,levels = c("linear-none","loess-none","gaussian-none","none-average","none-decompose","none-stl",
                                                        "linear-average", "loess-average","gaussian-average","linear-decompose","loess-decompose","gaussian-decompose",
                                                        "linear-stl","loess-stl","gaussian-stl"))) |>
  dplyr::arrange(Parameter) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

############ 
#Export Tables
############ 

write.csv(tableS1,"Results/supplementary_info/supplementary_tables/table_S1.csv",row.names = FALSE)
write.csv(tableS2,"Results/supplementary_info/supplementary_tables/table_S2.csv",row.names = FALSE)
write.csv(tableS3,"Results/supplementary_info/supplementary_tables/table_S3.csv",row.names = FALSE)
write.csv(tableS4,"Results/supplementary_info/supplementary_tables/table_S4.csv",row.names = FALSE)
write.csv(tableS5,"Results/supplementary_info/supplementary_tables/table_S5.csv",row.names = FALSE)
