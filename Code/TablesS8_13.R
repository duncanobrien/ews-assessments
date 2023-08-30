########################################################################################################################
# Supplementary Tables 6 - 11 #
########################################################################################################################

require(brms)
require(bayestestR)
require(tidyverse)

ews_mod_method_mth <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_mth.rds")
ews_mod_method_yr <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_yr.rds")
ews_mod_ind_mth_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_true.rds")
ews_mod_ind_yr_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_true.rds")
ews_mod_ind_mth_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_false.rds")
ews_mod_ind_yr_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_false.rds")

############ 
#Extract model coefficients
############ 

#drop third column (as is jsut stating the credible interval value - i.e. 95%)
tableS6 <- bayestestR::describe_posterior(ews_mod_method_mth, ci = 0.95, test="none")[,-3] |>
  dplyr::mutate(Parameter = gsub("b_method_code","",Parameter),
                Parameter = ifelse(Parameter == 'univariate_ML',"EWSNet",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS7 <- bayestestR::describe_posterior(ews_mod_method_yr, ci = 0.95, test="none")[,-3] |>
  dplyr::mutate(Parameter = gsub("b_method_code","",Parameter),
                Parameter = ifelse(Parameter == 'univariate_ML',"EWSNet",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS8 <- bayestestR::describe_posterior(ews_mod_ind_mth_true, ci = 0.95, test="none")[,-3] |>
  dplyr::mutate(Parameter = gsub("b_indicator","",Parameter),
                Parameter = gsub("P"," + ",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS9 <- bayestestR::describe_posterior(ews_mod_ind_yr_true, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_indicator","",Parameter),
                Parameter = gsub("P"," + ",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS10 <- bayestestR::describe_posterior(ews_mod_ind_mth_false, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_indicator","",Parameter),
                Parameter = gsub("P"," + ",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

tableS11 <- bayestestR::describe_posterior(ews_mod_ind_yr_false, ci = 0.95, test="none")[-1,-3] |>
  dplyr::mutate(Parameter = gsub("b_indicator","",Parameter),
                Parameter = gsub("P"," + ",Parameter)) |>
  dplyr::mutate(across(Median:CI_high,~round(.,digits =3))) |>
  dplyr::mutate(across(Rhat:ESS,~round(.,digits = 2)))

############ 
#Export Tables
############ 

write.csv(tableS6,"Results/supplementary_info/supplementary_tables/table_S6.csv",row.names = FALSE)
write.csv(tableS7,"Results/supplementary_info/supplementary_tables/table_S7.csv",row.names = FALSE)
write.csv(tableS8,"Results/supplementary_info/supplementary_tables/table_S8.csv",row.names = FALSE)
write.csv(tableS9,"Results/supplementary_info/supplementary_tables/table_S9.csv",row.names = FALSE)
write.csv(tableS10,"Results/supplementary_info/supplementary_tables/table_S10.csv",row.names = FALSE)
write.csv(tableS11,"Results/supplementary_info/supplementary_tables/table_S11.csv",row.names = FALSE)
