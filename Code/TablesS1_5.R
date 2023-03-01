########################################################################################################################
# Supplementary Tables #
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
tableS1 <- bayestestR::describe_posterior(ews_mod_detrend_mth_uni_roll, ci = 0.95, test="none")[-1,-3] 

tableS2 <- bayestestR::describe_posterior(ews_mod_detrend_mth_uni_exp, ci = 0.95, test="none")[-1,-3] 

tableS3 <- bayestestR::describe_posterior(ews_mod_detrend_mth_multi_roll, ci = 0.95, test="none")[-1,-3] 

tableS4 <- bayestestR::describe_posterior(ews_mod_detrend_mth_multi_exp, ci = 0.95, test="none")[-1,-3] 

tableS5 <- bayestestR::describe_posterior(ews_mod_detrend_mth_ml, ci = 0.95, test="none")[-1,-3] 

############ 
#Export Tables
############ 

write.csv(tableS1,"Results/supplementary_info/supplementary_tables/table_S1.csv")
write.csv(tableS2,"Results/supplementary_info/supplementary_tables/table_S2.csv")
write.csv(tableS3,"Results/supplementary_info/supplementary_tables/table_S3.csv")
write.csv(tableS4,"Results/supplementary_info/supplementary_tables/table_S4.csv")
write.csv(tableS5,"Results/supplementary_info/supplementary_tables/table_S5.csv")

