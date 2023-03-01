########################################################################################################################
# Supplementary Tables #
########################################################################################################################

require(brms)
require(bayestestR)
require(tidyverse)

load("Results/preprocessing_models.RData")

############ 
#Extract model coefficients
############ 

tableS1 <- bayestestR::describe_posterior(ews.mod.detrend.mth.uni.roll, ci = 0.95, test="none") |>
  dplyr::select(-CI)

tableS2 <- bayestestR::describe_posterior(ews.mod.detrend.mth.uni.exp, ci = 0.95, test="none")|>
  dplyr::select(-CI)

tableS3 <- bayestestR::describe_posterior(ews.mod.detrend.mth.multi.roll, ci = 0.95, test="none")|>
  dplyr::select(-CI)

tableS4 <- bayestestR::describe_posterior(ews.mod.detrend.mth.multi.exp, ci = 0.95, test="none")|>
  dplyr::select(-CI)

tableS5 <- bayestestR::describe_posterior(ews.mod.detrend.mth.ml, ci = 0.95, test="none")|>
  dplyr::select(-CI)

############ 
#Export Tables
############ 

write.csv(tableS1,"Results/table_S1.csv")
write.csv(tableS2,"Results/table_S2.csv")
write.csv(tableS3,"Results/table_S3.csv")
write.csv(tableS4,"Results/table_S4.csv")
write.csv(tableS4,"Results/table_S5.csv")

