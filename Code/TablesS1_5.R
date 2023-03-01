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
#drop first row (as is the reference none-none detrending and deseasoning method) and third column (as is jsut stating the credible interval value - i.e. 95%)
tableS1 <- bayestestR::describe_posterior(ews.mod.detrend.mth.uni.roll, ci = 0.95, test="none")[-1,-3] 

tableS2 <- bayestestR::describe_posterior(ews.mod.detrend.mth.uni.exp, ci = 0.95, test="none")[-1,-3] 

tableS3 <- bayestestR::describe_posterior(ews.mod.detrend.mth.multi.roll, ci = 0.95, test="none")[-1,-3] 

tableS4 <- bayestestR::describe_posterior(ews.mod.detrend.mth.multi.exp, ci = 0.95, test="none")[-1,-3] 

tableS5 <- bayestestR::describe_posterior(ews.mod.detrend.mth.ml, ci = 0.95, test="none")[-1,-3] 

############ 
#Export Tables
############ 

write.csv(tableS1,"Results/supplementary_tables/table_S1.csv")
write.csv(tableS2,"Results/supplementary_tables/table_S2.csv")
write.csv(tableS3,"Results/supplementary_tables/table_S3.csv")
write.csv(tableS4,"Results/supplementary_tables/table_S4.csv")
write.csv(tableS5,"Results/supplementary_tables/table_S5.csv")

