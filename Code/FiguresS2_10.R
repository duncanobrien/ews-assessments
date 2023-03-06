########################################################################################################################
# Supplementary Figures 2 - 10 #
########################################################################################################################
require(patchwork)
require(brms)

ews_mod_method_mth <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_mth.rds")
ews_mod_method_yr <- readRDS(file = "Results/ews_models/computation_models/ews_mod_method_yr.rds")

ews_mod_ind_mth_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_true.rds")
ews_mod_ind_yr_true <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_true.rds")
ews_mod_ind_mth_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_mth_false.rds")
ews_mod_ind_yr_false <- readRDS(file = "Results/ews_models/indicator_models/ews_mod_ind_yr_false.rds")

############ 
#Extract model trace plots
############ 

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S2.pdf",
    width = 5,height = 6)
  plot(ews_mod_method_mth,plot = TRUE,theme = theme_bw(),N=7,newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S3.pdf",
    width = 5,height = 6)
plot(ews_mod_method_yr,plot = TRUE,theme = theme_bw(),N=7,newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S4a.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_mth_true,plot = TRUE,theme = theme_bw(),variable = rownames(fixef(ews_mod_ind_mth_true))[1:16],N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S4b.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_mth_true,plot = TRUE,theme = theme_bw(),variable = c(rownames(fixef(ews_mod_ind_mth_true))[17:31],"sd_lake__Intercept"),N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S5a.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_yr_true,plot = TRUE,theme = theme_bw(),variable = rownames(fixef(ews_mod_ind_yr_true))[1:16],N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S5b.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_yr_true,plot = TRUE,theme = theme_bw(),variable = c(rownames(fixef(ews_mod_ind_yr_true))[17:31],"sd_lake__Intercept"),N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S6a.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_mth_false,plot = TRUE,theme = theme_bw(),variable = rownames(fixef(ews_mod_ind_mth_false))[1:16],N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S6b.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_mth_false,plot = TRUE,theme = theme_bw(),variable = c(rownames(fixef(ews_mod_ind_mth_false))[17:31],"sd_lake__Intercept"),N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S7a.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_yr_false,plot = TRUE,theme = theme_bw(),variable = rownames(fixef(ews_mod_ind_yr_false))[1:16],N=16,regex = T, newpage = FALSE)
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S7b.pdf",
    width = 5,height = 11)
plot(ews_mod_ind_yr_false,plot = TRUE,theme = theme_bw(),variable = c(rownames(fixef(ews_mod_ind_yr_false))[17:31],"sd_lake__Intercept"),N=16,regex = T, newpage = FALSE)
dev.off()

############ 
#Plot posterior predictive checks
############ 

ggsave(
  brms::pp_check(ews_mod_method_mth) + 
    ggtitle("Monthly computation\nmodel") +
    brms::pp_check(ews_mod_method_yr) + 
    ggtitle("Yearly computation\nmodel") +
    patchwork::plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'A') &
    theme_bw(),
  filename = "Results/supplementary_info/model_diagnoses/figure_S8.pdf",
  width = 6,height = 4)

ggsave(
  brms::pp_check(ews_mod_ind_mth_true) + 
    ggtitle("Monthly indicator\nmodel: true positive") +
    brms::pp_check(ews_mod_ind_yr_true) + 
    ggtitle("Yearly indicator\nmodel: true positive") +
    patchwork::plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'A') &
    theme_bw(),
  filename = "Results/supplementary_info/model_diagnoses/figure_S9.pdf",
  width = 6,height = 4)

ggsave(
  brms::pp_check(ews_mod_ind_mth_false) + 
    ggtitle("Monthly indicator\nmodel: true negative") +
    brms::pp_check(ews_mod_ind_yr_false) + 
    ggtitle("Yearly indicator\nmodel: true negative") +
    patchwork::plot_layout(guides = 'collect') + 
    plot_annotation(tag_levels = 'A') &
    theme_bw(),
  filename = "Results/supplementary_info/model_diagnoses/figure_S10.pdf",
  width = 6,height = 4)
