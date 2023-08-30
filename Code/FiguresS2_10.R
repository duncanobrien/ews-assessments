########################################################################################################################
# Supplementary Figures 2 - 10 #
########################################################################################################################

require(patchwork)
require(brms)
require(ggplot2)

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
    width = 5,height = 6,onefile = F)
ews_mod_method_mth$fit %>% 
  setNames(gsub("b_method_code|__Intercept", "", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:8],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S3.pdf",
    width = 5,height = 6,onefile = F)
ews_mod_method_yr$fit %>% 
  setNames(gsub("b_method_code|__Intercept", "", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:8],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S4a.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_mth_true$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:20],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S4b.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_mth_true$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[21:39],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S5a.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_yr_true$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:20],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S5b.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_yr_true$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[21:39],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S6a.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_mth_false$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:20],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S6b.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_mth_false$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[21:39],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S7a.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_yr_false$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[1:20],
                        gg_theme = theme_bw())
dev.off()

pdf(file = "Results/supplementary_info/model_diagnoses/figure_S7b.pdf",
    width = 5,height = 14,onefile = F)
ews_mod_ind_yr_false$fit %>% 
  setNames(gsub("b_indicator", "", names(.))) %>%
  setNames(gsub("P", " + ", names(.))) %>%
  bayesplot::mcmc_combo(pars = names(.)[21:39],
                        gg_theme = theme_bw())
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
