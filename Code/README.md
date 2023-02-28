# ews-assessments Code

## Recommended code workflow

[*lake_state_spaces.R*](lake_state_spaces.R) - from the raw plankton densities, fit threshold generalised additive models to both plankton time series and environmental state spaces. Coherence in break points between the two models indicates a critical transition.

[*plankton_ews_assessments.R*](plankton_ews_assessments.R) - from the processed plankton densities, assess each time series using the five forms of early warning signal (EWS): univariate rolling windows, univariate expanding windows, multivariate rolling windows, multivariate expanding windows and machine learning (represented by [EWSNet](https://doi.org/10.1098/rsos.211475)).

[*ews_success_analysis.R*](ews_success_analysis.R) - extract the success of each EWS using [*extract_ews_pred_fn.R*](extract_ews_pred_fn.R).

[*bayesian_success_models_ML.R*](bayesian_success_models_ML.R) - fit binomial models to EWS success rates and estimate probabilities of correct classification.

[*Figure2.R*](Figure2.R), [*Figure3.R*](Figure3.R), [*Figure4.R*](Figure4.R), [*supplementary_figures.R*](supplementary_figures.R) - generate figures.

## Supporting functions
[*threshold_gam.R*](threshold_gam.R) - collection of supporting functions to fit threshold generalised additive models.

[*perm_rollEWS_fn.R*](perm_rollEWS_fn.R) - a modfication of the `uniEWS(method = "rolling")` and `multiEWS(method = "rolling")` functions from the [EWSmethods R package](https://www.authorea.com/doi/full/10.22541/au.166801190.00303336) which performs permutations to assess the significance of rolling EWSs.

[*ewsnet_predict_impulse.R*](ewsnet_predict_impulse.R) - a modfication of the `ewsnet_predict()` function from the [EWSmethods R package](https://www.authorea.com/doi/full/10.22541/au.166801190.00303336) to use the Impulse version of the EWSNet machine learning model. These weights are provided in this [repository](https://github.com/duncanobrien/ews-assessments/tree/main/python/weights/Pretrained).

[*extract_ews_pred_fn.R*](extract_ews_pred_fn.R) - function to extract the prediction made by each EWS indicator for each EWS computation method.
