# ews-assessments Code

## Recommended code workflow

[*aggregated_plank_genus.R*](aggregated_plank_genus.R) - examples of how plankton data is wrangled prior to early warning signal assessment.

[*lake_state_spaces.R*](lake_state_spaces.R) - from the raw plankton densities, fit threshold generalised additive models to both plankton time series and environmental state spaces. Coherence in break points between the two models indicates a critical transition.

[*plankton_ews_assessment.R*](plankton_ews_assessment.R) - from the processed plankton densities, assess each time series using the five forms of early warning signal (EWS): univariate rolling windows, univariate expanding windows, multivariate rolling windows, multivariate expanding windows and machine learning (represented by [EWSNet](https://doi.org/10.1098/rsos.211475)).

[*ews_success_analysis.R*](ews_success_analysis.R) - extract the success of each EWS using [*extract_ews_pred_fn.R*](extract_ews_pred_fn.R).

[*bayesian_success_models.R*](bayesian_success_models.R) - fit binomial models to EWS success rates and estimate probabilities of correct classification. The saved models can be found in the [model folder](https://github.com/duncanobrien/ews-assessments/tree/main/Results/ews_models).

[*Figure4.R*](Figure4.R), [*Figure5.R*](Figure5.R), [*Figure6.R*](Figure6.R), [*FigureS1.R*](FigureS1.R), [*FiguresS2_10.R*](FiguresS2_10.R),[*TablesS3_7.R*](TablesS3_7.R), [*TablesS8_13.R*](TablesS8_13.R) - generate the main figures available [here](https://github.com/duncanobrien/ews-assessments/tree/main/Figures), the [supplementary figures](https://github.com/duncanobrien/ews-assessments/tree/main/Results/supplementary_info/model_diagnoses), and supplementary tables [here](https://github.com/duncanobrien/ews-assessments/tree/main/Results/supplementary_info/supplementary_tables).

## Supporting functions
[*threshold_gam.R*](threshold_gam.R) - collection of supporting functions to fit threshold generalised additive models.

[*ewsnet_predict_impulse.R*](ewsnet_predict_impulse.R) - a modification of the `ewsnet_predict()` function from the [EWSmethods R package](https://doi.org/10.1111/ecog.06674) to use the Impulse version of the EWSNet machine learning model. These weights are provided in this [repository](https://github.com/duncanobrien/ews-assessments/tree/main/python/weights/Pretrained).

[*extract_ews_pred_fn.R*](extract_ews_pred_fn.R) - function to extract the prediction made by each EWS indicator for each EWS computation method.

## Supplementary code
The [supplementary_code folder](https://github.com/duncanobrien/ews-assessments/tree/main/Code/supplementary_code) contains additional scripts to support statements made during the review process regarding classification of regime shifts in monthly time series and the use of Receiver Operator Curves (ROC) for assessing EWS ability.