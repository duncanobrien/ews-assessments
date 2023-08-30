# Supplementary code

Supporting code for the review process but not contributing to the manuscript.

[*descriptive_results.R*](descriptive_results.R) - script to generate the success probabilities reported in the Results section of the mansucript.

[*mth_lake_state_spaces.R*](mth_lake_state_spaces.R) - analagous script to [*lake_state_spaces.R*](https://github.com/duncanobrien/ews-assessments/tree/main/Code/lake_state_spaces.R) edited to fit threshold GAMs to monthly rather than yearly data. 

[*total_density_ews_assessment.R*](total_density_ews_assessment.R) - merged script combining [*plankton_ews_assessment.R*](https://github.com/duncanobrien/ews-assessments/tree/main/Code/plankton_ews_assessment.R) and [*bayesian_success_models.R*](https://github.com/duncanobrien/ews-assessments/tree/main/Code/bayesian_success_models.R) to repeat EWS assessments for trophic level aggregations rather than genus.

[*ROC_ability.R*](ROC_ability.R) - calculates Receiver Operator Curve and precision-recall based indices for the binary critical-no critical transtiion classification task performed by EWSs. Generates comparable results to [*bayesian_success_models.R*](https://github.com/duncanobrien/ews-assessments/tree/main/Code/bayesian_success_models.R).
