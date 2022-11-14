library(SamplingDesignTools)
# Load mini cohort
data("mini_cohort")
mini_cohort
# For illustrative purpose, assume subjects 5, 7 and 10 are newly collected
# controls, so that the resulting Kaplan-Meier type weights are identical to
# examples for compute_km_weights():
ncc_controls <- mini_cohort[c(5, 7, 10), c("id", "t")]
ncc_controls
# To use these newly sampled controls with existing cases, we need to have the 
# number of subjects at risk at event times in mini_cohort:
risk_table <- compute_risk_table(cohort = mini_cohort, t_name = "t", 
                                 y_name = "status")
# The following command computes the Kaplan-Meier type weights for the newly
# sampled controls:
ncc_controls <- compute_km_weights_controls(
  ncc_controls = ncc_controls, risk_table_manual = risk_table, 
  t_name = "t", n_per_case = 1
)
ncc_controls
