library(SamplingDesignTools)
# With common entry time:
data("mini_cohort")
mini_cohort
compute_risk_table(cohort = mini_cohort, t_name = "t", y_name = "status")
# With staggered entry time:
data("mini_cohort2")
mini_cohort2
compute_risk_table(cohort = mini_cohort2, t_start_name = "t_start", 
                   t_name = "t_end", y_name = "status")
