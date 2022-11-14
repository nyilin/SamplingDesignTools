library(SamplingDesignTools)
# Load cohort data
data("cohort_2") # Ignoring that we have the event time in this simulated cohort
head(cohort_2)
# Draw simple 1:2 case-control sample (i.e., with randomly selected controls):
dat_scc <- draw_mcc(cohort = cohort_2, y_name = "y", n_per_case = 2)
head(dat_scc)
# Draw simple 1:2 case-control sample, matched on age group and gender:
dat_mcc <- draw_mcc(cohort = cohort_2, y_name = "y", n_per_case = 2,
                    match_var_names = c("age_cat", "gender"),
                    weight_name = ".w")
head(dat_mcc)
