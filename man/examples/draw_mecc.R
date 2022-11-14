library(SamplingDesignTools)
# Load cohort data
data("cohort_1")
head(cohort_1)
# Draw simple 1:2 more extreme case-control sample, matched on gender.
# Let cases be subjects who had the event within 5 years, and controls be
# selected from those who did not have the event until the 15-th year.
dat_mecc <- draw_mecc(cohort = cohort_1, tau0 = 5, tau = 15,
                      id_name = "id", t_name = "t", delta_name = "y",
                      match_var_names = "gender", n_per_case = 2)
head(dat_mecc)
# Note that the new event indicator, `y_mecc`, is different from the original
# event, `y`, in the cohort:
identical(dat_mecc$y, dat_mecc$y_mecc) # Expect FALSE
