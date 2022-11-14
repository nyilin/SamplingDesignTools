library(SamplingDesignTools)
data("cohort_1")
head(cohort_1)
# Counter-match on binary indicator for age:
cohort_1$age_bin <- as.numeric(cohort_1$age < 50)
ncc_cm_bin <- draw_ncc_cm(cohort = cohort_1, y_name = "y", t_name = "t",
                          match_var_name = "age_bin",
                          include_var_name = c("age", "gender"), ml = 1)
head(ncc_cm_bin, 10)