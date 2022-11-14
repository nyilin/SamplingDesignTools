library(SamplingDesignTools)
# Load cohort data
data("cohort_1")
head(cohort_1)
# Draw simple 1:2 more extreme case-control sample, matched on gender.
# Let cases be subjects who had the event within 5 years, and controls be
# selected from those who did not have the event until the 15-th year.
set.seed(1)
dat_mecc <- draw_mecc(cohort = cohort_1, tau0 = 5, tau = 15,
                      id_name = "id", t_name = "t", delta_name = "y",
                      match_var_names = "gender", n_per_case = 2)
head(dat_mecc)
# To estimate the HR of age from MECC sample using the weighted approach,
# it is necessary to center age at the cohort average:
dat_mecc$age_c <- dat_mecc$age - mean(cohort_1$age)
result_mecc <- analyse_mecc_cond(
  y_name = "y_mecc", x_formula = ~ age_c, set_id_name = "set_id_mecc",
  surv = dat_mecc$surv, surv_tau = dat_mecc$surv_tau, mecc = dat_mecc,
  lower = -1, upper = 1
)
round(result_mecc$coef_mat[, -1], 3)
# Compare with the estimate from the full cohort:
library(survival)
result_cohort <- summary(coxph(Surv(t, y) ~ age + gender, data = cohort_1))$coef
round(result_cohort["age", ], 3)
# The MECC sample may also be analysed using a logistic regression to
# estimate the OR of age, which tends to overestimate the HR:
result_logit <- summary(glm(y_mecc ~ age + gender, data = dat_mecc))$coef
round(result_logit, 3)
