data("cohort_2")
data("ncc_2")
data("n_at_risk")

sample_stat <- numeric(nrow(cohort_2))
sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
ncc_nodup_old <- compute_km_weights_old(
  cohort = cohort_2, t_name = "t", y_name = "y",
  sample_stat = sample_stat, n_per_case = 5,
  match_var_names = c("age_cat", "gender")
)
ncc_nodup <- SamplingDesignTools::compute_km_weights(
  cohort = cohort_2, t_name = "t", 
  sample_stat = sample_stat, n_per_case = 5,
  match_var_names = c("age_cat", "gender")
)

ncc_nodup2 <- ncc_2 %>% select(-Set) %>% 
  SamplingDesignTools::compute_km_weights(
    ncc = ., id_name = "Map", risk_table_manual = risk_table, 
    t_name = "t", t_match_name = "Time", y_name = "Fail",
    match_var_names = c("age_cat", "gender"), n_per_case = 5
  )
test_that("compute_km_weights: cohort", {
  expect_equal(ncc_nodup_old$km_weight, ncc_nodup$.km_weight)
  expect_equal(ncc_nodup$.km_weight, ncc_nodup2$.km_weight)
})
ncc_nodup_old_coarse <- compute_km_weights_old(
  ncc = ncc_2, n_at_risk = n_at_risk,
  t_name = "t", y_name = "Fail",
  t_match_name = "Time",
  id_name = "Map", set_id_name = "Set",
  match_var_names = c("age_cat", "gender"),
  n_per_case = 5
)
risk_table_coarse <- data.frame(t_event = n_at_risk$Time, n_event = 1, 
                                n_at_risk = n_at_risk$n.risk, 
                                gender = n_at_risk$gender, 
                                age_cat = n_at_risk$age_cat)
ncc_nodup_coarse <- ncc_2 %>% select(-Set) %>% 
  SamplingDesignTools::compute_km_weights(
    ncc = ., id_name = "Map", risk_table_manual = risk_table_coarse, 
    t_name = "t", t_match_name = "Time", y_name = "Fail",
    match_var_names = c("age_cat", "gender"), n_per_case = 5
  )
test_that("compute_km_weights: ncc", {
  expect_equal(ncc_nodup_old_coarse$km_weight, ncc_nodup_coarse$.km_weight)
})
library(survival)
summary(coxph(Surv(t, y) ~ x * z + age + gender, data = cohort_2))
#             coef exp(coef)  se(coef)      z Pr(>|z|)    
# x       0.382501  1.465946  0.109906  3.480 0.000501 ***
# z       1.495078  4.459686  0.126331 11.835  < 2e-16 ***
# age     0.007139  1.007165  0.001898  3.762 0.000169 ***
# gender -0.086074  0.917526  0.038010 -2.265 0.023542 *  
# x:z     0.640698  1.897805  0.135081  4.743 2.11e-06 ***
summary(clogit(Fail ~ x * z + age + strata(Set), data = ncc_2))
#         coef exp(coef) se(coef)      z Pr(>|z|)    
# x   0.395587  1.485256 0.113481  3.486  0.00049 ***
# z   1.519171  4.568434 0.135113 11.244  < 2e-16 ***
# age 0.012406  1.012483 0.008178  1.517  0.12930    
# x:z 0.586911  1.798425 0.144971  4.048 5.16e-05 ***
summary(coxph(Surv(t, y) ~ x * z + age + gender, data = ncc_nodup, 
              weights = .km_weight))
#             coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
# x       0.391681  1.479466  0.109909  0.112464  3.483 0.000496 ***
# z       1.523390  4.587751  0.126334  0.132469 11.500  < 2e-16 ***
# age     0.007497  1.007525  0.001915  0.002291  3.272 0.001067 ** 
# gender -0.098255  0.906418  0.038017  0.044850 -2.191 0.028471 *  
# x:z     0.609781  1.840028  0.135081  0.142189  4.289  1.8e-05 ***
summary(coxph(Surv(t, y) ~ x * z + age + gender, data = ncc_nodup_old, 
              weights = km_weight))
#             coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
# x       0.391681  1.479466  0.109909  0.112464  3.483 0.000496 ***
# z       1.523390  4.587751  0.126334  0.132469 11.500  < 2e-16 ***
# age     0.007497  1.007525  0.001915  0.002291  3.272 0.001067 ** 
# gender -0.098255  0.906418  0.038017  0.044850 -2.191 0.028471 *  
# x:z     0.609781  1.840028  0.135081  0.142189  4.289  1.8e-05 ***
summary(coxph(Surv(t, Fail) ~ x * z + age + gender, data = ncc_nodup_coarse, 
              weights = .km_weight))
#             coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
# x       0.391681  1.479466  0.109909  0.112464  3.483 0.000496 ***
# z       1.523390  4.587751  0.126334  0.132469 11.500  < 2e-16 ***
# age     0.007497  1.007525  0.001915  0.002291  3.272 0.001067 ** 
# gender -0.098255  0.906418  0.038017  0.044850 -2.191 0.028471 *  
# x:z     0.609781  1.840028  0.135081  0.142189  4.289  1.8e-05 ***
summary(coxph(Surv(t, Fail) ~ x * z + age + gender, data = ncc_nodup_old_coarse, 
              weights = km_weight))
#             coef exp(coef)  se(coef) robust se      z Pr(>|z|)    
# x       0.389239  1.475857  0.109910  0.111693  3.485 0.000492 ***
# z       1.503388  4.496900  0.126335  0.130570 11.514  < 2e-16 ***
# age     0.005943  1.005961  0.001917  0.002165  2.745 0.006056 ** 
# gender -0.061402  0.940445  0.038018  0.042312 -1.451 0.146728    
# x:z     0.551678  1.736164  0.135085  0.139919  3.943 8.05e-05 ***


ncc_controls <- ncc_2 %>% filter(Fail == 0) %>% select(-Set, -Time) %>% distinct() %>%
  SamplingDesignTools::compute_km_weights_controls(
    ncc_controls = ., prob_table_manual = prob_table, t_name = "t", 
    match_var_names = c("age_cat", "gender")
  )
ncc_nodup <- SamplingDesignTools::compute_km_weights(
  cohort = cohort_2, t_name = "t", sample_stat = sample_stat, 
  match_var_names = c("age_cat", "gender"), n_per_case = 5, 
  attach_weight = TRUE
)
ncc_controls0 <- ncc_nodup %>% filter(y == 0)
ncc_controls2 <- ncc_controls %>% filter(Map %in% ncc_controls0$id) %>% 
  arrange(Map)
test_that("km_weights for new controls", {
  expect_equal(ncc_controls2$Map, ncc_controls0$id)
  expect_equal(ncc_controls2$t, ncc_controls0$t)
  expect_equal(ncc_controls2$.km_weight, ncc_controls0$.km_weight)
})

sample_stat <- numeric(nrow(cohort_2))
sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
#' # Then find the sampled subset and compute weights:
ncc_nodup <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y",
                                match_var_names = c("age_cat", "gender"),
                                sample_stat = sample_stat, n_per_case = 5)
head(ncc_nodup)
#' # Alternatively, if the cohort is not available, the weights can be computed 
#' # as long as number of subjects at risk at event times in each strata is 
#' # available elsewhere, and the actual time of event/censoring is available 
#' # for each subject in the NCC.
#' # Compute the number of subjects at risk from cohort_2:
risk_table <- compute_risk_table(cohort = cohort_2, t_name = "t", y_name = "y",
                                 match_var_names = c("age_cat", "gender"))
head(risk_table)
#' # he following command computes the same weights as in ncc_nodup:
ncc_nodup_v2 <- compute_km_weights(ncc = ncc_2[, -1], risk_table_manual = risk_table,
                                   id_name = "Map", t_match_name = "Time",
                                   t_name = "t", y_name = "Fail",
                                   match_var_names = c("age_cat", "gender"),
                                   n_per_case = 5)
head(ncc_nodup_v2)
expect_equal(ncc_nodup$km_weight, ncc_nodup_v2$km_weight[sort.list(ncc_nodup_v2$Map)])
