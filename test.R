library(survival)
library(dplyr)
data(cohort_1)
head(cohort_1)
data(ncc_1)
head(ncc_1)
data(cohort_2)
head(cohort_2)
data(ncc_2)
head(ncc_2)

# -----

cohort_1$x <- as.numeric(cohort_1$age > 60)
table(cohort_1$x)
ncc_cm <- draw_ncc_cm(cohort = cohort_1, y_name = "y", t_name = "t", 
                      match_var_name = "x", include_var_name = c("age", "gender"))
head(ncc_cm)
m_clogit <- clogit(y ~ age + gender + log(weight) + strata(set), data = ncc_cm)
summary(m_clogit)$coef

sample_stat <- numeric(nrow(cohort_1))
sample_stat[unique(ncc_cm$row_id[ncc_cm$y == 0])] <- 1
sample_stat[ncc_cm$row_id[ncc_cm$y == 1]] <- 2
table(sample_stat, x = cohort_1$x, useNA = "always")
table(sample_stat, y = cohort_1$y, useNA = "always")
ncc_cm_nodup <- compute_km_weights_cm(
  cohort = cohort_1, t_name = "t", y_name = "y", 
  sample_stat = sample_stat, match_var_name = "x", ml = 1
)
m_cox <- coxph(Surv(t, y) ~ age + gender, weights = km_weight, robust = TRUE, 
               data = ncc_cm_nodup)
summary(m_cox)$coef

# -----

ncc_2$t <- cohort_2$t[ncc_2$Map]
km <- survfit(Surv(t, y) ~ strata(gender, age_cat), data = cohort_2)
km_summ <- summary(km)
n_at_risk <- data.frame(strata = as.character(km_summ$strata), 
                        Time = km_summ$time, 
                        n.risk = km_summ$n.risk, stringsAsFactors = FALSE)
n_at_risk <- cbind(n_at_risk[, -1], do.call("rbind", lapply(n_at_risk$strata, function(s) {
  s_vec <- unlist(strsplit(s, split = "="))
  gender <- as.numeric(unlist(strsplit(s_vec[3], split = ","))[1])
  data.frame(gender = gender, age_cat = trimws(s_vec[4]), stringsAsFactors = FALSE)
})))

obj2 <- compute_kmw2(ncc = ncc_2, n_at_risk = n_at_risk, 
                     id_name = "Map", set_id_name = "Set",
                     t_match_name = "Time", y_name = "Fail", 
                     match_var_names = c("age_cat", "gender"), 
                     n_per_case = 5)
summary(obj2$km_tb) # same as obj1$km_tb

sample_stat <- numeric(nrow(cohort_2))
sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
obj1 <- compute_kmw1(cohort = cohort_2, t_name = "t", y_name = "y", 
                     sample_stat = sample_stat, 
                     match_var_names = c("age_cat", "gender"), n_per_case = 5)
summary(obj1$km_tb)

ncc_nodup1 <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y", 
                                 sample_stat = sample_stat, 
                                 match_var_names = c("age_cat", "gender"), 
                                 n_per_case = 5)
ncc_nodup2 <- compute_km_weights(ncc = ncc_2, n_at_risk = n_at_risk, 
                                 t_name = "t", y_name = "Fail", t_match_name = "Time",
                                 id_name = "Map", set_id_name = "Set", 
                                 match_var_names = c("age_cat", "gender"), 
                                 n_per_case = 5)
ncc_nodup1 <- ncc_nodup1 %>% arrange(id)
ncc_nodup2 <- ncc_nodup2 %>% arrange(Map)
all.equal(ncc_nodup1$id, ncc_nodup2$Map) # TRUE
all.equal(ncc_nodup1$t, ncc_nodup2$t) # TRUE
summary(ncc_nodup1$km_prob)
summary(ncc_nodup2$km_prob)
plot(x = ncc_nodup1$km_prob[ncc_nodup1$km_prob != 1], 
     y = ncc_nodup2$km_prob[ncc_nodup2$km_prob != 1])
abline(coef = c(0, 1), col = "red")
all.equal(ncc_nodup1$km_weight, ncc_nodup2$km_weight) # TRUE

n_at_risk <- prep_n_at_risk(ncc = ncc_2, t_match_name = "t_yr", y_name = "Fail", 
                            match_var_names = c("age_cat", "gender"))

# -----

ncc_2$t <- cohort_2$t[ncc_2$Map]
ncc_2$t_yr <- ceiling(ncc_2$t)
cohort_2$t_yr <- ceiling(cohort_2$t)
km <- survfit(Surv(t_yr, y) ~ strata(gender, age_cat), data = cohort_2)
km_summ <- summary(km)
n_at_risk <- data.frame(strata = as.character(km_summ$strata), 
                        t_yr = km_summ$time, 
                        n.risk = km_summ$n.risk, stringsAsFactors = FALSE)
n_at_risk <- cbind(n_at_risk[, -1], do.call("rbind", lapply(n_at_risk$strata, function(s) {
  s_vec <- unlist(strsplit(s, split = "="))
  gender <- as.numeric(unlist(strsplit(s_vec[3], split = ","))[1])
  data.frame(gender = gender, age_cat = trimws(s_vec[4]), stringsAsFactors = FALSE)
})))
ncc_nodup2 <- compute_km_weights(ncc = ncc_2, n_at_risk = n_at_risk, 
                                 t_name = "t", y_name = "Fail", t_match_name = "t_yr",
                                 id_name = "Map", set_id_name = "Set", 
                                 match_var_names = c("age_cat", "gender"), 
                                 n_per_case = 5)

m <- coxph(Surv(t, y) ~ x * z + age + gender, data = cohort_2)
round(summary(m)$coef, 3)
#          coef exp(coef) se(coef)      z Pr(>|z|)
# x       0.383     1.466    0.110  3.480    0.001
# z       1.495     4.460    0.126 11.835    0.000
# age     0.007     1.007    0.002  3.762    0.000
# gender -0.086     0.918    0.038 -2.265    0.024
# x:z     0.641     1.898    0.135  4.743    0.000
m_clogit <- clogit(Fail ~ x * z + strata(Set), data = ncc_2)
round(summary(m_clogit)$coef, 3)
#      coef exp(coef) se(coef)      z Pr(>|z|)
# x   0.397     1.488    0.113  3.502        0
# z   1.518     4.563    0.135 11.237        0
# x:z 0.586     1.797    0.145  4.045        0
m_cox <- coxph(Surv(Time, Fail) ~ x * z + age + gender, data = ncc_nodup2, 
               weights = km_weight, robust = TRUE)
round(summary(m_cox)$coef, 3)
#          coef exp(coef) se(coef) robust se      z Pr(>|z|)
# x       0.392     1.480    0.110     0.110  3.559    0.000
# z       1.438     4.214    0.126     0.126 11.388    0.000
# age     0.000     1.000    0.002     0.002  0.148    0.882
# gender -0.020     0.980    0.038     0.038 -0.526    0.599
# x:z     0.320     1.377    0.135     0.135  2.366    0.018

# -----

# First generate cohort_3 and ncc_3 in the data_generating script
load("../cohort_3.RData")
load("../ncc_3.RData")
ncc_3$t <- cohort_3$t[ncc_3$Map]
ncc_3$t_yr <- ceiling(ncc_3$t)
cohort_3$t_yr <- ceiling(cohort_3$t)
km <- survfit(Surv(t_yr, y) ~ strata(gender, age_cat), data = cohort_3)
km_summ <- summary(km)
n_at_risk <- data.frame(strata = as.character(km_summ$strata), 
                        t_yr = km_summ$time, 
                        n.risk = km_summ$n.risk, stringsAsFactors = FALSE)
n_at_risk <- cbind(n_at_risk[, -1], do.call("rbind", lapply(n_at_risk$strata, function(s) {
  s_vec <- unlist(strsplit(s, split = "="))
  gender <- as.numeric(unlist(strsplit(s_vec[3], split = ","))[1])
  data.frame(gender = gender, age_cat = trimws(s_vec[4]), stringsAsFactors = FALSE)
})))

sample_stat <- numeric(nrow(cohort_3))
sample_stat[unique(ncc_3$Map[ncc_3$Fail == 0])] <- 1
sample_stat[ncc_3$Map[ncc_3$Fail == 1]] <- 2
obj_auto <- SamplingDesignTools:::prep_km1(cohort = cohort_3, 
                                           t_name = "t",  y_name = "y",
                                           sample_stat = sample_stat,
                                           match_var_names = c("age_cat", "gender"), 
                                           n_per_case = 1)
km_tb <- obj$km_tb
km_tb$t_yr <- ceiling(km_tb$t)
unique(obj$match_var_ncc[obj$ncc_nodup$age_cat == "(-Inf,35]" & 
                           obj$ncc_nodup$gender == 0])
match_var_mat <- unique(data.frame(strata = as.character(obj$match_var_ncc), 
                                   gender = obj$ncc_nodup$gender, 
                                   age_cat = obj$ncc_nodup$age_cat, 
                                   stringsAsFactors = FALSE))
km_tb <- km_tb %>% left_join(match_var_mat) %>% rename(n.risk_manual = n.risk) %>% 
  left_join(n_at_risk)
plot(x = km_tb$n.risk_manual, y = km_tb$n.risk, pch = 20, cex = 0.5)
abline(coef = c(0, 1), col = "red")

km_tb %>% filter(t_yr == 1, strata == "match_var=1") %>% 
  select(n.risk) %>% summary()

ncc_nodup3 <- compute_km_weights(ncc = ncc_3, n_at_risk = n_at_risk, 
                                 t_name = "t", y_name = "Fail", t_match_name = "t_yr",
                                 id_name = "Map", set_id_name = "Set", 
                                 match_var_names = c("age_cat", "gender"), 
                                 n_per_case = 1)
ncc_nodup3_auto <- compute_km_weights(
  cohort = cohort_3, t_name = "t", y_name = "y", sample_stat = sample_stat, 
  match_var_names = c("gender", "age_cat"), n_per_case = 1
)
df <- data.frame(w = ncc_nodup3$km_weight, w_auto = ncc_nodup3_auto$km_weight)

plot(x = df$w_auto, y = df$w, pch = 20, cex = 0.5)
abline(coef = c(0, 1), col = "red")

summary(df$w_auto)
summary(df$w)

m_cox <- coxph(Surv(Time, Fail) ~ x * z + age + gender, data = ncc_nodup3, 
               weights = km_weight, robust = TRUE)
summary(m_cox)
#        exp(coef) exp(-coef) lower .95 upper .95
# x         1.5571     0.6422    1.4412    1.6824
# z         3.3051     0.3026    3.0061    3.6339
# age       1.0037     0.9963    1.0021    1.0052
# gender    0.9675     1.0336    0.9386    0.9972
# x:z       1.3455     0.7432    1.2163    1.4884
m_clogit <- clogit(Fail ~ x * z + strata(Set), data = ncc_3)
summary(m_clogit)
#     exp(coef) exp(-coef) lower .95 upper .95
# x       1.606     0.6226     1.476     1.748
# z       3.966     0.2521     3.551     4.430
# x:z     2.016     0.4960     1.789     2.271
m <- coxph(Surv(t, y) ~ x * z + age + gender, data = cohort_3)
summary(m)
#        exp(coef) exp(-coef) lower .95 upper .95
# x          1.622     0.6164     1.514     1.738
# z          4.073     0.2455     3.758     4.414
# age        1.010     0.9905     1.008     1.011
# gender     1.034     0.9676     1.010     1.058
# x:z        1.952     0.5122     1.792     2.126


ncc <- ncc_3
obj <- SamplingDesignTools:::prep_km2(ncc = ncc, n_at_risk = n_at_risk, 
                                      id_name = "Map", set_id_name = "Set",
                                      t_match_name = "t_yr", y_name = "Fail", 
                                      match_var_names = c("age_cat", "gender"), 
                                      n_per_case = 1)
obj_auto <- SamplingDesignTools:::prep_km1(cohort = cohort_3, 
                                           t_name = "t",  y_name = "y",
                                           sample_stat = sample_stat,
                                           match_var_names = c("age_cat", "gender"), 
                                           n_per_case = 1)
km_tb_auto <- obj_auto$km_tb %>%
  group_by(strata) %>% 
  mutate(prob_not_sampled = 1 - (n_per_case / Rj),
         cumulative_product = cumprod(prob_not_sampled),
         sampling_prob = 1 - cumulative_product)
p_ncc_auto <- unlist(lapply(1:nrow(obj_auto$ncc_nodup), function(j) {
  if (obj_auto$ncc_nodup[j, y_name] == 1) {
    1
  } else {
    km_tb_i_auto <- km_tb_auto[km_tb_auto$t < obj_auto$ncc_nodup[j, t_name] & 
                                 km_tb_auto$strata == obj_auto$match_var_ncc[j], ]
    r <- nrow(km_tb_i_auto)
    if (r == 0) {
      0
    } else {
      km_tb_i_auto$sampling_prob[r]
    }
  }
}))

ncc_nodup3 <- cbind(obj$ncc_nodup, km_prob = p_ncc, km_weight = 1 / p_ncc)
