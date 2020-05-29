# 0 Preparation

Install the SamplingDesignTools package from GitHub (package devtools
needed):

``` r
# devtools::install_github("nyilin/SamplingDesignTools")
```

Load packages used:

``` r
library(SamplingDesignTools)
library(survival)
library(dplyr)
## 
## Attaching package: 'dplyr'
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
library(Epi) # To draw (non-counter-matched) nested case-control sample
```

# 1 Drawing counter-matched sample

``` r
data(cohort_1)
head(cohort_1)
##   id y        t age gender
## 1  1 0 25.00000  47      1
## 2  2 0 10.65152  58      0
## 3  3 0 25.00000  46      0
## 4  4 0 15.84131  52      0
## 5  5 0 22.57659  49      0
## 6  6 0 25.00000  63      0
m_cox_cohort_1 <- coxph(Surv(t, y) ~ age + gender, data = cohort_1)
```

In this example, I am using a simulated cohort with 10,000 subjects,
where a continuous exposure (x ~ N(55, 10^2)) represents their age and a
binary covariate indicates their gender (c = 1 for male, Pr(c=1)=0.5).

The event time t was generated from a Weibull distribution with hazard
function:

log{h(t; x, c)} = log{h\_0} + log(1.1) \* x + log(2) \* c,

where the baseline hazard h\_0=10^{-5} and the scale parameter was 1.
All individuals entered the cohort at time t=0 and were followed until
the event of interest or a censoring time (generated from an exponential
distribution with a rate of 0.05) with a maximum follow-up time of t=25.
There were 582 cases (5.8%) in this cohort.

## Counter-match on binary surrogate

In this analysis, I dichotomised the continuous exposure at 50 years of
age, and drew a 1:1 time-matched sample counter matched on this binary
variable that is considered as a surrogate for the exposure. This sample
was analysed using a weighted conditional logistic model, where log of
sampling weights were specified as offset in the model.

``` r
cohort_1$age_bin <- as.numeric(cohort_1$age < 50)
table(cohort_1$age_bin)
## 
##    0    1 
## 7181 2819
ncc_cm_bin <- draw_ncc_cm(cohort = cohort_1, y_name = "y", t_name = "t", 
                          match_var_name = "age_bin", 
                          include_var_name = c("age", "gender"), ml = 1)
head(ncc_cm_bin, 10)
##      set row_id         t n_at_risk n_sampled weight y age_bin age gender
## 11    11     11  6.952504      4866         1   4866 1       0  60      0
## 6682  11   6682  6.952504      2002         1   2002 0       1  41      1
## 29    29     29 11.157250      3802         1   3802 1       0  64      0
## 6543  29   6543 11.157250      1594         1   1594 0       1  40      0
## 58    58     58 24.578458      1860         1   1860 1       0  63      0
## 1434  58   1434 24.578458       810         1    810 0       1  43      1
## 67    67     67  1.434794      6613         1   6613 1       0  65      1
## 3810  67   3810  1.434794      2642         1   2642 0       1  48      1
## 89    89     89 19.254871      2457         1   2457 1       0  64      1
## 6004  89   6004 19.254871      1062         1   1062 0       1  46      1
table(ncc_cm_bin$age_bin, ncc_cm_bin$y)
##    
##       0   1
##   0  32 550
##   1 550  32
m_clogit_bin <- clogit(y ~ age + gender + strata(set) + offset(log(weight)), 
                       data = ncc_cm_bin)
```

## Counter-match on categorical surrogate

In this analysis, I categorised the continuous exposure into 4
categories at 40, 50 and 60 years of age, and drew a 1:3 time-matched
sample counter matched on this categorical variable that is considered
as a surrogate for the exposure. This sample was again analysed using a
weighted conditional logistic
model.

``` r
cohort_1$age_quart <- cut(cohort_1$age, breaks = c(-Inf, 40, 50, 60, Inf), 
                          labels = 1:4, include.lowest = TRUE)
table(cohort_1$age_quart)
## 
##    1    2    3    4 
##  712 2453 3908 2927
ncc_cm_quart <- draw_ncc_cm(cohort = cohort_1, y_name = "y", t_name = "t", 
                            match_var_name = "age_quart", 
                            include_var_name = c("age", "gender"), ml = 1)
head(ncc_cm_quart, 20)
##      set row_id         t n_at_risk n_sampled weight y age_quart age gender
## 5123  11   5123  6.952504       512         1    512 0         1  39      1
## 1050  11   1050  6.952504      1756         1   1756 0         2  42      0
## 11    11     11  6.952504      2672         1   2672 1         3  60      0
## 9857  11   9857  6.952504      1928         1   1928 0         4  63      0
## 4445  29   4445 11.157250       413         1    413 0         1  28      1
## 3770  29   3770 11.157250      1392         1   1392 0         2  42      1
## 4913  29   4913 11.157250      2116         1   2116 0         3  57      0
## 29    29     29 11.157250      1475         1   1475 1         4  64      0
## 9234  58   9234 24.578458       211         1    211 0         1  36      1
## 7715  58   7715 24.578458       704         1    704 0         2  50      1
## 3727  58   3727 24.578458      1066         1   1066 0         3  51      0
## 58    58     58 24.578458       689         1    689 1         4  63      0
## 7369  67   7369  1.434794       665         1    665 0         1  25      1
## 9614  67   9614  1.434794      2300         1   2300 0         2  43      1
## 1404  67   1404  1.434794      3607         1   3607 0         3  60      0
## 67    67     67  1.434794      2683         1   2683 1         4  65      1
## 170   89    170 19.254871       269         1    269 0         1  39      1
## 3947  89   3947 19.254871       933         1    933 0         2  50      0
## 8190  89   8190 19.254871      1395         1   1395 0         3  53      0
## 89    89     89 19.254871       922         1    922 1         4  64      1
table(ncc_cm_quart$age_quart)
## 
##   1   2   3   4 
## 582 582 582 582
m_clogit_quart <- clogit(y ~ age + gender + strata(set) + offset(log(weight)), 
                         data = ncc_cm_quart)
```

## Compare results

``` r
results_1 <- rbind(summary(m_cox_cohort_1)$coef, 
                   summary(m_clogit_bin)$coef, 
                   summary(m_clogit_quart)$coef)
rownames(results_1) <- NULL
kable(data.frame(
  Data = c("Full cohort", "", "1:1 NCC-CM", "", "1:3 NCC-CM", ""), 
  Variable = rep(c("Age", "Male"), 3), 
  `True HR` = rep(c(1.1, 2), 3),
  `Estimated HR` = results_1[, "exp(coef)"], 
  `SE of log(HR)` = results_1[, "se(coef)"], 
  `p-value` = results_1[, "Pr(>|z|)"], check.names = FALSE
), digits = c(0, 0, 1, 2, 3, 3))
```

| Data        | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
| :---------- | :------- | ------: | -----------: | ------------: | ------: |
| Full cohort | Age      |     1.1 |         1.11 |         0.004 |       0 |
|             | Male     |     2.0 |         2.18 |         0.088 |       0 |
| 1:1 NCC-CM  | Age      |     1.1 |         1.13 |         0.013 |       0 |
|             | Male     |     2.0 |         3.21 |         0.303 |       0 |
| 1:3 NCC-CM  | Age      |     1.1 |         1.11 |         0.006 |       0 |
|             | Male     |     2.0 |         2.56 |         0.133 |       0 |

**Notes:** 1:1 NCC-CM indicates the analysis of the time-matched sample
counter-matched on a binary surrogate. 1:3 NCC-CM indicates the analysis
of the time-matched sample counter-matched on a surrogate with 4
categories.

# 2 Dropping controls in NCC

``` r
data(cohort_2)
head(cohort_2)
##   id y         t x age   age_cat gender z
## 1  1 0 25.000000 1  -2   (45,55]      0 0
## 2  2 0 19.819801 1  -4   (45,55]      1 0
## 3  3 0 25.000000 1  -5   (45,55]      0 0
## 4  4 0 12.414616 1  20 (75, Inf]      1 0
## 5  5 0 25.000000 1  -2   (45,55]      0 1
## 6  6 0  1.019023 0 -15   (35,45]      1 0
m_cox_cohort_2 <- coxph(Surv(t, y) ~ x * z + age + gender, data = cohort_2)
```

In this example, I use a simulated cohort with 100,000 subjects with
age~N(55, 10^2), Pr(male=1)=0.5, a binary exposure (x) that depends on
gender (Pr(x = 1)=0.8 for female and 0.5 for male) and a binary effect
modifier (Pr(z = 1)=0.3). The event time t was generated from a Weibull
distribution with hazard function:

log{h(t;x, age, gender, z)} = log{h\_0} + log(1.5) \* x + log(1.1) \*
gender + log(4) \* z + log(2) \* x \* z + log(1.02) \* age,

where the baseline hazard h\_0=0.0005, and the scale parameter was 1.
All individuals entered the cohort at time t=0 and were followed until
the event of interest or a censoring time (generated from an exponential
distribution with a rate of 0.05) with a maximum follow-up time of t=25.
There were 2773 cases (2.8%) in this cohort.

## 1:5 NCC matched on age and gender

Subjects were divided into 6 age groups: \<35, 36-45, 46-55, 56-65,
66-75 and \>75. A 1:5 NCC sample was draw from the cohort, matched on
age group and gender.

``` r
n_per_case <- 5
ncc_2 <- ccwc(exit = t, fail = y, controls = n_per_case, 
              match = list(age_cat, gender), include = list(x, age, z), 
              data = cohort_2)
## 
## Sampling risk sets: .....................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................
names(ncc_2)[-(1:4)] <- c("age_cat", "gender", "x", "age", "z")
head(ncc_2, 12)
##    Set   Map      Time Fail age_cat gender x age z
## 1    1    36 0.2016047    1 (35,45]      1 1 -16 0
## 2    1 78786 0.2016047    0 (35,45]      1 0 -14 0
## 3    1 32614 0.2016047    0 (35,45]      1 1 -11 1
## 4    1 26970 0.2016047    0 (35,45]      1 1 -12 0
## 5    1 23770 0.2016047    0 (35,45]      1 1 -14 0
## 6    1 85830 0.2016047    0 (35,45]      1 1 -15 0
## 7    2   888 4.9218685    1 (35,45]      1 1 -20 1
## 8    2 65447 4.9218685    0 (35,45]      1 1 -11 0
## 9    2 85760 4.9218685    0 (35,45]      1 0 -17 1
## 10   2 56587 4.9218685    0 (35,45]      1 0 -12 0
## 11   2 41962 4.9218685    0 (35,45]      1 1 -13 0
## 12   2 80273 4.9218685    0 (35,45]      1 1 -12 0
# Create the sampling and status indicator
sample_stat <- numeric(nrow(cohort_2))
sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
table(sample_stat)
## sample_stat
##     0     1     2 
## 84916 12311  2773
ncc_2_nodup <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y",
                                  sample_stat = sample_stat, 
                                  match_var_names = c("age_cat", "gender"), 
                                  n_per_case = n_per_case)
head(ncc_2_nodup, 12)
##    id y          t x age age_cat gender z    km_prob km_weight
## 1   1 0 25.0000000 1  -2 (45,55]      0 0 0.21851027  4.576444
## 16 16 0 25.0000000 1  10 (65,75]      0 1 0.24202210  4.131854
## 35 35 0 16.5202891 1  -5 (45,55]      1 1 0.14312411  6.986943
## 36 36 1  0.2016047 1 -16 (35,45]      1 0 1.00000000  1.000000
## 54 54 0 25.0000000 0   3 (55,65]      1 1 0.21862668  4.574007
## 61 61 1 14.8098399 1   7 (55,65]      0 0 1.00000000  1.000000
## 69 69 0 25.0000000 1  -1 (45,55]      1 0 0.19666881  5.084690
## 72 72 1  4.1087132 1  -1 (45,55]      0 1 1.00000000  1.000000
## 76 76 1 20.1285244 1  15 (65,75]      0 1 1.00000000  1.000000
## 82 82 0  5.7436608 0   0 (45,55]      0 1 0.06188087 16.160084
## 84 84 1  5.9032539 1  -9 (45,55]      1 1 1.00000000  1.000000
## 86 86 1  6.0871929 1  -3 (45,55]      1 0 1.00000000  1.000000
summary(ncc_2_nodup$km_weight)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.000   4.218   4.805   6.547   6.436 482.572
m_cox_ncc_2 <- coxph(Surv(t, y) ~ x * z + gender + age, data = ncc_2_nodup,
                     weights = km_weight, robust = TRUE)
m_clogit_ncc_2 <- clogit(Fail ~ x * z + strata(Set), data = ncc_2)
```

## 1:5 NCC with 3 controls dropped

I dropped 3 controls from each sampled set, and only kept the
information on the exposure for the cases and controls left.

``` r
n_kept <- 2
# Only keep the case and the first 2 controls
i_case <- which(ncc_2$Fail == 1)
i_kept <- rep(i_case, each = n_kept + 1) + rep(0:n_kept, length(i_case))
ncc_2_drop3 <- ncc_2[i_kept, ]
keep_stat <- numeric(nrow(cohort_2))
keep_stat[ncc_2_drop3$Map] <- 1
table(sample_stat, keep_stat)
##            keep_stat
## sample_stat     0     1
##           0 84916     0
##           1  7092  5219
##           2     0  2773
ncc_2_nodup_drop3 <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y",
                                        sample_stat = sample_stat, keep_stat = keep_stat,
                                        match_var_names = c("age_cat", "gender"), 
                                        n_per_case = n_per_case, n_kept = n_kept)
head(ncc_2_nodup_drop3, 12)
##      id y          t x age age_cat gender z    km_prob km_weight
## 36   36 1  0.2016047 1 -16 (35,45]      1 0 1.00000000   1.00000
## 61   61 1 14.8098399 1   7 (55,65]      0 0 1.00000000   1.00000
## 72   72 1  4.1087132 1  -1 (45,55]      0 1 1.00000000   1.00000
## 76   76 1 20.1285244 1  15 (65,75]      0 1 1.00000000   1.00000
## 84   84 1  5.9032539 1  -9 (45,55]      1 1 1.00000000   1.00000
## 86   86 1  6.0871929 1  -3 (45,55]      1 0 1.00000000   1.00000
## 116 116 0 11.6328597 1   1 (55,65]      1 0 0.04928111  20.29175
## 121 121 0 12.5638565 1  -2 (45,55]      1 0 0.04349269  22.99237
## 145 145 0 25.0000000 1  -2 (45,55]      1 0 0.08385424  11.92546
## 151 151 0 23.4817749 1  -1 (45,55]      0 0 0.09013661  11.09427
## 176 176 1 16.2877312 1   3 (55,65]      0 1 1.00000000   1.00000
## 187 187 0 25.0000000 1  -3 (45,55]      0 0 0.09389862  10.64978
summary(ncc_2_nodup_drop3$km_weight)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1.00    1.00   10.64   12.25   13.32 1205.79
m_cox_ncc_2_drop3 <- coxph(Surv(t, y) ~ x * z + gender + age, data = ncc_2_nodup_drop3,
                           weights = km_weight, robust = TRUE)
m_clogit_ncc_2_drop3 <- clogit(Fail ~ x * z + strata(Set), data = ncc_2_drop3)
```

## Compare results

``` r
results_2 <- rbind(summary(m_cox_cohort_2)$coef, 
                   summary(m_clogit_ncc_2)$coef, 
                   summary(m_cox_ncc_2)$coef[, -3], 
                   summary(m_clogit_ncc_2_drop3)$coef,
                   summary(m_cox_ncc_2_drop3)$coef[, -3])
results_2 <- data.frame(Variable = rownames(results_2), results_2, 
                        check.names = FALSE)
rownames(results_2) <- NULL
kable(data.frame(
  Data = c("Full cohort", rep("", 4), 
           "1:5 NCC (clogit)", rep("", 2), "1:5 NCC (weighted Cox)", rep("", 4),
           "1:5 NCC keep 2 controls (clogit)", rep("", 2),
           "1:5 NCC keep 2 controls (weighted Cox)", rep("", 4)), 
  Variable = results_2$Variable, 
  `True HR` = c(c(1.5, 4, 1.01, 1.01, 2), 
                c(1.5, 4, 2), c(1.5, 4, 1.01, 1.01, 2), 
                c(1.5, 4, 2), c(1.5, 4, 1.01, 1.01, 2)),
  `Estimated HR` = results_2[, "exp(coef)"], 
  `SE of log(HR)` = results_2[, "se(coef)"], 
  `p-value` = results_2[, "Pr(>|z|)"], check.names = FALSE
), digits = c(0, 0, 2, 2, 3, 3))
```

| Data                                   | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
| :------------------------------------- | :------- | ------: | -----------: | ------------: | ------: |
| Full cohort                            | x        |    1.50 |         1.47 |         0.110 |   0.001 |
|                                        | z        |    4.00 |         4.46 |         0.126 |   0.000 |
|                                        | age      |    1.01 |         1.01 |         0.002 |   0.000 |
|                                        | gender   |    1.01 |         0.92 |         0.038 |   0.024 |
|                                        | x:z      |    2.00 |         1.90 |         0.135 |   0.000 |
| 1:5 NCC (clogit)                       | x        |    1.50 |         1.51 |         0.114 |   0.000 |
|                                        | z        |    4.00 |         4.54 |         0.135 |   0.000 |
|                                        | x:z      |    2.00 |         1.92 |         0.145 |   0.000 |
| 1:5 NCC (weighted Cox)                 | x        |    1.50 |         1.50 |         0.112 |   0.000 |
|                                        | z        |    4.00 |         4.41 |         0.132 |   0.000 |
|                                        | gender   |    1.01 |         0.90 |         0.045 |   0.017 |
|                                        | age      |    1.01 |         1.01 |         0.002 |   0.001 |
|                                        | x:z      |    2.00 |         2.00 |         0.142 |   0.000 |
| 1:5 NCC keep 2 controls (clogit)       | x        |    1.50 |         1.50 |         0.120 |   0.001 |
|                                        | z        |    4.00 |         4.27 |         0.148 |   0.000 |
|                                        | x:z      |    2.00 |         1.90 |         0.160 |   0.000 |
| 1:5 NCC keep 2 controls (weighted Cox) | x        |    1.50 |         1.52 |         0.116 |   0.000 |
|                                        | z        |    4.00 |         4.47 |         0.141 |   0.000 |
|                                        | gender   |    1.01 |         0.95 |         0.055 |   0.313 |
|                                        | age      |    1.01 |         1.01 |         0.003 |   0.007 |
|                                        | x:z      |    2.00 |         1.98 |         0.153 |   0.000 |

**Notes:** 1:5 NCC indicates the weighted analysis of the NCC without
dropping any controls. 1:5 NCC keep 2 controls indicates the weighted
analysis of the NCC data where only the cases and 2 controls from each
set are used. Each NCC data was first analysed using conditional
logistic regression, and then weighted Cox model after breaking the
matching.
