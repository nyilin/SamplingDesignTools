0 Preparation
=============

Install the SamplingDesignTools package from GitHub (package devtools
needed):

    # devtools::install_github("nyilin/SamplingDesignTools")

Load packages used:

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

1 Drawing counter-matched sample \[\#S1\]
=========================================

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

In this example, I am using a simulated cohort with 10,000 subjects,
where a continuous exposure (x \~ N(55, 10^2)) represents their age and
a binary covariate indicates their gender (c = 1 for male, Pr(c=1)=0.5).

The event time t was generated from a Weibull distribution with hazard
function:

log{h(t; x, c)} = log{h\_0} + log(1.1) \* x + log(2) \* c,

where the baseline hazard h\_0=10^{-5} and the scale parameter was 1.
All individuals entered the cohort at time t=0 and were followed until
the event of interest or a censoring time (generated from an exponential
distribution with a rate of 0.05) with a maximum follow-up time of t=25.
There were 582 cases (5.8%) in this cohort.

Counter-match on binary surrogate
---------------------------------

In this analysis, I dichotomised the continuous exposure at 50 years of
age, and drew a 1:1 time-matched sample counter matched on this binary
variable that is considered as a surrogate for the exposure. This sample
was analysed using a weighted conditional logistic model, where log of
sampling weights were specified as offset in the model.

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

Counter-match on categorical surrogate
--------------------------------------

In this analysis, I categorised the continuous exposure into 4
categories at 40, 50 and 60 years of age, and drew a 1:3 time-matched
sample counter matched on this categorical variable that is considered
as a surrogate for the exposure. This sample was again analysed using a
weighted conditional logistic model.

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

Compare results
---------------

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

| Data        | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
|:------------|:---------|--------:|-------------:|--------------:|--------:|
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

2 Compute KM-type weights for NCC sample
========================================

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

In this example, I use a simulated cohort with 100,000 subjects with
age\~N(55, 10^2), Pr(male=1)=0.5, a binary exposure (x) that depends on
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

1:5 NCC matched on age and gender
---------------------------------

Subjects were divided into 6 age groups: &lt;35, 36-45, 46-55, 56-65,
66-75 and &gt;75. A 1:5 NCC sample was draw from the cohort, matched on
age group and gender.

    n_per_case <- 5
    ncc_2 <- ccwc(exit = t, fail = y, controls = n_per_case, 
                  match = list(age_cat, gender), include = list(x, age, z), 
                  data = cohort_2, silent = TRUE)
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
    ## Start time is 0 for all subjects. Event/censoring time is given by variable t.
    ## Joining, by = c("age_cat", "gender")
    ## Joining, by = "strata"
    head(ncc_2_nodup, 12)
    ##    id y          t x age age_cat gender z    km_prob km_weight
    ## 1   1 0 25.0000000 1  -2 (45,55]      0 0 0.21851027  4.576444
    ## 2  16 0 25.0000000 1  10 (65,75]      0 1 0.24202210  4.131854
    ## 3  35 0 16.5202891 1  -5 (45,55]      1 1 0.14312411  6.986943
    ## 4  36 1  0.2016047 1 -16 (35,45]      1 0 1.00000000  1.000000
    ## 5  54 0 25.0000000 0   3 (55,65]      1 1 0.21862668  4.574007
    ## 6  61 1 14.8098399 1   7 (55,65]      0 0 1.00000000  1.000000
    ## 7  69 0 25.0000000 1  -1 (45,55]      1 0 0.19666881  5.084690
    ## 8  72 1  4.1087132 1  -1 (45,55]      0 1 1.00000000  1.000000
    ## 9  76 1 20.1285244 1  15 (65,75]      0 1 1.00000000  1.000000
    ## 10 82 0  5.7436608 0   0 (45,55]      0 1 0.06188087 16.160084
    ## 11 84 1  5.9032539 1  -9 (45,55]      1 1 1.00000000  1.000000
    ## 12 86 1  6.0871929 1  -3 (45,55]      1 0 1.00000000  1.000000
    summary(ncc_2_nodup$km_weight)
    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.000   4.218   4.805   6.547   6.436 482.572
    m_cox_ncc_2 <- coxph(Surv(t, y) ~ x * z + gender + age, data = ncc_2_nodup,
                         weights = km_weight, robust = TRUE)
    m_clogit_ncc_2 <- clogit(Fail ~ x * z + strata(Set), data = ncc_2)

Compare results
---------------

    results_2 <- rbind(summary(m_cox_cohort_2)$coef, 
                       summary(m_clogit_ncc_2)$coef, 
                       summary(m_cox_ncc_2)$coef[, -3])
    results_2 <- data.frame(Variable = rownames(results_2), results_2, 
                            check.names = FALSE)
    rownames(results_2) <- NULL
    kable(data.frame(
      Data = c("Full cohort", rep("", 4), 
               "1:5 NCC (clogit)", rep("", 2), 
               "1:5 NCC (weighted Cox)", rep("", 4)), 
      Variable = results_2$Variable, 
      `True HR` = c(c(1.5, 4, 1.01, 1.01, 2), 
                    c(1.5, 4, 2), c(1.5, 4, 1.01, 1.01, 2)),
      `Estimated HR` = results_2[, "exp(coef)"], 
      `SE of log(HR)` = results_2[, "se(coef)"], 
      `p-value` = results_2[, "Pr(>|z|)"], check.names = FALSE
    ), digits = c(0, 0, 2, 2, 3, 3))

| Data                   | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
|:-----------------------|:---------|--------:|-------------:|--------------:|--------:|
| Full cohort            | x        |    1.50 |         1.47 |         0.110 |   0.001 |
|                        | z        |    4.00 |         4.46 |         0.126 |   0.000 |
|                        | age      |    1.01 |         1.01 |         0.002 |   0.000 |
|                        | gender   |    1.01 |         0.92 |         0.038 |   0.024 |
|                        | x:z      |    2.00 |         1.90 |         0.135 |   0.000 |
| 1:5 NCC (clogit)       | x        |    1.50 |         1.51 |         0.114 |   0.000 |
|                        | z        |    4.00 |         4.54 |         0.135 |   0.000 |
|                        | x:z      |    2.00 |         1.92 |         0.145 |   0.000 |
| 1:5 NCC (weighted Cox) | x        |    1.50 |         1.50 |         0.112 |   0.000 |
|                        | z        |    4.00 |         4.41 |         0.132 |   0.000 |
|                        | gender   |    1.01 |         0.90 |         0.045 |   0.017 |
|                        | age      |    1.01 |         1.01 |         0.002 |   0.001 |
|                        | x:z      |    2.00 |         2.00 |         0.142 |   0.000 |

**Notes:** 1:5 NCC indicates the weighted analysis of the NCC without
dropping any controls. 1:5 NCC keep 2 controls indicates the weighted
analysis of the NCC data where only the cases and 2 controls from each
set are used. Each NCC data was first analysed using conditional
logistic regression, and then weighted Cox model after breaking the
matching.

3 Compute weights for NCC without full cohort
=============================================

In the previous section, we assumed the full cohort was available when
computing the KM-type weight for each subject in the NCC. However, in
reality this may not always be the case. When the cohort is not
available, the KM-type weights can be computed for the NCC sample as
long the time of event/censoring for each subject is available, and the
number of subjects at risk can be obtained (or approximated) elsewhere.

To illustrate this scenario, I reuse the simulated cohort from the
previous section:

    head(cohort_2)
    ##   id y         t x age   age_cat gender z
    ## 1  1 0 25.000000 1  -2   (45,55]      0 0
    ## 2  2 0 19.819801 1  -4   (45,55]      1 0
    ## 3  3 0 25.000000 1  -5   (45,55]      0 0
    ## 4  4 0 12.414616 1  20 (75, Inf]      1 0
    ## 5  5 0 25.000000 1  -2   (45,55]      0 1
    ## 6  6 0  1.019023 0 -15   (35,45]      1 0

Note that `Time` is the event time of each case in a matched set, while
`t` is the actual event/censoring time of each subject. When the full
cohort is not available, it is often difficult to obtain the exact
number of subject at risk at each event time if the cohort is not
available. I simulate such situation by rounding the event/censoring
time (`t`) to the next integer, and use this coarsened time to compute
the number at risk:

    ncc_2$t <- cohort_2$t[ncc_2$Map]
    ncc_2$t_yr <- ceiling(ncc_2$t)
    cohort_2$t_yr <- ceiling(cohort_2$t)
    risk_table_coarse <- compute_risk_table(cohort = cohort_2, t_name = "t_yr", 
                                            y_name = "y", 
                                            match_var_names = c("age_cat", "gender"))
    ## Start time is 0 for all subjects. Event/censoring time is given by variable t_yr.
    ## Joining, by = c("age_cat", "gender")
    ## Joining, by = "strata"
    head(risk_table_coarse)
    ##   t_event n_event n_at_risk   age_cat gender
    ## 1       1       3      1165 (-Inf,35]      0
    ## 2       3       2      1045 (-Inf,35]      0
    ## 3       4       2       983 (-Inf,35]      0
    ## 4       5       1       928 (-Inf,35]      0
    ## 5       8       1       806 (-Inf,35]      0
    ## 6       9       3       772 (-Inf,35]      0

Now I use this information to compute the KM-type weights and
subsequently fit the weighted Cox model:

    ncc_nodup2 <- compute_km_weights(ncc = ncc_2[, -1], 
                                     risk_table_manual = risk_table_coarse, 
                                     t_name = "t", y_name = "Fail", 
                                     t_match_name = "t_yr",
                                     id_name = "Map", 
                                     match_var_names = c("age_cat", "gender"), 
                                     n_per_case = 5)
    ## Make sure input ncc does not include ID of matched sets.
    ## Joining, by = c("age_cat", "gender")
    ## Joining, by = c(".t_event", "age_cat", "gender", "n_event")
    ## Start time is 0 for all subjects. Event/censoring time is given by variable t.
    ## Joining, by = c("age_cat", "gender")
    ## Returned data contains 16359 rows for the 16359 unique subjects in the input ncc (identified by Map).
    m_cox_ncc_2_v2 <- coxph(Surv(t, Fail) ~ x * z + age + gender, 
                            data = ncc_nodup2, weights = km_weight, robust = TRUE)

Compare with results when the full cohort is available:

    results_3 <- rbind(summary(m_cox_cohort_2)$coef, 
                       summary(m_cox_ncc_2)$coef[, -3], 
                       summary(m_cox_ncc_2_v2)$coef[, -3])
    results_3 <- data.frame(Variable = rownames(results_3), results_3, 
                            check.names = FALSE)
    rownames(results_3) <- NULL
    kable(data.frame(
      Data = c("Full cohort", rep("", 4), 
               "1:5 NCC (weighted Cox)", rep("", 4),
               "1:5 NCC (weighted Cox, v2)*", rep("", 4)), 
      Variable = results_3$Variable, 
      `True HR` = rep(c(1.5, 4, 1.01, 1.01, 2), 3),
      `Estimated HR` = results_3[, "exp(coef)"], 
      `SE of log(HR)` = results_3[, "se(coef)"], 
      `p-value` = results_3[, "Pr(>|z|)"], check.names = FALSE
    ), digits = c(0, 0, 2, 2, 3, 3))

| Data                         | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
|:-----------------------------|:---------|--------:|-------------:|--------------:|--------:|
| Full cohort                  | x        |    1.50 |         1.47 |         0.110 |   0.001 |
|                              | z        |    4.00 |         4.46 |         0.126 |   0.000 |
|                              | age      |    1.01 |         1.01 |         0.002 |   0.000 |
|                              | gender   |    1.01 |         0.92 |         0.038 |   0.024 |
|                              | x:z      |    2.00 |         1.90 |         0.135 |   0.000 |
| 1:5 NCC (weighted Cox)       | x        |    1.50 |         1.50 |         0.112 |   0.000 |
|                              | z        |    4.00 |         4.41 |         0.132 |   0.000 |
|                              | gender   |    1.01 |         0.90 |         0.045 |   0.017 |
|                              | age      |    1.01 |         1.01 |         0.002 |   0.001 |
|                              | x:z      |    2.00 |         2.00 |         0.142 |   0.000 |
| 1:5 NCC (weighted Cox, v2)\* | x        |    1.50 |         1.50 |         0.110 |   0.000 |
|                              | z        |    4.00 |         4.15 |         0.127 |   0.000 |
|                              | age      |    1.01 |         1.00 |         0.002 |   0.823 |
|                              | gender   |    1.01 |         0.97 |         0.038 |   0.494 |
|                              | x:z      |    2.00 |         1.74 |         0.135 |   0.000 |

\*: Using the KM-type weights computed from the approximate number at
risk obtained in this section.

In reality, number at risk at each event time may be approximated by,
e.g., size of the relevant sub-population at mid-year. In such case,
user may use the following function to generate a template for
`risk_table_coarse` to fill in:

    risk_table_coarse <- prepare_risk_table(ncc = ncc_2, t_match_name = "t_yr", 
                                            y_name = "Fail", 
                                            match_var_names = c("gender", "age_cat"), 
                                            csv_file = NULL)
    head(risk_table_coarse)
    ##   t_event gender   age_cat n_at_risk
    ## 1       1      0 (-Inf,35]        NA
    ## 2       1      0   (35,45]        NA
    ## 3       1      0   (45,55]        NA
    ## 4       1      0   (55,65]        NA
    ## 5       1      0   (65,75]        NA
    ## 6       1      0 (75, Inf]        NA

This template will be written to a `csv` if specified by `csv_file`,
making it easier to supply information regarding the cohort that is
required for computing the KM-type weights.

4 Weighted analysis of (more) extreme case-control samples
==========================================================

In this section, I illustrate the weighted analysis of (more) extreme
case-control (MECC) samples using the conditional approach described in
Section 2.2 of Salim et al (2014), using the simulated cohort described
in [Section 1](#S1) of this document:

    # Load cohort data
    data(cohort_1)
    head(cohort_1)
    ##   id y        t age gender
    ## 1  1 0 25.00000  47      1
    ## 2  2 0 10.65152  58      0
    ## 3  3 0 25.00000  46      0
    ## 4  4 0 15.84131  52      0
    ## 5  5 0 22.57659  49      0
    ## 6  6 0 25.00000  63      0

For illustrative purpose, I draw a matched MECC sample, with 2 controls
frequency-matched to each case on gender. Unlike a NCC study where any
subject who had the event during the entire follow-up is sampled as a
case, in this MECC study I define cases as subjects who had the event
within *τ*<sub>0</sub> = 5 years of follow-up, and define eligible
controls as subjects who did not have the event in the first *τ* = 15
years of follow-up. Note that an event in a MECC study (`y_mecc`) has a
different definition as an event in the original cohort study.

    set.seed(1)
    dat_mecc <- draw_mecc(cohort = cohort_1, tau0 = 5, tau = 15,
                          id_name = "id", t_name = "t", delta_name = "y",
                          match_var_names = "gender", n_per_case = 2)
    ## Joining, by = "gender"
    ## Joining, by = c(".t", ".strata")
    ## Joining, by = ".strata"
    kable(head(dat_mecc))

|   id |   y |         t | age | gender | y\_mecc |      surv | surv\_tau | set\_id\_mecc |
|-----:|----:|----------:|----:|-------:|--------:|----------:|----------:|:--------------|
| 8756 |   1 | 2.1351357 |  76 |      0 |       1 | 0.9943501 |  0.957786 | set\_8756     |
| 8069 |   1 | 3.6774974 |  74 |      0 |       1 | 0.9898270 |  0.957786 | set\_8069     |
| 5939 |   1 | 0.8364486 |  49 |      0 |       1 | 0.9979733 |  0.957786 | set\_5939     |
| 1756 |   1 | 4.3993170 |  63 |      0 |       1 | 0.9878869 |  0.957786 | set\_1756     |
| 7955 |   1 | 4.4904526 |  64 |      0 |       1 | 0.9876400 |  0.957786 | set\_7955     |
| 1531 |   1 | 2.5664238 |  69 |      0 |       1 | 0.9918995 |  0.957786 | set\_1531     |

    table(y = dat_mecc$y, y_mecc = dat_mecc$y_mecc)
    ##    y_mecc
    ## y     0   1
    ##   0 397   0
    ##   1  11 204

As the conditional analysis of the MECC sample assumes individual
matching, the function `draw_mecc` randomly matches controls to each
case to form matched sets, indicated by the variable `set_id_mecc`. In
addition, this function provides the estimated baseline survival
probabilities (i.e., *Ŝ*(*t*<sub>*i*</sub>) and *Ŝ*(*τ*)) that are
needed to compute the weighted likelihood.

The weighted approach estimates the HR, where all covariates need to be
centred at the cohort average:

    dat_mecc$age_c <- dat_mecc$age - mean(cohort_1$age)
    # When there is only a single coefficient in, users should provide reasonable 
    # lower and upper limits for the estimate:
    m_mecc <- analyse_mecc_cond(
      y_name = "y_mecc", x_formula = ~ age_c, set_id_name = "set_id_mecc",
      surv = dat_mecc$surv, surv_tau = dat_mecc$surv_tau, mecc = dat_mecc,
      lower = -1, upper = 1 
    )
    ## Warning in optimize(function(par) fn(par, ...)/con$fnscale, lower = lower, : NA/
    ## Inf replaced by maximum positive value

    ## Warning in optimize(function(par) fn(par, ...)/con$fnscale, lower = lower, : NA/
    ## Inf replaced by maximum positive value
    round(m_mecc$coef_mat[, -1], 3)
    ##        est exp_est   se pval
    ## age_c 0.11   1.116 0.01    0
    coef_cox_cohort_1 <- summary(m_cox_cohort_1)$coef
    kable(data.frame(
      Data = c("Full cohort", "1:2 MECC, weighted"), 
      Variable = "Age", 
      `True HR` = 1.1,
      `Estimated HR` = c(coef_cox_cohort_1["age", "exp(coef)"], 
                         m_mecc$coef_mat["age_c", "exp_est"]), 
      `SE of log(HR)` = c(coef_cox_cohort_1["age", "se(coef)"], 
                         m_mecc$coef_mat["age_c", "se"]), 
      `p-value` = c(coef_cox_cohort_1["age", "Pr(>|z|)"], 
                         m_mecc$coef_mat["age_c", "pval"]), 
      check.names = FALSE
    ), digits = c(0, 0, 1, 2, 3, 3))

| Data               | Variable | True HR | Estimated HR | SE of log(HR) | p-value |
|:-------------------|:---------|--------:|-------------:|--------------:|--------:|
| Full cohort        | Age      |     1.1 |         1.11 |         0.004 |       0 |
| 1:2 MECC, weighted | Age      |     1.1 |         1.12 |         0.010 |       0 |

**References:**

-   Salim A, Ma X, Fall K, et al. Analysis of incidence and prognosis
    from ‘extreme’ case–control designs. Stat Med 2014; 33: 5388–5398.
-   Støer NC, Salim A, Bokenberger K, et al. Is the matched extreme
    case–control design more powerful than the nested case–control
    design?. Stat Methods Med Res 2019; 28(6): 1911-1923.
