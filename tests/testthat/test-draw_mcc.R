library(dplyr)
library(testthat)
cohort <- data.frame(id = 1:40, 
                     y = c(rep(1, 10), rep(0, 30)), 
                     a = rep(letters[1:2], 20), 
                     b = rep(c(1, 1, 2, 2), 10))
count_cohort <- cohort %>% 
  group_by(a, b) %>% 
  summarise(n_cases = sum(y == 1), n_coltrols = sum(y == 0))
table(cohort$a, cohort$b, cohort$y)
test_that("draw scc (private)", {
  dat_scc <- SamplingDesignTools:::draw_scc0(cohort = cohort, y_name = "y", 
                                             n_per_case = 2)
  expect_true(!any(duplicated(dat_scc$id)))
  expect_true(sum(dat_scc$y == 1) * 2 == sum(dat_scc$y == 0))
  expect_equal(unique(dat_scc$.w[dat_scc$y == 1]), 1)
  expect_equal(unique(dat_scc$.w[dat_scc$y == 0]), 
               sum(cohort$y == 0) / (sum(cohort$y == 1) * 2))
})
test_that("draw mcc (private)", {
  dat_mcc <- SamplingDesignTools:::draw_mcc0(cohort = cohort, y_name = "y", 
                                             n_per_case = 2, 
                                             match_var_names = "a")
  expect_true(!any(duplicated(dat_mcc$id)))
  count <- matrix(table(dat_mcc$y, dat_mcc$a), nrow = 2)
  # 1st row for y=0, 2nd row for y=1
  expect_equivalent(count[1, ], count[2, ] * 2)
  expect_equal(unique(dat_mcc$.w[dat_mcc$y == 1]), 1)
  expect_equal(unique(dat_mcc$.w[dat_mcc$y == 0 & dat_mcc$a == "a"]), 
               sum(cohort$y == 0 & cohort$a == "a") / 
                 (sum(cohort$y == 1 & cohort$a == "a") * 2))
  expect_equal(unique(dat_mcc$.w[dat_mcc$y == 0 & dat_mcc$a == "b"]), 
               sum(cohort$y == 0 & cohort$a == "b") / 
                 (sum(cohort$y == 1 & cohort$a == "b") * 2))
})
test_that("draw scc", {
  dat_scc <- SamplingDesignTools::draw_mcc(cohort = cohort, y_name = "y", 
                                           n_per_case = 2, 
                                           weight_name = "weight")
  expect_true(!any(duplicated(dat_scc$id)))
  expect_equivalent(sum(dat_scc$y == 1), sum(cohort$y == 1))
  expect_true(sum(dat_scc$y == 1) * 2 == sum(dat_scc$y == 0))
  expect_equivalent(names(dat_scc), c(names(cohort), "weight"))
})
test_that("draw mcc", {
  dat_mcc <- SamplingDesignTools::draw_mcc(cohort = cohort, y_name = "y", 
                                           n_per_case = 2, 
                                           match_var_names = c("a", "b"), 
                                           weight_name = "weight")
  expect_true(!any(duplicated(dat_mcc$id)))
  count <- dat_mcc %>% 
    group_by(a, b) %>% 
    summarise(n_cases = sum(y == 1), n_coltrols = sum(y == 0))
  expect_equivalent(count$n_cases, count_cohort$n_cases)
  expect_equivalent(count$n_coltrols, count$n_cases * 2)
  expect_equivalent(names(dat_mcc), c(names(cohort), "weight"))
})
test_that("draw mcc part of cases", {
  dat_mcc <- SamplingDesignTools::draw_mcc(cohort = cohort, y_name = "y", 
                                           n_per_case = 2, n_cases = 5,
                                           match_var_names = c("a", "b"), 
                                           weight_name = "weight")
  expect_true(!any(duplicated(dat_mcc$id)))
  count <- dat_mcc %>% 
    group_by(a, b) %>% 
    summarise(n_cases = sum(y == 1), n_coltrols = sum(y == 0))
  expect_equivalent(sum(dat_mcc$y == 1), 5)
  expect_equivalent(count$n_coltrols, count$n_cases * 2)
  expect_equivalent(names(dat_mcc), c(names(cohort), "weight"))
  expect_equivalent(unique(dat_mcc$weight[dat_mcc$y == 1]), sum(cohort$y == 1) / 5)
})
