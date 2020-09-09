library(testthat)
library(survival)
library(dplyr)
cohort <- data.frame(id = 1:10, 
                     t1 = c(0, 1, 0, 1, 3, 5, 2, 2, 6, 2), 
                     t2 = c(2, 10, 4, 5, 5, 7, 4, 8, 9, 8), 
                     status = c(1, 0, 0, 1, 0, 1, 1, 0, 1, 1), 
                     a = c(1, 2, 1, 2, 1, 1, 1, 2, 2, 2), 
                     b = c("b", "a", "b", "b", "a", "a", "a", "b", "b", "a"))
# Let everyone start at t=0 and end at t2:
# plot(x = c(0, 10), y = c(1, 10), type = "n", xlab = "time", ylab = "id")
# apply(cohort, 1, function(row) {
#   lines(x = c(0, row[3]), y = rep(row[1], 2), 
#         col = ifelse(row[5] == 1, "black", "red"), 
#         lty = ifelse(row[6] == "a", 1, 2))
#   if (row[4] == 1) points(x = row[3], y = row[1], pch = "x")
# })

test_that("prepare_cohort", {
  cohort2 <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = "a"
  )
  expect_equal(cohort2$id, cohort$id)
  expect_equal(cohort2$.y, cohort$status)
  expect_equal(cohort2$.t, cohort$t2)
  expect_equal(cohort2$.strata0, factor(cohort$a))
  expect_equal(cohort2$.strata, paste0("match_var=", cohort$a))
})

test_that("compute_risk_table: no matching", {
  tb_manual <- data.frame(t = c(2, 4, 5, 7, 8, 9), 
                          n = c(9, 8, 6, 4, 3, 1))
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = NULL
  ) %>% 
    SamplingDesignTools:::compute_risk_tb(cohort = ., match_var_names = NULL, 
                                          staggered = FALSE)
  expect_equal(tb_manual$n + 1, tb$n_at_risk)
})
test_that("compute_risk_table: matching on a", {
  tb_manual <- data.frame(t = c(2, 4, 5, 7, 8, 9), 
                          a = c(1, 1, 2, 1, 2, 2), 
                          n = c(4, 3, 4, 0, 3, 1))
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = "a"
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = "a", 
                                              staggered = FALSE) %>% 
    arrange(t_event)
  expect_equal(tb_manual$n + 1, tb$n_at_risk)
  expect_equal(tb_manual$a, tb$a)
})
test_that("compute_risk_table: matching on a&b", {
  tb_manual <- data.frame(t = c(2, 4, 5, 7, 8, 9), 
                          a = c(1, 1, 2, 1, 2, 2), 
                          b = c("b", "a", "b", "a", "a", "b"),
                          n = c(1, 2, 2, 0, 1, 0))
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = FALSE) %>% 
    arrange(t_event)
  expect_equal(tb_manual$n + 1, tb$n_at_risk)
  expect_equal(tb_manual$a, tb$a)
  expect_equal(tb_manual$b, tb$b)
})

test_that("compute_p: normal", {
  p <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 6, 
                                           n_per_case = 2, n_kept = 2)
  expect_identical(p, 1 - 2 / 5)
})
test_that("compute_p: not enough non-cases", {
  p <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 3, 
                                           n_per_case = 5, n_kept = 5)
  expect_identical(p, 0)
})
test_that("compute_p: normal, drop control", {
  p <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 6, 
                                           n_per_case = 2, n_kept = 1)
  p_selected <- 2 / 5
  p_kept <- 1 / 2
  expect_identical(p, (1 - p_selected) + p_selected * (1 - p_kept))
})
test_that("compute_p: normal, keep too many", {
  expect_warning({
    p <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 6, 
                                             n_per_case = 2, n_kept = 3)
  })
  expect_identical(p, 1 - 2 / 5)
})
test_that("compute_p: not enough non-cases, but drop some", {
  p <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 6, 
                                           n_per_case = 5, n_kept = 2)
  p2 <- SamplingDesignTools:::p_not_sampled(n_event = 1, n_at_risk = 6, 
                                            n_per_case = 2, n_kept = 2)
  expect_identical(p, p2)
})

test_that("km_weight: no matching", {
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = NULL
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., match_var_names = NULL, 
                                              staggered = FALSE) %>% 
    mutate(p = SamplingDesignTools:::p_not_sampled(n_event = n_event, 
                                                   n_at_risk = n_at_risk, 
                                                   n_per_case = 1, 
                                                   n_kept = 1))
  tb_kmw <- SamplingDesignTools:::compute_kmw0(risk_table = tb)
  tb_kmw2 <- data.frame(
    t_event = c(2, 4, 5, 7, 8, 9), 
    kmw = c(NA, 9 / 2, 1 / (1 - 35 / 54), NA, 1 / (1 - 35 / 108), 1)
  )
  tb <- left_join(tb_kmw, tb_kmw2) %>% filter(complete.cases(.))
  expect_equal(tb$km_weight, tb$kmw)
})
test_that("km_weight: matched on a&b", {
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = NULL, t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = FALSE) %>% 
    mutate(p = SamplingDesignTools:::p_not_sampled(n_event = n_event, 
                                                   n_at_risk = n_at_risk, 
                                                   n_per_case = 1, 
                                                   n_kept = 1))
  tb_kmw <- SamplingDesignTools:::compute_kmw0(risk_table = tb)
  tb_kmw2 <- data.frame(t_event = c(2, 4, 5, 8), kmw = c(1, 2, 2, 1))
  tb <- left_join(tb_kmw, tb_kmw2)
  expect_equal(tb$km_weight, tb$kmw)
})

test_that("compute_kmw_cohort: no matching", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  sample_stat <- numeric(nrow(cohort))
  sample_stat[unique(ncc$id[ncc$status == 0])] <- 1
  sample_stat[ncc$id[ncc$status == 1]] <- 2
  ncc_nodup <- SamplingDesignTools:::compute_kmw_cohort(
    cohort = cohort, t_name = "t2", sample_stat = sample_stat, 
    n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(unique(ncc$id)))
  expect_equal(names(ncc_nodup), c(names(cohort), c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, 
               # 1  2  3      4  6  7  8                   9  10
               c(1, 1, 9 / 2, 1, 1, 1, 1 / (1 - 35 / 108), 1, 1))
})
test_that("compute_kmw_cohort: matched on a&b", {
  ncc <- data.frame(set_id = rep(1:4, each = 2), 
                    id = c(1, 3, 
                           7, 5, 
                           4, 8, 
                           10, 2), 
                    t = rep(c(2, 4, 5, 8), each = 2), 
                    y = rep(c(1, 0), 4))
  sample_stat <- numeric(nrow(cohort))
  sample_stat[unique(ncc$id[ncc$y == 0])] <- 1
  sample_stat[ncc$id[ncc$y == 1]] <- 2
  ncc_nodup <- SamplingDesignTools:::compute_kmw_cohort(
    cohort = cohort, t_name = "t2", sample_stat = sample_stat, 
    match_var_names = c("a", "b"), n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(ncc$id))
  expect_equal(names(ncc_nodup), c(names(cohort), c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, c(1, 1, 1, 1, 2, 1, 2, 1))
})

# Staggered entry:
# plot(x = c(0, 10), y = c(1, 10), type = "n", xlab = "time", ylab = "id")
# apply(cohort, 1, function(row) {
#   lines(x = c(row[2], row[3]), y = rep(row[1], 2), 
#         col = ifelse(row[5] == 1, "black", "red"), 
#         lty = ifelse(row[6] == "a", 1, 2))
#   if (row[4] == 1) points(x = row[3], y = row[1], pch = "x")
# })
test_that("prepare_cohort: staggered entry", {
  cohort2 <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = "a"
  )
  expect_equal(cohort2$id, cohort$id)
  expect_equal(cohort2$.y, cohort$status)
  expect_equal(cohort2$.t, cohort$t2)
  expect_equal(cohort2$.t_start, cohort$t1)
  expect_equal(cohort2$.strata0, factor(cohort$a))
  expect_equal(cohort2$.strata, paste0("match_var=", cohort$a))
})

test_that("compute_risk_table: no matching, staggered entry", {
  tb_manual <- data.frame(t = c(2, 4, 5, 7, 8, 9), 
                          n = c(3, 6, 4, 4, 3, 1))
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = NULL
  ) %>% 
    SamplingDesignTools:::compute_risk_tb(cohort = ., match_var_names = NULL, 
                                          staggered = TRUE)
  expect_equal(tb_manual$n + 1, tb$n_at_risk)
})
test_that("compute_risk_table: matching on a&b, staggered entry", {
  tb_manual <- data.frame(t = c(2, 4, 5, 7, 8, 9), 
                          a = c(1, 1, 2, 1, 2, 2), 
                          b = c("b", "a", "b", "a", "a", "b"),
                          n = c(1, 1, 1, 0, 1, 0))
  tb <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = TRUE) %>% 
    arrange(t_event)
  expect_equal(tb_manual$n + 1, tb$n_at_risk)
  expect_equal(tb_manual$a, tb$a)
  expect_equal(tb_manual$b, tb$b)
})

test_that("compute_kmw_cohort: no matching, staggered entry", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  sample_stat <- numeric(nrow(cohort))
  sample_stat[unique(ncc$id[ncc$status == 0])] <- 1
  sample_stat[ncc$id[ncc$status == 1]] <- 2
  ncc_nodup <- SamplingDesignTools:::compute_kmw_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", sample_stat = sample_stat, 
    n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(unique(ncc$id)))
  expect_equal(names(ncc_nodup), c(names(cohort), c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, 
               # 1  2  3      4  6  7  8        9  10
               c(1, 1, 9 / 4, 1, 1, 1, 16 / 11, 1, 1))
})
test_that("compute_kmw_cohort: matched on a&b, staggered entry", {
  ncc <- data.frame(set_id = rep(1:4, each = 2), 
                    id = c(1, 3, 
                           7, 5, 
                           4, 8, 
                           10, 2), 
                    t = rep(c(2, 4, 5, 8), each = 2), 
                    y = rep(c(1, 0), 4))
  sample_stat <- numeric(nrow(cohort))
  sample_stat[unique(ncc$id[ncc$y == 0])] <- 1
  sample_stat[ncc$id[ncc$y == 1]] <- 2
  ncc_nodup <- SamplingDesignTools:::compute_kmw_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", sample_stat = sample_stat, 
    match_var_names = c("a", "b"), n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(ncc$id))
  expect_equal(names(ncc_nodup), c(names(cohort), c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, c(1, 1, 1, 1, 1, 1, 1, 1))
})

test_that("prepare_ncc_cases: no matching", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  ncc_cases <- ncc %>% filter(fail == 1) %>% 
    SamplingDesignTools:::prepare_ncc_cases(ncc_cases = ., t_name = "t", 
                                            match_var_names = NULL)
  expect_identical(ncc_cases$t, ncc_cases$.t_event)
})

test_that("match_risk_table: no matching", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    s = 1,
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  risk_table_manual <- data.frame(t_event = c(2, 4, 5, 7, 8, 9), 
                                  n_at_risk = c(3, 6, 4, 4, 3, 1) + 1, 
                                  s = 1)
  risk_table <- ncc %>% filter(fail == 1) %>% 
    SamplingDesignTools::match_risk_table(
      ncc_cases = ., risk_table_manual = risk_table_manual, t_coarse_name = "t", 
      t_name = "t", match_var_names = "s"
    )
  expect_true("s" %in% names(risk_table))
  expect_equal(risk_table$n_at_risk, risk_table_manual$n_at_risk)
})
test_that("match_risk_table: matched on a&b, staggered entry", {
  ncc <- data.frame(set_id = rep(1:3, each = 2), 
                    id = c(1, 3, 
                           7, 5, 
                           10, 2), 
                    t = rep(c(2, 4, 8), each = 2), 
                    y = rep(c(1, 0), 3)) %>% 
    left_join(cohort) %>% 
    arrange(set_id, status)
  risk_table_manual <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = TRUE) %>% 
    filter(t_event %in% ncc$t) %>% 
    select(t_event, n_at_risk, a, b)
  risk_table <- ncc %>% filter(y == 1) %>% 
    SamplingDesignTools::match_risk_table(
      ncc_cases = ., risk_table_manual = risk_table_manual, t_coarse_name = "t", 
      t_name = "t", match_var_names = c("a", "b")
    )
  expect_equal(risk_table$n_at_risk, risk_table_manual$n_at_risk)
})
test_that("attach km_weight to ncc controls: no matching", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  risk_table_manual <- data.frame(t_event = c(2, 4, 5, 7, 8, 9), 
                                  n_at_risk = c(3, 6, 4, 4, 3, 1) + 1)
  risk_table <- ncc %>% filter(fail == 1) %>% 
    SamplingDesignTools::match_risk_table(
      ncc_cases = ., risk_table_manual = risk_table_manual, t_coarse_name = "t", 
      t_name = "t", match_var_names = NULL
    ) %>% 
    mutate(p = SamplingDesignTools:::p_not_sampled(
      n_event = n_event, n_at_risk = n_at_risk, n_per_case = 1, n_kept = 1
    ))
  ncc_controls <- cohort %>% 
    filter(id %in% setdiff(ncc$id[ncc$fail == 0], ncc$id[ncc$fail == 1])) %>% 
    SamplingDesignTools:::prepare_cohort(
      cohort = ., t_start_name = "t1", t_name = "t2", 
      y_name = "status", match_var_names = NULL
    )
  ncc_controls <- SamplingDesignTools:::assign_kmw0(ncc_nodup = ncc_controls, 
                                                    risk_table = risk_table)
  expect_equal(ncc_controls$.km_weight, 
               # 2  3      8      
               c(1, 9 / 4, 16 / 11))
})

test_that("compute km_weight given ncc: no matching, staggered entry", {
  ncc <- data.frame(set = rep(1:6, each = 2), 
                    id = c(1, 3, 7, 10, 4, 8, 6, 10, 10, 2, 9, 2), 
                    t = rep(c(2, 4, 5, 7, 8, 9), each = 2),
                    fail = rep(c(1, 0), 6)) %>% 
    left_join(cohort) %>% 
    arrange(set, status)
  risk_table_manual <- data.frame(t_event = c(2, 4, 5, 7, 8, 9), 
                                  n_at_risk = c(3, 6, 4, 4, 3, 1) + 1)
  ncc_nodup <- SamplingDesignTools:::compute_kmw_ncc(
    ncc = ncc[, -1], risk_table_manual = risk_table_manual, t_match_name = "t", 
    y_name = "fail", n_per_case = 1, id_name = "id", 
    t_start_name = "t1", t_name = "t2"
  )
  expect_equal(sort(ncc_nodup$id), sort(unique(ncc$id)))
  expect_equal(names(ncc_nodup), c(names(ncc)[-c(1, 3)], c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, 
               # 1  2  3      4  6  7  8        9  10
               c(1, 1, 9 / 4, 1, 1, 1, 16 / 11, 1, 1))
  expect_equal(ncc_nodup$fail, ncc_nodup$status)
})
test_that("compute km_weight given ncc: matched on a&b, staggered entry", {
  ncc <- data.frame(set_id = rep(1:4, each = 2), 
                    id = c(1, 3, 
                           7, 5, 
                           4, 8, 
                           10, 2), 
                    t = rep(c(2, 4, 5, 8), each = 2), 
                    y = rep(c(1, 0), 4)) %>% 
    left_join(cohort) %>% 
    arrange(set_id, status)
  risk_table_manual <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = TRUE) %>% 
    filter(t_event %in% ncc$t) %>% 
    select(t_event, n_at_risk, a, b)
  ncc_nodup <- SamplingDesignTools:::compute_kmw_ncc(
    ncc = ncc[, -1], risk_table_manual = risk_table_manual, id_name = "id",
    t_start_name = "t1", t_name = "t2", t_match_name = "t", y_name = "y",
    match_var_names = c("a", "b"), n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(ncc$id))
  expect_equal(names(ncc_nodup), c(names(ncc)[-c(1, 3)], c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, c(1, 1, 1, 1, 1, 1, 1, 1))
})
test_that("compute km_weight given ncc: matched on a&b, staggered entry, cases only", {
  ncc <- data.frame(set_id = rep(1:4, each = 2), 
                    id = c(1, 3, 
                           7, 5, 
                           4, 8, 
                           10, 2), 
                    t = rep(c(2, 4, 5, 8), each = 2), 
                    y = rep(c(1, 0), 4)) %>% 
    left_join(cohort) %>% 
    arrange(set_id, status)
  risk_table_manual <- SamplingDesignTools:::prepare_cohort(
    cohort = cohort, t_start_name = "t1", t_name = "t2", y_name = "status", 
    match_var_names = c("a", "b")
  ) %>% SamplingDesignTools:::compute_risk_tb(cohort = ., 
                                              match_var_names = c("a", "b"), 
                                              staggered = TRUE) %>% 
    filter(t_event %in% ncc$t) %>% 
    select(t_event, n_at_risk, a, b)
  ncc_nodup <- SamplingDesignTools:::compute_kmw_ncc(
    ncc = ncc[ncc$y == 1, -1], risk_table_manual = risk_table_manual, id_name = "id",
    t_start_name = "t1", t_name = "t2", t_match_name = "t", y_name = "y",
    match_var_names = c("a", "b"), n_per_case = 1
  )
  expect_equal(sort(ncc_nodup$id), sort(ncc$id[ncc$y == 1]))
  expect_equal(names(ncc_nodup), c(names(ncc)[-c(1, 3)], c(".km_prob", ".km_weight")))
  expect_equal(ncc_nodup$.km_weight, c(1, 1, 1, 1))
})

data("cohort_2")
data("ncc_2")

test_that("prepare_risk_table", {
  risk_table <- SamplingDesignTools::compute_risk_table(
    cohort = cohort_2, t_name = "t", y_name = "y",
    match_var_names = c("age_cat", "gender")
  ) %>% arrange(t_event)
  risk_table_temp <- SamplingDesignTools::prepare_risk_table(
    ncc = ncc_2, t_match_name = "Time", y_name = "Fail", 
    match_var_names = c("age_cat", "gender")
  )
  expect_equal(risk_table$t_event, risk_table_temp$t_event)
  expect_equal(risk_table$age_cat, risk_table_temp$age_cat)
})

sample_stat <- numeric(nrow(cohort_2))
sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
output <- SamplingDesignTools::compute_km_weights(
  cohort = cohort_2, t_name = "t", sample_stat = sample_stat, 
  match_var_names = c("age_cat", "gender"), n_per_case = 5, 
  return_risk_table = TRUE
)
risk_table <- SamplingDesignTools::compute_risk_table(
  cohort = cohort_2, t_name = "t", y_name = "y",
  match_var_names = c("age_cat", "gender")
)
test_that("compute_km_weights: cohort, check risk_table", {
  expect_equal(output$risk_table$t_event, risk_table$t_event)
  expect_equal(output$risk_table$n_at_risk, risk_table$n_at_risk)
})
output2 <- ncc_2 %>% select(-Set) %>% 
  SamplingDesignTools::compute_km_weights(
    ncc = ., risk_table_manual = risk_table, 
    id_name = "Map", t_match_name = "Time", t_name = "t", y_name = "Fail", 
    match_var_names = c("age_cat", "gender"), n_per_case = 5, 
    return_risk_table = TRUE
  ) 
test_that("compute_km_weights: ncc, check risk_table", {
  expect_equal(output2$risk_table$t_event, risk_table$t_event)
  expect_equal(output2$risk_table$n_at_risk, risk_table$n_at_risk)
})
km_table_ncc <- ncc_2 %>% filter(Fail == 1) %>% select(-Set) %>% 
  SamplingDesignTools::match_risk_table(
    ncc_cases = ., risk_table_manual = risk_table, 
    t_coarse_name = "Time", t_name = "Time",
    match_var_names = c("age_cat", "gender")
  ) %>% 
  mutate(p = SamplingDesignTools:::p_not_sampled(
    n_event = n_event, n_at_risk = n_at_risk, n_per_case = 5, n_kept = 5
  )) %>% 
  SamplingDesignTools:::compute_kmw0(risk_table = .)
test_that("match_risk_table", {
  expect_equal(km_table_ncc$t_event, risk_table$t_event)
})

output4 <- ncc_2 %>% filter(Fail == 1) %>% select(-Set) %>% 
  SamplingDesignTools::compute_km_weights(
    ncc = ., risk_table_manual = risk_table, 
    id_name = "Map", t_name = "t", y_name = "Fail", 
    match_var_names = c("age_cat", "gender"), n_per_case = 5, 
    return_risk_table = TRUE
  ) 
test_that("compute_km_weights: ncc cases, check risk_table", {
  expect_equal(output4$risk_table$t_event, risk_table$t_event)
  expect_equal(output4$risk_table$n_at_risk, risk_table$n_at_risk)
})

test_that("km_weights", {
  expect_equal(output$dat$km_weight, output2$dat$km_weight)
  expect_true(all(output$dat$km_weight[output$dat$y == 0] %in% 
                    km_table_ncc$km_weight))
})

# ncc_nodup <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y",
#                                 match_var_names = c("age_cat", "gender"),
#                                 sample_stat = sample_stat, n_per_case = 5)
# risk_table <- compute_risk_table(cohort = cohort_2, t_name = "t", y_name = "y",
#                                  match_var_names = c("age_cat", "gender"))
# ncc_nodup_v2 <- compute_km_weights(ncc = ncc_2[, -1], risk_table_manual = risk_table,
#                                    id_name = "Map", t_match_name = "Time",
#                                    t_name = "t", y_name = "Fail",
#                                    match_var_names = c("age_cat", "gender"),
#                                    n_per_case = 5)
# all.equal(ncc_nodup$km_weight, ncc_nodup_v2$km_weight)
