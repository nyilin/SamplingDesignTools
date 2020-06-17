
#' Compute Kaplan-Meier type weights for counter-matched NCC
#' @inheritParams compute_km_weights
#' @param match_var_name Name of the categorical variable in \code{cohort} to
#'   count-match on, which can be the exposure or the surrogate for the
#'   exposure. A \code{string}. If a vector is supplied, only the first element
#'   will be used.
#' @param ml Number of subjects to draw from each strata defined by the match
#'   variable (including the case). 
#' @import survival
#' @import dplyr
#' @export
compute_km_weights_cm <- function(cohort, t_name = NULL, y_name = NULL,
                                  sample_stat, match_var_name = NULL,
                                  ml) {
  cohort <- as.data.frame(cohort)
  if (length(sample_stat) != nrow(cohort)) {
    stop(simpleError("Length of sample_stat must be same as number of rows in cohort."))
  }
  if (is.null(y_name)) {
    stop(simpleError("Please sapply name of event status."))
  } else {
    y_name <- y_name[1]
    if (!(y_name %in% names(cohort))) {
      stop(simpleError(paste(y_name, "not found in cohort.")))
    }
  }
  if (is.null(t_name)) {
    stop(simpleError("Please sapply name of event/censoring time."))
  } else {
    t_name <- t_name[1]
    if (!(t_name %in% names(cohort))) {
      stop(simpleError(paste(t_name, "not found in cohort.")))
    }
  }
  if (is.null(match_var_name)) {
    stop(simpleError("Please sapply name of variable to counter-match on."))
  } else {
    match_var_name <- match_var_name[1]
    if (!(match_var_name %in% names(cohort))) {
      stop(simpleError(paste(match_var_name, "not found in cohort.")))
    }
  }
  t <- cohort[, t_name]
  y <- as.numeric(sample_stat >= 2)
  match_var <- cohort[, match_var_name]
  # Make sure the levels are integers starting from 1
  match_var <- factor(as.numeric(factor(match_var)))
  km <- survfit(Surv(t, y) ~ match_var)
  km_summ <- summary(km)
  df <- data.frame(t = t[which(y == 1)], strata_case = match_var[which(y == 1)])
  km_tb <- data.frame(t = km_summ$time, Rj = km_summ$n.risk - km_summ$n.event,
                      strata = as.character(km_summ$strata), 
                      stringsAsFactors = FALSE)
  km_tb <- merge(df, km_tb, by = "t", all = TRUE) %>%
    arrange(strata, t)
  km_tb <- km_tb %>%
    group_by(strata) %>% 
    mutate(prob_not_sampled = 1 - ((ml - as.numeric(strata == strata_case)) / Rj),
           cumulative_product = cumprod(prob_not_sampled),
           sampling_prob = 1 - cumulative_product)
  in_ncc <- sample_stat > 0
  ncc_nodup <- cohort[in_ncc, ]
  match_var_ncc <- paste0("match_var=", match_var[in_ncc])
  t_ncc <- t[in_ncc]
  p_ncc <- unlist(lapply(1:nrow(ncc_nodup), function(j) {
    if (ncc_nodup[j, y_name] == 1) {
      1
    } else {
      km_tb_i <- km_tb[km_tb$t < t_ncc[j] & km_tb$strata == match_var_ncc[j], ]
      r <- nrow(km_tb_i)
      if (r == 0) {
        0
      } else {
        km_tb_i$sampling_prob[r]
      }
    }
  }))
  cbind(ncc_nodup, km_prob = p_ncc, km_weight = 1 / p_ncc)
}