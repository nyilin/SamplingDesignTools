#' Draw counter-matched NCC from cohort
#' @param cohort Cohort data. A \code{data.frame} or a matrix with column names.
#' @param y_name Name of the variable in \code{cohort} indicating the
#'   event-censoring status of each subject, where 1 indicates events. A
#'   \code{string}. If a vector is supplied, only the first element will be
#'   used.
#' @param t_name Name of the variable in \code{cohort} for the time of event or
#'   censoring. A \code{string}. Currently do not support scenarios where
#'   subjects enter the cohort at different time. If a vector is supplied, only
#'   the first element will be used.
#' @param match_var_name Name of the categorical variable in \code{cohort} to
#'   count-match on, which can be the exposure or the surrogate for the
#'   exposure. A \code{string}. If a vector is supplied, only the first element
#'   will be used.
#' @param include_var_name A string vector containing additional variables in
#'   \code{cohort} to include in the counter-matched NCC sample. Default is
#'   \code{NULL}.
#' @param ml Number of subjects to draw from each strata defined by the match
#'   variable (including the case). Default is 1.
#' @return Returns a \code{data.frame} of the counter-matched NCC sample.
#' @export
#' @examples 
#' data(cohort_1)
#' head(cohort_1)
#' # Counter-match on binary indicator for age:
#' cohort_1$age_bin <- as.numeric(cohort_1$age < 50)
#' ncc_cm_bin <- draw_ncc_cm(cohort = cohort_1, y_name = "y", t_name = "t", 
#'                           match_var_name = "age_bin", 
#'                           include_var_name = c("age", "gender"), ml = 1)
#' head(ncc_cm_bin, 10)
draw_ncc_cm <- function(cohort, y_name = NULL, t_name = NULL,
                        match_var_name = NULL, include_var_name = NULL, ml = 1) {
  cohort <- as.data.frame(cohort)
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
  cohort$row_id <- 1:nrow(cohort)
  cohort[, match_var_name] <- factor(cohort[, match_var_name])
  strata <- levels(cohort[, match_var_name])
  row_ids_case <- cohort$row_id[cohort[, y_name] == 1]
  df_ncc_cm <- do.call("rbind", lapply(row_ids_case, function(j) {
    # Subjects at risk at tj
    tj <- cohort[j, t_name]
    cohort_j <- cohort[cohort[, t_name] >= tj, ]
    match_var_j <- cohort[j, match_var_name]
    do.call("rbind", lapply(strata, function(l) {
      n_at_risk <- sum(cohort_j[, match_var_name] == l)
      row_ids_l <- cohort_j$row_id[cohort_j[, match_var_name] == l]
      if (match_var_j == l) {
        # Case is in strata l
        row_ids_l <- setdiff(row_ids_l, j)
        data.frame(set = j, row_id = c(j, sample(row_ids_l, size = ml - 1)),
                   t = tj, y = c(1, numeric(ml - 1)), match_var = l,
                   n_at_risk = n_at_risk, n_sampled = ml)
      } else {
        # Case is not in strata l
        data.frame(set = j, row_id = sample(row_ids_l, size = ml),
                   t = tj, y = numeric(ml), match_var = l,
                   n_at_risk = n_at_risk, n_sampled = ml)
      }
    }))
  }))
  ncc_cm <- data.frame(set = df_ncc_cm$set, row_id = df_ncc_cm$row_id,
                       t = df_ncc_cm$t,
                       n_at_risk = df_ncc_cm$n_at_risk,
                       n_sampled = df_ncc_cm$n_sampled,
                       weight = df_ncc_cm$n_at_risk / df_ncc_cm$n_sampled,
                       y = df_ncc_cm$y, match_var = df_ncc_cm$match_var)
  nc <- ncol(ncc_cm)
  names(ncc_cm)[c(3, 7, 8)] <- c(t_name, y_name, match_var_name)
  if (!is.null(include_var_name)) {
    ncc_cm <- cbind(ncc_cm, cohort[ncc_cm$row_id, include_var_name])
    names(ncc_cm)[(nc + 1):(nc + length(include_var_name))] <- include_var_name
    ncc_cm
  } else {
    ncc_cm
  }
}
#' Compute Kaplan-Meier type weights for (matched) NCC, possibly with dropped
#' controls
#' @inheritParams draw_ncc_cm
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing the NCC. A \code{string} vector. Default is \code{NULL}, i.e.,
#'   the NCC was only time-matched.
#' @param sample_stat A numeric vector containing sampling and status
#'   information for each subject in \code{cohort}: use 0 for non-sampled
#'   controls, 1 for sampled controls, and integers >=2 for events. The length 
#'   of this vector must be the same as the number of rows in \code{cohort}.
#' @param keep_stat A numeric vector indicating whether each subject in
#'   \code{cohort} are kept in the final NCC: use 1 for subjects who were in the
#'   final NCC, and 0 for subjects who were kept or never selected. The length
#'   of this vector must be the same as the number of rows in \code{cohort}.
#'   When unspecified, the function assumes all subjects are kept in the final
#'   NCC.
#' @param n_per_case Number of controls matched to each case.
#' @param n_kept Number of sampled controls in each set that were kept in the
#'   final NCC. When unspecified, the function assumes all subjects are kept in
#'   the final NCC.
#' @import survival
#' @import dplyr
#' @export
#' @examples 
#' # Load cohort data
#' data(cohort_1)
#' head(cohort_1)
#' # Load an NCC sampled drawn from cohort_1 using the following code
#' # ncc_1 <- Epi::ccwc(exit = t, fail = y, controls = 1, 
#' #                    include = list(age, gender), data = cohort_1)
#' data(ncc_1)
#' head(ncc_1)
#' # Map the NCC sample to the original cohort, break the matching, identify the 
#' # subjects selected into the NCC, and return this subset with KM type weights 
#' # computed for them.
#' # First create the sampling and status indicator:
#' sample_stat <- numeric(nrow(cohort_1))
#' sample_stat[unique(ncc_1$Map[ncc_1$Fail == 0])] <- 1
#' sample_stat[ncc_1$Map[ncc_1$Fail == 1]] <- 2
#' # Then find the sampled subset and compute weights:
#' ncc_nodup <- compute_km_weights(cohort = cohort_1, t_name = "t", y_name = "y",
#'                                 sample_stat = sample_stat, n_per_case = 1)
#' head(ncc_nodup)
compute_km_weights <- function(cohort, t_name = NULL, y_name = NULL,
                               sample_stat, keep_stat = NULL, 
                               match_var_names = NULL,
                               n_per_case, n_kept = NULL) {
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
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(cohort))) {
      stop(simpleError("Make sure all match variables are in cohort."))
    }
  }
  if (is.null(n_kept)) {
    keep_stat <- NULL
  }
  if (is.null(keep_stat)) {
    n_kept <- NULL
    keep_stat <- rep(1, nrow(cohort))
  } else {
    if (length(keep_stat) != nrow(cohort)) {
      stop(simpleError("Length of keep_stat must be same as number of rows in cohort."))
    }
  }
  t <- cohort[, t_name]
  y <- as.numeric(sample_stat >= 2)
  if (is.null(match_var_names)) {
    match_var <- rep(1, nrow(cohort))
    km <- survfit(Surv(t, y) ~ 1)
    km_summ <- summary(km)
    km_tb <- data.frame(t = km_summ$time, Rj = km_summ$n.risk - km_summ$n.event,
                        strata = "match_var=1", 
                        stringsAsFactors = FALSE) %>%
      arrange(strata, t)
  } else {
    mat <- cohort[, match_var_names]
    if (length(match_var_names) == 1) {
      mat <- matrix(mat, ncol = length(match_var_names))
    }
    match_var <- apply(mat, 1, function(row) paste(row, collapse = "-"))
    # Make sure the levels are integers starting from 1
    match_var <- factor(as.numeric(factor(match_var)))
    km <- survfit(Surv(t, y) ~ match_var)
    km_summ <- summary(km)
    km_tb <- data.frame(t = km_summ$time, Rj = km_summ$n.risk - km_summ$n.event,
                        strata = as.character(km_summ$strata), 
                        stringsAsFactors = FALSE) %>%
      arrange(strata, t)
  }
  if (is.null(n_kept)) { # All are kept
    km_tb <- km_tb %>%
      group_by(strata) %>% 
      mutate(prob_not_sampled = 1 - (n_per_case / Rj),
             cumulative_product = cumprod(prob_not_sampled),
             sampling_prob = 1 - cumulative_product)
  } else {
    km_tb <- km_tb %>%
      group_by(strata) %>% 
      # p_not_sampled = p(not selected) + p(selected but then dropped)
      mutate(prob_not_sampled = 1 - (n_per_case / Rj),
             prob_not_sampled2 = prob_not_sampled + # not selected
               (1 - prob_not_sampled) * 
               (choose(n = n_per_case - 1, k = n_kept) / 
                  choose(n = n_per_case, k = n_kept)), # selected, but not kept
             cumulative_product = cumprod(prob_not_sampled2),
             sampling_prob = 1 - cumulative_product)
  }
  in_ncc <- sample_stat > 0 & keep_stat == 1
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

#' Simulated cohort 1
#'
#' @format A data frame with 10000 rows and 5 variables:
#' \describe{
#'   \item{id}{Subject ID of each subject.}
#'   \item{y}{Event/censoring indicator. Events are indicated by y=1.}
#'   \item{t}{Event/censoring time (in years). Maximum follow-up time is 25 years.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#' }
"cohort_1"

#' Nested case-control (NCC) data 1
#' 
#' @description Time-matched data drawn from \code{\link{cohort_1}}, with one  
#' controls matched to each case.
#' 
#' @format A data frame with 1164 rows and 6 variables:
#' \describe{
#'   \item{Set}{Set ID.}
#'   \item{Map}{Row numbers in \code{cohort_1}.}
#'   \item{Time}{Event time (in years) of the case in each \code{Set}.}
#'   \item{Fail}{Case-control indicator.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#' }
"ncc_1"

#' Simulated cohort 2
#'
#' @format A data frame with 100000 rows and 8 variables:
#' \describe{
#'   \item{id}{Subject ID of each subject.}
#'   \item{y}{Event/censoring indicator. Events are indicated by y=1.}
#'   \item{t}{Event/censoring time (in years). Maximum follow-up time is 25 years.}
#'   \item{x}{Binary exposure.}
#'   \item{age}{Age of subjects (rounded to integers and mean-centred).}
#'   \item{age_cat}{Age category of subjects.}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#'   \item{z}{Binary effect modifier.}
#' }
"cohort_2"

#' Nested case-control (NCC) data 2
#' 
#' @description NCC data drawn from \code{\link{cohort_2}} matched on , with one  
#' controls matched to each case.
#' 
#' @format A data frame with 16638 rows and 9 variables:
#' \describe{
#'   \item{Set}{Set ID.}
#'   \item{Map}{Row numbers in \code{cohort_2}.}
#'   \item{Time}{Event time (in years) of the case in each \code{Set}.}
#'   \item{Fail}{Case-control indicator.}
#'   \item{age_cat}{Age category of subjects.}
#'   \item{gender}{Gender of subjects (1 for male and 0 for female).}
#'   \item{x}{Binary exposure.}
#'   \item{age}{Age of subjects (rounded to integers).}
#'   \item{z}{Binary effect modifier.}
#' }
"ncc_2"