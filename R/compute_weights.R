#' Assemble KM table and NCC data without duplicates when cohort data is
#' available
#' @inheritParams compute_km_weights
#' @import survival
#' @import dplyr
prep_km1 <- function(cohort, t_start_name = NULL, t_name, y_name, sample_stat, 
                     keep_stat = NULL, match_var_names = NULL, n_per_case) {
  cohort <- as.data.frame(cohort)
  if (!(y_name %in% names(cohort))) {
    stop(simpleError(paste(y_name, "not found in cohort.")))
  }
  if (!is.null(t_start_name)) {
    if (!(t_start_name %in% names(cohort))) {
      stop(simpleError(paste(t_start_name, "not found in cohort.")))
    }
    t_start <- cohort[, t_start_name] # Different subject start follow-up at different time
  } else {
    t_start <- rep(0, nrow(cohort)) # All subjects started follow-up at time 0
  }
  if (!(t_name %in% names(cohort))) {
    stop(simpleError(paste(t_name, "not found in cohort.")))
  }
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(cohort))) {
      stop(simpleError("Make sure all match variables are in cohort."))
    }
  }
  if (is.null(keep_stat)) {
    keep_stat <- rep(1, nrow(cohort))
  } else {
    if (length(keep_stat) != nrow(cohort)) {
      stop(simpleError("Length of keep_stat must be same as number of rows in cohort."))
    }
  }
  if (is.null(sample_stat)) {
    stop(simpleError("Please supply an indicator for cases and sampled controls."))
  } else if (length(sample_stat) != nrow(cohort)) {
    stop(simpleError("Length of sample_stat must be same as number of rows in cohort."))
  }
  t <- cohort[, t_name]
  y <- as.numeric(sample_stat >= 2)
  if (is.null(match_var_names)) {
    match_var <- rep(1, nrow(cohort))
    km <- survfit(Surv(t_start, t, y) ~ 1)
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
    km <- survfit(Surv(t_start, t, y) ~ match_var)
    km_summ <- summary(km)
    km_tb <- data.frame(t = km_summ$time, n.risk = km_summ$n.risk, 
                        n.event = km_summ$n.event, 
                        Rj = km_summ$n.risk - km_summ$n.event,
                        strata = as.character(km_summ$strata), 
                        stringsAsFactors = FALSE) %>%
      arrange(strata, t)
  }
  in_ncc <- sample_stat > 0 & keep_stat == 1
  list(km_tb = km_tb, match_var_ncc = paste0("match_var=", match_var[in_ncc]), 
       t_start = t_start[in_ncc], ncc_nodup = cohort[in_ncc, ])
}
#' Assemble KM table and NCC data without duplicates when cohort data is not
#' available
#' @inheritParams compute_km_weights
#' @import dplyr
prep_km2 <- function(ncc, n_at_risk, id_name, set_id_name, 
                     t_name, t_match_name, y_name, 
                     match_var_names = NULL, n_per_case) {
  ncc <- as.data.frame(ncc, stringsAsFactors = FALSE)
  n_at_risk <- as.data.frame(n_at_risk, stringsAsFactors = FALSE)
  if (!(id_name %in% names(ncc))) {
    stop(simpleError(paste(id_name, "not found in ncc")))
  }
  if (!(set_id_name %in% names(ncc))) {
    stop(simpleError(paste(set_id_name, "not found in ncc")))
  }
  if (!(y_name %in% names(ncc))) {
    stop(simpleError(paste(y_name, "not found in ncc")))
  }
  if (!(t_match_name %in% names(ncc))) {
    stop(simpleError(paste(t_match_name, "not found in ncc")))
  }
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(ncc))) {
      stop(simpleError("Make sure all match variables are in ncc."))
    }
  }
  if (!t_match_name %in% names(n_at_risk)) {
    stop(simpleError("n_at_risk must include a column for unique event times with the same column name as in ncc."))
  }
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(n_at_risk))) {
      stop(simpleError("n_at_risk must include one column for each matching variable and with the same column name as in ncc."))
    }
  }
  if (!"n.risk" %in% names(n_at_risk)) {
    stop(simpleError("n_at_risk must include one column for the number of subjects at risk at each event time with column name 'n.risk'."))
  }
  strata_symbol <- as.symbol("strata")
  t_symbol <- as.symbol("t")
  # Compile KM table and prepare match variable for ncc -----
  vars <- c(t_match_name, match_var_names)
  y_name_symbol <- as.symbol(y_name)
  km_summ <- ncc %>% 
    filter({{y_name_symbol}} == 1) %>% 
    unique() %>%
    group_by(across({{vars}})) %>% 
    summarise(n.event = n()) %>%
    left_join(n_at_risk)
  km_summ <- merge(km_summ, 
                   data.frame(t_match = ncc[, t_match_name], t = ncc[, t_name]), 
                   by.x = t_match_name, by.y = "t_match", all = TRUE)
  if (is.null(match_var_names)) {
    km_tb <- data.frame(t = km_summ$t, 
                        n.risk = km_summ$n.risk, n.event = km_summ$n.event, 
                        Rj = km_summ$n.risk - km_summ$n.event,
                        strata = "match_var=1", stringsAsFactors = FALSE) %>%
      arrange({{strata_symbol}}, {{t_symbol}})
    match_var_ncc <- rep("match_var=1", nrow(ncc))
  } else {
    # Create a strata variable for km_summ. Do not turn into factor
    mat <- km_summ[, match_var_names]
    if (length(match_var_names) == 1) {
      mat <- matrix(mat, ncol = length(match_var_names))
    }
    km_tb <- data.frame(
      t = km_summ$t, 
      n.risk = km_summ$n.risk, n.event = km_summ$n.event, 
      Rj = km_summ$n.risk - km_summ$n.event,
      strata = apply(mat, 1, function(row) paste(row, collapse = "-")), 
      stringsAsFactors = FALSE
    ) %>%
      arrange({{strata_symbol}}, {{t_symbol}})
  }
  rownames(km_tb) <- NULL
  # Obtain ncc_nodup from ncc -----
  # Break matching and keep only one row for each subject
  # Need to figure out from ncc who are cases and who are true controls
  id_name_symbol <- as.symbol(id_name)
  t_match_name_symbol <- as.symbol(t_match_name)
  set_id_name_symbol <- as.symbol(set_id_name)
  # Straightforward for subjects with y=1: they must be cases and each row
  # corresponds to a unique subject.
  ncc_nodup_cases <- ncc %>% 
    filter({{y_name_symbol}} == 1) %>% 
    select(-{{set_id_name_symbol}})
  # For subjects with y=0, only take them if they are not already taken as cases.
  # If there are multiple rows for the same subject, only keep the latest entry.
  ncc_nodup_controls <- ncc %>% 
    filter({{y_name_symbol}} == 0, 
           !({{id_name_symbol}} %in% ncc_nodup_cases[, id_name])) %>% 
    arrange(-{{t_match_name_symbol}}) %>% 
    select(-{{set_id_name_symbol}}) %>%
    filter(!duplicated({{id_name_symbol}}))
  ncc_nodup <- rbind(ncc_nodup_cases, ncc_nodup_controls) %>% 
    arrange({{id_name_symbol}})
  if (!is.null(match_var_names)) {# NULL case is already handled previously
    # Create a strata variable for ncc. Do not turn into factor
    mat_ncc <- ncc_nodup[, match_var_names]
    if (length(match_var_names) == 1) {
      mat_ncc <- matrix(mat_ncc, ncol = length(match_var_names))
    }
    match_var_ncc <- apply(mat_ncc, 1, function(row) paste(row, collapse = "-"))
  }
  list(km_tb = km_tb, match_var_ncc = match_var_ncc, ncc_nodup = ncc_nodup)
}
#' Prepare a template file for number of subjects at risk at each event time in
#' a nested case-control data
#' @description Given a nested case-control (NCC) data, create a template for
#'   user to fill in the number of subjects at risk in the cohort at each event
#'   time within each stratum (if applicable).
#' @param ncc Nested case-control data. A \code{data.frame} or a matrix with
#'   column names.
#' @param t_match_name Name of the column of event time in each matched set in
#'   \code{ncc}. A \code{string}. Time should be recorded in a reasonable
#'   resolution of measurement (e.g., in months or years) such that the number
#'   of subjects at risk in the cohort is available.
#' @param csv_file Name of CSV file to write the template to. If left
#'   unspecified, the template will be returned as a \code{data.frame} instead.
#' @return If \code{csv_file} is not supplied, this function returns a
#'   \code{data.frame} that contains the unique event times, additional matching
#'   variables (if \code{match_var_names} is supplied), and an empty column
#'   named \code{n.risk}. Column names in this \code{data.frame} should not be
#'   modified. If \code{csv_file} is supplied, this function write the
#'   \code{data.frame} to the designated file instead. Users should fill this
#'   column with the number of subjects at risk at each time point in the
#'   stratum, and later use this \code{data.frame} or file to compute the
#'   KM-type weights for the subjects in the NCC data.
#' @import dplyr
#' @export
#' @examples 
#' data(ncc_2)
#' n_at_risk <- prep_n_at_risk(ncc = ncc_2, t_match_name = "Time", y_name = "Fail", 
#'                             match_var_names = c("age_cat", "gender"))
prep_n_at_risk <- function(ncc, t_match_name, y_name, match_var_names = NULL, 
                           csv_file = NULL) {
  y_name <- y_name[1]
  if (!(y_name %in% names(ncc))) {
    stop(simpleError(paste(y_name, "not found in ncc.")))
  }
  t_match_name <- t_match_name[1]
  if (!(t_match_name %in% names(ncc))) {
    stop(simpleError(paste(t_match_name, "not found in ncc.")))
  }
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(ncc))) {
      stop(simpleError("Make sure all match variables are in ncc."))
    }
  }
  vars <- c(t_match_name, match_var_names)
  y_name_symbol <- as.symbol(y_name)
  n_at_risk <- ncc %>% 
    filter({{y_name_symbol}} == 1) %>% 
    group_by(across({{vars}})) %>% 
    summarise(n.risk = n())
  n_at_risk$n.risk <- NA
  if (is.null(csv_file)) {
    as.data.frame(n_at_risk)
  } else {
    write.csv(n_at_risk, file = csv_file, na = "", row.names = FALSE)
    message(paste0("Cohort summary template is saved as '", csv_file, "'."))
    NULL
  }
}
#' Compute Kaplan-Meier type weights for (matched) NCC, possibly with dropped
#' controls
#' @inheritParams draw_ncc_cm
#' @param ncc (Matched) NCC data, if \code{cohort} is not available. A
#'   \code{data.frame} or a matrix with column names. In addition to the time of 
#'   event in each matched set, \code{ncc} must also include another column for 
#'   the actual time of event or censoring of each subject. Weights computed 
#'   based on the time of event in each matched set (and subsequent weighted 
#'   analyses of this time) is generally inappropriate.
#' @param n_at_risk Number of subjects at risk at time of each cases in the NCC,
#'   if \code{cohort} is not available. A \code{data.frame} or a matrix with
#'   column names. See Details.
#' @param t_start_name Name of the variable in \code{cohort} for the start time
#'   of follow-up. A \code{string}. Default is \code{NULL}, i.e., every subject
#'   started the follow-up at time 0. This is not yet implemented when the
#'   cohort data is not available.
#' @param t_name Name of the variable in \code{cohort} for the time of event or
#'   censoring. A \code{string}. Note that if \code{ncc} is supplied, in order
#'   to correctly compute the weight for each sampled control this should be the
#'   actual time of censoring, not the time of event of the case in the same
#'   matched set.
#' @param t_match_name Name of the column of event time in each matched set in
#'   \code{ncc}. A \code{string}. Time should be recorded in a reasonable
#'   resolution of measurement (e.g., in months or years) such that the number
#'   of subjects at risk in the cohort is available. The corresponding column in
#'   \code{n_at_risk} must have the same name. See Details.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing the NCC. A \code{string} vector. Default is \code{NULL}, i.e.,
#'   the NCC was only time-matched. The corresponding column in \code{n_at_risk}
#'   must have the same name.
#' @param id_name Name of the column indicating subject ID in \code{ncc}, if
#'   \code{cohort} is not available.
#' @param set_id_name Name of the column indicating the matched sets in
#'   \code{ncc}, if \code{cohort} is not available.
#' @param sample_stat A numeric vector containing sampling and status
#'   information for each subject in \code{cohort}: use 0 for non-sampled
#'   controls, 1 for sampled controls, and integers >=2 for events. The length 
#'   of this vector must be the same as the number of rows in \code{cohort}.
#' @param keep_stat A numeric vector indicating whether each subject in
#'   \code{cohort} are kept in the final NCC: use 1 for subjects who were in the
#'   final NCC, and 0 for subjects who were kept or never selected. The length
#'   of this vector must be the same as the number of rows in \code{cohort}.
#'   When unspecified, the function keeps all subjects in the final NCC.
#' @param n_per_case Number of controls matched to each case.
#' @param n_kept Number of sampled controls in each set that were kept in the
#'   final NCC. When unspecified, the function assumes all subjects are kept in
#'   the final NCC.
#' @param attach_weight If \code{TRUE} (default), returns a \code{data.frame}
#'   containing all the unique subjects selected in the NCC sample, with a
#'   column for the KM-type weight associated with each subject. If \code{FALSE}, 
#'   returns a \code{data.frame} containing the KM-type weight (with column name
#'   \code{sampling_prob}) for controls that were available to be sampled up
#'   until each event time. 
#' @details When the full cohort is not available, in order to compute the
#'   correct weights for each sampled control in the NCC sample, it is important
#'   to keep the actual time of event or censoring of each subject in the NCC
#'   sample, which should be specified as \code{t_name} in the input. Since the 
#'   number of subjects in each risk set will be supplied separately (i.e., as 
#'   \code{n_at_risk}) in such scenario, \code{t_match_name} is required to 
#'   map each control to the appropriate risk set. \code{t_match_name} may be 
#'   the same as \code{t_name} if the exact risk set is available in 
#'   \code{n_at_risk}, but when the full cohort is not available the risk set is 
#'   usually approximated by using a coarsened version of \code{t_name}. For 
#'   example, when controls were drawn from a population registry by matching on
#'   the exact date of death of cases, birth cohort and gender, the number at
#'   risk may be approximated by using the population size in the year of event
#'   in the same birth cohort of the same gender. In this scenario
#'   \code{t_match_name} would be the year of \code{t_name}.
#' @import dplyr
#' @export
#' @seealso \code{\link{prep_n_at_risk}}
#' @examples 
#' # Load cohort data
#' data(cohort_2)
#' head(cohort_2)
#' # Load an NCC sampled drawn from cohort_2 using the following code
#' # ncc_2 <- Epi::ccwc(exit = t, fail = y, controls = 5, 
#' #                    match = list(age_cat, gender), include = list(x, age, z), 
#' #                    data = cohort_2)
#' data(ncc_2)
#' head(ncc_2)
#' # Map the NCC sample to the original cohort, break the matching, identify the 
#' # subjects selected into the NCC, and return this subset with KM type weights 
#' # computed for them.
#' # First create the sampling and status indicator:
#' sample_stat <- numeric(nrow(cohort_2))
#' sample_stat[unique(ncc_2$Map[ncc_2$Fail == 0])] <- 1
#' sample_stat[ncc_2$Map[ncc_2$Fail == 1]] <- 2
#' # Then find the sampled subset and compute weights:
#' ncc_nodup <- compute_km_weights(cohort = cohort_2, t_name = "t", y_name = "y",
#'                                 sample_stat = sample_stat, n_per_case = 5)
#' head(ncc_nodup)
#' # Alternatively, if the cohort is not available, the weights can be computed 
#' # as long as number of subjects at risk at event times in each strata is 
#' # available elsewhere, and the actual time of event/censoring is available 
#' # for each subject in the NCC.
#' # Load the number of subjects at risk from cohort_2:
#' data(n_at_risk)
#' head(n_at_risk)
#' # Note that column names in `n_at_risk` must match column names in `ncc_2`.
#' # In ncc_2, `Time` is the time of event within each matched set, and `t` is 
#' # the actual time of event/censoring. The following command computes the same 
#' # weights as in ncc_nodup:
#' ncc_nodup_v2 <- compute_km_weights(ncc = ncc_2, n_at_risk = n_at_risk, 
#'                                    t_name = "t", y_name = "Fail", 
#'                                    t_match_name = "Time",
#'                                    id_name = "Map", set_id_name = "Set", 
#'                                    match_var_names = c("age_cat", "gender"), 
#'                                    n_per_case = 5)
#' head(ncc_nodup_v2)
compute_km_weights <- function(cohort = NULL, ncc = NULL, n_at_risk = NULL, 
                               t_start_name = NULL, t_name = NULL, y_name = NULL, 
                               t_match_name = t_name, 
                               id_name = NULL, set_id_name = NULL,
                               sample_stat = NULL, keep_stat = NULL, 
                               match_var_names = NULL,
                               n_per_case, n_kept = NULL, 
                               attach_weight = TRUE) {
  if (is.null(y_name)) {
    stop(simpleError("Please sapply name of event status."))
  } else {
    y_name <- y_name[1]
  }
  if (is.null(t_name)) {
    stop(simpleError("Please sapply name of event/censoring time."))
  } else {
    t_name <- t_name[1]
  }
  if (is.null(n_kept)) {
    keep_stat <- NULL
  }
  if (!is.null(cohort)) {
    # Full cohort is available
    if (!is.null(t_start_name)) {
      message(simpleMessage(sprintf(
        "Start time is given by variable %s. Event/censoring time is given by variable %s.\n", 
        t_start_name, t_name
      )))
    } else {
      message(simpleMessage(sprintf(
        "Start time is 0 for all subjects. Event/censoring time is given by variable %s.\n", 
        t_name
      )))
    }
    obj <- prep_km1(cohort = cohort, t_start_name = t_start_name, t_name = t_name, 
                    y_name = y_name, sample_stat = sample_stat, 
                    match_var_names = match_var_names, n_per_case = n_per_case)
    if (is.null(n_kept)) { # All are kept
      km_tb <- obj$km_tb %>%
        group_by(strata) %>% 
        mutate(prob_not_sampled = 1 - (n_per_case / Rj))
    } else {
      km_tb <- obj$km_tb %>%
        group_by(strata) %>% 
        # p_not_sampled = p(not selected) + p(selected but then dropped)
        mutate(prob_not_sampled0 = 1 - (n_per_case / Rj),
               prob_not_sampled = prob_not_sampled0 + # not selected
                 (1 - prob_not_sampled0) * 
                 (choose(n = n_per_case - 1, k = n_kept) / 
                    choose(n = n_per_case, k = n_kept))) # selected, but not kept
    }
    if (!attach_weight) {
      if (is.null(t_start_name)) {
        km_tb <- km_tb %>%
          ungroup() %>%
          group_by(strata) %>% 
          mutate(cumulative_product = cumprod(prob_not_sampled),
                 sampling_prob = 1 - cumulative_product, 
                 km_weight = 1 / sampling_prob) %>% 
          as.data.frame(stringsAsFactors = FALSE)
        return(km_tb)
      } else {
        warning(simpleWarning("No common KM-type weight associated with each event time when subjects do not share a common start time."))
        return(NULL)
      }
    }
    do.call("rbind", lapply(1:nrow(obj$ncc_nodup), function(j) {
      ncc_nodup_j <- obj$ncc_nodup[j, ]
      if (ncc_nodup_j[, y_name] == 1) {
        km_prob <- 1
      } else {
        # A subject is in the risk set at time t if this subject is still under
        # observation at t-, i.e., right before t. Hence a subject censored
        # exactly at time t is still in the risk set for an event at t.
        km_tb_i <- km_tb %>% 
          filter(strata == obj$match_var_ncc[j], 
                 t > obj$t_start[j], t <= ncc_nodup_j[, t_name])
        if (nrow(km_tb_i) == 0) {
          km_prob <- 0
        } else {
          km_prob <- 1 - prod(km_tb_i$prob_not_sampled)
        }
      }
      ncc_nodup_j$km_prob <- km_prob
      ncc_nodup_j$km_weight <- 1 / km_prob
      ncc_nodup_j
    }))
  } else {
    # Full cohort is not available
    if (!is.null(t_start_name)) {
      stop(simpleError("Staggered time not yet implemented when full cohort is not available."))
    }
    if (is.null(ncc)) {
      stop(simpleError("If full cohort is not available, please supply the ncc data."))
    }
    if (is.null(n_at_risk)) {
      stop(simpleError("If full cohort is not available, please supply the number of subjects at risk prepared using 'prep_n_at_risk' function."))
    }
    if (is.null(id_name)) {
      stop(simpleError("If full cohort is not available, please supply the column name of subject ID in the ncc data."))
    } else {
      id_name <- id_name[1]
    }
    if (is.null(set_id_name)) {
      stop(simpleError("If full cohort is not available, please supply the column name of ID of matched sets in the ncc data."))
    } else {
      set_id_name <- set_id_name[1]
    }
    obj <- prep_km2(ncc = ncc, n_at_risk = n_at_risk, 
                    id_name = id_name, set_id_name = set_id_name,
                    t_name = t_name, t_match_name = t_match_name, y_name = y_name, 
                    match_var_names = match_var_names, 
                    n_per_case = n_per_case)
    if (is.null(n_kept)) { # All are kept
      km_tb <- obj$km_tb %>%
        group_by(strata) %>% 
        mutate(prob_not_sampled = 1 - (n_per_case / Rj),
               cumulative_product = cumprod(prob_not_sampled),
               sampling_prob = 1 - cumulative_product)
    } else {
      km_tb <- obj$km_tb %>%
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
    if (!attach_weight) {
      km_tb <- km_tb %>% mutate(km_weight = 1 / sampling_prob) %>% 
        as.data.frame(stringsAsFactors = FALSE)
      return(km_tb)
    }
    p_ncc <- unlist(lapply(1:nrow(obj$ncc_nodup), function(j) {
      if (obj$ncc_nodup[j, y_name] == 1) {
        1
      } else {
        # A subject is in the risk set at time t if this subject is still under
        # observation at t-, i.e., right before t. Hence a subject censored
        # exactly at time t is still in the risk set for an event at t.
        km_tb_i <- km_tb[km_tb$t <= obj$ncc_nodup[j, t_name] &
                           km_tb$strata == obj$match_var_ncc[j], ]
        r <- nrow(km_tb_i)
        if (r == 0) {
          0
        } else {
          km_tb_i$sampling_prob[r]
        }
      }
    }))
    cbind(obj$ncc_nodup, km_prob = p_ncc, km_weight = 1 / p_ncc)
  }
}
