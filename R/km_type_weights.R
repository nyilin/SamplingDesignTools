#' <Private function> Prepare the skeleton cohort for subsequent steps.
#' @param cohort Cohort data with at least the following information on each
#'   subject: start time (if not 0 for all subjects) and end time of follow-up,
#'   censoring status and matching variables (if any). A \code{data.frame} or a
#'   matrix with column names.
#' @param t_start_name Name of the variable in \code{cohort} for the start time
#'   of follow-up. A \code{string}.
#' @param t_name Name of the variable in \code{cohort} for the event or
#'   censoring time. A \code{string}.
#' @param y_name Name of the column of censoring status in each matched set in
#'   \code{cohort}, with 1 for event and 0 for censoring. A \code{string}.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing the NCC. A \code{string} vector.
#' @param print_message Whether to print the message regarding start and exit
#'   time. Default is \code{TRUE}.
#' @import dplyr
prepare_cohort <- function(cohort, t_start_name, t_name, y_name,
                           match_var_names, print_message = TRUE) {
  cohort <- unique(as.data.frame(cohort, stringsAsFactors = FALSE))
  if (!(y_name %in% names(cohort))) {
    stop(simpleError(paste(y_name, "not found in cohort.")))
  } else {
    cohort$.y <- cohort[, y_name[1]]
  }
  if (!is.null(t_start_name)) {
    if (!(t_start_name %in% names(cohort))) {
      stop(simpleError(paste(t_start_name, "not found in cohort.")))
    }
    if (print_message) {
      message(simpleMessage(sprintf(
        "Start time is given by variable %s. Event/censoring time is given by variable %s.\n", 
        t_start_name, t_name
      )))
    }
    cohort$.t_start <- cohort[, t_start_name]
  } else {
    if (print_message) {
      message(simpleMessage(sprintf(
        "Start time is 0 for all subjects. Event/censoring time is given by variable %s.\n", 
        t_name
      )))
    }
    cohort$.t_start <- -Inf # Not used in Surv()
  }
  if (!(t_name %in% names(cohort))) {
    stop(simpleError(paste(t_name, "not found in cohort.")))
  } else {
    cohort$.t <- cohort[, t_name[1]]
  }
  if (is.null(match_var_names)) {
    cohort %>% mutate(.strata0 = 1, .strata = "match_var=1")
  } else {
    if (!all(match_var_names %in% names(cohort))) {
      stop(simpleError("Make sure all match variables are in cohort."))
    }
    mat <- unique(as.data.frame(cohort[, match_var_names]))
    names(mat) <- match_var_names
    match_var <- apply(mat, 1, function(row) paste(row, collapse = "-"))
    # Make sure the levels are integers starting from 1
    mat$.strata0 <- factor(as.numeric(factor(match_var)))
    mat$.strata <- paste0("match_var=", mat$.strata0)
    cohort %>% left_join(mat)
  }
}
#' <Private function> Prepare the NCC cases for subsequent steps.
#' @param ncc_cases Cases in NCC data. A \code{data.frame} or a matrix with
#'   column names. This data should not include the ID of each matched set.
#' @param t_name Name of the variable in \code{ncc_cases} for event time. A
#'   \code{string}.
#' @param match_var_names Name(s) of the match variable(s) in
#'   \code{ncc_cases} used when drawing the NCC. A \code{string} vector.
#' @import dplyr
prepare_ncc_cases <- function(ncc_cases, t_name, match_var_names) {
  ncc <- unique(as.data.frame(ncc_cases, stringsAsFactors = FALSE))
  if (!(t_name %in% names(ncc))) {
    stop(simpleError(paste(t_name, "not found in ncc.")))
  } else {
    ncc$.t_event <- ncc[, t_name[1]]
  }
  if (is.null(match_var_names)) {
    ncc %>% mutate(.strata0 = 1, .strata = "match_var=1")
  } else {
    if (!all(match_var_names %in% names(ncc))) {
      stop(simpleError("Make sure all match variables are in ncc."))
    }
    mat <- unique(as.data.frame(ncc[, match_var_names]))
    names(mat) <- match_var_names
    match_var <- apply(mat, 1, function(row) paste(row, collapse = "-"))
    # Make sure the levels are integers starting from 1
    mat$.strata0 <- factor(as.numeric(factor(match_var)))
    mat$.strata <- paste0("match_var=", mat$.strata0)
    ncc %>% left_join(mat)
  }
}
#' <Private function> Compute number at risk at each event time from cohort
#' @param cohort Cohort data prepared by \code{\link{prepare_cohort}}.
#' @param match_var_names Name(s) of the match variable(s).
#' @param staggered Whether time is staggered (i.e., t_start is specified in
#'   input).
#' @return Returns a data.frame with the following columns:
#' \describe{
#'   \item{t_event}{Unique event times in \code{cohort_skeleton}.}
#'   \item{n_event}{Number of events at each event time.}
#'   \item{n_at_risk}{Number of subjects at risk at each event time.}
#'   \item{strata}{Strata defined by matching variables. \code{match_var=1} if the
#' NCC is only time-matched.}
#'   \item{<each matching variable>}{If the NCC is matched on additional
#' confounders, each matching variable will be included as a column to the right
#' of \code{strata}.}
#' }
#' @import survival
#' @import dplyr
compute_risk_tb <- function(cohort, match_var_names, staggered) {
  match_var <- cohort$.strata0
  if (staggered) {
    km <- survfit(Surv(.t_start, .t, .y) ~ match_var, data = cohort)
  } else {
    km <- survfit(Surv(.t, .y) ~ match_var, data = cohort)
  }
  km_summ <- summary(km)
  if (is.null(km_summ$strata)) {
    risk_table <- data.frame(t_event = km_summ$time, n_event = km_summ$n.event, 
                             n_at_risk = km_summ$n.risk,
                             strata = "match_var=1",
                             stringsAsFactors = FALSE) %>%
      arrange(strata, t_event)
  } else {
    mat <- unique(cohort[, c(".strata", match_var_names)])
    names(mat) <- c("strata", match_var_names)
    risk_table <- data.frame(t_event = km_summ$time, n_event = km_summ$n.event, 
                             n_at_risk = km_summ$n.risk,
                             strata = as.character(km_summ$strata), 
                             stringsAsFactors = FALSE) %>%
      arrange(strata, t_event) %>% 
      left_join(unique(mat))
  }
  n0 <- sum(risk_table$n_at_risk == 0)
  if (n0 > 0) {
    warning(simpleWarning(sprintf(
      "No subject at risk for %d of the %d event times.", n0, nrow(risk_table)
    )))
  }
  risk_table
}
#' <Private function> Compute the probability that an eligible non-case is not
#' sampled at each event time
#' @param n_event A vector of number of event at each event time (1 if no tie in
#'   event time, or >1 in the presence of ties). Non-integers are rounded down.
#' @param n_at_risk A vector of the number of subject at risk at each event
#'   time. Non-integers are rounded down.
#' @param n_per_case Number of controls to select per case. A positive integer.
#'   Non-integers are rounded down.
#' @return Returns a vector of probabilities.
#' @import dplyr
p_not_sampled <- function(n_event, n_at_risk, n_per_case) {
  n_event <- as.integer(as.vector(n_event))
  n_at_risk <- as.integer(as.vector(n_at_risk))
  if (anyNA(n_event)) stop(simpleError("Some entries in n_event are not numeric."))
  if (anyNA(n_at_risk)) stop(simpleError("Some entries in n_at_risk are not numeric."))
  if (length(n_event) != length(n_at_risk)) {
    stop(simpleError("n_event, n_at_risk and strata should be vectors of the same length."))
  }
  if (any(n_at_risk < n_event)) {
    stop(simpleError("Some entries in n_at_risk are smaller than n_event. Number at risk must be larger than or equal to number of events."))
  }
  n_per_case <- as.integer(n_per_case)
  if (is.na(n_per_case) || n_per_case < 1) {
    stop(simpleError("n_per_case must be a positive integer."))
  }
  df <- data.frame(n_event = n_event, Rj = n_at_risk - n_event) %>% 
    mutate(n_target = n_event * n_per_case, 
           n_sampled = ifelse(Rj > n_target, n_target, Rj), 
           p0 = 1 - n_sampled / Rj, # not sampled
           p = p0)
  as.numeric(df$p)
}
#' <Private function> Compute KM-type weight from risk table, provided t=0 for all
#' @param risk_table Output from \code{\link{compute_risk_tb}}, with an additional
#'   column \code{p}.
#' @return Returns a \code{data.frame} with \code{t_event}, \code{km_weight} and
#'   matching variables (if any).
#' @import dplyr
compute_kmw0 <- function(risk_table) {
  km_table <- risk_table %>%
    # filter(!is.na(p)) %>%
    group_by(strata) %>% 
    mutate(cumulative_product = cumprod(p),
           sampling_prob = 1 - cumulative_product, 
           km_weight = 1 / sampling_prob) 
  if (any(is.na(km_table$km_weight))) {
    warning(simpleWarning(sprintf("%d entries in KM-type weight are NA.", 
                                  sum(is.na(km_table$km_weight)))))
  }
  if (any(is.infinite(km_table$km_weight))) {
    warning(simpleWarning(sprintf("%d entries in KM-type weight are infinite.", 
                                  sum(is.infinite(km_table$km_weight)))))
  }
  match_var_names <- setdiff(names(risk_table), 
                             c("t_event", "n_event", "n_at_risk", "strata", "p"))
  km_table[, c("t_event", "n_event", "strata", match_var_names, "km_weight")] %>% 
    arrange(strata, t_event) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
#' <Private function> Assign appropriate KM-type weight to unique subjects in NCC
#' @inheritParams compute_kmw0
#' @param ncc_nodup NCC data prepared using \code{\link{prepare_cohort}}
#'   function. This data should not include the ID of each matched set or the
#'   time of event in each set.
#' @return Returns \code{ncc_nodup} with two additional columns: \code{.km_prob}
#'   for KM-type probability, and \code{.km_weight} for KM-type weight.
assign_kmw0 <- function(ncc_nodup, risk_table) {
  ncc_nodup <- do.call("rbind", lapply(1:nrow(ncc_nodup), function(j) {
    ncc_nodup_j <- ncc_nodup[j, ]
    if (ncc_nodup_j$.y == 1) {
      km_prob <- 1
    } else {
      # A subject is in the risk set at time t if this subject is still under
      # observation at t-, i.e., right before t. Hence a subject censored
      # exactly at time t is still in the risk set for an event at t.
      risk_table_i <- risk_table %>% 
        filter(strata == ncc_nodup_j$.strata, 
               t_event > ncc_nodup_j$.t_start, t_event <= ncc_nodup_j$.t)
      if (nrow(risk_table_i) == 0) {
        km_prob <- 0
      } else {
        km_prob <- 1 - prod(risk_table_i$p)
      }
    }
    ncc_nodup_j %>% mutate(.km_prob = km_prob, .km_weight = 1 / km_prob) %>% 
      select(-.strata0, -.strata, -.y, -.t_start, -.t) %>% 
      as.data.frame(stringsAsFactors = FALSE)
  }))
  if (any(is.na(ncc_nodup$.km_weight))) {
    warning(simpleWarning(sprintf("%d entries in KM-type weight are NA.", 
                                  sum(is.na(ncc_nodup$.km_weight)))))
  }
  if (any(is.infinite(ncc_nodup$.km_weight))) {
    warning(simpleWarning(sprintf("%d entries in KM-type weight are infinite.", 
                                  sum(is.infinite(ncc_nodup$.km_weight)))))
  }
  ncc_nodup
}
#' <Private function> Compute KM-type weights for NCC sample given full cohort
#' @inheritParams prepare_cohort
#' @param t_start_name Name of the variable in \code{cohort_skeleton} for the
#'   start time of follow-up. A \code{string}. Default is \code{NULL}, where all 
#'   subjects started the follow-up at time 0.
#' @param sample_stat A numeric vector containing sampling and status
#'   information for each subject in \code{cohort}: use 0 for non-sampled
#'   controls, 1 for sampled (and kept) controls, and integers >=2 for events.
#'   The length of this vector must be the same as the number of rows in
#'   \code{cohort}.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing the NCC. A \code{string} vector. Default is \code{NULL}, i.e.,
#'   the NCC was only time-matched. 
#' @param n_per_case Number of controls matched to each case.
#' @param return_risk_table Whether the risk table should be returned. Default
#'   is \code{FALSE}.
#' @param km_names Column names for the KM-type probability (the first element)
#'   and weight (the second element) computed, if these two columns are to be 
#'   attached to each subject in the input data. Default is
#'   \code{c(".km_prob", ".km_weight")}.
#' @return If \code{return_risk_table = FALSE} (the default), returns the
#'   subcohort of sampled subjects with the appropriate KM-type probability and
#'   weight attached to each subject. If \code{return_risk_table = TRUE},
#'   returns a list containing this subcohort (\code{dat}) and the risk table
#'   (\code{risk_table}), which is a \code{data.frame} containing the distinct
#'   event time (\code{t_event}), matching variables (if any), and the number of
#'   subject at risk at each event time in each strata defined by matching
#'   variables (\code{n_at_risk}).
#' @import dplyr
#' @import purrr
compute_kmw_cohort <- function(cohort, t_start_name = NULL, t_name, sample_stat, 
                               match_var_names = NULL, n_per_case, 
                               return_risk_table = FALSE, 
                               km_names = c(".km_prob", ".km_weight")) {
  cohort <- unique(as.data.frame(cohort, stringsAsFactors = FALSE))
  if (length(sample_stat) != nrow(cohort)) {
    stop(simpleError("The length of sample_stat must be the same as the number of rows in cohort."))
  }
  cohort$.y <- as.numeric(sample_stat >= 2)
  cohort <- prepare_cohort(cohort = cohort, t_start_name = t_start_name, 
                           t_name = t_name, y_name = ".y", 
                           match_var_names = match_var_names)
  risk_table <- compute_risk_tb(cohort = cohort, match_var_names = match_var_names, 
                                staggered = !is.null(t_start_name)) %>% 
    mutate(p = p_not_sampled(n_event = n_event, n_at_risk = n_at_risk, 
                             n_per_case = n_per_case))
  # Assign km_weight to each subject
  dat <- assign_kmw0(ncc_nodup = cohort[sample_stat > 0, ], 
                     risk_table = risk_table)
  dat <- change_km_names(dat = dat, km_names = km_names)
  if (return_risk_table) { 
    # Return time of event and number of subject at risk at each event time
    list(dat = dat, 
         risk_table = risk_table[, c("t_event", "n_event", match_var_names, "n_at_risk")])
  } else { 
    dat
  }
}
#' Match number at risk at coarsened time to NCC cases
#' @inheritParams prepare_ncc_cases
#' @param risk_table_manual A \code{data.frame} with columns \code{t_event}
#'   (unique event times, possibly coarsened), \code{n_at_risk} (number of
#'   subjects at risk at each \code{t_event} in the underlying cohort), and
#'   additional columns for matching variables, if any. Make sure the matching
#'   variables have the same column names as in \code{ncc_cases} and
#'   \code{match_var_names}.
#' @param t_coarse_name Name of the column of event time in each matched set in
#'   \code{ncc}, possibly coarsened to the same level as \code{t_event} in
#'   \code{risk_table_manual}. A \code{string}. 
#' @param t_name Name of the column of the exact event time. A \code{string}. 
#' @return Returns a data.frame with the following columns:
#' \describe{
#'   \item{t_event}{Unique event times (exact, not coarsened).}
#'   \item{n_event}{Number of events at each event time.}
#'   \item{n_at_risk}{Number of subjects at risk at each event time.}
#'   \item{strata}{Strata defined by matching variables. \code{match_var=1} if the
#' NCC is only time-matched.}
#'   \item{<each matching variable>}{If the NCC is matched on additional
#' confounders, each matching variable will be included as a column to the right
#' of \code{strata}.}
#' }
#' @import dplyr
#' @export
match_risk_table <- function(ncc_cases, risk_table_manual, t_coarse_name, t_name, 
                             match_var_names = NULL) {
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(risk_table_manual))) {
      stop(simpleError("Make sure all match variables are in risk_table_manual."))
    }
  }
  risk_table_manual <- risk_table_manual %>% rename(.t_event = t_event)
  vars <- c(".t_event", match_var_names)
  risk_table <- prepare_ncc_cases(ncc_cases = ncc_cases, 
                                  t_name = t_coarse_name, 
                                  match_var_names = match_var_names) %>% 
    group_by(across({{vars}})) %>% 
    summarise(n_event = n(), strata = unique(.strata)) %>%
    ungroup() %>%
    left_join(risk_table_manual) %>% 
    rename(t_coarse = .t_event)
  risk_table[, c("t_coarse", "n_event", "n_at_risk", "strata", match_var_names)] %>% 
    left_join(unique(data.frame(t_coarse = ncc_cases[, t_coarse_name], 
                                t_event = ncc_cases[, t_name])), 
              by = "t_coarse") %>% 
    arrange(strata, t_event) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
#' <Private function> Compute KM-type weight for NCC sample given information on 
#' underlying cohort
#' @inheritParams match_risk_table
#' @param ncc NCC data. A \code{data.frame} or a matrix with column names. This
#'   data should not include the ID of each matched set, but should include the
#'   actual event/censoring time of each subject.
#' @param id_name Name of the column of the subject ID. A \code{string}. 
#' @param t_name Name of the column of the exact event/censoring time of each
#'   subject. Note that for controls, this should not be the time of the event
#'   in the same matched set. A \code{string}.
#' @param t_match_name Name of the column of event time in each matched set in
#'   \code{ncc}, possibly coarsened to the same level as \code{t_event} in
#'   \code{risk_table_manual}. A \code{string}. Default is \code{t_name}, i.e., 
#'   not coarsened.
#' @param y_name Name of the column of censoring status in each matched set in
#'   \code{ncc}, with 1 for event and 0 for censoring. A \code{string}.
#' @param n_per_case Number of controls matched to each case.
#' @param return_risk_table Whether the risk table should be returned. Default
#'   is \code{FALSE}.
#' @param km_names Column names for the KM-type probability (the first element)
#'   and weight (the second element) computed, if these two columns are to be 
#'   attached to each subject in the input data. Default is
#'   \code{c(".km_prob", ".km_weight")}.
#' @return If \code{return_risk_table = FALSE} (the default), returns a
#'   \code{data.frame} containing all the unique subjects selected in the NCC
#'   sample, with a column for the KM-type weight associated with each subject.
#'   If \code{return_risk_table = TRUE}, returns a list containing this
#'   subcohort (\code{dat}) and the risk table (\code{risk_table}), which is a
#'   \code{data.frame} containing the distinct (and exact) event time
#'   (\code{t_event}), matching variables (if any), and the number of subject at
#'   risk at each event time in each strata defined by matching variables
#'   (\code{n_at_risk}).
#' @import dplyr
compute_kmw_ncc <- function(ncc, id_name, risk_table_manual, 
                            t_start_name = NULL, t_name, t_match_name = t_name, 
                            y_name, match_var_names = NULL, n_per_case, 
                            return_risk_table = FALSE, 
                            km_names = c(".km_prob", ".km_weight")) {
  message(simpleMessage("Make sure input ncc does not include ID of matched sets. Failing to do so results in Errors.\n"))
  ncc <- unique(as.data.frame(ncc, stringsAsFactors = FALSE))
  risk_table_manual <- as.data.frame(risk_table_manual, stringsAsFactors = FALSE)
  if (!(y_name %in% names(ncc))) {
    stop(simpleError(paste(y_name, "not found in ncc.")))
  } else {
    y_name <- y_name[1]
  }
  if (is.null(id_name)) {
    stop(simpleError("Cannot attach weight to each subject if subject id is unknown."))
  } else if (!(id_name %in% names(ncc))) {
    stop(simpleError(paste(id_name, "not found in ncc.")))
  }
  if (is.null(t_name)) {
    stop(simpleError("Cannot attach weight to each subject if the event/censoring time is unknown."))
  } else if (!(t_name %in% names(ncc))) {
    stop(simpleError(paste(t_name, "not found in ncc.")))
  }
  if (is.null(t_match_name)) {
    stop(simpleError("Please specify the variable for event time in each matched set."))
  } else if (!(t_match_name %in% names(ncc))) {
    stop(simpleError(paste(t_match_name, "not found in ncc.")))
  }
  if (!is.null(t_start_name)) {
    if (!(t_start_name %in% names(ncc))) {
      stop(simpleError(paste(t_start_name, "not found in ncc.")))
    }
  }
  risk_tb_vars <- c("t_event", "n_event", "n_at_risk", match_var_names)
  if (!all(risk_tb_vars %in% names(risk_table_manual))) {
    stop(simpleError(
      paste("risk_table_manual must include the following columns:", 
            toString(risk_tb_vars))
    ))
  }
  if (!all(ncc[, t_match_name] %in% risk_table_manual$t_event)) {
    stop(simpleError(sprintf(
      "There is mismatch in event time in ncc ('%s') and risk_table_manual$t_event.", 
      t_match_name
    )))
  }
  # Make sure each row in ncc_cases uniquely corresponds to a case, and if a 
  # case is also selected as a control to other cases, we make sure we keep only 
  # the row where a case is selected as a case (s.t. if time is not coarsened, 
  # then t = t_mathc):
  t_match_name_symbol <- as.symbol(t_match_name)
  ncc_cases <- unique(ncc[ncc[, y_name] == 1, ]) %>% 
    arrange(desc({{t_match_name_symbol}}))
  row_id_cases <- which(!duplicated(ncc_cases[, id_name]))
  ncc_cases <- ncc_cases[row_id_cases, ]
  risk_table <- match_risk_table(ncc_cases = ncc_cases, 
                                 risk_table_manual = risk_table_manual,
                                 t_coarse_name = t_match_name, t_name = t_name, 
                                 match_var_names = match_var_names)
  if (anyNA(risk_table)) {
    stop(simpleError(sprintf(
      "There is mismatch in event time in ncc ('%s') and risk_table_manual$t_event.", 
      t_match_name
    )))
  }
  risk_table <- risk_table %>% 
    mutate(p = p_not_sampled(n_event = n_event, n_at_risk = n_at_risk, 
                             n_per_case = n_per_case))
  ncc_cases <- ncc_cases[, -which(names(ncc_cases) == t_match_name)] 
  if (sum(ncc[, y_name] != 1) == 0) {
    dat <- ncc_cases %>% mutate(.km_prob = 1, .km_weight = 1 / .km_prob) %>% 
      change_km_names(km_names = km_names) %>%
      as.data.frame(stringsAsFactors = FALSE)
  } else {
    id_cases <- unique(ncc_cases[, id_name])
    id_controls <- setdiff(ncc[ncc[, y_name] != 1, id_name], id_cases)
    # For each control, only take the row where it first appeared in the NCC
    row_id_controls <- which(ncc[, id_name] %in% id_controls)
    ids <- ncc[row_id_controls, id_name]
    row_id_controls <- row_id_controls[which(!duplicated(ids))]
    ncc_controls <- unique(ncc[row_id_controls, -which(names(ncc) == t_match_name)])
    dat <- rbind(ncc_cases, ncc_controls)
    dat <- dat[sort.list(dat[, id_name]), ] %>% 
      prepare_cohort(cohort = ., t_start_name = t_start_name, t_name = t_name, 
                     y_name = y_name, match_var_names = match_var_names, 
                     print_message = FALSE) %>% 
      assign_kmw0(ncc_nodup = ., risk_table = risk_table) %>% 
      change_km_names(km_names = km_names) %>%
      as.data.frame(stringsAsFactors = FALSE)
  }
  message(simpleMessage(sprintf("There are %d unique subjects (identified by %s) in the input ncc with %d rows, therefore the return data only has %d rows.\n", 
                                length(unique(ncc[, id_name])), id_name, nrow(ncc), nrow(dat))))
  if (return_risk_table) {
    list(dat = dat, 
         risk_table = risk_table[, c("t_event", "n_event", match_var_names, "n_at_risk")])
  } else {
    dat
  }
}
#' Prepare a template risk table given a nested case-control sample
#' @description Given a nested case-control (NCC) data, create a template for
#'   users to fill in the number of subjects at risk in the cohort at each event
#'   time within each sampling stratum (if applicable).
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
#' risk_table <- prepare_risk_table(ncc = ncc_2, t_match_name = "Time", y_name = "Fail", 
#'                                  match_var_names = c("age_cat", "gender"))
prepare_risk_table <- function(ncc, t_match_name, y_name, match_var_names = NULL, 
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
  t_name_symbol <- as.symbol(t_match_name)
  risk_table <- ncc %>% 
    filter({{y_name_symbol}} == 1) %>% 
    group_by(across({{vars}})) %>% 
    summarise(n_at_risk = n()) %>% 
    rename(t_event = {{t_name_symbol}})
  risk_table$n_at_risk <- NA
  if (is.null(csv_file)) {
    as.data.frame(risk_table)
  } else {
    write.csv(risk_table, file = csv_file, na = "", row.names = FALSE)
    message(paste0("Cohort summary template is saved as '", csv_file, "'."))
    NULL
  }
}
#' Compute number at risk at each event time from a "skeleton cohort"
#' @inheritParams prepare_cohort
#' @param t_start_name Name of the variable in \code{cohort} for the start time
#'   of follow-up. A \code{string}. Default is \code{NULL}, i.e., every subject
#'   started the follow-up at time 0.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing the NCC. A \code{string} vector. Default is \code{NULL}, i.e.,
#'   the NCC was only time-matched.
#' @return Returns a data.frame with the following columns:
#' \describe{
#'   \item{t_event}{Unique event times in \code{cohort_skeleton}.}
#'   \item{n_event}{Number of events at each event time.}
#'   \item{n_at_risk}{Number of subjects at risk at each event time.}
#'   \item{<each matching variable>}{If the NCC is matched on additional
#' confounders, each matching variable will be included as a column to the right
#' of \code{strata}.}
#' }
#' @export
compute_risk_table <- function(cohort, t_start_name = NULL, t_name, y_name, 
                               match_var_names = NULL) {
  cohort <- prepare_cohort(cohort = cohort, t_start_name = t_start_name, 
                           t_name = t_name, y_name = y_name, 
                           match_var_names = match_var_names)
  compute_risk_tb(cohort = cohort, match_var_names = match_var_names, 
                  staggered = !is.null(t_start_name)) %>% 
    select(-strata) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
#' <private function> Change variable names corresponding to KM-type probability
#' and weight
#' @param dat Cohort or NCC data generated from \code{\link{compute_km_weights}}
#'   or \code{\link{compute_km_weights_controls}}.
#' @param km_names Names for KM-type probability and weight, corresponding to the 
#' last two columns in \code{dat}.
#' @return Returns \code{dat} with names changed for the last two columns.
change_km_names <- function(dat, km_names) {
  nc <- ncol(dat)
  names(dat)[c(nc - 1, nc)] <- km_names
  dat
}
#' Compute Kaplan-Meier type weights for (matched) nested case-control (NCC)
#' sample
#' @inheritParams compute_kmw_cohort
#' @param ncc (Matched) NCC data, if \code{cohort} is not available.
#'   \strong{This data should not include the ID of each matched set, but should
#'   include the actual event/censoring time of each subject.} A
#'   \code{data.frame} or a matrix with column names.
#' @param risk_table_manual Number of subjects at risk at time of each cases in
#'   the NCC, if \code{cohort} is not available. A \code{data.frame} or a matrix
#'   with column names. See Details.
#' @param t_start_name Name of the variable in \code{cohort} or \code{ncc} for
#'   the start time of follow-up. A \code{string}. Default is \code{NULL}, i.e.,
#'   every subject started the follow-up at time 0.
#' @param t_name Name of the variable in \code{cohort} or \code{ncc} for the
#'   time of event or censoring. A \code{string}. Note that if \code{ncc} is
#'   supplied, in order to correctly compute the weight for each sampled control
#'   this should be the actual time of censoring, not the time of event of the
#'   case in the same matched set.
#' @param t_match_name Name of the column of event time in each matched set in
#'   \code{ncc}, possibly coarsened to the same level as \code{t_event} in
#'   \code{risk_table_manual}. A \code{string}. Default is \code{t_name}, i.e., 
#'   not coarsened.
#' @param y_name Name of the column of censoring status in \code{cohort} or
#'   \code{ncc}, with 1 for event and 0 for censoring. A \code{string}.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} or
#'   \code{ncc} used when drawing the NCC. A \code{string} vector. Default is
#'   \code{NULL}, i.e., the NCC was only time-matched. The corresponding column
#'   in \code{risk_table_manual} must have the same name.
#' @param id_name Name of the column indicating subject ID in \code{ncc}, if
#'   \code{cohort} is not available.
#' @param km_names Column names for the KM-type probability (the first element)
#'   and weight (the second element) computed, if these two columns are to be 
#'   attached to each subject in the input data. Default is
#'   \code{c("km_prob", "km_weight")}.
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
#' @seealso \code{\link{compute_risk_table}}
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
#'                                 match_var_names = c("age_cat", "gender"), 
#'                                 sample_stat = sample_stat, n_per_case = 5)
#' head(ncc_nodup)
#' # Alternatively, if the cohort is not available, the weights can be computed 
#' # as long as number of subjects at risk at event times in each strata is 
#' # available elsewhere, and the actual time of event/censoring is available 
#' # for each subject in the NCC.
#' # Compute the number of subjects at risk from cohort_2:
#' risk_table <- compute_risk_table(cohort = cohort_2, t_name = "t", y_name = "y",
#'                                  match_var_names = c("age_cat", "gender"))
#' head(risk_table)
#' # he following command computes the same weights as in ncc_nodup:
#' ncc_nodup_v2 <- compute_km_weights(ncc = ncc_2[, -1], risk_table_manual = risk_table, 
#'                                    id_name = "Map", t_match_name = "Time", 
#'                                    t_name = "t", y_name = "Fail", 
#'                                    match_var_names = c("age_cat", "gender"),
#'                                    n_per_case = 5)
#' head(ncc_nodup_v2)
compute_km_weights <- function(cohort = NULL, ncc = NULL, id_name = NULL, 
                               risk_table_manual = NULL, 
                               t_start_name = NULL, t_name = NULL, 
                               sample_stat = NULL, 
                               t_match_name = t_name, y_name = NULL, 
                               match_var_names = NULL, n_per_case, 
                               return_risk_table = FALSE, 
                               km_names = c("km_prob", "km_weight")) {
  if (!is.null(cohort)) {
    if (is.null(sample_stat)) {
      stop(simpleError("sample_stat is required if cohort is available."))
    }
    sample_stat <- as.vector(sample_stat)
    # Full cohort is available
    compute_kmw_cohort(cohort = cohort, 
                       t_start_name = t_start_name, t_name = t_name, 
                       sample_stat = sample_stat, 
                       match_var_names = match_var_names, 
                       n_per_case = n_per_case, 
                       return_risk_table = return_risk_table, 
                       km_names = km_names)
  } else {
    # Full cohort is not available
    if (is.null(risk_table_manual)) {
      stop(simpleError("risk_table_manual is required if cohort is not available."))
    }
    compute_kmw_ncc(ncc = ncc, id_name = id_name, 
                    risk_table_manual = risk_table_manual, 
                    t_start_name = t_start_name, t_name = t_name, 
                    t_match_name = t_match_name, y_name = y_name, 
                    match_var_names = match_var_names, 
                    n_per_case = n_per_case, 
                    return_risk_table = return_risk_table, 
                    km_names = km_names)
  }
}
#' Compute Kaplan-Meier type weights for newly collected NCC controls
#' @param ncc_controls Newly collected NCC controls, where each row corresponds
#'   to a unique subject. Make sure this dataset does not include any subject
#'   that later became cases. This data should include the actual
#'   event/censoring time of each subject. A \code{data.frame} or a matrix with
#'   column names.
#' @param risk_table_manual Number of subjects at risk at time of each cases in
#'   the NCC, prepared using function \code{\link{match_risk_table}}.
#' @param t_start_name Name of the variable in \code{ncc_controls} for the start
#'   time of follow-up. A \code{string}. Default is \code{NULL}, i.e., every
#'   subject started the follow-up at time 0.
#' @param t_name Name of the variable in \code{ncc_controls} for the time of
#'   event or censoring. A \code{string}.
#' @param match_var_names Name(s) of the match variable(s) in
#'   \code{ncc_controls} used when drawing the NCC. A \code{string} vector.
#'   Default is \code{NULL}, i.e., the NCC was only time-matched. The
#'   corresponding column in \code{prob_table_manual} must have the same name.
#' @param km_names Column names for the KM-type probability (the first element)
#'   and weight (the second element) computed. Default is
#'   \code{c("km_prob", "km_weight")}.
#' @import dplyr
#' @export
#' @seealso \code{\link{match_risk_table}}
compute_km_weights_controls <- function(ncc_controls, risk_table_manual, 
                                        t_start_name = NULL, t_name, 
                                        match_var_names = NULL, n_per_case,
                                        km_names = c("km_prob", "km_weight")) {
  ncc_controls <- as.data.frame(ncc_controls, stringsAsFactors = FALSE)
  ncc_controls$.y <- 0
  risk_table <- risk_table_manual %>% 
    mutate(p = p_not_sampled(n_event = n_event, n_at_risk = n_at_risk, 
                             n_per_case = n_per_case))
  ncc_controls %>% 
    prepare_cohort(cohort = ., t_start_name = t_start_name, t_name = t_name, 
                   y_name = ".y", match_var_names = match_var_names) %>%
    assign_kmw0(ncc_nodup = ., risk_table = risk_table) %>% 
    change_km_names(km_names = km_names) %>%
    as.data.frame(stringsAsFactors = FALSE)
}