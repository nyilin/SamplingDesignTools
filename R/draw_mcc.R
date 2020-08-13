#' Draw a (matched) case-control (CC) sample from a cohort (or cross-sectional data)
#' @author Yilin Ning, Anastasia Lam, Marie Reilly
#' @param cohort Cohort (or cross-sectional) data. A \code{data.frame} or a
#'   matrix with column names.
#' @param y_name Name of the outcome variable in \code{cohort}. Cases are
#'   required to be indicated by integer \code{1}. It is recommended to code the 
#'   outcome variable as an integer vector, where 1=case and 0=non-case.
#' @param n_per_case Number of controls to draw for each case.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing a matched CC. A \code{string} vector. Default is \code{NULL},
#'   i.e., the controls will be sampled randomly.
#' @param weight_name Name to assign to the column for the inverse probability
#'   weight, which will be created as the last column in the (matched) CC
#'   data returned. Default is \code{"ip_weight"}.
#' @return Returns a \code{data.frame} of the (matched) CC sample, with an 
#'   additional column for the inverse probability weight for each subject, i.e., 
#'   1 for cases and the inverse probability of the sampling fraction for each 
#'   control.
#' @import dplyr
#' @export
#' @examples 
#' # Load cohort data
#' data(cohort_2) # Ignoring that we have the event time in this simulated cohort
#' head(cohort_2)
#' # Draw simple 1:2 case-control sample (i.e., with randomly selected controls):
#' dat_scc <- draw_mcc(cohort = cohort_2, y_name = "y", n_per_case = 2)
#' head(dat_scc)
#' # Draw simple 1:2 case-control sample, matched on age group and gender:
#' dat_mcc <- draw_mcc(cohort = cohort_2, y_name = "y", n_per_case = 2,
#'                     match_var_names = c("age_cat", "gender"),
#'                     weight_name = ".w")
#' head(dat_mcc)
draw_mcc <- function(cohort, y_name, n_per_case, match_var_names = NULL, 
                     weight_name = "ip_weight") {
  cohort <- as.data.frame(cohort)
  y_name <- y_name[1]
  if (!(y_name %in% names(cohort))) {
    stop(simpleError(paste(y_name, "not found in cohort.")))
  }
  if (sum(cohort[, y_name] == 1) == 0) {
    stop(simpleError(sprintf("There are no cases (i.e., %s=1) in the input cohort.", 
                             y_name)))
  }
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(cohort))) {
      stop(simpleError("Make sure all match variables are in cohort."))
    }
  }
  n_per_case <- as.integer(n_per_case)
  if (n_per_case <= 0) stop(simpleError("n_per_case should be a positive integer."))
  message(simpleMessage(sprintf(
    "%d of the %d subjects in the input cohort are cases (i.e., with %s=1).\n", 
    sum(cohort[, y_name] == 1), nrow(cohort), y_name
  )))
  message(simpleMessage(sprintf("%d controls are sampled per case.\n", n_per_case)))
  if (is.null(match_var_names)) {
    n_cases <- sum(cohort[, y_name] == 1)
    n_sampled <- n_cases * n_per_case
    if (n_sampled >= sum(cohort[, y_name] != 1)) {
      warnings(simpleWarning(
        "Not enough non-cases to sample from. All subjects in the input cohort are returned."
      ))
      dat_mcc <- cbind(cohort, .w = 1)
    }
    i_controls <- sample(x = which(cohort[, y_name] != 1), size = n_sampled, replace = FALSE)
    dat_mcc <- rbind(cohort[cohort[, y_name] == 1, ], cohort[i_controls, ])
    # Sampling weights:
    dat_mcc$.w <- ifelse(dat_mcc[, y_name] == 1, 1, 
                         sum(cohort[, y_name] != 1) / sum(dat_mcc[, y_name] != 1))
  } else {
    y_name_symbol <- as.symbol(y_name)
    vars1 <- c(match_var_names, ".u")
    vars2 <- c(match_var_names, y_name)
    dat_mcc <- cohort %>% 
      mutate(.u = runif(n = n())) %>% 
      arrange(across({{vars1}})) %>%
      group_by(across({{match_var_names}})) %>% 
      mutate(.n_cases = sum({{y_name_symbol}} == 1)) %>%
      ungroup() %>%
      group_by(across({{vars2}})) %>%
      mutate(.tmp_id = 1:n(), 
             .select = ifelse(.tmp_id <= n_per_case * .n_cases, 1, 0), 
             .w = n() / sum(.select == 1)) %>%
      filter(.select == 1 | {{y_name_symbol}} == 1)
    check_n_control <- dat_mcc %>%
      summarise(enough_control = max(.tmp_id) < n_per_case * unique(.n_cases))
    dat_mcc <- dat_mcc %>%
      select(-.tmp_id, -.select, -.u, -.n_cases)
    dat_mcc$.w[dat_mcc[, y_name] == 1] <- 1
    if (any(check_n_control$enough_control)) {
      warnings(simpleWarning(
        "Not enough non-cases to sample from in at least one matched set. All subjects in the input cohort are returned for those sets."
      ))
    }
  }
  names(dat_mcc)[which(names(dat_mcc) == ".w")] <- weight_name
  dat_mcc %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
