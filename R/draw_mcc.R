#' Private function: draw simple case-control data from cohort, taking all 
#' cases available
draw_scc0 <- function(cohort, y_name, n_per_case) {
  n_cases <- sum(cohort[, y_name] == 1)
  n_controls <- n_cases * n_per_case
  if (n_controls >= sum(cohort[, y_name] != 1)) {
    warnings(simpleWarning(
      "Not enough non-cases to sample from. All non-cases are selected."
    ))
    dat <- cbind(cohort, .w = 1)
  } else {
    i_controls <- sample(x = which(cohort[, y_name] != 1), size = n_controls, 
                         replace = FALSE)
    dat <- rbind(cohort[cohort[, y_name] == 1, ], cohort[i_controls, ])
    # Sampling weights:
    dat$.w <- ifelse(dat[, y_name] == 1, 1, sum(cohort[, y_name] != 1) / n_controls)
  }
  dat
}
#' Private function: draw matched case-control data from cohort, taking all 
#' cases available
#' @import dplyr
draw_mcc0 <- function(cohort, y_name, n_per_case, match_var_names) {
  y_name_symbol <- as.symbol(y_name)
  vars1 <- c(match_var_names, ".u")
  vars2 <- c(match_var_names, y_name)
  dat <- cohort %>% 
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
  check_n_control <- dat %>%
    summarise(enough_control = max(.tmp_id) < n_per_case * unique(.n_cases))
  dat <- dat %>%
    mutate(.w = ifelse({{y_name_symbol}} == 1, 1, .w)) %>%
    select(-.tmp_id, -.select, -.u, -.n_cases) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  if (any(check_n_control$enough_control)) {
    warnings(simpleWarning(
      "Not enough non-cases to sample from in at least one matched set. All non-cases in those sets are selected."
    ))
  }
  dat
}
#' Draw a (matched) case-control (CC) sample from a cohort (or cross-sectional data)
#' @author Yilin Ning, Anastasia Lam, Marie Reilly
#' @param cohort Cohort (or cross-sectional) data. A \code{data.frame} or a
#'   matrix with column names.
#' @param y_name Name of the outcome variable in \code{cohort}. Cases are
#'   required to be indicated by integer \code{1}. It is recommended to code the 
#'   outcome variable as an integer vector, where 1=case and 0=non-case.
#' @param n_cases Number of cases to draw. Default is \code{NULL}, in which case 
#'   all cases will be selected.
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
draw_mcc <- function(cohort, y_name, n_cases = NULL, n_per_case, 
                     match_var_names = NULL, weight_name = "ip_weight") {
  cohort <- as.data.frame(cohort)
  y_name <- y_name[1]
  # Check input -----
  # Outcome
  if (!(y_name %in% names(cohort))) {
    stop(simpleError(paste(y_name, "not found in cohort.")))
  }
  # Matching variables
  if (!is.null(match_var_names)) {
    if (!all(match_var_names %in% names(cohort))) {
      stop(simpleError("Make sure all match variables are in cohort."))
    }
  }
  # Cases
  if (sum(cohort[, y_name] == 1) == 0) {
    stop(simpleError(sprintf("There are no cases (i.e., %s=1) in the input cohort.", 
                             y_name)))
  } else {
    n_cases0 <- sum(cohort[, y_name] == 1) # Number of cases available
    message(simpleMessage(sprintf(
      "%d of the %d subjects in the input cohort are cases (i.e., with %s=1). ", 
      n_cases0, nrow(cohort), y_name
    )))
  }
  # Controls
  n_per_case <- as.integer(n_per_case)
  if (is.na(n_per_case) || n_per_case <= 0) {
    stop(simpleError("n_per_case should be a positive integer."))
  }
  # Draw sample -----
  if (is.null(n_cases)) {
    # By default, select all cases in cohort
    message(simpleMessage(sprintf("All cases are sampled.\n")))
    n_cases <- n_cases0
  } else {
    # If specified by user, select the number of cases as requested
    # Remove unwanted cases and reuse the function that draws all cases
    n_cases <- as.integer(n_cases)
    if (is.na(n_cases) || n_cases <= 0) {
      stop(simpleError("n_cases should be a positive integer."))
    }
    if (n_cases > n_cases0) {
      n_cases <- n_cases0
      message(simpleMessage(sprintf("Number of cases requested is larger than number of cases available. All cases are sampled.\n")))
    } else {
      i_cases_remove <- sample(x = which(cohort[, y_name] == 1), 
                               size = n_cases0 - n_cases, replace = FALSE)
      message(simpleMessage(
        sprintf("%d cases are sampled, as requested.\n", n_cases)
      ))
      cohort <- cohort[-i_cases_remove, ]
    }
  }
  message(simpleMessage(sprintf("%d controls are sampled per case.\n", n_per_case)))
  if (is.null(match_var_names)) {
    dat <- draw_scc0(cohort = cohort, y_name = y_name, n_per_case = n_per_case)
  } else {
    dat <- draw_mcc0(cohort = cohort, y_name = y_name, n_per_case = n_per_case, 
                     match_var_names = match_var_names)
  }
  # Correct weight for cases: no longer 1 if not all cases are sampled
  dat$.w[dat[, y_name] == 1] <- n_cases0 / n_cases
  # Last column is sampling weight. Change the name of this column
  nc <- ncol(dat)
  names(dat)[nc] <- weight_name
  dat
}
