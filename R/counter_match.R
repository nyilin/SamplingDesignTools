#' Draw counter-matched NCC from cohort
#' @param cohort Cohort data. A \code{data.frame} or a matrix with column names.
#' @param y_name Name of the variable in \code{cohort} indicating the
#'   event-censoring status of each subject, where 1 indicates events. A
#'   \code{string}. 
#' @param t_name Name of the variable in \code{cohort} for the time of event or
#'   censoring. A \code{string}. 
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
