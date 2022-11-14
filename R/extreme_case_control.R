#' <Private function> Recode subjects as cases and potential controls for (M)ECC
#' sampling
#' @param t Event/censoring time.
#' @param delta Event/censoring status.
#' @inheritParams draw_mecc
recode_extreme <- function(t, delta, tau0, tau) {
  y <- rep(2, length(delta))
  y[t <= tau0 & delta == 1] <- 1
  y[t > tau] <- 0
  y
}
#' <private function> Change variable names in MECC sample
#' @param dat MECC data generated from \code{\link{draw_mecc}}.
#' @param mecc_names Names to assign to the following four variables (in
#'   sequence): the indicator for extreme events, the estimated survival
#'   probability at each time point, the estimated survival probability at
#'   \code{tau}, and the ID for matched sets (see Details). These four variables 
#'   correspond to the last four columns of the MECC sample.
#' @return Returns \code{dat} with names changed for the outcome and set id
#'   columns.
change_mecc_names <- function(dat, mecc_names) {
  nc <- ncol(dat)
  names(dat)[(nc - 3):nc] <- mecc_names
  dat
}
#' Draw matched (more) extreme case-control samples
#' @param cohort Cohort data.
#' @param tau0 Only subjects who had event before \code{tau0} are considered as
#'   cases.
#' @param tau Only subjects who did not have event at or before \code{tau} are
#'   eligible to be selected as controls. Extreme case-control (ECC) has 
#'   \code{tau = tau0}, and more extreme case-control (MECC) has 
#'   \code{tau > tau0}.
#' @param id_name Name of subject ID.
#' @param t_name Name of time variable.
#' @param delta_name Name of event/censoring indicator.
#' @param match_var_names Name(s) of the match variable(s) in \code{cohort} used
#'   when drawing a matched MECC sample. A \code{string} vector. Default is
#'   \code{NULL}, i.e., the controls will be sampled randomly.
#' @param n_per_case Number of controls to match to each case.
#' @param mecc_names Names to assign to the following four variables (in
#'   sequence): the indicator for extreme events, the estimated survival
#'   probability at each time point, the estimated survival probability at
#'   \code{tau}, and the ID for matched sets (see Details). Default value is
#'   \code{c("y_mecc", "surv", "surv_tau", "set_id_mecc")}.
#' @details This function draws a (more) extreme case-control sample from a
#'   cohort, with or without matching on confounder(s), with \code{n_per_case}
#'   controls drawn for each case. Eligible cases are subjects who had the event
#'   before time \code{tau0} (and all of them are sampled), and eligible
#'   controls are subjects who do not have the event until time \code{tau}. This
#'   function throws an error if there is not enough controls to select from,
#'   otherwise it uses the \code{\link{draw_mcc}} function to draw a
#'   case-control sample. 
#'   
#'   This package implements the conditional approach described in Section 2.2
#'   of Salim et al 2014, where cases and controls are "individually matched".
#'   Therefore, although controls are selected randomly or by frequency
#'   matching, \code{n_per_case} controls would be randomly assigned to a case
#'   to form a "matched set" after drawing the sample.
#' @return Returns a \code{data.frame} of the (more) extreme case-control
#'   sample, with four additional columns as described in parameter
#'   \code{mecc_names}.
#' @importFrom survival survfit Surv
#' @importFrom dplyr rename group_by summarise mutate filter left_join arrange desc
#' @importFrom rlang .data
#' @importFrom stats approx
#' @export
#' @example man/examples/draw_mecc.R
#' @references
#' \itemize{
#'   \item{Salim A, Ma X, Fall K, et al. Analysis of incidence and prognosis 
#'   from 'extreme' case–control designs. Stat Med 2014; 33: 5388–5398.}
#'   \item{Støer NC, Salim A, Bokenberger K, et al. Is the matched extreme
#'   case–control design more powerful than the nested case–control design?.
#'   Stat Methods Med Res 2019; 28(6): 1911-1923.}
#' }
#' @author Yilin Ning, Nathalie C Støer
#' @seealso \code{\link{llh_mecc_cond}}, \code{\link{analyse_mecc_cond}}
draw_mecc <- function(cohort, tau0, tau, id_name, t_name, delta_name, 
                      match_var_names = NULL, n_per_case, 
                      mecc_names = c("y_mecc", "surv", "surv_tau", "set_id_mecc")) {
  cohort <- prepare_cohort(cohort = cohort, t_start_name = NULL, t_name = t_name, 
                           y_name = delta_name, match_var_names = match_var_names, 
                           print_message = FALSE) %>% 
    dplyr::rename(.delta = .data$.y)
  match_var <- cohort$.strata0
  km <- survfit(as.formula("Surv(.t, .delta) ~ match_var"), data = cohort)
  km_summ <- summary(km)
  if (is.null(km_summ$strata)) {
    surv_df <- data.frame(.t = km_summ$time, .surv = km_summ$surv, 
                          .strata = "match_var=1", 
                          stringsAsFactors = FALSE)
    # Estimated survival probabilities that will be used to compute weights for 
    # controls in unmatched MECC:
  } else {
    surv_df <- data.frame(.t = km_summ$time, .surv = km_summ$surv, 
                          .strata = as.character(km_summ$strata), 
                          stringsAsFactors = FALSE)
  }
  # Estimated survival probabilities that will be used to compute weights for 
  # controls in matched MECC:
  surv_tau <- surv_df %>% 
    dplyr::group_by(.data$.strata) %>% 
    dplyr::summarise(.surv_tau = approx(x = .data$.t, y = .data$.surv, 
                                        xout = tau)$y)
  cohort <- cohort %>% 
    dplyr::mutate(.y = recode_extreme(t = .data$.t, delta = .data$.delta, 
                                      tau0 = tau0, tau = tau)) %>% 
    dplyr::filter(.data$.y %in% c(0, 1)) %>%
    dplyr::left_join(surv_df) %>%
    dplyr::left_join(surv_tau) %>% 
    dplyr::arrange(.data$.strata, .data$.t)
  check_controls <- cohort %>% 
    dplyr::group_by(.data$.strata) %>% 
    dplyr::summarise(n_cases = sum(.data$.y == 1), 
                     n_controls = sum(.data$.y == 0), 
                     enough = .data$n_controls >= .data$n_cases * n_per_case)
  if (any(!check_controls$enough)) {
    stop(simpleError("Not enough control to sample from."))
  }
  i_na <- which(is.na(cohort$.surv))
  surv_na <- unlist(lapply(i_na, function(i) {
    t_i <- cohort$.t[i]
    t_strata <- surv_df$.t[surv_df$.strata == cohort$.strata[i]]
    if (t_i < min(t_strata)) {
      1
    } else {
      surv_strata <- surv_df$.surv[surv_df$.strata == cohort$.strata[i]]
      t_max <- max(t_strata[t_strata < t_i])
      surv_strata[which(t_strata == t_max)]
    }
  }))
  cohort$.surv[i_na] <- surv_na
  mecc <- suppressMessages(draw_mcc(cohort = cohort, y_name = ".y", 
                                    n_per_case = n_per_case,
                                    match_var_names = ".strata", 
                                    weight_name = ".w_unused")) %>%
    dplyr::arrange(dplyr::desc(.data$.y), .data$.strata) %>%
    dplyr::select(-.data$.w_unused, -.data$.strata0, -.data$.strata, 
                  -.data$.t_start, -.data$.t, -.data$.delta) %>%
    as.data.frame(stringsAsFactors = FALSE)
  set_id <- paste0("set_", mecc[mecc$.y == 1, id_name])
  mecc$.set_id <- c(set_id, rep(set_id, each = n_per_case))
  change_mecc_names(dat = mecc, mecc_names = mecc_names)
}
#' <Private function> Negative log-likelihood in weighted analysis of (M)ECC sample
#' @param beta A vector of regression coefficients corresponding to
#'   \code{x_formula}.
#' @param x_mat A \code{model.matrix} (without the first column that has value 
#'   1 for all subjects).
#' @param y_name A string indicating cases and controls in \code{mecc}. Note
#'   that this is not the original indicator for event/censoring.
#' @param set_id_name A string indicating the ID of matched sets in the (M)ECC
#'   data. See Details.
#' @param mecc (M)ECC data. A \code{data.frame}. Make sure the covariates in
#'   \code{x_formula} are all centred at cohort average.
#' @param surv Estimated baseline survival probability for each subject in
#'   \code{mecc}, based on the underlying cohort. The length of this variable 
#'   must be the same as the number of rows in \code{mecc}. See Details.
#' @param surv_tau Estimated baseline survival probability evaluate at time
#'   \code{tau}, based on the underlying cohort. The length of this variable 
#'   must be the same as the number of rows in \code{mecc}. See Details.
#' @details This function implements the conditional approach described in
#'   Section 2.2.1 of Salim et al 2014. This requires controls to be
#'   "individually matched" to each control. If users used the function
#'   \code{\link{draw_mecc}} to draw the MECC sample, the sample would include a
#'   variable that indicates the matched sets, and two variables corresponding
#'   to \code{surv} and \code{surv_tau}. Otherwise users should randomly match
#'   each case to controls (within strata if the MECC sample is matched on
#'   confounder(s)) and compute them from the full cohort as described in
#'   Section 2.2 of Salim et al 2014.
#' @return Returns the negative log-likelihood in the weighted analysis of a
#'   (more) extreme case-control sample.
#' @importFrom dplyr rename mutate group_by summarise
#' @importFrom rlang .data
#' @author Yilin Ning, Nathalie C Støer
#' @references
#' \itemize{
#'   \item{Salim A, Ma X, Fall K, et al. Analysis of incidence and prognosis 
#'   from 'extreme' case–control designs. Stat Med 2014; 33: 5388–5398.}
#'   \item{Støer NC, Salim A, Bokenberger K, et al. Is the matched extreme
#'   case–control design more powerful than the nested case–control design?.
#'   Stat Methods Med Res 2019; 28(6): 1911-1923.}
#' }
#' @seealso \code{\link{draw_mecc}}, \code{\link{analyse_mecc_cond}}
llh_mecc_cond <- function(beta, x_mat, y_name, set_id_name, surv, surv_tau, mecc) {
  if (y_name != ".y") {
    y_name_symbol <- as.symbol(y_name)
    mecc <- mecc %>% dplyr::rename(.y = {{y_name_symbol}})
  }
  if (set_id_name != ".set_id") {
    set_id_name_symbol <- as.symbol(set_id_name)
    mecc <- mecc %>% dplyr::rename(.set_id = {{set_id_name_symbol}})
  }
  lh <- mecc %>% 
    dplyr::mutate(lp = as.matrix(x_mat) %*% beta, 
                  .surv = .data$surv, .surv_tau = .data$surv_tau) %>%
    dplyr::group_by(.data$.set_id) %>%
    dplyr::mutate(w = (.data$.surv[.data$.y == 1] / 
                         .data$.surv_tau[.data$.y == 1]) ^ 
                    exp(.data$lp)) %>% 
    dplyr::summarise(lh_i = (.data$w[.data$.y == 1] * 
                               exp(.data$lp[.data$.y == 1])) / 
                       sum(.data$w * exp(.data$lp)))
  - sum(log(lh$lh_i))
}
#' Perform weighted analysis of (M)ECC data using the conditional approach
#' @inheritParams llh_mecc_cond
#' @param x_formula A \code{formula} object specifying the model assumed for the
#'   covariates, starting with \code{~} (but does not include the outcome
#'   variable).
#' @param lower,upper Parameters of the function \code{\link{optim}}, which
#'   specifies the lower and upper bounds of the estimated coefficient for
#'   method "Brent" if there is only one coefficient.
#' @details This function uses the function \code{\link{optim}} to estimate the
#'   regression coefficients and the Hessian matrix. The default method of
#'   \code{\link{optim}} (i.e., "Nelder-Mead") is used if there are more than 1
#'   coefficient to estimate, otherwise the "Brent" method is used, with the
#'   search range given by \code{lower} and \code{upper}.
#' @return Returns a list consisting of two components:
#' \describe{
#'   \item{coef_mecc}{A \code{data.frame} containing the name of covariates in
#'   the regression model (\code{var}), estimated regression coefficients
#'   (\code{est}), estimated hazard ratios (\code{exp_est = exp(est)}) standard
#'   error of the estimates (\code{se}), and the p-values (\code{pval}).}
#'   \item{optim_obj}{The \code{optim} output from maximising the log-likelihood
#'   in the weighted analysis.}
#' }
#' @importFrom stats model.matrix optim pnorm
#' @export
#' @example man/examples/analyse_mecc_cond.R
#' @references
#' \itemize{
#'   \item{Salim A, Ma X, Fall K, et al. Analysis of incidence and prognosis 
#'   from 'extreme' case–control designs. Stat Med 2014; 33: 5388–5398.}
#'   \item{Støer NC, Salim A, Bokenberger K, et al. Is the matched extreme
#'   case–control design more powerful than the nested case–control design?.
#'   Stat Methods Med Res 2019; 28(6): 1911-1923.}
#' }
#' @author Yilin Ning, Nathalie C Støer
#' @seealso \code{\link{draw_mecc}}, \code{\link{llh_mecc_cond}}
analyse_mecc_cond <- function(y_name, x_formula, set_id_name, surv, surv_tau, 
                              mecc, lower = -10, upper = 10) {
  x_mat <- model.matrix(x_formula, data = mecc)
  x_names <- colnames(x_mat)[-1]
  if (ncol(x_mat) == 2) {
    x_mat <- matrix(x_mat[, -1], ncol = 1)
  } else {
    x_mat <- x_mat[, -1]
  }
  if (ncol(x_mat) == 1) {
    opt_obj <- optim(par = numeric(ncol(x_mat)), fn = llh_mecc_cond, 
                     x_mat = x_mat, y_name = y_name, set_id_name = set_id_name, 
                     mecc = mecc, surv = surv, surv_tau = surv_tau, 
                     hessian = TRUE, method = "Brent", lower = lower, upper = upper)
  } else {
    opt_obj <- optim(par = numeric(ncol(x_mat)), fn = llh_mecc_cond, 
                     x_mat = x_mat, y_name = y_name, set_id_name = set_id_name, 
                     mecc = mecc, surv = surv, surv_tau = surv_tau, 
                     hessian = TRUE)
  }
  coef_mecc <- data.frame(var = x_names, est = opt_obj$par, 
                          exp_est = exp(opt_obj$par),
                          se = 1 / sqrt(diag(opt_obj$hessian)), 
                          stringsAsFactors = FALSE)
  coef_mecc$pval <- 2 * pnorm(q = abs(coef_mecc$est / coef_mecc$se), 
                              lower.tail = FALSE)
  rownames(coef_mecc) <- x_names
  list(coef_mat = coef_mecc, optim_obj = opt_obj)
}
