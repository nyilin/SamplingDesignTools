#' Recode subjects as cases and potential controls for (M)ECC sampling
#' @param t Event/censoring time.
#' @param delta Event/censoring status.
#' @inheritParams draw_mecc_matched
recode_extreme <- function(t, delta, tau0, tau) {
  y <- rep(2, length(delta))
  y[t <= tau0 & delta == 1] <- 1
  y[t > tau] <- 0
  y
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
#' @param match_var_names Confounder(s) to match on.
#' @param n_per_case Number of controls to match to each case.
#' @import survival
#' @import dplyr
#' @export
draw_mecc_matched <- function(cohort, tau0, tau, id_name, t_name, delta_name, 
                              match_var_names, n_per_case) {
  cohort <- prepare_cohort(cohort = cohort, t_start_name = NULL, t_name = t_name, 
                           y_name = delta_name, match_var_names = match_var_names, 
                           print_message = FALSE) %>% 
    rename(.delta = .y)
  match_var <- cohort$.strata0
  km <- survfit(Surv(.t, .delta) ~ match_var, data = cohort)
  km_summ <- summary(km)
  surv_df <- data.frame(.t = km_summ$time, surv = km_summ$surv, 
                        .strata = as.character(km_summ$strata), 
                        stringsAsFactors = FALSE)
  surv_tau <- surv_df %>% 
    group_by(.strata) %>% 
    summarise(surv_tau = approx(x = .t, y = surv, xout = tau)$y)
  cohort <- cohort %>% 
    mutate(.y = recode_extreme(t = .t, delta = .delta, tau0 = tau0, tau = tau)) %>% 
    filter(.y %in% c(0, 1)) %>%
    left_join(surv_df) %>%
    left_join(surv_tau) %>% 
    arrange(.strata, .t)
  i_na <- which(is.na(cohort$surv))
  surv_na <- unlist(lapply(i_na, function(i) {
    t_i <- cohort$.t[i]
    t_strata <- surv_df$.t[surv_df$.strata == cohort$.strata[i]]
    if (t_i < min(t_strata)) {
      1
    } else {
      surv_strata <- surv_df$surv[surv_df$.strata == cohort$.strata[i]]
      t_max <- max(t_strata[t_strata < t_i])
      surv_strata[which(t_strata == t_max)]
    }
  }))
  cohort$surv[i_na] <- surv_na
  id <- cohort[, id_name]
  do.call("rbind", lapply(cohort[cohort$.y == 1, id_name], function(id_case) {
    id_controls <- sample(x = id[cohort$.y == 0 & 
                                   cohort$.strata == cohort$.strata[id == id_case]], 
                          size = n_per_case, replace = FALSE)
    cohort[id %in% c(id_case, id_controls), ] %>% 
      mutate(.set_id = paste0("set_", id_case))
  })) %>% 
    select(-.strata0, -.strata, -.t_start, -.t, -.delta) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}
#' Negative log-likelihood in weighted analysis of matched (M)ECC sample
#' @param beta A vector of regression coefficients corresponding to
#'   \code{x_formula}.
#' @param x_mat A \code{model.matrix} (without the first column that has value 
#'   1 for all subjects).
#' @param y_name A string indicating cases and controls in \code{mecc}. Note
#'   that this is not the original indicator for event/censoring.
#' @param set_id_name A string indicating the ID of matched sets in the (M)ECC
#'   data.
#' @param mecc (M)ECC data. A \code{data.frame}. Make sure the covariates in
#'   \code{x_formula} are all centred at cohort average.
#' @param surv Estimated baseline survival probability for each subject in
#'   \code{mecc}, based on the underlying cohort. The length of this variable 
#'   must be the same as the number of rows in \code{mecc}.
#' @param surv_tau Estimated baseline survival probability evaluate at time
#'   \code{tau}, based on the underlying cohort. The length of this variable 
#'   must be the same as the number of rows in \code{mecc}.
#' @import dplyr
llh_matched <- function(beta, x_mat, y_name, set_id_name, surv, surv_tau, mecc) {
  if (y_name != ".y") {
    y_name_symbol <- as.symbol(y_name)
    mecc <- mecc %>% rename(.y = {{y_name_symbol}})
  }
  if (set_id_name != ".set_id") {
    set_id_name_symbol <- as.symbol(set_id_name)
    mecc <- mecc %>% rename(.set_id = {{set_id_name_symbol}})
  }
  lh <- mecc %>% 
    mutate(lp = as.matrix(x_mat) %*% beta, .surv = surv, .surv_tau = surv_tau) %>%
    group_by(.set_id) %>%
    mutate(w = (.surv[.y == 1] / .surv_tau[.y == 1]) ^ exp(lp)) %>% 
    summarise(lh_i = (w[.y == 1] * exp(lp[.y == 1])) / sum(w * exp(lp)))
  - sum(log(lh$lh_i))
}
#' Perform weighted analysis of matched (M)ECC data
#' @inheritParams llh_matched
#' @param x_formula A \code{formula} object specifying the model assumed for the
#'   covariates, starting with \code{~} (but does not include the outcome).
#' @import dplyr
#' @export
analyse_mecc_matched <- function(y_name, x_formula, set_id_name, surv, surv_tau, 
                                 mecc) {
  x_mat <- model.matrix(x_formula, data = mecc)
  x_names <- colnames(x_mat)[-1]
  if (ncol(x_mat) == 2) {
    x_mat <- matrix(x_mat[, -1], ncol = 1)
  } else {
    x_mat <- x_mat[, -1]
  }
  opt_obj <- optim(par = numeric(ncol(x_mat)), fn = llh_matched, 
                   x_mat = x_mat, y_name = y_name, set_id_name = set_id_name, 
                   mecc = mecc, surv = surv, surv_tau = surv_tau, hessian = TRUE)
  coef_mecc <- data.frame(var = x_names, est = opt_obj$par, 
                          exp_est = exp(opt_obj$par),
                          se = 1 / sqrt(diag(opt_obj$hessian))) %>% 
    mutate(pval = 2 * pnorm(q = abs(est / se), lower.tail = FALSE)) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  rownames(coef_mecc) <- x_names
  list(coef_mat = coef_mecc, optim_obj = opt_obj)
}
