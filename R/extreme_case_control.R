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
#' @details This function draws a matched (more) extreme case-control samples
#'   from a cohort, with \code{n_per_case} controls drawn for each case.
#'   Eligible cases are subjects who had the event before time \code{tau0} (and
#'   all of them are sampled), and eligible controls are subjects who do not
#'   have the event until time \code{tau}. This function throws an error if
#'   there is not enough controls to select from, otherwise it uses the
#'   \code{\link{draw_mcc}} function to draw a matched case-control sample.
#' @return Returns a \code{data.frame} of the (more) extreme case-control
#'   sample, with column \code{.set_id} indicating the matched sets in the
#'   sample, \code{.y} indicating the new case-control indicator (which is
#'   likely to be different from the event indicator in the cohort), and the
#'   estimated baseline survival probabilities (\code{surv} and \code{surv_tau})
#'   that will be used in the weighted analysis of this sample.
#' @import survival
#' @import dplyr
#' @export
#' @examples 
#' # Load cohort data
#' data(cohort_2) 
#' head(cohort_2)
#' # Draw simple 1:2 more extreme case-control sample, matched on gender. 
#' # Let cases be subjects who had the event within 5 years, and controls be 
#' # selected from those who did not have the event until the 15-th year.
#' dat_mecc <- draw_mecc_matched(cohort = cohort_2, tau0 = 5, tau = 15, 
#'                               id_name = "id", t_name = "t", delta_name = "y", 
#'                               match_var_names = "gender", n_per_case = 2)
#' head(dat_mecc)
#' # Note that the new event indicator, `.y`, is different from the original one, 
#' # `y`, in the cohort:
#' identical(dat_mecc$y, dat_mecc$.y)
#' @references
#' \itemize{
#'   \item{Salim A, Ma X, Fall K, et al. Analysis of incidence and prognosis 
#'   from 'extreme' case–control designs. Stat Med 2014; 33: 5388–5398.}
#'   \item{Støer NC, Salim A, Bokenberger K, et al. Is the matched extreme
#'   case–control design more powerful than the nested case–control design?.
#'   Stat Methods Med Res 2019; 28(6): 1911-1923.}
#' }
#' @author Yilin Ning, Nathalie C Støer
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
  check_controls <- cohort %>% group_by(.strata) %>% 
    summarise(n_cases = sum(.y == 1), n_controls = sum(.y == 0), 
              enough = n_controls >= n_cases * n_per_case)
  if (any(!check_controls$enough)) {
    stop(simpleError("Not enough control to sample from."))
  }
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
  mecc <- suppressMessages(draw_mcc(cohort = cohort, y_name = ".y", 
                                    n_per_case = n_per_case,
                                    match_var_names = ".strata", 
                                    weight_name = ".w")) %>%
    arrange(desc(.y), .strata) %>%
    select(-.w, -.strata0, -.strata, -.t_start, -.t, -.delta, -.y) %>%
    as.data.frame(stringsAsFactors = FALSE)
  set_id <- paste0("set_", mecc[mecc$.y == 1, id_name])
  mecc$.set_id <- c(set_id, rep(set_id, each = n_per_case))
  mecc
  # id <- cohort[, id_name]
  # do.call("rbind", lapply(cohort[cohort$.y == 1, id_name], function(id_case) {
  #   id_controls <- sample(x = id[cohort$.y == 0 &
  #                                  cohort$.strata == cohort$.strata[id == id_case]],
  #                         size = n_per_case, replace = FALSE)
  #   cohort[id %in% c(id_case, id_controls), ] %>%
  #     mutate(.set_id = paste0("set_", id_case))
  # })) %>%
  #   select(-.strata0, -.strata, -.t_start, -.t, -.delta) %>%
  #   as.data.frame(stringsAsFactors = FALSE)
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
#' @return The negative log-likelihood in the weighted analysis of a matched 
#'   (more) extreme case-control sample.
#' @import dplyr
#' @author Yilin Ning, Nathalie C Støer
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
#' @return Returns a list consisting of two components:
#' \describe{
#'   \item{coef_mecc}{A \code{data.frame} containing the name of covariates in
#'   the regression model (\code{var}), estimated regression coefficients
#'   (\code{est}), estimated hazard ratios (\code{exp_est = exp(est)}) standard
#'   error of the estimates (\code{se}), and the p-values (\code{pval}).}
#'   \item{optim_obj}{The \code{optim} output from maximising the log-likelihood
#'   in the weighted analysis.}
#' }
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
