utils::globalVariables(c("..w_names", "A", "M", "Y", "R", "v_prime"))


#' Title
#'
#' @param fold .
#' @param data_in .
#' @param contrast .
#' @param g_learners .
#' @param h_learners .
#' @param b_learners .
#' @param q_learners .
#' @param r_learners .
#' @param u_learners .
#' @param v_learners .
#' @param d_learners .
#' @param l_learners .
#' @param effect_type .
#' @param l_names .
#' @param w_names .
#' @param z_names .
#' @param g_bounds .
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table data.table copy
#' @importFrom origami training validation fold_index
#' @importFrom sl3 Lrnr_mean
cv_eif_RCT = function(fold,
                      data_in,
                      contrast,
                      g_learners,
                      h_learners,
                      b_learners,
                      q_learners,
                      r_learners,
                      u_learners,
                      v_learners,
                      d_learners,
                      l_learners,
                      effect_type = c("shift_k", "shift_k_order", "Y_1"),
                      w_names,
                      z_names,
                      l_names,
                      g_bounds = c(0.005, 0.995)){
  # make training and validation data
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)

  # 1) fit regression for propensity score regression
  g_out <- fit_treat_mech_RCT(train_data = train_data,
                              valid_data = valid_data,
                              contrast = contrast,
                              learners = g_learners,
                              w_names = w_names,
                              type = "g",
                              bounds = g_bounds)

  # 2) fit clever regression for treatment, conditional on mediators
  if(effect_type == "Y_1"){
    h_learners <- sl3::Lrnr_mean$new()
  }
  h_out <- fit_treat_mech_RCT(train_data = train_data,
                              valid_data = valid_data,
                              contrast = contrast,
                              learners = h_learners,
                              w_names = w_names,
                              type = "h",
                              bounds = g_bounds)


  b_out <- fit_out_mech_RCT(train_data = train_data,
                            valid_data = valid_data,
                            contrast = contrast,
                            learners = b_learners,
                            w_names = w_names,
                            z_names = z_names,
                            l_names = l_names)
  if (effect_type == "Y_1") {
    q_learners <- sl3::Lrnr_mean$new()
  }
  q_out <- fit_moc_mech_RCT(train_data = train_data,
                            valid_data = valid_data,
                            contrast = contrast,
                            learners = q_learners,
                            z_names = z_names,
                            w_names = w_names,
                            type = "q")

  if (effect_type == "Y_1") {
    r_learners <- sl3::Lrnr_mean$new()
  }
  r_out <- fit_moc_mech_RCT(train_data = train_data,
                            valid_data = valid_data,
                            contrast = contrast,
                            learners = r_learners,
                            z_names = z_names,
                            w_names = w_names,
                            type = "r")

  b_prime <- b_out$b_est_valid$b_pred_A_prime
  h_star <- h_out$treat_est_valid$treat_pred_A_star
  g_star <- g_out$treat_est_valid$treat_pred_A_star[valid_data$R == 1]
  h_prime <- h_out$treat_est_valid$treat_pred_A_prime
  g_prime <- g_out$treat_est_valid$treat_pred_A_prime[valid_data$R == 1]

  q_star_M_one <- q_out$moc_est_valid_M_one$moc_pred_A_star[valid_data$R == 1]
  r_prime_M_one <- r_out$moc_est_valid_M_one$moc_pred_A_prime
  q_star_M_natural <- q_out$moc_est_valid_M_natural$moc_pred_A_star[valid_data$R == 1]
  r_prime_M_natural <- r_out$moc_est_valid_M_natural$moc_pred_A_prime

  # set A in the validation data.
  valid_data_a_prime <- data.table::copy(valid_data)[, `:=`(A, contrast[1])]
  valid_data_a_star <- data.table::copy(valid_data)[, `:=`(A, contrast[2])]

  u_out <- fit_nuisance_u_RCT(train_data = train_data,
                              valid_data = valid_data_a_prime,
                              learners = u_learners,
                              b_out = b_out,
                              q_out = q_out,
                              r_out = r_out,
                              g_out = g_out,
                              h_out = h_out,
                              w_names = w_names)
  u_prime <- u_out$u_pred

  v_out <- fit_nuisance_v_RCT(train_data = train_data,
                              valid_data = valid_data_a_prime,
                              contrast = contrast,
                              learners = v_learners,
                              b_out = b_out,
                              q_out = q_out,
                              r_out = r_out,
                              z_names = z_names,
                              w_names = w_names,
                              effect_type = effect_type)
  v_prime <- v_out$v_pred

  # NOTE: assuming M in {0,1}; other cases not supported yet
  u_int_eif <- lapply(c(1, 0), function(m_val) {

    valid_data_m_interv <- data.table::copy(valid_data[R == 1, ])

    valid_data_m_interv[, `:=`(M = m_val,
                               A = contrast[1],
                               U_pseudo = u_prime)]

    u_task_valid_m_interv <- sl3::sl3_Task$new(data = valid_data_m_interv,
                                               weights = "obs_weights",
                                               covariates = c("M", "A", w_names),
                                               outcome = "U_pseudo",
                                               outcome_type = "continuous")

    out_valid <- u_out[["u_fit"]]$predict(u_task_valid_m_interv)
    return(out_valid)
  })
  u_int_eif <- do.call(`-`, u_int_eif)

  # IPW
  ipw_a_prime <- as.numeric(valid_data[R == 1, A] == contrast[1])/g_prime
  ipw_a_star <- as.numeric(valid_data[R == 1, A] == contrast[2])/g_star

  # q(m_k|a^star)/r(m_k|a^prime)
  c_star <- (q_star_M_natural/r_prime_M_natural)

  eif_y <- ipw_a_prime * c_star / mean(ipw_a_prime * c_star) *
    (valid_data[R == 1, Y] - b_prime)
  eif_u <- ipw_a_star/mean(ipw_a_star) * u_int_eif *
    (valid_data[R == 1, M] - q_star_M_one)
  if (effect_type == "shift_k") {
    eif_v <- ipw_a_prime/mean(ipw_a_prime) * (v_out$v_pseudo - v_prime)
  } else if (effect_type == "shift_k_order") {
    eif_v <- ipw_a_prime/mean(ipw_a_prime) * (v_out$s_pseudo - v_prime)
  }

  if (effect_type == "Y_1") {
    assertthat::assert_that(all(eif_u == 0))
  }

  if(effect_type == "Y_1") {
    eif <- ipw_a_prime/mean(ipw_a_prime) * (valid_data[R == 1, Y] - b_prime) + b_prime
    eif_save = NULL
  } else if(effect_type == "shift_k") {
    eif <- eif_y + eif_u + eif_v + v_prime
    eif_save = NULL
  } else if(effect_type == "shift_k_order") {
    l_out <- fit_nuisance_l_RCT(train_data = train_data,
                              valid_data = valid_data_a_prime,
                              contrast = contrast,
                              learners = l_learners,
                              b_out = b_out,
                              z_names = z_names,
                              w_names = w_names,
                              effect_type = effect_type)
    l_prime <- l_out$l_pred
    eif_ad = ipw_a_prime * c_star / mean(ipw_a_prime * c_star) *
      (b_prime - l_prime)
    eif_save = eif_v + eif_ad + v_prime
    eif <- eif_y + eif_u + eif_v + eif_ad + v_prime
  }


  if (!all(data_in$R == 1) || !(all(data_in$two_phase_weights == 1))) {
    plugin_est <- est_plugin(v_pred = v_prime)
    centered_eif <- eif - plugin_est
    d_out <- medoutcon:::fit_nuisance_d(train_data = train_data, valid_data = valid_data,
                                        contrast = contrast, learners = d_learners, b_out = b_out,
                                        g_out = g_out, h_out = h_out, q_out = q_out, r_out = r_out,
                                        u_out = u_out, v_out = v_out, z_names = z_names,
                                        w_names = w_names)
    centered_eif_pred <- d_out$d_pred
    full_eif <- two_phase_eif(R = valid_data$R, two_phase_weights = valid_data$two_phase_weights,
                              eif = centered_eif, eif_predictions = centered_eif_pred,
                              plugin_est = plugin_est)
  } else {
    full_eif <- eif
    centered_eif_pred <- NA
  }
  out <- list(tmle_components = data.table::data.table(g_prime = g_prime,
                                                       g_star = g_star, h_prime = h_prime, h_star = h_star,
                                                       q_star_M_natural = q_star_M_natural, q_star_M_one = q_star_M_one,
                                                       r_prime_M_natural = r_prime_M_natural, r_prime_M_one = r_prime_M_one,
                                                       v_prime = v_prime, u_int_diff = u_int_eif, b_prime = b_prime,
                                                       b_prime_M_zero = v_out$b_A_prime_M_zero, b_prime_M_one = v_out$b_A_prime_M_one,
                                                       D_star = eif, fold = origami::fold_index(), eif_save = eif_save),
              D_star = full_eif,
              D_pred = centered_eif_pred)
  return(out)
}







#' Title
#'
#' @param data .
#' @param contrast .
#' @param g_learners .
#' @param h_learners .
#' @param b_learners .
#' @param q_learners .
#' @param r_learners .
#' @param u_learners .
#' @param v_learners .
#' @param l_learners .
#' @param d_learners .
#' @param w_names .
#' @param z_names .
#' @param l_names .
#' @param y_bounds .
#' @param g_bounds .
#' @param effect_type .
#' @param svy_weights .
#' @param cv_folds .
#'
#' @importFrom assertthat assert_that
#' @importFrom stats var weighted.mean
#' @importFrom origami make_folds cross_validate folds_vfold
est_onestep_RCT = function(data,
                           contrast,
                           g_learners,
                           h_learners,
                           b_learners,
                           q_learners,
                           r_learners,
                           u_learners,
                           v_learners,
                           d_learners,
                           l_learners,
                           w_names,
                           z_names,
                           l_names,
                           y_bounds,
                           g_bounds = c(0.005, 0.995),
                           effect_type = c("shift_k", "shift_k_order", "Y_1"),
                           svy_weights = NULL,
                           cv_folds = 5L,
                           cv_strat = FALSE,
                           strat_pmin = 0.1){
  # make sure that more than one fold is specified
  assertthat::assert_that(cv_folds > 1L)

  # create cross-validation folds
  if (cv_strat && data[, mean(Y) <= strat_pmin]) {
    # check that outcome is binary for stratified V-fold cross-validation
    assertthat::assert_that(data[, all(unique(Y) %in% c(0, 1))])
    folds <- origami::make_folds(
      data,
      fold_fun = origami::folds_vfold,
      V = cv_folds,
      strata_ids = data$Y)
  } else {
    folds <- origami::make_folds(data,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds)
  }

  # estimate the EIF on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = cv_eif_RCT,
    folds = folds,
    data_in = data,
    contrast = contrast,
    g_learners = g_learners,
    h_learners = h_learners,
    b_learners = b_learners,
    q_learners = q_learners,
    r_learners = r_learners,
    u_learners = u_learners,
    v_learners = v_learners,
    d_learners = d_learners,
    l_learners = l_learners,
    effect_type = effect_type,
    w_names = w_names,
    z_names = z_names,
    l_names = l_names,
    g_bounds = g_bounds,
    use_future = FALSE,
    .combine = FALSE)

  # get estimated efficient influence function
  v_prime <- do.call(rbind, cv_eif_results[[1]])$v_prime
  obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
  cv_eif_est <- unlist(cv_eif_results$D_star)[order(obs_valid_idx)]

  # re-scale efficient influence function
  eif_est_rescaled <- cv_eif_est %>%
    scale_from_unit(y_bounds[2], y_bounds[1])

  # calculate one-step estimate and variance from efficient influence function
  if (is.null(svy_weights)) {
    os_est <- mean(eif_est_rescaled)
    eif_est_out <- eif_est_rescaled
  } else {
    os_est <- stats::weighted.mean(eif_est_rescaled, svy_weights)
    eif_est_out <- eif_est_rescaled * svy_weights
  }
  os_var <- stats::var(eif_est_out) / length(eif_est_out)

  # return results
  os_est_out <- list(theta = os_est,
                     theta_plugin = ifelse(effect_type == "Y_1", NA, est_plugin(v_prime)),
                     var = os_var,
                     eif = (eif_est_out - os_est),
                     type = "onestep")
  return(os_est_out)
}


est_plugin <- function(v_pred) {
  mean(v_pred)
}


#' Title
#'
#' @param data .
#' @param contrast .
#' @param g_learners .
#' @param h_learners .
#' @param b_learners .
#' @param q_learners .
#' @param r_learners .
#' @param u_learners .
#' @param v_learners .
#' @param l_learners .
#' @param d_learners .
#' @param w_names .
#' @param z_names .
#' @param l_names .
#' @param y_bounds .
#' @param g_bounds .
#' @param effect_type .
#' @param svy_weights .
#' @param cv_folds .
#' @param max_iter .
#' @param tiltmod_tol .
#' @importFrom dplyr "%>%"
#' @importFrom assertthat assert_that
#' @importFrom origami make_folds cross_validate folds_vfold
#' @importFrom stats var as.formula plogis qlogis coef predict weighted.mean
#'   binomial
#' @importFrom glm2 glm2
est_tml_RCT <- function(data,
                        contrast,
                        g_learners,
                        h_learners,
                        b_learners,
                        q_learners,
                        r_learners,
                        u_learners,
                        v_learners,
                        l_learners,
                        d_learners,
                        w_names,
                        z_names,
                        l_names,
                        y_bounds,
                        g_bounds = c(0.005, 0.995),
                        effect_type = c("shift_k", "shift_k_order", "Y_1"),
                        svy_weights = NULL,
                        cv_folds = 5L,
                        cv_strat = FALSE,
                        strat_pmin = 0.1,
                        max_iter = 5L,
                        tiltmod_tol = 5) {
  # make sure that more than one fold is specified
  assertthat::assert_that(cv_folds > 1L)

  # create cross-validation folds
  if (cv_strat && data[, mean(Y) <= strat_pmin]) {
    # check that outcome is binary for stratified V-fold cross-validation
    assertthat::assert_that(data[, all(unique(Y) %in% c(0, 1))])
    folds <- origami::make_folds(
      data,
      fold_fun = origami::folds_vfold,
      V = cv_folds,
      strata_ids = data$Y)
  } else {
    folds <- origami::make_folds(data,
                                 fold_fun = origami::folds_vfold,
                                 V = cv_folds)
  }

  # perform the cv_eif procedure on a per-fold basis
  cv_eif_results <- origami::cross_validate(
    cv_fun = cv_eif_RCT,
    folds = folds,
    data_in = data,
    contrast = contrast,
    g_learners = g_learners,
    h_learners = h_learners,
    b_learners = b_learners,
    q_learners = q_learners,
    r_learners = r_learners,
    u_learners = u_learners,
    v_learners = v_learners,
    l_learners = l_learners,
    d_learners = d_learners,
    effect_type = effect_type,
    w_names = w_names,
    z_names = z_names,
    l_names = l_names,
    g_bounds = g_bounds,
    use_future = FALSE,
    .combine = FALSE
  )

  # concatenate nuisance function function estimates
  # make sure that data is in the same order as the concatenated validations
  # sets
  cv_eif_est <- do.call(rbind, cv_eif_results[[1]])
  obs_valid_idx <- do.call(c, lapply(folds, `[[`, "validation_set"))
  data <- data[obs_valid_idx]

  # extract nuisance function estimates and auxiliary quantities
  g_prime <- cv_eif_est$g_prime
  h_prime <- cv_eif_est$h_prime
  g_star <- cv_eif_est$g_star
  h_star <- cv_eif_est$h_star
  r_prime_M_natural <- cv_eif_est$r_prime_M_natural
  r_prime_M_one <- cv_eif_est$r_prime_M_one
  q_star_M_one <- cv_eif_est$q_star_M_one
  q_star_M_natural <- cv_eif_est$q_star_M_natural


  b_prime_M_one <- cv_eif_est$b_prime_M_one
  b_prime_M_zero <- cv_eif_est$b_prime_M_zero
  b_prime_M_natural <- cv_eif_est$b_prime

  # generate inverse weights and multiplier for auxiliary covariates
  ipw_prime <- as.numeric(data[R == 1, A] == contrast[1]) / g_prime
  ipw_star <- as.numeric(data[R == 1, A] == contrast[2]) / g_star

  # prepare for iterative targeting
  eif_stop_crit <- FALSE
  n_iter <- 0
  n_obs <- nrow(data)
  se_eif <- sqrt(var(cv_eif_est$D_star) / n_obs)
  tilt_stop_crit <- se_eif / log(n_obs)
  b_score <- q_score <- Inf
  tilt_two_phase_weights <- sum(data$R) != nrow(data)
  d_pred <- unlist(cv_eif_results$D_pred)[order(obs_valid_idx)]

  # perform iterative targeting
  if(effect_type != "Y_1"){
    while (!eif_stop_crit && n_iter <= max_iter) {
      # NOTE: check convergence condition for outcome regression
      if (mean(b_score) > tilt_stop_crit) {
        # compute auxiliary covariates from updated estimates
        c_star_M_natural <- q_star_M_natural/r_prime_M_natural
        c_star_M_one <- q_star_M_one / r_prime_M_one
        c_star_M_zero <- (1 - q_star_M_one) / (1 - r_prime_M_one)

        # bound and transform nuisance estimates for tilting regressions
        b_prime_M_natural_logit <- b_prime_M_natural %>%
          bound_precision() %>%
          stats::qlogis()
        b_prime_M_one_logit <- b_prime_M_one %>%
          bound_precision() %>%
          stats::qlogis()
        b_prime_M_zero_logit <- b_prime_M_zero %>%
          bound_precision() %>%
          stats::qlogis()

        # fit tilting model for the outcome mechanism
        c_star_b_tilt <- c_star_M_natural
        if (tilt_two_phase_weights) {
          weights_b_tilt <- as.numeric(data[R == 1, A] == contrast[1]) /
            g_prime * as.numeric(data[R == 1, two_phase_weights])
        } else {
          weights_b_tilt <- data$obs_weights * (data$A == contrast[1]) / g_prime
        }

        suppressWarnings(
          b_tilt_fit <- glm2::glm2(
            stats::as.formula("y_scaled ~ -1 + offset(b_prime_logit) + c_star"),
            data = data.table::as.data.table(list(
              y_scaled = data[R == 1, Y],
              b_prime_logit = b_prime_M_natural_logit,
              c_star = c_star_b_tilt
            )),
            weights = weights_b_tilt,
            family = stats::binomial(),
            start = 0
          )
        )
        if (is.na(stats::coef(b_tilt_fit))) {
          b_tilt_fit$coefficients <- 0
        } else if (!b_tilt_fit$converged || abs(max(stats::coef(b_tilt_fit))) >
                   tiltmod_tol) {
          b_tilt_fit$coefficients <- 0
        }
        b_tilt_coef <- unname(stats::coef(b_tilt_fit))

        # update nuisance estimates via tilting regressions for outcome
        b_prime_M_natural <- stats::plogis(b_prime_M_natural_logit +
                                             b_tilt_coef * c_star_M_natural)
        b_prime_M_one <- stats::plogis(b_prime_M_one_logit +
                                         b_tilt_coef * c_star_M_one)
        b_prime_M_zero <- stats::plogis(b_prime_M_zero_logit +
                                          b_tilt_coef * c_star_M_zero)

        # compute efficient score for outcome regression component
        b_score <- data[R == 1, two_phase_weights] *
          ipw_prime * c_star_M_natural * (data[R == 1, Y] - b_prime_M_natural)
      } else {
        b_score <- 0
      }
      # NOTE: check convergence condition for intermediate confounding
      if (mean(q_score) > tilt_stop_crit) {
        # perform iterative targeting for intermediate confounding
        q_star_M_one_logit <- q_star_M_one %>%
          bound_precision() %>%
          stats::qlogis()

        # fit tilting regressions for intermediate confounding
        u_prime_diff_q_tilt <- cv_eif_est$u_int_diff
        if (tilt_two_phase_weights) {
          weights_q_tilt <- as.numeric(data[R == 1, A] == contrast[2]) /
            g_star * as.numeric(data[R == 1, two_phase_weights])
        } else {
          weights_q_tilt <- data$obs_weights * (data$A == contrast[2]) / g_star
        }
        suppressWarnings(
          q_tilt_fit <- glm2::glm2(
            stats::as.formula("M ~ -1 + offset(q_star_logit) + u_prime_diff"),
            data = data.table::as.data.table(list(
              M = data[R == 1, M],
              q_star_logit = q_star_M_one_logit,
              u_prime_diff = u_prime_diff_q_tilt
            )),
            weights = weights_q_tilt,
            family = stats::binomial(),
            start = 0
          )
        )

        if (is.na(stats::coef(q_tilt_fit))) {
          q_tilt_fit$coefficients <- 0
        } else if (!q_tilt_fit$converged || abs(max(stats::coef(q_tilt_fit))) >
                   tiltmod_tol) {
          q_tilt_fit$coefficients <- 0
        }
        q_tilt_coef <- unname(stats::coef(q_tilt_fit))

        # update nuisance estimates via tilting of intermediate confounder
        q_star_M_one <- stats::plogis(q_star_M_one_logit + q_tilt_coef *
                                        cv_eif_est$u_int_diff)
        q_star_M_natural <- (data[R == 1, M] * q_star_M_one) +
          ((1 - data[R == 1, M]) * (1 - q_star_M_one))

        # compute efficient score for intermediate confounding component
        q_score <- ipw_star * cv_eif_est$u_int_diff *
          (data[R == 1, M] - q_star_M_one) *
          (data[R == 1, two_phase_weights])
      } else {
        q_score <- 0
      }

      # check convergence and iterate the counter
      eif_stop_crit <- all(
        abs(c(mean(b_score), mean(q_score))) < tilt_stop_crit
      )
      n_iter <- n_iter + 1
    }
    # update auxiliary covariates after completion of iterative targeting
    c_star_M_natural <- q_star_M_natural / r_prime_M_natural
    c_star_M_one <- q_star_M_one / r_prime_M_one
    c_star_M_zero <- (1 - q_star_M_one) / (1 - r_prime_M_one)

    if(effect_type == "shift_k") {
      # compute updated substitution estimator and prepare for tilting regression
      v_pseudo <- ((b_prime_M_one * q_star_M_one) +
                     (b_prime_M_zero * (1 - q_star_M_one))) %>%
        bound_precision()
      v_prime_logit <- cv_eif_est$v_prime %>%
        bound_precision() %>%
        stats::qlogis()

      # fit tilting regression for substitution estimator
      if (tilt_two_phase_weights) {
        weights_v_tilt <- (as.numeric(data[R == 1, A]) == contrast[1]) / g_prime *
          (as.numeric(data[R == 1, two_phase_weights]))
      } else {
        weights_v_tilt <- data$obs_weights * (data$A == contrast[1]) / g_prime
      }
      suppressWarnings(
        v_tilt_fit <- glm2::glm2(
          stats::as.formula("v_pseudo ~ offset(v_prime_logit)"),
          data = data.table::as.data.table(list(
            v_pseudo = v_pseudo,
            v_prime_logit = v_prime_logit
          )),
          weights = weights_v_tilt,
          family = stats::binomial(),
          start = 0
        )
      )
      v_prime_tmle <- unname(stats::predict(v_tilt_fit, type = "response"))
    } else if (effect_type == "shift_k_order") {
      v_prime_tmle <- b_score + q_score + cv_eif_est$eif_save
    }
  } else if (effect_type == "Y_1") {

    # bound and transform nuisance estimates for tilting regressions
    b_prime_M_natural_logit <- b_prime_M_natural %>%
      bound_precision() %>%
      stats::qlogis()

    # fit tilting model for the outcome mechanism
    if (tilt_two_phase_weights) {
      weights_b_tilt <- as.numeric(data[R == 1, A] == contrast[1]) /
        g_prime * as.numeric(data[R == 1, two_phase_weights])
    } else {
      weights_b_tilt <- data$obs_weights * (data$A == contrast[1]) / g_prime
    }

    suppressWarnings(
      b_tilt_fit <- glm2::glm2(
        stats::as.formula("y_scaled ~ offset(b_prime_logit)"),
        data = data.table::as.data.table(list(
          y_scaled = data[R == 1, Y],
          b_prime_logit = b_prime_M_natural_logit
        )),
        weights = weights_b_tilt,
        family = stats::binomial(),
        start = 0
      )
    )
    v_prime_tmle = unname(stats::predict(b_tilt_fit, type = "response"))
  }

  # compute influence function with centering at the TML estimate
  # make sure that it's in the same order as the original data
  eif_est <- unlist(cv_eif_results$D_star)[order(obs_valid_idx)]

  # re-scale efficient influence function
  v_prime_tmle_rescaled <- v_prime_tmle %>%
    scale_from_unit(y_bounds[2], y_bounds[1])
  eif_est_rescaled <- eif_est %>%
    scale_from_unit(y_bounds[2], y_bounds[1])

  # compute TML estimator and variance from efficient influence function
  if (is.null(svy_weights)) {
    tml_est <- mean(v_prime_tmle_rescaled)
    eif_est_out <- eif_est_rescaled
  } else {
    # compute a re-weighted TMLE, with re-weighted influence function
    # NOTE: make sure that survey weights are ordered like the concatenated
    #       validation sets
    svy_weights <- svy_weights[obs_valid_idx]
    tml_est <- stats::weighted.mean(v_prime_tmle_rescaled, svy_weights)
    eif_est_out <- eif_est_rescaled * svy_weights
  }
  tmle_var <- stats::var(eif_est_out) / length(eif_est_out)

  # output
  tmle_out <- list(
    theta = tml_est,
    theta_plugin = ifelse(effect_type == "Y_1", NA, est_plugin(cv_eif_est$v_prime)),
    var = tmle_var,
    eif = (eif_est_out - tml_est),
    n_iter = n_iter,
    type = "tmle"
  )

  return(tmle_out)
}


