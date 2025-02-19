utils::globalVariables(c(
  "..w_names", "A", "M", "R", "V_pseudo", "obs_weights", "two_phase_weights",
  "eif"
))


#' Fit propensity scores for treatment contrasts
#'
#' @param train_data .
#' @param valid_data .
#' @param contrast .
#' @param learners .
#' @param w_names .
#' @param type .
#' @param bounds .
#'
#' @importFrom data.table as.data.table copy setnames ":="
#' @importFrom sl3 sl3_Task
fit_treat_mech_RCT = function(train_data, valid_data = NULL, contrast, learners,
                              w_names, type = c("g", "h"), bounds = c(0.01, 0.99)) {
  if (type == "g") {
    cov_names <- w_names
  } else if (type == "h") {
    cov_names <- c("M", w_names)
    train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    train_data <- train_data[R == 1,]
    valid_data <- valid_data[R == 1,]
  }
  treat_task <- sl3::sl3_Task$new(data = train_data,
                                  weights = "obs_weights",
                                  covariates = cov_names,
                                  outcome = "A",
                                  outcome_type = "binomial")

  treat_fit <- learners$train(treat_task)
  treat_pred <- treat_fit$predict()

  if (is.null(valid_data)) {
    treat_pred_A_prime <- contrast[1] * treat_pred + (1 - contrast[1]) * (1 - treat_pred)
    treat_pred_A_star <- contrast[2] * treat_pred + (1 - contrast[2]) * (1 - treat_pred)
    out_treat_mat <- cbind(treat_pred_A_prime, treat_pred_A_star)
    out_treat_est <- apply(out_treat_mat, 2, function(x) {
      x_precise <- bound_precision(x)
      x_bounded <- bound_propensity(x_precise, bounds = bounds)
      return(x_bounded)
    })

    out_treat_est <- data.table::as.data.table(out_treat_est)
    data.table::setnames(out_treat_est, c("treat_pred_A_prime",
                                          "treat_pred_A_star"))
    out <- list(treat_est = out_treat_est, treat_fit = treat_fit)
  } else {
    out_treat_est <- lapply(list(train_data, valid_data),
                            function(data) {
                              treat_task <- sl3::sl3_Task$new(data = data,
                                                              weights = "obs_weights",
                                                              covariates = cov_names,
                                                              outcome = "A",
                                                              outcome_type = "binomial")

                              treat_pred <- treat_fit$predict(treat_task)

                              treat_pred_A_prime <- contrast[1] * treat_pred +
                                (1 - contrast[1]) * (1 - treat_pred)

                              treat_pred_A_star <- contrast[2] * treat_pred +
                                (1 - contrast[2]) * (1 - treat_pred)

                              out_treat_mat <- cbind(treat_pred_A_prime, treat_pred_A_star)
                              out_treat_est <- apply(out_treat_mat, 2, function(x) {
                                x_precise <- bound_precision(x)
                                x_bounded <- bound_propensity(x_precise, bounds = bounds)
                                return(x_bounded)
                              })
                              out_treat_est <- data.table::as.data.table(out_treat_est)
                              data.table::setnames(out_treat_est, c("treat_pred_A_prime",
                                                                    "treat_pred_A_star"))
                            })
    out <- list(treat_est_train = out_treat_est[[1]], treat_est_valid = out_treat_est[[2]],
                treat_fit = treat_fit)
  }
  return(out)
}






#' Title Fit outcome regression
#'
#' @param train_data .
#' @param valid_data .
#' @param contrast .
#' @param learners .
#' @param z_names .
#' @param w_names .
#' @param l_names .
#'
#' @importFrom data.table as.data.table copy setnames ":="
#' @importFrom sl3 sl3_Task
fit_out_mech_RCT = function(train_data, valid_data = NULL, contrast, learners,
                            z_names, w_names, l_names) {
  train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  train_data <- train_data[R == 1, ]
  valid_data <- valid_data[R == 1, ]

  b_natural_task <- sl3::sl3_Task$new(data = train_data,
                                      weights = "obs_weights",
                                      covariates = c(z_names, "M", "A", l_names, w_names),
                                      outcome = "Y")

  b_natural_fit <- learners$train(b_natural_task)
  b_natural_pred <- b_natural_fit$predict()

  if (is.null(valid_data)) {
    train_data_intervene <- data.table::copy(train_data)
    train_data_intervene[, `:=`(A, contrast[1])]
    b_natural_pred <- b_natural_fit$predict()
    b_intervened_prime_task <- sl3::sl3_Task$new(data = train_data_intervene,
                                                 weights = "obs_weights",
                                                 covariates = c(z_names, "M", "A", l_names, w_names),
                                                 outcome = "Y")
    b_intervened_pred_A_prime <- b_natural_fit$predict(b_intervened_prime_task)
    train_data_intervene[, `:=`(A, contrast[2])]
    b_intervened_star_task <- sl3::sl3_Task$new(data = train_data_intervene,
                                                weights = "obs_weights",
                                                covariates = c(z_names, "M", "A", l_names, w_names),
                                                outcome = "Y")
    b_intervened_pred_A_star <- b_natural_fit$predict(b_intervened_star_task)

    out_b_est <- data.table::as.data.table(cbind(b_natural_pred,
                                                 b_intervened_pred_A_prime, b_intervened_pred_A_star))
    data.table::setnames(out_b_est, c("b_pred_A_natural",
                                      "b_pred_A_prime", "b_pred_A_star"))
    out <- list(b_est = out_b_est, b_fit = b_natural_fit)
  } else {
    train_data_intervene <- data.table::copy(train_data)
    valid_data_intervene <- data.table::copy(valid_data)
    b_natural_pred_train <- b_natural_fit$predict()
    b_natural_task_valid <- sl3::sl3_Task$new(data = valid_data,
                                              weights = "obs_weights",
                                              covariates = c(z_names, "M", "A", l_names, w_names),
                                              outcome = "Y")
    b_natural_pred_valid <- b_natural_fit$predict(b_natural_task_valid)
    out_b_est <- lapply(list(train_data_intervene, valid_data_intervene),
                        function(data_intervene) {
                          data_intervene[, `:=`(A, contrast[1])]
                          b_intervened_prime_task <- sl3::sl3_Task$new(data = data_intervene,
                                                                       weights = "obs_weights",
                                                                       covariates = c(z_names, "M", "A", l_names, w_names),
                                                                       outcome = "Y")
                          b_intervened_pred_A_prime <- b_natural_fit$predict(b_intervened_prime_task)
                          data_intervene[, `:=`(A, contrast[2])]
                          b_intervened_star_task <- sl3::sl3_Task$new(data = data_intervene,
                                                                      weights = "obs_weights",
                                                                      covariates = c(z_names, "M", "A", l_names, w_names),
                                                                      outcome = "Y")
                          b_intervened_pred_A_star <- b_natural_fit$predict(b_intervened_star_task)

                          out_b_est <- data.table::as.data.table(cbind(b_intervened_pred_A_prime,
                                                                       b_intervened_pred_A_star))
                          return(out_b_est)
                        })
    out_b_est[[1]] <- cbind(b_natural_pred_train, out_b_est[[1]])
    out_b_est[[2]] <- cbind(b_natural_pred_valid, out_b_est[[2]])
    lapply(out_b_est, function(x) {
      data.table::setnames(x, c("b_pred_A_natural", "b_pred_A_prime",
                                "b_pred_A_star"))
    })
    out <- list(b_est_train = out_b_est[[1]], b_est_valid = out_b_est[[2]],
                b_fit = b_natural_fit)
  }
  return(out)
}




#' Title Fit q(M_k|a,w) & r(M_k|a,z,w)
#'
#' @param train_data .
#' @param valid_data .
#' @param contrast .
#' @param learners .
#' @param z_names .
#' @param w_names .
#' @param type .
#'
#' @importFrom data.table as.data.table copy setnames ":="
#' @importFrom sl3 sl3_Task
fit_moc_mech_RCT = function(train_data, valid_data = NULL, contrast, learners, z_names, w_names, type = c("q", "r")){
  if (type == "q") {
    cov_names <- w_names
    train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  } else if (type == "r") {
    cov_names <- c(z_names, w_names)
    train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    train_data <- train_data[R == 1, ]
    valid_data <- valid_data[R == 1, ]
  }
  suppressWarnings(moc_task <- sl3::sl3_Task$new(data = train_data,
                               weights = "obs_weights",
                               covariates = c("A", cov_names),
                               outcome = "M",
                               outcome_type = "binomial"))
  # fit model on observed data
  moc_fit <- learners$train(moc_task)

  if (is.null(valid_data)) {
    train_data_intervene <- data.table::copy(train_data)

    train_data_intervene[, `:=`(A, contrast[1])]
    moc_pred_A_natural <- moc_fit$predict()
    moc_prime_task <- sl3::sl3_Task$new(data = train_data_intervene,
                                        weights = "obs_weights",
                                        covariates = c("A", cov_names),
                                        outcome = "M",
                                        outcome_type = "binomial")
    moc_pred_A_prime <- moc_fit$predict(moc_prime_task)

    train_data_intervene[, `:=`(A, contrast[2])]
    moc_star_task <- sl3::sl3_Task$new(data = train_data_intervene,
                                       weights = "obs_weights",
                                       covariates = c("A", cov_names),
                                       outcome = "M",
                                       outcome_type = "binomial")
    moc_pred_A_star <- moc_fit$predict(moc_star_task)

    out_moc_est <- data.table::as.data.table(cbind(moc_pred_A_natural,
                                                   moc_pred_A_prime, moc_pred_A_star))
    data.table::setnames(out_moc_est, c("moc_pred_A_natural",
                                        "moc_pred_A_prime", "moc_pred_A_star"))
    out <- list(moc_est = out_moc_est,
                moc_fit = moc_fit)
  } else {
    train_data_intervene <- data.table::copy(train_data)
    valid_data_intervene <- data.table::copy(valid_data)
    moc_pred_A_natural_train <- moc_fit$predict()

    suppressWarnings(moc_task_valid <- sl3::sl3_Task$new(data = valid_data,
                                                         weights = "obs_weights",
                                                         covariates = c("A", cov_names),
                                                         outcome = "M",
                                                         outcome_type = "binomial"))

    moc_pred_A_natural_valid <- moc_fit$predict(moc_task_valid)

    out_moc_est <- lapply(list(train_data_intervene, valid_data_intervene),
                          function(data_intervene) {
                            data_intervene[, `:=`(A, contrast[1])]
                            suppressWarnings(moc_prime_task <- sl3::sl3_Task$new(data = data_intervene,
                                                                weights = "obs_weights",
                                                                covariates = c("A", cov_names),
                                                                outcome = "M",
                                                                outcome_type = "binomial"))
                            moc_pred_A_prime <- moc_fit$predict(moc_prime_task)

                            data_intervene[, `:=`(A, contrast[2])]
                            suppressWarnings(moc_star_task <- sl3::sl3_Task$new(data = data_intervene,
                                                               weights = "obs_weights",
                                                               covariates = c("A", cov_names),
                                                               outcome = "M",
                                                               outcome_type = "binomial"))
                            moc_pred_A_star <- moc_fit$predict(moc_star_task)
                            out_moc_est <- data.table::as.data.table(cbind(moc_pred_A_prime,
                                                                           moc_pred_A_star))
                          })
    out_moc_est[[1]] <- cbind(moc_pred_A_natural_train, out_moc_est[[1]])
    out_moc_est[[2]] <- cbind(moc_pred_A_natural_valid, out_moc_est[[2]])
    lapply(out_moc_est, function(x) {
      data.table::setnames(x, c("moc_pred_A_natural", "moc_pred_A_prime",
                                "moc_pred_A_star"))
    })
    out <- list(moc_est_train_M_one = out_moc_est[[1]], moc_est_valid_M_one = out_moc_est[[2]],
                moc_est_train_M_natural = out_moc_est[[1]] * train_data$M +
                  (1 - out_moc_est[[1]]) * (1 - train_data$M),
                moc_est_valid_M_natural = out_moc_est[[2]] * valid_data$M +
                  (1 - out_moc_est[[2]]) * (1 - valid_data$M),
                moc_fit = moc_fit)
  }
  return(out)
}




fit_nuisance_u_RCT = function(train_data,
                              valid_data,
                              learners,
                              b_out,
                              q_out,
                              r_out,
                              g_out,
                              h_out,
                              w_names){
  train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  b_prime <- b_out$b_est_train$b_pred_A_prime
  h_star <- h_out$treat_est_train$treat_pred_A_star
  g_star <- g_out$treat_est_train$treat_pred_A_star[train_data$R == 1]
  h_prime <- h_out$treat_est_train$treat_pred_A_prime
  g_prime <- g_out$treat_est_train$treat_pred_A_prime[train_data$R == 1]
  q_star_M_natural <- q_out$moc_est_train_M_natural$moc_pred_A_star[train_data$R == 1]
  r_prime_M_natural <- r_out$moc_est_train_M_natural$moc_pred_A_prime
  train_data <- train_data[R == 1, ]
  valid_data <- valid_data[R == 1, ]
  d_star <- (g_star/g_prime) * (q_star_M_natural/r_prime_M_natural) *
    (h_prime/h_star)
  u_pseudo_train <- b_prime * d_star
  if (stats::sd(u_pseudo_train) < .Machine$double.eps) {
    warning("U: constant pseudo-outcome, using intercept model.")
    learners <- sl3::Lrnr_mean$new()
  }
  u_data_train <- data.table::as.data.table(cbind(train_data[, ..w_names],
                                                  train_data$A,
                                                  train_data$M,
                                                  u_pseudo_train,
                                                  train_data$obs_weights))
  data.table::setnames(u_data_train, c(w_names,
                                       "A",
                                       "M",
                                       "U_pseudo",
                                       "obs_weights"))
  u_task_train <- sl3::sl3_Task$new(data = u_data_train,
                                    weights = "obs_weights",
                                    covariates = c("M", "A", w_names),
                                    outcome = "U_pseudo",
                                    outcome_type = "continuous")
  u_param_fit <- learners$train(u_task_train)
  u_data_valid <- data.table::as.data.table(cbind(valid_data[, ..w_names],
                                                  valid_data$A,
                                                  valid_data$M,
                                                  rep(0, nrow(valid_data)),
                                                  valid_data$obs_weights))
  data.table::setnames(u_data_valid, c(w_names,
                                       "A",
                                       "M",
                                       "U_pseudo",
                                       "obs_weights"))
  suppressWarnings(u_task_valid <- sl3::sl3_Task$new(data = u_data_valid,
                                                     weights = "obs_weights",
                                                     covariates = c("M", "A", w_names),
                                                     outcome = "U_pseudo",
                                                     outcome_type = "continuous"))
  u_valid_pred <- u_param_fit$predict(u_task_valid)
  u_train_pred <- u_param_fit$predict(u_task_train)
  return(list(u_fit = u_param_fit,
              u_pred = as.numeric(u_valid_pred),
              u_train_pred = as.numeric(u_train_pred)))
}




fit_nuisance_v_RCT = function(train_data,
                              valid_data,
                              contrast,
                              learners,
                              b_out,
                              q_out,
                              r_out,
                              z_names,
                              w_names,
                              effect_type){
  if (effect_type == "shift_k"){
    q_train_star_M_one <- q_out$moc_est_train_M_one$moc_pred_A_star[train_data$R == 1]
    q_valid_star_M_one <- q_out$moc_est_valid_M_one$moc_pred_A_star[valid_data$R == 1]
    train_data <- train_data[R == 1, ]
    valid_data <- valid_data[R == 1, ]

    v_pseudo <- lapply(unique(train_data$M), function(m_val) {
      train_data_m_interv <- data.table::copy(train_data)
      train_data_m_interv[, `:=`(obs_weights, two_phase_weights * obs_weights)]
      train_data_m_interv[, `:=`(
        M = m_val,
        A = contrast[1]
      )]

      b_reg_train_v_subtask <- sl3::sl3_Task$new(
        data = train_data_m_interv,
        weights = "obs_weights",
        covariates = c(z_names, "M", "A", w_names, "L"),
        outcome = "Y"
      )
      b_pred_train_m_interv <- b_out$b_fit$predict(b_reg_train_v_subtask)
      q_train_star_m_val <- (m_val * q_train_star_M_one) +
        (1 - m_val) * (1 - q_train_star_M_one)


      valid_data_m_interv <- data.table::copy(valid_data)
      valid_data_m_interv[, `:=`(obs_weights, two_phase_weights * obs_weights)]
      valid_data_m_interv[, `:=`(
        M = m_val,
        A = contrast[1]
      )]

      b_reg_valid_v_subtask <- sl3::sl3_Task$new(data = valid_data_m_interv,
                                                 weights = "obs_weights",
                                                 covariates = c(z_names, "M", "A", w_names, "L"),
                                                 outcome = "Y")
      b_pred_valid_m_interv <- b_out$b_fit$predict(b_reg_valid_v_subtask)
      q_valid_star_m_val <- (m_val * q_valid_star_M_one) +
        (1 - m_val) * (1 - q_valid_star_M_one)
      out_train <- b_pred_train_m_interv * q_train_star_m_val
      out_valid <- b_pred_valid_m_interv * q_valid_star_m_val
      out <- list(training = out_train, validation = out_valid,
                  b_train = b_pred_train_m_interv, b_valid = b_pred_valid_m_interv)
      return(out)
    })

    if (length(unique(train_data$M)) > 1) {
      v_pseudo_train <- v_pseudo[[1]]$training + v_pseudo[[2]]$training
      v_pseudo_valid <- v_pseudo[[1]]$validation + v_pseudo[[2]]$validation
    } else {
      v_pseudo_train <- v_pseudo[[1]]$training
      v_pseudo_valid <- v_pseudo[[1]]$validation
    }


    ## extract outcome model predictions with intervened M for TMLE fluctuation
    if (length(unique(train_data$M)) > 1) {
      b_pred_A_prime_M_zero <- v_pseudo[[1]]$b_valid
      b_pred_A_prime_M_one <- v_pseudo[[2]]$b_valid
    } else {
      b_pred_A_prime_M_zero <- rep(0, nrow(valid_data))
      b_pred_A_prime_M_one <- v_pseudo[[1]]$b_valid
    }

    if (stats::sd(v_pseudo_train) < .Machine$double.eps) {
      warning("V: constant pseudo-outcome, using intercept model.")
      learners <- sl3::Lrnr_mean$new()
    }

    train_data[, `:=`(V_pseudo, v_pseudo_train)]
    v_task_train <- sl3::sl3_Task$new(data = train_data,
                                      weights = "obs_weights",
                                      covariates = c("A", w_names),
                                      outcome = "V_pseudo",
                                      outcome_type = "continuous")

    valid_data[, `:=`(
      V_pseudo = v_pseudo_valid,
      A = contrast[1]
    )]

    v_task_valid <- sl3::sl3_Task$new(
      data = valid_data,
      weights = "obs_weights",
      covariates = c("A", w_names),
      outcome = "V_pseudo",
      outcome_type = "continuous"
    )
    v_param_fit <- learners$train(v_task_train)
    v_valid_pred <- v_param_fit$predict(v_task_valid)
    v_train_pred <- v_param_fit$predict(v_task_train)
    return(list(v_fit = v_param_fit,
                v_pred = as.numeric(v_valid_pred),
                v_train_pred = as.numeric(v_train_pred),
                v_pseudo = as.numeric(v_pseudo_valid),
                v_pseudo_train = as.numeric(v_pseudo_train),
                b_A_prime_M_zero = as.numeric(b_pred_A_prime_M_zero),
                b_A_prime_M_one = as.numeric(b_pred_A_prime_M_one)))
  } else if (effect_type == "shift_k_order") {
    train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
    b_prime <- b_out$b_est_train$b_pred_A_prime
    q_star_M_natural <- q_out$moc_est_train_M_natural$moc_pred_A_star[train_data$R == 1]
    r_prime_M_natural <- r_out$moc_est_train_M_natural$moc_pred_A_prime
    train_data <- train_data[R == 1, ]
    valid_data <- valid_data[R == 1, ]
    v_pseudo_train <- b_prime * q_star_M_natural / r_prime_M_natural
    if (stats::sd(v_pseudo_train) < .Machine$double.eps) {
      warning("V: constant pseudo-outcome, using intercept model.")
      learners <- sl3::Lrnr_mean$new()
    }
    var_names = c(w_names, z_names)
    v_data_train <- data.table::as.data.table(cbind(train_data[, ..var_names],
                                                    train_data$A,
                                                    train_data$M,
                                                    v_pseudo_train,
                                                    train_data$obs_weights))
    data.table::setnames(v_data_train, c(w_names, z_names,
                                         "A", "M",
                                         "V_pseudo", "obs_weights"))

    v_task_train <- sl3::sl3_Task$new(data = v_data_train,
                                      weights = "obs_weights",
                                      covariates = c("A", w_names),
                                      outcome = "V_pseudo",
                                      outcome_type = "continuous")

    v_param_fit <- learners$train(v_task_train)
    v_data_valid <- data.table::as.data.table(cbind(valid_data[, ..var_names],
                                                    valid_data$A,
                                                    valid_data$M,
                                                    rep(0, nrow(valid_data)),
                                                    valid_data$obs_weights))
    data.table::setnames(v_data_valid, c(w_names,
                                         z_names,
                                         "A",
                                         "M",
                                         "V_pseudo",
                                         "obs_weights"))
    suppressWarnings(v_task_valid <- sl3::sl3_Task$new(data = v_data_valid,
                                                       weights = "obs_weights",
                                                       covariates = c("A", w_names),
                                                       outcome = "V_pseudo",
                                                       outcome_type = "continuous"))
    v_valid_pred <- v_param_fit$predict(v_task_valid)
    v_train_pred <- v_param_fit$predict(v_task_train)

    s_task_train_pseudo <- sl3::sl3_Task$new(data = v_data_train,
                                             weights = "obs_weights",
                                             covariates = c("A", w_names, z_names),
                                             outcome = "V_pseudo",
                                             outcome_type = "continuous")
    s_param_fit_pseudo <- learners$train(s_task_train_pseudo)
    suppressWarnings(s_task_valid_pseudo <- sl3::sl3_Task$new(data = v_data_valid,
                                                              weights = "obs_weights",
                                                              covariates = c("A", w_names, z_names),
                                                              outcome = "V_pseudo",
                                                              outcome_type = "continuous"))
    s_valid_pred_pseudo <- s_param_fit_pseudo$predict(s_task_valid_pseudo)
    s_train_pred_pesudo <- s_param_fit_pseudo$predict(s_task_train_pseudo)

    return(list(v_fit = v_param_fit, v_pred = as.numeric(v_valid_pred),
                v_train_pred = as.numeric(v_train_pred),
                s_pseudo = as.numeric(s_valid_pred_pseudo),
                s_pseudo_train = as.numeric(s_train_pred_pesudo),
                b_prime_M_zero = NULL,
                b_prime_M_one = NULL))
  }
}



fit_nuisance_l_RCT = function(train_data,
                              valid_data,
                              contrast,
                              learners,
                              b_out,
                              z_names,
                              w_names,
                              effect_type){
  train_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  valid_data[, `:=`(obs_weights, two_phase_weights * obs_weights)]
  b_prime <- b_out$b_est_train$b_pred_A_prime
  train_data <- train_data[R == 1, ]
  valid_data <- valid_data[R == 1, ]
  if (stats::sd(b_prime) < .Machine$double.eps) {
    warning("B: constant pseudo-outcome, using intercept model.")
    learners <- sl3::Lrnr_mean$new()
  }
  var_names = c(w_names, z_names)
  l_data_train <- data.table::as.data.table(cbind(train_data[, ..var_names],
                                                  train_data$A,
                                                  train_data$M,
                                                  b_prime,
                                                  train_data$obs_weights))
  data.table::setnames(l_data_train, c(w_names,
                                       z_names,
                                       "A",
                                       "M",
                                       "b_prime",
                                       "obs_weights"))

  l_task_train <- sl3::sl3_Task$new(data = l_data_train,
                                    weights = "obs_weights",
                                    covariates = c("A", w_names, z_names, 'M'),
                                    outcome = "b_prime",
                                    outcome_type = "continuous")

  l_param_fit <- learners$train(l_task_train)
  l_data_valid <- data.table::as.data.table(cbind(valid_data[, ..var_names],
                                                  valid_data$A,
                                                  valid_data$M,
                                                  rep(0, nrow(valid_data)),
                                                  valid_data$obs_weights))
  data.table::setnames(l_data_valid, c(w_names,
                                       z_names,
                                       "A",
                                       "M",
                                       "b_prime",
                                       "obs_weights"))
  suppressWarnings(l_task_valid <- sl3::sl3_Task$new(data = l_data_valid,
                                                     weights = "obs_weights",
                                                     covariates = c("A", w_names, z_names, "M"),
                                                     outcome = "b_prime",
                                                     outcome_type = "continuous"))
  l_valid_pred <- l_param_fit$predict(l_task_valid)
  l_train_pred <- l_param_fit$predict(l_task_train)

  return(list(l_fit = l_param_fit, l_pred = as.numeric(l_valid_pred),
              l_train_pred = as.numeric(l_train_pred)))
}

