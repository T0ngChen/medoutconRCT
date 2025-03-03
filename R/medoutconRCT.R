#' Title .
#'
#' @param W .
#' @param A .
#' @param Z .
#' @param M .
#' @param L .
#' @param Y .
#' @param R .
#' @param obs_weights .
#' @param svy_weights .
#' @param two_phase_weights .
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
#' @param estimator .
#' @param estimator_args .
#' @param g_bounds .
#'
#'
#' @importFrom data.table as.data.table setnames set
#' @importFrom sl3 Lrnr_glm_fast Lrnr_hal9001
#' @importFrom stats var
medoutconRCT = function(
  W,
  A,
  Z,
  M,
  Y,
  L,
  R = rep(1, length(Y)),
  obs_weights = rep(1, length(Y)),
  svy_weights = NULL,
  two_phase_weights = rep(1, length(Y)),
  contrast = c(1, 0),
  g_learners = sl3::Lrnr_glm_fast$new(),
  h_learners = sl3::Lrnr_glm_fast$new(),
  b_learners = sl3::Lrnr_glm_fast$new(),
  q_learners = sl3::Lrnr_glm_fast$new(),
  r_learners = sl3::Lrnr_glm_fast$new(),
  u_learners = sl3::Lrnr_hal9001$new(),
  v_learners = sl3::Lrnr_hal9001$new(),
  l_learners = sl3::Lrnr_hal9001$new(),
  d_learners = sl3::Lrnr_glm_fast$new(),
  estimator = c("tmle", "onestep"),
  estimator_args = list(
    cv_folds = 5L,
    cv_strat = FALSE,
    strat_pmin = 0.1,
    max_iter = 5L,
    tiltmod_tol = 5
  ),
  g_bounds = c(0.005, 0.995)
) {
  # set defaults
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)
  est_args_os <- estimator_args[
    names(estimator_args) %in%
      names(formals(est_onestep_RCT))
  ]
  est_args_tmle <- estimator_args[
    names(estimator_args) %in%
      names(formals(est_tml_RCT))
  ]

  # define effect type shift_k, shift_k_order
  if (is.null(Z)) {
    Z = rep(1, length(Y))
  }

  w_names <- paste(
    "W",
    seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  z_names <- paste(
    "Z",
    seq_len(dim(data.table::as.data.table(Z))[2]),
    sep = "_"
  )

  if (is.null(L)) {
    L <- rep(1, length(Y))
    l_names = "L"
    effect_type <- "shift_k"
  } else {
    effect_type <- "shift_k_order"
    l_names <- paste(
      "L",
      seq_len(dim(data.table::as.data.table(L))[2]),
      sep = "_"
    )
  }

  # construct input data structure
  data <- data.table::as.data.table(cbind(
    Y,
    M,
    R,
    Z,
    L,
    A,
    W,
    obs_weights,
    two_phase_weights
  ))

  data.table::setnames(
    data,
    c(
      "Y",
      "M",
      "R",
      z_names,
      l_names,
      "A",
      w_names,
      "obs_weights",
      "two_phase_weights"
    )
  )

  # bound outcome Y in unit interval
  min_y <- min(data[["Y"]])
  max_y <- max(data[["Y"]])
  data.table::set(data, j = "Y", value = medoutcon:::scale_to_unit(data[["Y"]]))

  contrast_grid = list(contrast)
  est_params <- lapply(contrast_grid, function(contrast) {
    if (estimator == "onestep") {
      # EFFICIENT ONE-STEP ESTIMATOR
      onestep_est_args <- list(
        data = data,
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
        w_names = w_names,
        z_names = z_names,
        l_names = l_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )

      onestep_est_args <- unlist(
        list(
          onestep_est_args,
          est_args_os
        ),
        recursive = FALSE
      )
      est_out <- do.call(est_onestep_RCT, onestep_est_args)
    } else if (estimator == "tmle") {
      # TARGETED MINIMUM LOSS ESTIMATOR
      tmle_est_args <- list(
        data = data,
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
        w_names = w_names,
        z_names = z_names,
        l_names = l_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      tmle_est_args <- unlist(
        list(tmle_est_args, est_args_tmle),
        recursive = FALSE
      )
      est_out <- do.call(est_tml_RCT, tmle_est_args)
    }
    est_out$outcome <- as.numeric(Y)
    class(est_out) <- "medoutconRCT"
    est_params <- est_out
  })

  # calculate E[Y_1]
  est_params_out <- lapply(contrast_grid, function(contrast) {
    if (estimator == "onestep") {
      # EFFICIENT ONE-STEP ESTIMATOR
      data_out <- data.table::copy(data)
      data_out[, `:=`(M = 1)]
      data_out[, (z_names) := lapply(.SD, function(x) 1), .SDcols = z_names]
      data_out[, (l_names) := lapply(.SD, function(x) 1), .SDcols = l_names]
      onestep_est_args_out <- list(
        data = data_out,
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
        w_names = w_names,
        z_names = z_names,
        l_names = l_names,
        y_bounds = c(min_y, max_y),
        effect_type = "Y_1",
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      onestep_est_args_out <- unlist(
        list(onestep_est_args_out, est_args_os),
        recursive = FALSE
      )
      est_out <- do.call(est_onestep_RCT, onestep_est_args_out)
    } else if (estimator == "tmle") {
      # TARGETED MINIMUM LOSS ESTIMATOR
      data_out <- data.table::copy(data)
      data_out[, `:=`(M = 1)]
      data_out[, (z_names) := lapply(.SD, function(x) 1), .SDcols = z_names]
      data_out[, (l_names) := lapply(.SD, function(x) 1), .SDcols = l_names]
      tmle_est_args_out <- list(
        data = data_out,
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
        w_names = w_names,
        z_names = z_names,
        l_names = l_names,
        y_bounds = c(min_y, max_y),
        effect_type = "Y_1",
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      tmle_est_args_out <- unlist(
        list(tmle_est_args_out, est_args_tmle),
        recursive = FALSE
      )
      est_out <- do.call(est_tml_RCT, tmle_est_args_out)
    }
    est_out$outcome <- as.numeric(Y)
    class(est_out) <- "medoutconRCT"
    est_params <- est_out
  })

  de_theta_est <- est_params_out[[1]]$theta - est_params[[1]]$theta
  de_eif_est <- est_params_out[[1]]$eif - est_params[[1]]$eif
  de_eif_org <- est_params_out[[1]]$eif_est_rescaled -
    est_params[[1]]$eif_est_rescaled
  de_var_est <- stats::var(de_eif_est) / nrow(data)
  de_est_out <- list(
    theta = de_theta_est,
    var = de_var_est,
    eif = de_eif_est,
    type = estimator,
    param = paste(effect_type),
    outcome = as.numeric(Y)
  )
  class(de_est_out) <- "medoutconRCT"
  return(de_est_out)
}


medoutconPall = function(
  W,
  A,
  Z = NULL,
  M,
  Y,
  L = NULL,
  R = rep(1, length(Y)),
  obs_weights = rep(1, length(Y)),
  svy_weights = NULL,
  two_phase_weights = rep(1, length(Y)),
  contrast = NULL,
  g_learners = sl3::Lrnr_glm_fast$new(),
  h_learners = sl3::Lrnr_glm_fast$new(),
  b_learners = sl3::Lrnr_glm_fast$new(),
  q_learners = sl3::Lrnr_glm_fast$new(),
  r_learners = sl3::Lrnr_glm_fast$new(),
  u_learners = sl3::Lrnr_hal9001$new(),
  v_learners = sl3::Lrnr_hal9001$new(),
  d_learners = sl3::Lrnr_glm_fast$new(),
  estimator = c("tmle", "onestep"),
  estimator_args = list(
    cv_folds = 5L,
    cv_strat = FALSE,
    strat_pmin = 0.1,
    max_iter = 5L,
    tiltmod_tol = 5
  ),
  g_bounds = c(0.005, 0.995)
) {
  # set defaults
  estimator <- match.arg(estimator)
  estimator_args <- unlist(estimator_args, recursive = FALSE)
  est_args_os <- estimator_args[
    names(estimator_args) %in%
      names(formals(medoutcon:::est_onestep))
  ]
  est_args_tmle <- estimator_args[
    names(estimator_args) %in%
      names(formals(medoutcon:::est_tml))
  ]

  effect_type <- "natural"
  Z = rep(1, length(Y))
  L = rep(1, length(Y))
  # construct input data structure
  data <- data.table::as.data.table(cbind(
    Y,
    M,
    R,
    Z,
    L,
    A,
    W,
    obs_weights,
    two_phase_weights
  ))
  w_names <- paste(
    "W",
    seq_len(dim(data.table::as.data.table(W))[2]),
    sep = "_"
  )
  m_names <- paste(
    "M",
    seq_len(dim(data.table::as.data.table(M))[2]),
    sep = "_"
  )
  data.table::setnames(
    data,
    c(
      "Y",
      m_names,
      "R",
      "Z",
      "L",
      "A",
      w_names,
      "obs_weights",
      "two_phase_weights"
    )
  )

  # bound outcome Y in unit interval
  min_y <- min(data[["Y"]])
  max_y <- max(data[["Y"]])
  data.table::set(data, j = "Y", value = scale_to_unit(data[["Y"]]))

  contrast_grid = list(contrast)
  est_params <- lapply(contrast_grid, function(contrast) {
    if (estimator == "onestep") {
      # EFFICIENT ONE-STEP ESTIMATOR
      onestep_est_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        h_learners = h_learners,
        b_learners = b_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        d_learners = d_learners,
        w_names = w_names,
        m_names = m_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      onestep_est_args <- unlist(
        list(onestep_est_args, est_args_os),
        recursive = FALSE
      )
      est_out <- do.call(medoutcon:::est_onestep, onestep_est_args)
    } else if (estimator == "tmle") {
      # TARGETED MINIMUM LOSS ESTIMATOR
      tmle_est_args <- list(
        data = data,
        contrast = contrast,
        g_learners = g_learners,
        h_learners = h_learners,
        b_learners = b_learners,
        q_learners = q_learners,
        r_learners = r_learners,
        u_learners = u_learners,
        v_learners = v_learners,
        d_learners = d_learners,
        w_names = w_names,
        m_names = m_names,
        y_bounds = c(min_y, max_y),
        effect_type = effect_type,
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      tmle_est_args <- unlist(
        list(tmle_est_args, est_args_tmle),
        recursive = FALSE
      )
      est_out <- do.call(medoutcon:::est_tml, tmle_est_args)
    }

    # lazily create output as classed list
    est_out$outcome <- as.numeric(Y)
    class(est_out) <- "medoutcon"
    return(est_out)
  })

  # calculate E[Y_1]
  est_params_out <- lapply(contrast_grid, function(contrast) {
    if (estimator == "onestep") {
      # EFFICIENT ONE-STEP ESTIMATOR
      data_out <- data.table::copy(data)
      data_out[, `:=`(M = 1)]
      onestep_est_args <- list(
        data = data_out,
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
        w_names = w_names,
        z_names = "Z",
        l_names = "L",
        y_bounds = c(min_y, max_y),
        effect_type = "Y_1",
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      onestep_est_args <- unlist(
        list(onestep_est_args, est_args_os),
        recursive = FALSE
      )
      est_out <- do.call(est_onestep_RCT, onestep_est_args)
    } else if (estimator == "tmle") {
      # TARGETED MINIMUM LOSS ESTIMATOR
      data_out <- data.table::copy(data)
      data_out[, `:=`(M = 1)]
      tmle_est_args <- list(
        data = data_out,
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
        w_names = w_names,
        z_names = "Z",
        l_names = "L",
        y_bounds = c(min_y, max_y),
        effect_type = "Y_1",
        svy_weights = svy_weights,
        g_bounds = g_bounds
      )
      tmle_est_args <- unlist(
        list(tmle_est_args, est_args_tmle),
        recursive = FALSE
      )
      est_out <- do.call(est_tml_RCT, tmle_est_args)
    }
    est_out$outcome <- as.numeric(Y)
    class(est_out) <- "medoutconRCT"
    est_params <- est_out
  })

  de_theta_est <- est_params_out[[1]]$theta - est_params[[1]]$theta
  de_eif_est <- est_params_out[[1]]$eif - est_params[[1]]$eif
  de_eif_org <- est_params_out[[1]]$eif_est_rescaled -
    est_params[[1]]$eif_est_rescaled
  de_var_est <- stats::var(de_eif_est) / nrow(data)
  de_est_out <- list(
    theta = de_theta_est,
    var = de_var_est,
    eif = de_eif_est,
    type = estimator,
    param = paste(effect_type),
    outcome = as.numeric(Y)
  )
  class(de_est_out) <- "medoutconRCT"
  return(de_est_out)
}
