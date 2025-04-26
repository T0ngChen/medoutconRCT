#' Interventional (in)direct effects for estimating estimands \eqn{\theta_k} and \eqn{\theta'_k} using causal ML
#'
#' @param W A data.frame or matrix of baseline covariates (confounders).
#' @param A A vector of the binary exposure/treatment assignments.
#' @param Z A numeric vector, matrix, or data frame of intermediate confounders. When estimating \eqn{\theta'_k}, they are the causal ancestors of mediators.
#' @param M A vector of binary mediator of interest.
#' @param L A numeric vector, matrix, or data frame containing variables that are causal descendants of the mediators of interest. This arugment is only relevant when estimating \eqn{\theta'_k}.
#' @param Y Outcome vector
#' @param contrast Numeric vector of length 2 giving treatment contrast (e.g. c(1,0)).
#' @param g_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters g
#' @param h_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' @param b_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating outcome regression b()
#' @param q_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters q
#' @param r_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters r
#' @param u_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters u
#' @param v_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters v
#' @param l_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' @param d_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' @param estimator Character string, either "tmle" or "onestep", specifying the estimator.
#' @param estimator_args estimator_args List of control parameters for the estimator:
#'   - cv_folds: integer, number of cross-validation folds
#'   - cv_strat: logical, whether to stratify folds by outcome
#'   - strat_pmin: numeric, minimum stratum proportion
#'   - max_iter: integer, maximum TMLE updating iterations
#'   - tiltmod_tol: numeric, convergence tolerance for tilting
#' @param g_bounds A numeric vector specifying the upper and lower bound for the estimated propensity scores.
#'
#'
#' @importFrom data.table as.data.table setnames set
#' @importFrom sl3 Lrnr_glm_fast Lrnr_hal9001
#' @importFrom stats var
#'
#' @export
medoutconRCT = function(
  W,
  A,
  Z,
  M,
  Y,
  L,
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
  # set default for two-phase sampling parameters in medoutcon (unused)
  R = rep(1, length(Y))
  obs_weights = rep(1, length(Y))
  svy_weights = NULL
  two_phase_weights = rep(1, length(Y))
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


#' Interventional (in)direct effects for estimating estimands \eqn{\theta_{all}} using causal ML
#'
#' @param W A data.frame or matrix of baseline covariates (confounders).
#' @param A A vector of the binary exposure/treatment assignments.
#' @param M A vector of binary mediator of interest.
#' @param Y Outcome vector
#' @param contrast Numeric vector of length 2 giving treatment contrast (e.g. c(1,0)).
#' @param g_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters g
#' @param h_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' @param b_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating outcome regression b()
#' @param q_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters q
#' @param r_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters r
#' @param u_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters u
#' @param v_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' This learner is used for estimating nuisance parameters v
#' @param d_learners A learner object. It is either an \code{\link[sl3]{Stack}} or any class inheriting from \code{\link[sl3]{Lrnr_base}} that contains one or more instantiated \pkg{sl3} learners.
#' @param estimator Character string, either "tmle" or "onestep", specifying the estimator.
#' @param estimator_args estimator_args List of control parameters for the estimator:
#'   - cv_folds: integer, number of cross-validation folds
#'   - cv_strat: logical, whether to stratify folds by outcome
#'   - strat_pmin: numeric, minimum stratum proportion
#'   - max_iter: integer, maximum TMLE updating iterations
#'   - tiltmod_tol: numeric, convergence tolerance for tilting
#' @param g_bounds A numeric vector specifying the upper and lower bound for the estimated propensity scores.
#'
#'
#' @importFrom data.table as.data.table setnames set
#' @importFrom sl3 Lrnr_glm_fast Lrnr_hal9001
#' @importFrom stats var
#'
#' @export
medoutconPall = function(
  W,
  A,
  M,
  Y,
  contrast = c(1, 0),
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
  # set default for two-phase sampling parameters in medoutcon (unused)
  R = rep(1, length(Y))
  obs_weights = rep(1, length(Y))
  svy_weights = NULL
  two_phase_weights = rep(1, length(Y))
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
  l_learners = sl3::Lrnr_hal9001$new()

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
