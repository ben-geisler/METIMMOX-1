# ===============================================================================
# EVIDENCE UPDATING OF THE JOINT COEFFICIENT DISTRIBUTION (issue #202, paper 3 step 3)
# ===============================================================================
# The survival coefficient draws are multivariate normal on flexsurv's
# optimisation scale (`mvn_joint_v2`, #156/#159): theta ~ N(m, S), with m the
# OS and PFS maximum-likelihood estimates and S the recoloured joint
# covariance. Evidence that is normal in a linear function of theta updates
# (m, S) in closed form (conjugate_update(), the Kalman form):
#
#   K = S A' (A S A' + V)^-1,  m' = m + K (y - A m),  S' = S - K (A S A' + V) K'.
#
# Two kinds of evidence use it:
#   * Sceptical priors (Spiegelhalter, Freedman & Parmar 1994) on single
#     coefficients: a zero-centred normal pseudo-observation of the
#     coefficient, SD log(2) / qnorm(0.95) so that P(rate ratio < 0.5) = 5% on
#     the gamma rate scale (sceptical_prior_sd()); the ridge variant keeps the
#     centre at 0 and calibrates the SD to a target shrinkage
#     (calibrate_sceptical_sd()), never centring on the ridge estimate, which
#     would use METIMMOX's data twice.
#   * External trials: trial k reports a Cox log hazard ratio y_k with
#     standard error s_k. The model is not proportional hazards, so the
#     summary a draw implies is the limit of the misspecified Cox estimator
#     (Struthers & Kalbfleisch 1986), g_k(theta): the root of the limiting
#     partial-likelihood score for trial k's two arms, allocation, accrual and
#     follow-up (cox_limit_from_curves()), with both arms averaged over
#     METIMMOX's observed joint biomarker mix held fixed (per-draw Dirichlet
#     weights would break the conjugate form) and, for a maintenance design,
#     conditioned on being progression-free at the randomisation landmark
#     (make_trial_summary_map()). The likelihood is
#     y ~ N(g(theta), V), V = diag(s_k^2) + tau_shared^2 (same endpoint, all
#     trials; a floor that more trials cannot average away) + diag(tau_k^2)
#     (evidence_covariance()). g is linearised at the current mean by central
#     differences, the conjugate update applied, and the linearisation
#     repeated at the updated mean until the mean moves by less than `tol`
#     (iterated_conjugate_update()).
#
# Importance weighting is the check, not the engine (decisions.md, 10 October
# 2026): weights proportional to N(y; g(theta_i), V) with the exact g on the
# original draws give a self-normalised estimate of every posterior mean
# (importance_weight_check()), used in its control-variate form (weighted
# mean minus the draws' own deviation from the exact prior mean, which the
# conjugate update does not see); it must agree with the conjugate mean within
# `tol_mcse` Monte Carlo SEs, else the step stops. With p coefficients a
# per-coefficient two-SE rule stops a correct update with probability
# 1 - 0.95^p (60% for the 18 coefficients), so the threshold is
# Bonferroni-adjusted by default (adjust = "bonferroni": two-sided level
# 2 * pnorm(-tol_mcse) shared over the p coefficients); adjust = "none" applies
# the literal two-SE rule. The effective sample size is reported.
#
# Which external rows enter (ruling (2) of 11 October 2026) is decided by
# select_evidence_rows(), which build_evidence_state() applies by default: one
# row per trial and endpoint, the pMMR/MSS population over the ITT row (the
# ITT row's dMMR share recorded otherwise), the latest data cut, the 95%
# interval where the same analysis has two levels; single-arm or
# non-randomised designs (#214), interaction P values, non-HR measures, zero
# lower limits and secondary duplicates excluded;
# biomarker subgroups kept apart (subgroups = TRUE) and matched against
# METIMMOX's definitions (evidence_subgroup_match(): none matches today).
#
# build_evidence_state() combines the steps and regenerates the draws with the
# existing sampler from the posterior (regenerate_evidence_component():
# build_sampled_flexsurv_model() for the posterior-centred fits,
# draw_joint_survival_coefficients() for the draws), so the result has the
# sampling-cache layout (method "mvn_joint_v2") and sampled_survival_models()
# rebuilds models from it. Its fingerprint and cache path come from
# R/voi_cache.R; the base sampling cache is never overwritten.
#
# Requires R/pfs_os_violation_diagnostics.R (flexsurv_design_matrix()),
# R/joint_survival_sampling.R (draw_joint_survival_coefficients()),
# R/prediction_functions.R (build_sampled_flexsurv_model()) and
# R/design_information.R (censoring_survival(), design_cell_codes()).
# ===============================================================================

#' Weeks per month, the convention of R/survival_validation.R (12 / 52)
WEEKS_PER_MONTH <- 52 / 12

#' Required columns of the external-evidence file
EXTERNAL_EVIDENCE_COLUMNS <- c("trial", "population", "endpoint", "estimate", "ci_lower", "ci_upper",
                               "ci_level", "allocation_ratio", "randomisation_landmark",
                               "median_follow_up_months")

#' Location of the external-evidence table (#201)
EXTERNAL_EVIDENCE_PATH <- file.path("data", "external", "voi_external_evidence.csv")

#' Closed-form update of a multivariate normal by a normal linear likelihood
#'
#' Rows whose variance is infinite carry no information and are dropped
#' (they must be uncorrelated with the other rows), so an infinitely vague
#' likelihood returns the prior exactly.
#'
#' @param mean Prior mean (named vector, length p).
#' @param cov Prior covariance (p x p).
#' @param A Observation matrix (m x p) or a vector for m = 1.
#' @param y Observed values of A theta (length m).
#' @param V Observation covariance (m x m) or a vector of variances.
#' @return List with `mean`, `cov` (symmetrised), `gain` (K) and `innovation`
#'   (y - A mean), or the prior with `gain = NULL` when nothing is informative.
conjugate_update <- function(mean, cov, A, y, V) {
  A <- if (is.matrix(A)) A else matrix(A, nrow = 1)
  V <- if (is.matrix(V)) V else diag(V, nrow = length(V))
  if (ncol(A) != length(mean) || nrow(cov) != length(mean) || ncol(cov) != length(mean) ||
      length(y) != nrow(A) || !identical(dim(V), c(nrow(A), nrow(A)))) {
    stop("conjugate_update(): incompatible dimensions")
  }
  if (any(is.na(V)) || any(diag(V) < 0)) stop("Observation variances must be non-negative")
  v_inf <- is.infinite(diag(V))
  if (any(v_inf)) {
    if (any(is.infinite(V[v_inf, !v_inf])) || any(V[v_inf, !v_inf] != 0)) {
      stop("Rows with infinite variance must be uncorrelated with the other rows")
    }
    A <- A[!v_inf, , drop = FALSE]; y <- y[!v_inf]; V <- V[!v_inf, !v_inf, drop = FALSE]
  }
  if (nrow(A) == 0) return(list(mean = mean, cov = cov, gain = NULL, innovation = numeric(0)))
  if (any(!is.finite(A)) || any(!is.finite(y)) || any(!is.finite(V))) stop("Non-finite evidence")
  m_s <- A %*% cov %*% t(A) + V
  m_s <- (m_s + t(m_s)) / 2
  m_gain <- cov %*% t(A) %*% solve(m_s)
  v_innov <- as.vector(y - A %*% mean)
  v_mean <- as.vector(mean + m_gain %*% v_innov)
  names(v_mean) <- names(mean)
  m_cov <- cov - m_gain %*% m_s %*% t(m_gain)
  m_cov <- (m_cov + t(m_cov)) / 2
  dimnames(m_cov) <- dimnames(cov)
  list(mean = v_mean, cov = m_cov, gain = m_gain, innovation = v_innov)
}

#' SD of a sceptical normal prior on a log rate ratio
#'
#' @param rate_ratio Rate ratio whose lower tail probability is fixed (0.5).
#' @param probability Prior probability of a rate ratio below `rate_ratio`.
#' @return -log(rate_ratio) / qnorm(1 - probability); log(2) / qnorm(0.95)
#'   for the defaults.
sceptical_prior_sd <- function(rate_ratio = 0.5, probability = 0.05) {
  if (rate_ratio <= 0 || rate_ratio >= 1 || probability <= 0 || probability >= 0.5) {
    stop("rate_ratio must lie in (0, 1) and probability in (0, 0.5)")
  }
  -log(rate_ratio) / stats::qnorm(1 - probability)
}

#' Zero-centred sceptical prior terms as pseudo-observations
#'
#' @param coef_names Joint coefficient names receiving the prior.
#' @param all_names All joint coefficient names (column order).
#' @param sd Prior SD (scalar or one per coefficient).
#' @return List with `A` (one unit row per coefficient), `y` (zeros), `V`
#'   (variances) and `label`.
sceptical_prior_terms <- function(coef_names, all_names, sd = sceptical_prior_sd()) {
  v_missing <- setdiff(coef_names, all_names)
  if (length(v_missing)) stop("Unknown coefficients for the sceptical prior: ", paste(v_missing, collapse = ", "))
  sd <- rep_len(sd, length(coef_names))
  if (any(!is.finite(sd)) || any(sd <= 0)) stop("Sceptical prior SDs must be positive and finite")
  m_a <- matrix(0, length(coef_names), length(all_names), dimnames = list(coef_names, all_names))
  m_a[cbind(seq_along(coef_names), match(coef_names, all_names))] <- 1
  list(A = m_a, y = rep(0, length(coef_names)), V = sd^2, label = paste0("sceptical::", coef_names))
}

#' Apply a zero-centred sceptical prior to a normal coefficient distribution
#'
#' @param mean,cov Joint mean and covariance (names "os::<coef>", ...).
#' @param coef_names Coefficients receiving the prior.
#' @param sd Prior SD.
#' @return conjugate_update() result.
apply_sceptical_prior <- function(mean, cov, coef_names, sd = sceptical_prior_sd()) {
  terms <- sceptical_prior_terms(coef_names, names(mean), sd)
  conjugate_update(mean, cov, terms$A, terms$y, terms$V)
}

#' Calibrate a zero-centred prior SD to reproduce a target shrinkage
#'
#' Solves for the common SD of zero-centred priors on `prior_coefs` such that
#' the posterior mean of `target_coef` equals `target_shrinkage` times its
#' prior mean (shrinkage on the coefficient, i.e. log rate ratio, scale). The
#' prior stays centred at 0.
#'
#' @param mean,cov Joint mean and covariance.
#' @param prior_coefs Coefficients receiving the prior.
#' @param target_coef Coefficient whose shrinkage is matched.
#' @param target_shrinkage Ratio posterior / prior mean, in (0, 1).
#' @param interval Search interval for log(SD).
#' @return List with `sd`, `shrinkage` (achieved), `probability_below_half`
#'   (P(rate ratio < 0.5) implied by the SD) and the update.
calibrate_sceptical_sd <- function(mean, cov, prior_coefs, target_coef, target_shrinkage,
                                   interval = log(c(1e-4, 1e4))) {
  if (!is.finite(target_shrinkage) || target_shrinkage <= 0 || target_shrinkage >= 1) {
    stop("target_shrinkage must lie in (0, 1)")
  }
  if (!target_coef %in% names(mean) || mean[[target_coef]] == 0) stop("target_coef must have a non-zero mean")
  shrink <- function(log_sd) {
    apply_sceptical_prior(mean, cov, prior_coefs, exp(log_sd))$mean[[target_coef]] / mean[[target_coef]]
  }
  root <- stats::uniroot(function(x) shrink(x) - target_shrinkage, interval, tol = 1e-10)
  sd <- exp(root$root)
  list(sd = sd, shrinkage = shrink(root$root), probability_below_half = stats::pnorm(log(0.5) / sd),
       update = apply_sceptical_prior(mean, cov, prior_coefs, sd))
}

#' Limit of a Cox estimate under non-proportional hazards (two arms)
#'
#' Root in beta of the limiting partial-likelihood score
#'   U(beta) = sum_j G_j [p1 dF1_j - (p0 dF0_j + p1 dF1_j) E_j(beta)],
#'   E_j(beta) = p1 S1_j e^beta / (p0 S0_j + p1 S1_j e^beta),
#' over quadrature intervals j with exact interval masses dF, survival at the
#' interval midpoints and censoring survival G at the midpoints.
#'
#' @param nodes Quadrature nodes (time since randomisation).
#' @param s_control,s_experimental Population survival at the nodes.
#' @param censor_survival Censoring survival at the interval midpoints.
#' @param p_experimental Allocation probability of the experimental arm.
#' @return The limiting log hazard ratio (experimental vs control).
cox_limit_from_curves <- function(nodes, s_control, s_experimental, censor_survival, p_experimental) {
  n_nodes <- length(nodes)
  if (length(s_control) != n_nodes || length(s_experimental) != n_nodes ||
      length(censor_survival) != n_nodes - 1) stop("cox_limit_from_curves(): length mismatch")
  p0 <- 1 - p_experimental; p1 <- p_experimental
  v_df0 <- s_control[-n_nodes] - s_control[-1]
  v_df1 <- s_experimental[-n_nodes] - s_experimental[-1]
  v_s0 <- (s_control[-n_nodes] + s_control[-1]) / 2
  v_s1 <- (s_experimental[-n_nodes] + s_experimental[-1]) / 2
  v_keep <- censor_survival > 0 & (v_s0 + v_s1) > 0
  v_g <- censor_survival[v_keep]; v_df0 <- v_df0[v_keep]; v_df1 <- v_df1[v_keep]
  v_s0 <- v_s0[v_keep]; v_s1 <- v_s1[v_keep]
  if (sum(v_g * v_df0) <= 0 || sum(v_g * v_df1) <= 0) stop("No expected events in one arm")
  score <- function(beta) {
    v_e <- p1 * v_s1 * exp(beta) / (p0 * v_s0 + p1 * v_s1 * exp(beta))
    sum(v_g * (p1 * v_df1 - (p0 * v_df0 + p1 * v_df1) * v_e))
  }
  stats::uniroot(score, c(-20, 20), tol = 1e-12)$root
}

#' Specification of the summary a trial reports (the map g_k)
#'
#' @param label Row label (e.g. "CheckMate 9X8 PFS").
#' @param endpoint "os" or "pfs".
#' @param p_experimental Allocation probability of the experimental arm.
#' @param min_follow_up_weeks Follow-up after the last patient (F).
#' @param accrual_weeks Uniform accrual duration (A).
#' @param landmark_weeks Randomisation landmark of a maintenance design (0 for
#'   randomisation at the start of first-line treatment).
#' @param cells Cell codes the trial's population maps to (NULL: all).
#' @param step Quadrature spacing (weeks; default the weekly grid).
#' @return List of class "voi_trial_summary_spec".
trial_summary_spec <- function(label, endpoint, p_experimental, min_follow_up_weeks, accrual_weeks = 0,
                               landmark_weeks = 0, cells = NULL, step = 1) {
  if (!endpoint %in% c("os", "pfs")) stop("endpoint must be 'os' or 'pfs'")
  if (!is.finite(p_experimental) || p_experimental <= 0 || p_experimental >= 1) stop("p_experimental must lie in (0, 1)")
  if (!is.finite(landmark_weeks) || landmark_weeks < 0) stop("landmark_weeks must be non-negative")
  structure(list(label = label, endpoint = endpoint, p_experimental = p_experimental,
                 min_follow_up_weeks = min_follow_up_weeks, accrual_weeks = accrual_weeks,
                 landmark_weeks = landmark_weeks, cells = cells, step = step),
            class = "voi_trial_summary_spec")
}

#' Joint coefficient names of a sampling component
#'
#' @param component Joint sampling component.
#' @return Character vector "os::<par>", "pfs::<par>" over the optimised
#'   parameters, in the order of the joint covariance.
joint_coefficient_names <- function(component) {
  c(paste0("os::", rownames(component$original_os$res)[component$original_os$optpars]),
    paste0("pfs::", rownames(component$original_pfs$res)[component$original_pfs$optpars]))
}

#' Joint mean and covariance of a sampling component
#'
#' @param component Joint sampling component with `joint_covariance$covariance`.
#' @return List with `mean` (the fitted estimates) and `cov`, both named by
#'   joint_coefficient_names().
joint_prior_moments <- function(component) {
  v_names <- joint_coefficient_names(component)
  m_cov <- component$joint_covariance$covariance
  if (is.null(m_cov)) stop("The sampling component has no joint covariance (regenerate it with mvn_joint_v2)")
  if (!identical(dim(m_cov), rep(length(v_names), 2))) stop("Joint covariance has the wrong dimension")
  dimnames(m_cov) <- list(v_names, v_names)
  v_mean <- c(component$original_os$opt$par, component$original_pfs$opt$par)
  names(v_mean) <- v_names
  list(mean = v_mean, cov = m_cov)
}

#' Joint draw matrix of a sampling component
#'
#' @param component Joint sampling component with `draws$os`, `draws$pfs`.
#' @return Matrix, one row per draw, columns joint_coefficient_names().
joint_draw_matrix <- function(component) {
  m <- cbind(component$draws$os[, rownames(component$original_os$res)[component$original_os$optpars], drop = FALSE],
             component$draws$pfs[, rownames(component$original_pfs$res)[component$original_pfs$optpars], drop = FALSE])
  colnames(m) <- joint_coefficient_names(component)
  m
}

#' Build the map g_k from joint coefficients to a trial's Cox log hazard ratio
#'
#' @param spec trial_summary_spec().
#' @param component Joint sampling component (fits supply formulas and levels).
#' @param reference_population Reference patients (Age, sex, Rx, biomarkers);
#'   METIMMOX's prediction population in live runs.
#' @param weights Fixed reference weights (METIMMOX's observed joint mix; not
#'   per-draw Dirichlet weights).
#' @param biomarkers Biomarker columns defining the cells.
#' @return Function of a joint coefficient vector (names "os::<par>",
#'   "pfs::<par>") returning g_k(theta).
make_trial_summary_map <- function(spec, component, reference_population, weights = NULL,
                                   biomarkers = design_biomarkers()) {
  if (!inherits(spec, "voi_trial_summary_spec")) stop("spec must come from trial_summary_spec()")
  if (is.null(weights)) weights <- rep(1, nrow(reference_population))
  if (!is.null(spec$cells)) {
    v_keep <- design_cell_codes(reference_population, biomarkers) %in% spec$cells
    if (!any(v_keep)) stop("No reference patients in the cells of ", spec$label)
    reference_population <- reference_population[v_keep, , drop = FALSE]
    weights <- weights[v_keep]
  }
  v_w <- weights / sum(weights)
  fit <- component[[paste0("original_", spec$endpoint)]]
  fit_pfs <- component$original_pfs
  v_levels <- model_rx_levels(fit)
  design_for <- function(model, level) {
    df <- reference_population
    df$Rx <- factor(level, levels = model_rx_levels(model))
    flexsurv_design_matrix(model, df)
  }
  m_x_control <- design_for(fit, v_levels[1])
  m_x_exp <- design_for(fit, v_levels[2])
  m_x_sel <- if (spec$landmark_weeks > 0) design_for(fit_pfs, model_rx_levels(fit_pfs)[1]) else NULL
  v_nodes <- design_time_nodes(spec$min_follow_up_weeks, spec$accrual_weeks, spec$step)
  v_mid <- (v_nodes[-1] + v_nodes[-length(v_nodes)]) / 2
  v_g_mid <- censoring_survival(v_mid, spec$min_follow_up_weeks, spec$accrual_weeks)
  v_names_end <- rownames(fit$res)
  v_names_pfs <- rownames(fit_pfs$res)
  landmark <- spec$landmark_weeks
  endpoint_par <- function(theta, model, outcome) {
    v_full <- stats::setNames(model$res.t[, "est"], rownames(model$res))
    v_key <- paste0(outcome, "::", names(v_full))
    v_have <- v_key %in% names(theta)
    v_full[v_have] <- theta[v_key[v_have]]
    v_full
  }
  population_curve <- function(par, m_x, v_sel_w) {
    k <- exp(par[["shape"]])
    v_rate <- exp(as.vector(par[["rate"]] + m_x %*% par[colnames(m_x)]))
    m_s <- outer(v_rate, landmark + v_nodes, function(r, t) stats::pgamma(t, shape = k, rate = r, lower.tail = FALSE))
    if (landmark > 0) m_s <- m_s / stats::pgamma(landmark, shape = k, rate = v_rate, lower.tail = FALSE)
    as.vector(crossprod(v_sel_w, m_s))
  }
  function(theta) {
    par <- endpoint_par(theta, fit, spec$endpoint)
    v_sel <- v_w
    if (landmark > 0) {
      par_pfs <- endpoint_par(theta, fit_pfs, "pfs")
      v_rate_sel <- exp(as.vector(par_pfs[["rate"]] + m_x_sel %*% par_pfs[colnames(m_x_sel)]))
      v_sel <- v_w * stats::pgamma(landmark, shape = exp(par_pfs[["shape"]]), rate = v_rate_sel, lower.tail = FALSE)
      v_sel <- v_sel / sum(v_sel)
    }
    cox_limit_from_curves(v_nodes, population_curve(par, m_x_control, v_sel),
                          population_curve(par, m_x_exp, v_sel), v_g_mid, spec$p_experimental)
  }
}

#' Central-difference Jacobian
#'
#' @param fn Function of a named numeric vector returning a numeric vector.
#' @param x Point (named).
#' @param rel_step Step relative to max(1, |x_j|).
#' @return Matrix length(fn(x)) x length(x) with column names names(x).
numerical_jacobian <- function(fn, x, rel_step = 1e-5) {
  v_f0 <- fn(x)
  m_j <- matrix(NA_real_, length(v_f0), length(x), dimnames = list(names(v_f0), names(x)))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[[j]]))
    x_p <- x; x_m <- x
    x_p[j] <- x[[j]] + h; x_m[j] <- x[[j]] - h
    m_j[, j] <- (fn(x_p) - fn(x_m)) / (2 * h)
  }
  m_j
}

#' Combine evidence observations into one vector-valued map
#'
#' @param observations List of observations, each list(label, y, fn) for a
#'   nonlinear map or list(label, y, a) for a linear row `a` (named by joint
#'   coefficient).
#' @param all_names Joint coefficient names.
#' @return List with `fn` (theta -> vector), `y`, `labels` and `linear`
#'   (logical per observation) plus `a_linear` (rows of the linear ones).
combine_observations <- function(observations, all_names) {
  if (!length(observations)) stop("No evidence observations")
  v_labels <- vapply(observations, function(o) o$label, character(1))
  if (anyDuplicated(v_labels)) stop("Duplicate evidence labels: ", paste(v_labels[duplicated(v_labels)], collapse = ", "))
  v_linear <- vapply(observations, function(o) !is.null(o$a), logical(1))
  rows <- lapply(observations[v_linear], function(o) {
    v <- stats::setNames(numeric(length(all_names)), all_names)
    if (length(setdiff(names(o$a), all_names))) stop("Linear evidence row on unknown coefficients")
    v[names(o$a)] <- o$a
    v
  })
  m_a <- if (length(rows)) do.call(rbind, rows) else NULL
  list(fn = function(theta) {
         v <- numeric(length(observations))
         for (i in seq_along(observations)) {
           o <- observations[[i]]
           v[i] <- if (is.null(o$a)) o$fn(theta) else sum(o$a * theta[names(o$a)])
         }
         stats::setNames(v, v_labels)
       },
       y = stats::setNames(vapply(observations, function(o) o$y, numeric(1)), v_labels),
       labels = v_labels, linear = v_linear, a_linear = m_a)
}

#' Iterated linearised conjugate update (relinearised at the updated mean)
#'
#' @param mean,cov Prior mean and covariance (joint names).
#' @param observations Evidence (see combine_observations()).
#' @param V Observation covariance (evidence_covariance()).
#' @param tol Convergence tolerance on the maximum absolute change of the mean.
#' @param max_iter Maximum number of linearisations.
#' @param rel_step Jacobian step.
#' @return List with `mean`, `cov`, `jacobian` (at the final linearisation
#'   point), `iterations`, `max_change` (per iteration) and `predicted`
#'   (g at the posterior mean).
iterated_conjugate_update <- function(mean, cov, observations, V, tol = 1e-6, max_iter = 50,
                                      rel_step = 1e-5) {
  obs <- combine_observations(observations, names(mean))
  v_point <- mean
  v_change <- numeric(0)
  for (it in seq_len(max_iter)) {
    v_g <- obs$fn(v_point)
    m_j <- numerical_jacobian(obs$fn, v_point, rel_step)
    if (any(obs$linear)) m_j[obs$linear, ] <- obs$a_linear
    v_y_lin <- obs$y - v_g + as.vector(m_j %*% v_point)
    post <- conjugate_update(mean, cov, m_j, v_y_lin, V)
    v_change[it] <- max(abs(post$mean - v_point))
    v_point <- post$mean
    if (v_change[it] < tol) break
  }
  if (v_change[length(v_change)] >= tol) {
    stop("Linearised update did not converge in ", max_iter, " iterations (last change ",
         signif(v_change[length(v_change)], 3), ")")
  }
  list(mean = post$mean, cov = post$cov, jacobian = m_j, iterations = length(v_change),
       max_change = v_change, predicted = obs$fn(post$mean), observed = obs$y)
}

#' Observation covariance with shared and trial-specific heterogeneity
#'
#' @param df_obs Data frame with `label`, `trial`, `endpoint`, `se` and `type`
#'   ("external" or "prior"); prior rows get variance se^2 only.
#' @param tau_shared Shared heterogeneity SD (common to every external trial
#'   of an endpoint: covariance tau_shared^2 between them).
#' @param tau_trial Trial-specific heterogeneity SD (scalar or per row).
#' @param within_trial_correlation Correlation of the sampling errors of the
#'   OS and PFS results of one trial (paired bootstrap value when available).
#' @return Covariance matrix with dimnames `label`.
evidence_covariance <- function(df_obs, tau_shared = 0, tau_trial = 0, within_trial_correlation = 0) {
  v_ext <- df_obs$type == "external"
  tau_trial <- rep_len(tau_trial, nrow(df_obs))
  m_v <- diag(df_obs$se^2 + ifelse(v_ext, tau_trial^2, 0), nrow = nrow(df_obs))
  for (i in seq_len(nrow(df_obs))) for (j in seq_len(nrow(df_obs))) {
    if (i == j || !v_ext[i] || !v_ext[j]) next
    if (df_obs$endpoint[i] == df_obs$endpoint[j]) m_v[i, j] <- m_v[i, j] + tau_shared^2
    if (df_obs$trial[i] == df_obs$trial[j] && df_obs$endpoint[i] != df_obs$endpoint[j]) {
      m_v[i, j] <- m_v[i, j] + within_trial_correlation * df_obs$se[i] * df_obs$se[j]
    }
    if (df_obs$trial[i] == df_obs$trial[j] && df_obs$endpoint[i] == df_obs$endpoint[j]) {
      stop("Two results of one trial and endpoint: ", df_obs$label[i], ", ", df_obs$label[j])
    }
  }
  diag(m_v) <- diag(m_v) + ifelse(v_ext, tau_shared^2, 0)
  dimnames(m_v) <- list(df_obs$label, df_obs$label)
  m_v
}

#' Importance-weight check of a conjugate (linearised) update
#'
#' Weights w_i proportional to N(y; g(theta_i), V) with the exact maps on the
#' prior draws. The weighted posterior mean is compared with the conjugate
#' mean in its control-variate form,
#'   mu_cv = sum_i w_i theta_i - (mean_i theta_i - m),
#' which removes the sampling error that the prior draws carry before any
#' update (m is the exact prior mean); its MCSE is sd(psi) / sqrt(n) with
#' psi_i = n w_i (theta_i - mu_w) - (theta_i - mean theta), the delta-method
#' influence of the self-normalised estimator. The plain weighted mean and its
#' MCSE sqrt(sum w_i^2 (theta_i - mu_w)^2) are reported alongside.
#'
#' @param draws Joint draw matrix (joint_draw_matrix()) from the prior.
#' @param observations Evidence (exact maps; see combine_observations()).
#' @param V Observation covariance.
#' @param conjugate_mean Posterior mean of the conjugate update.
#' @param prior_mean Exact prior mean the draws were generated from.
#' @param conjugate_cov Posterior covariance (for the SD comparison).
#' @param tol_mcse Tolerance in Monte Carlo SEs (2).
#' @param adjust "bonferroni" (default) or "none" (see the file header).
#' @param stop_on_disagreement Stop when any coefficient disagrees.
#' @return List with `table` (coefficient, prior_mean, conjugate_mean,
#'   weighted_mean, weighted_mcse, cv_mean, cv_mcse, z, conjugate_sd,
#'   weighted_sd, agree), `ess`, `n_draws`, `threshold`, `adjust` and
#'   `log_weights`.
importance_weight_check <- function(draws, observations, V, conjugate_mean, prior_mean, conjugate_cov = NULL,
                                    tol_mcse = 2, adjust = c("bonferroni", "none"),
                                    stop_on_disagreement = TRUE) {
  adjust <- match.arg(adjust)
  v_names <- colnames(draws)
  obs <- combine_observations(observations, v_names)
  m_g <- t(apply(draws, 1, function(theta) obs$fn(stats::setNames(theta, v_names))))
  if (length(obs$y) == 1) m_g <- t(m_g)
  m_r <- sweep(m_g, 2, obs$y, "-")
  v_inf <- is.infinite(diag(V))
  m_vinv <- matrix(0, nrow(V), ncol(V))
  if (any(!v_inf)) m_vinv[!v_inf, !v_inf] <- solve(V[!v_inf, !v_inf, drop = FALSE])
  v_logw <- -0.5 * rowSums((m_r %*% m_vinv) * m_r)
  v_w <- exp(v_logw - max(v_logw))
  v_w <- v_w / sum(v_w)
  n <- nrow(draws)
  v_wmean <- colSums(v_w * draws)
  v_dmean <- colMeans(draws)
  m_c <- sweep(draws, 2, v_wmean, "-")
  v_wmcse <- sqrt(colSums(v_w^2 * m_c^2))
  v_wsd <- sqrt(colSums(v_w * m_c^2))
  m_psi <- n * v_w * m_c - sweep(draws, 2, v_dmean, "-")
  v_cv_mcse <- apply(m_psi, 2, stats::sd) / sqrt(n)
  v_cv_mean <- v_wmean - (v_dmean - prior_mean[v_names])
  v_cmean <- conjugate_mean[v_names]
  p <- ncol(draws)
  threshold <- if (adjust == "bonferroni") stats::qnorm(1 - stats::pnorm(-tol_mcse) / p) else tol_mcse
  # A coefficient untouched by the evidence has zero shift and zero cv MCSE.
  v_z <- ifelse(v_cv_mcse > 0, abs(v_cv_mean - v_cmean) / v_cv_mcse,
                ifelse(abs(v_cv_mean - v_cmean) <= 1e-12 * pmax(1, abs(v_cmean)), 0, Inf))
  df <- data.frame(coefficient = v_names, prior_mean = prior_mean[v_names], conjugate_mean = v_cmean,
                   weighted_mean = v_wmean, weighted_mcse = v_wmcse, cv_mean = v_cv_mean,
                   cv_mcse = v_cv_mcse, z = v_z,
                   conjugate_sd = if (is.null(conjugate_cov)) NA_real_ else sqrt(diag(conjugate_cov))[v_names],
                   weighted_sd = v_wsd, agree = v_z <= threshold, row.names = NULL, stringsAsFactors = FALSE)
  ess <- 1 / sum(v_w^2)
  if (stop_on_disagreement && !all(df$agree)) {
    df_bad <- df[!df$agree, , drop = FALSE]
    stop("Conjugate and importance-weighted posterior means disagree (threshold ", signif(threshold, 3),
         " MCSE; ESS ", round(ess), "): ",
         paste(sprintf("%s conjugate %.5f, weighted %.5f (z = %.2f)", df_bad$coefficient,
                       df_bad$conjugate_mean, df_bad$cv_mean, df_bad$z), collapse = "; "))
  }
  list(table = df, ess = ess, n_draws = n, threshold = threshold, adjust = adjust, log_weights = v_logw)
}


#' First a:b ratio of a text, as the experimental share
#'
#' The external-evidence contract records allocation as experimental:control
#' (e.g. "2:1 (experimental:control)"); "control:experimental" written
#' directly after the ratio swaps the order.
#'
#' @param text Allocation text.
#' @return Experimental allocation probability, NA when no ratio is found.
parse_allocation_ratio <- function(text) {
  vapply(text, function(x) {
    if (is.na(x)) return(NA_real_)
    m <- regmatches(x, regexec("([0-9.]+)\\s*:\\s*([0-9.]+)", x))[[1]]
    if (!length(m)) return(NA_real_)
    a <- as.numeric(m[2]); b <- as.numeric(m[3])
    swap <- grepl("^[^(]*[0-9.]+\\s*:\\s*[0-9.]+\\s*\\(\\s*control\\s*:\\s*experimental", x, ignore.case = TRUE)
    if (swap) b / (a + b) else a / (a + b)
  }, numeric(1), USE.NAMES = FALSE)
}

#' First number of a free-text field
#'
#' @param text Character vector.
#' @return Numeric vector (NA when no number is present).
parse_first_number <- function(text) {
  v <- regmatches(text, regexpr("[0-9]+(\\.[0-9]+)?", text))
  out <- rep(NA_real_, length(text))
  out[!is.na(text) & grepl("[0-9]", text)] <- as.numeric(v)
  out
}

#' Check the external-evidence table contract
#'
#' @param df Data frame read from the evidence file.
#' @return `df` with numeric `estimate`, `ci_lower`, `ci_upper`, `ci_level`;
#'   stops listing every violation (missing columns, non-positive estimates,
#'   an estimate outside its own CI, a CI level outside (0, 100)).
check_external_evidence <- function(df) {
  v_missing <- setdiff(EXTERNAL_EVIDENCE_COLUMNS, names(df))
  if (length(v_missing)) stop("External evidence lacks required column(s): ", paste(v_missing, collapse = ", "))
  for (col in c("estimate", "ci_lower", "ci_upper", "ci_level")) {
    v <- suppressWarnings(as.numeric(df[[col]]))
    if (any(is.na(v) & !is.na(df[[col]]))) stop("Column ", col, " is not numeric")
    df[[col]] <- v
  }
  v_problems <- character(0)
  v_has_ci <- !is.na(df$estimate) & !is.na(df$ci_lower) & !is.na(df$ci_upper)
  v_bad_ci <- v_has_ci & (df$ci_lower > df$estimate | df$ci_upper < df$estimate)
  if (any(v_bad_ci)) v_problems <- c(v_problems, paste("estimate outside its CI in rows", paste(which(v_bad_ci), collapse = ", ")))
  v_bad_level <- !is.na(df$ci_level) & (df$ci_level <= 0 | df$ci_level >= 100)
  if (any(v_bad_level)) v_problems <- c(v_problems, paste("CI level outside (0, 100) in rows", paste(which(v_bad_level), collapse = ", ")))
  if (any(is.na(df$trial) | !nzchar(df$trial))) v_problems <- c(v_problems, "rows without a trial name")
  if (length(v_problems)) stop("External evidence contract violated: ", paste(v_problems, collapse = "; "))
  df
}

#' Read the external-evidence table with its contract checked
#'
#' @param path CSV path (default data/external/voi_external_evidence.csv).
#' @return Checked data frame (check_external_evidence()); attribute `md5`.
read_external_evidence <- function(path = here::here(EXTERNAL_EVIDENCE_PATH)) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("NA", ""), check.names = FALSE)
  df <- check_external_evidence(df)
  attr(df, "md5") <- unname(tools::md5sum(path))
  df
}

#' Hazard-ratio rows of the evidence table as normal likelihood terms
#'
#' Keeps hazard ratios on OS or PFS with an estimate and a CI; converts the CI
#' at its own level (z = qnorm(1 - (1 - level) / 2); 1.2816 for 80%) to the
#' SE of the log hazard ratio. Rows whose lower limit prints as 0 (infinite
#' SE) are skipped too, as are rows whose `design` (when the column exists) is
#' single-arm or non-randomised (EVIDENCE_NONRANDOMISED_PATTERN, #214).
#' Skipped rows are reported with a message. It returns
#' every usable row, including overlapping ones; select_evidence_rows()
#' chooses one per trial and endpoint.
#'
#' @param df check_external_evidence() output.
#' @param effect_measures Accepted `effect_measure` values (when the column
#'   exists).
#' @return Data frame with `label`, `trial`, `population`, `endpoint` ("os",
#'   "pfs"), `log_hr`, `se`, `z`, `p_experimental`, `follow_up_weeks`,
#'   `maintenance` (randomisation landmark recorded) and the source row.
external_evidence_likelihood <- function(df, effect_measures = "HR") {
  v_row <- seq_len(nrow(df))
  v_skip <- character(0)
  v_endpoint <- ifelse(grepl("^\\s*OS", df$endpoint), "os", ifelse(grepl("^\\s*PFS", df$endpoint), "pfs", NA))
  v_keep <- rep(TRUE, nrow(df))
  if ("effect_measure" %in% names(df)) {
    v_measure <- !is.na(df$effect_measure) & df$effect_measure %in% effect_measures
    v_skip <- c(v_skip, sprintf("row %d (%s): effect measure '%s'", v_row[!v_measure], df$trial[!v_measure],
                                df$effect_measure[!v_measure]))
    v_keep <- v_keep & v_measure
  }
  if ("design" %in% names(df)) {
    v_nonrandom <- v_keep & !is.na(df$design) & grepl(EVIDENCE_NONRANDOMISED_PATTERN, df$design, ignore.case = TRUE)
    v_skip <- c(v_skip, sprintf("row %d (%s, %s): single-arm or non-randomised design", v_row[v_nonrandom],
                                df$trial[v_nonrandom], df$endpoint[v_nonrandom]))
    v_keep <- v_keep & !v_nonrandom
  }
  v_na <- v_keep & is.na(df$estimate)
  v_skip <- c(v_skip, sprintf("row %d (%s, %s): no estimate", v_row[v_na], df$trial[v_na], df$endpoint[v_na]))
  v_keep <- v_keep & !v_na
  v_noci <- v_keep & (is.na(df$ci_lower) | is.na(df$ci_upper) | is.na(df$ci_level))
  v_skip <- c(v_skip, sprintf("row %d (%s, %s): no confidence interval", v_row[v_noci], df$trial[v_noci], df$endpoint[v_noci]))
  v_keep <- v_keep & !v_noci
  v_zero <- v_keep & df$ci_lower <= 0
  v_skip <- c(v_skip, sprintf("row %d (%s, %s): lower limit prints as 0 (SE not finite)", v_row[v_zero],
                              df$trial[v_zero], df$endpoint[v_zero]))
  v_keep <- v_keep & !v_zero
  v_noend <- v_keep & is.na(v_endpoint)
  v_skip <- c(v_skip, sprintf("row %d (%s): endpoint '%s' is neither OS nor PFS", v_row[v_noend], df$trial[v_noend], df$endpoint[v_noend]))
  v_keep <- v_keep & !v_noend
  if (length(v_skip)) message("External evidence rows skipped:\n  ", paste(v_skip, collapse = "\n  "))
  d <- df[v_keep, , drop = FALSE]
  v_level <- ifelse(d$ci_level > 1, d$ci_level / 100, d$ci_level)
  v_z <- stats::qnorm(1 - (1 - v_level) / 2)
  data.frame(label = paste(d$trial, d$population, d$endpoint, sep = " | "), trial = d$trial,
             population = d$population, endpoint = v_endpoint[v_keep], log_hr = log(d$estimate),
             se = (log(d$ci_upper) - log(d$ci_lower)) / (2 * v_z), z = v_z,
             p_experimental = parse_allocation_ratio(d$allocation_ratio),
             follow_up_weeks = parse_first_number(d$median_follow_up_months) * WEEKS_PER_MONTH,
             maintenance = !is.na(d$randomisation_landmark), source_row = v_row[v_keep],
             stringsAsFactors = FALSE)
}

# -------------------------------------------------------------------------------
# Selection of the evidence rows (ruling (2) of 11 October 2026)
# -------------------------------------------------------------------------------

#' TMB cut-off of METIMMOX's TMB/BRAF biomarker (mutations per megabase,
#' baseline NGS panel; TMB >= 9 or any BRAF mutation is positive)
METIMMOX_TMB_CUTOFF <- 9
#' Default randomisation landmark of maintenance trials in the evidence table
#' (MODUL cohort 2: 8 two-weekly induction cycles, 16 weeks). ASSUMPTION; the
#' VoI driver passes E2_MAINTENANCE_LANDMARK_WEEKS (R/voi_research_value.R).
EVIDENCE_MAINTENANCE_LANDMARK_WEEKS <- 16
#' Markers of a biomarker subgroup in `subgroup_definition` (any of them makes
#' the row a subgroup row) and in `population` (RAS and BRAF are left out
#' there because a population label can describe trial eligibility, e.g.
#' MODUL cohort 2's "BRAF V600 wild-type" ITT population)
EVIDENCE_SUBGROUP_MARKERS <- c(tmb = "TMB", immunoscore = "Immunoscore", determaio = "DetermaIO",
                               braf = "BRAF", ras = "\\bK?RAS\\b", pdl1 = "PD-L1", cms = "\\bCMS",
                               liver = "liver")
#' Population markers of a pMMR/MSS population
EVIDENCE_PMMR_PATTERN <- "pMMR|\\bMSS\\b|proficient mismatch|microsatellite stable"
#' Markers of a design without a randomised contrast in the `design` column
#' (single-arm studies and non-randomised or cross-trial comparisons, #214);
#' "randomised" alone does not match
EVIDENCE_NONRANDOMISED_PATTERN <- "single[- ]arm|non-?randomi[sz]ed|uncontrolled|cross-trial"

#' Biomarker-subgroup kind of evidence rows
#'
#' @param subgroup_definition,population Character vectors.
#' @return Character vector: the first matching name of
#'   EVIDENCE_SUBGROUP_MARKERS ("tmb", "immunoscore", "determaio", "braf",
#'   "ras", "pdl1", "cms", "liver"), NA for none.
evidence_subgroup_kind <- function(subgroup_definition, population) {
  v_pop_markers <- EVIDENCE_SUBGROUP_MARKERS[setdiff(names(EVIDENCE_SUBGROUP_MARKERS), c("braf", "ras"))]
  vapply(seq_along(population), function(i) {
    for (k in names(EVIDENCE_SUBGROUP_MARKERS)) {
      in_def <- !is.na(subgroup_definition[i]) && grepl(EVIDENCE_SUBGROUP_MARKERS[[k]], subgroup_definition[i])
      in_pop <- k %in% names(v_pop_markers) && !is.na(population[i]) &&
        grepl(v_pop_markers[[k]], population[i])
      if (in_def || in_pop) return(k)
    }
    NA_character_
  }, character(1))
}

#' Population class of evidence rows
#'
#' @param subgroup_definition,population Character vectors.
#' @return "pmmr" (pMMR/MSS population without a further biomarker
#'   restriction), "itt" (ITT or mITT without a subgroup definition),
#'   "subgroup" (a biomarker subgroup, evidence_subgroup_kind()) or "other".
evidence_population_class <- function(subgroup_definition, population) {
  v_kind <- evidence_subgroup_kind(subgroup_definition, population)
  v_pmmr <- grepl(EVIDENCE_PMMR_PATTERN, population, ignore.case = TRUE)
  v_itt <- grepl("^\\s*m?ITT", population) & is.na(subgroup_definition)
  ifelse(!is.na(v_kind), "subgroup", ifelse(v_pmmr, "pmmr", ifelse(v_itt, "itt", "other")))
}

#' Data cut of evidence rows (a leading ISO date only)
#'
#' A date later in the text ("not stated (data collected after primary
#' completion 2021-02-01)") is not the cut and gives NA.
#'
#' @param text `data_cut` column.
#' @return Date vector.
parse_data_cut <- function(text) {
  v <- ifelse(!is.na(text) & grepl("^\\s*[0-9]{4}-[0-9]{2}-[0-9]{2}", text),
              sub("^\\s*([0-9]{4}-[0-9]{2}-[0-9]{2}).*$", "\\1", text), NA_character_)
  as.Date(v)
}

#' Numeric dMMR (MSI) share from the free-text dMMR field
#'
#' Reads the earliest "a/b" (no spaces) or "a <word> of b" count; a count
#' split across arms ("4 (7.0%) serplulimab / 0 placebo") gives NA.
#'
#' @param text `dmmr_share_in_population` column.
#' @return Numeric share, NA when no count is given.
parse_dmmr_share <- function(text) {
  vapply(text, function(x) {
    if (is.na(x)) return(NA_real_)
    v_patterns <- c("([0-9]+)/([0-9]+)", "([0-9]+)\\s+[A-Za-z-]+\\s+of\\s+([0-9]+)")
    v_pos <- vapply(v_patterns, function(p) regexpr(p, x)[1], numeric(1))
    if (all(v_pos < 0)) return(NA_real_)
    m <- regmatches(x, regexec(v_patterns[which(v_pos == min(v_pos[v_pos > 0]))[1]], x))[[1]]
    as.numeric(m[2]) / as.numeric(m[3])
  }, numeric(1), USE.NAMES = FALSE)
}

#' Does a biomarker-subgroup row match METIMMOX's biomarker definition?
#'
#' TMB-high and BRAF results enter the biomarker-specific mapping only where
#' the definition matches METIMMOX's (ruling (2) of 11 October 2026): the
#' TMB/BRAF-positive group is TMB >= METIMMOX_TMB_CUTOFF mut/Mb (baseline NGS
#' panel) or any BRAF mutation, and week-4 CRP is the other biomarker. A TMB
#' row matches only when it states a ">= 9 mut/Mb" cut-off; a BRAF row only
#' when it reports the combined TMB-high-or-BRAF group. With the current
#' table no row matches (AtezoTRIBE uses 10 mut/Mb on FoundationOne CDx,
#' CheckMate 9X8 whole-exome tertiles; BRAF-mutated alone is a subset of the
#' TMB/BRAF-positive group; no trial reports CRP).
#'
#' @param kind evidence_subgroup_kind() of the rows.
#' @param subgroup_definition,population Character vectors.
#' @param tmb_cutoff METIMMOX's TMB cut-off (mut/Mb).
#' @return Data frame with `matches_metimmox` (logical) and `match_reason`.
evidence_subgroup_match <- function(kind, subgroup_definition, population, tmb_cutoff = METIMMOX_TMB_CUTOFF) {
  v_text <- ifelse(is.na(subgroup_definition), population, subgroup_definition)
  rows <- lapply(seq_along(kind), function(i) {
    x <- v_text[i]
    k <- kind[i]
    if (is.na(k)) return(list(FALSE, "not a biomarker subgroup"))
    if (k == "tmb") {
      if (grepl("exome", x, ignore.case = TRUE)) {
        return(list(FALSE, paste0("TMB by whole-exome study-specific tertiles (mutations/exome), not METIMMOX's ",
                                  tmb_cutoff, " mut/Mb panel cut-off")))
      }
      v_cut <- suppressWarnings(as.numeric(sub("^.*?([0-9]+(\\.[0-9]+)?)\\s*(mut|mutations)/Mb.*$", "\\1", x,
                                             perl = TRUE)))
      if (!grepl("mut(ations)?/Mb", x) || !is.finite(v_cut)) {
        return(list(FALSE, "TMB cut-off not stated"))
      }
      if (v_cut != tmb_cutoff) {
        panel <- if (grepl("FoundationOne", x)) " (FoundationOne CDx)" else ""
        return(list(FALSE, sprintf("TMB cut-off %s mut/Mb%s, not METIMMOX's %s mut/Mb", format(v_cut), panel,
                                   format(tmb_cutoff))))
      }
      if (grepl("low|<", x)) return(list(FALSE, "TMB-low row; METIMMOX's TMB/BRAF-negative group also requires BRAF wild type"))
      if (!grepl("BRAF", x)) {
        return(list(FALSE, "TMB-high alone is a subset of METIMMOX's TMB/BRAF-positive group (TMB-high or any BRAF mutation)"))
      }
      return(list(TRUE, "TMB-high or BRAF-mutated at METIMMOX's cut-off"))
    }
    if (k == "braf") {
      if (grepl("TMB", x) && grepl(paste0(">=\\s*", tmb_cutoff, "\\s*mut"), x)) {
        return(list(TRUE, "TMB-high or BRAF-mutated at METIMMOX's cut-off"))
      }
      if (grepl("wild", x, ignore.case = TRUE)) {
        return(list(FALSE, "BRAF wild type mixes TMB/BRAF-positive (TMB-high) and TMB/BRAF-negative patients"))
      }
      return(list(FALSE, paste0("BRAF-mutated alone is a subset of METIMMOX's TMB/BRAF-positive group (TMB >= ",
                                tmb_cutoff, " mut/Mb or any BRAF mutation); the combined group is not reported")))
    }
    list(FALSE, paste0(c(immunoscore = "Immunoscore IC", determaio = "DetermaIO signature", ras = "RAS/KRAS status",
                         pdl1 = "PD-L1 expression", cms = "consensus molecular subtype",
                         liver = "liver metastases")[[k]], " is not a METIMMOX biomarker"))
  })
  data.frame(matches_metimmox = vapply(rows, `[[`, logical(1), 1),
             match_reason = vapply(rows, `[[`, character(1), 2), stringsAsFactors = FALSE)
}

#' Select one external result per trial, population and endpoint
#'
#' The named rule of ruling (2) of 11 October 2026, applied in this order:
#'
#' 1. Rows that cannot enter the likelihood are excluded: rows whose `design`
#'    (when the column exists) is single-arm or non-randomised
#'    (EVIDENCE_NONRANDOMISED_PATTERN; no randomised contrast for g_k, so this
#'    reason comes first even when the row carries a hazard ratio);
#'    interaction P-value rows; effect measures other than "HR" (ORR, medians, the stratification-
#'    adjusted "HR (adjusted)" secondary analysis); endpoints other than OS or
#'    PFS; no estimate or no confidence interval; a lower limit that prints
#'    as 0 or an SE that is not finite; rows marked as secondary duplicates
#'    (`read_from` "secondary ..." or a note saying "do not use" or "does not
#'    match"); rows without the median follow-up or allocation the Cox-limit
#'    map g_k needs; maintenance rows when no landmark is given.
#' 2. Biomarker subgroups (TMB, Immunoscore, RAS/KRAS, BRAF, PD-L1, DetermaIO,
#'    CMS, liver metastases) are excluded from the all-comers update and, with
#'    `subgroups = TRUE`, returned separately with evidence_subgroup_match()
#'    (by default none matches METIMMOX's definitions).
#' 3. Per trial and endpoint the pMMR/MSS population is preferred over the
#'    ITT population; the ITT rows of that trial and endpoint are then not
#'    used. When an ITT row is used its dMMR share is recorded.
#' 4. Within the chosen population: the latest data cut (rows whose cut is
#'    not a date rank after dated rows), then the 95% interval where the same
#'    analysis is reported at two levels, then the full text over a
#'    registry, abstract or supplement, then table order.
#'
#' @param df Raw evidence table (read_external_evidence()) or
#'   external_evidence_likelihood() rows (columns the rule needs but the
#'   input lacks are taken as missing: effect measure "HR", cut and CI level
#'   unknown).
#' @param subgroups Also return the biomarker-subgroup rows.
#' @param maintenance_landmark_weeks Landmark of maintenance designs
#'   (`landmark_weeks` of the selected rows; NA excludes maintenance rows).
#' @param tmb_cutoff METIMMOX's TMB cut-off (mut/Mb).
#' @return List with `selected` (external_evidence_likelihood() columns plus
#'   `population_class`, `data_cut`, `ci_level`, `landmark_weeks`,
#'   `dmmr_share`, `dmmr_share_text` (ITT rows) and `selection_reason`; one
#'   row per trial and endpoint), `excluded` (every other row with
#'   `source_row`, `trial`, `population`, `endpoint`, `exclusion_reason`),
#'   `subgroups` (NULL, or the subgroup rows with `subgroup_kind`,
#'   `matches_metimmox`, `match_reason`, `usable`) and `diagnostics` (the
#'   selected ITT rows with their dMMR share).
select_evidence_rows <- function(df, subgroups = FALSE,
                                 maintenance_landmark_weeks = EVIDENCE_MAINTENANCE_LANDMARK_WEEKS,
                                 tmb_cutoff = METIMMOX_TMB_CUTOFF) {
  likelihood_input <- "log_hr" %in% names(df)
  column_or <- function(name, default = NA) if (name %in% names(df)) df[[name]] else rep(default, nrow(df))
  v_row <- if (likelihood_input && "source_row" %in% names(df)) df$source_row else seq_len(nrow(df))
  d <- data.frame(row = seq_len(nrow(df)), source_row = v_row, trial = df$trial, population = df$population,
                  subgroup_definition = as.character(column_or("subgroup_definition")),
                  endpoint_text = as.character(df$endpoint), stringsAsFactors = FALSE)
  d$effect_measure <- as.character(column_or("effect_measure", "HR"))
  d$read_from <- as.character(column_or("read_from"))
  d$note <- as.character(column_or("note"))
  d$design <- as.character(column_or("design"))
  d$data_cut <- parse_data_cut(as.character(column_or("data_cut")))
  d$ci_level <- suppressWarnings(as.numeric(column_or("ci_level")))
  d$dmmr_share_text <- as.character(column_or("dmmr_share_in_population"))
  if (likelihood_input) {
    d$endpoint <- df$endpoint
    d$log_hr <- df$log_hr
    d$se <- df$se
    d$follow_up_weeks <- df$follow_up_weeks
    d$p_experimental <- df$p_experimental
    d$maintenance <- df$maintenance
    d$has_estimate <- is.finite(d$log_hr)
    d$has_ci <- !is.na(d$se)
  } else {
    d$endpoint <- ifelse(grepl("^\\s*OS", df$endpoint), "os", ifelse(grepl("^\\s*PFS", df$endpoint), "pfs", NA))
    v_level <- ifelse(d$ci_level > 1, d$ci_level / 100, d$ci_level)
    d$log_hr <- suppressWarnings(log(df$estimate))
    d$se <- suppressWarnings((log(df$ci_upper) - log(df$ci_lower)) / (2 * stats::qnorm(1 - (1 - v_level) / 2)))
    d$follow_up_weeks <- parse_first_number(df$median_follow_up_months) * WEEKS_PER_MONTH
    d$p_experimental <- parse_allocation_ratio(df$allocation_ratio)
    d$maintenance <- !is.na(df$randomisation_landmark)
    d$has_estimate <- !is.na(df$estimate)
    d$has_ci <- !is.na(df$ci_lower) & !is.na(df$ci_upper) & !is.na(df$ci_level)
  }
  d$population_class <- evidence_population_class(d$subgroup_definition, d$population)
  d$subgroup_kind <- evidence_subgroup_kind(d$subgroup_definition, d$population)

  # Step 1: rows that cannot enter the likelihood.
  v_reason <- rep(NA_character_, nrow(d))
  set_reason <- function(v_cond, text) {
    v_new <- is.na(v_reason) & v_cond
    v_reason[v_new] <<- if (length(text) == 1) text else text[v_new]
  }
  set_reason(!is.na(d$design) & grepl(EVIDENCE_NONRANDOMISED_PATTERN, d$design, ignore.case = TRUE),
             paste0("single-arm or non-randomised design (", d$design,
                    "): no randomised contrast for the Cox-limit map g_k"))
  set_reason(grepl("interaction", d$effect_measure, ignore.case = TRUE), "interaction P value, not an effect estimate")
  set_reason(d$effect_measure != "HR" & grepl("^HR", d$effect_measure),
             paste0("effect measure '", d$effect_measure, "': secondary analysis of the HR row"))
  set_reason(is.na(d$effect_measure) | d$effect_measure != "HR",
             paste0("effect measure '", d$effect_measure, "' is not a hazard ratio"))
  set_reason(is.na(d$endpoint), paste0("endpoint '", d$endpoint_text, "' is neither OS nor PFS"))
  set_reason(!d$has_estimate, "no estimate")
  set_reason(!d$has_ci, "no confidence interval")
  set_reason(!is.finite(d$se) | d$se <= 0, "lower limit prints as 0 (SE not finite)")
  v_secondary <- (!is.na(d$read_from) & grepl("^secondary", d$read_from)) |
    (!is.na(d$note) & grepl("do not use|does not match|duplicate of the primary", d$note, ignore.case = TRUE))
  set_reason(v_secondary, "marked as a secondary duplicate in the evidence table")
  set_reason(!is.finite(d$follow_up_weeks), "no median follow-up (needed for the Cox-limit map g_k)")
  set_reason(!is.finite(d$p_experimental) | d$p_experimental <= 0 | d$p_experimental >= 1, "no allocation ratio")
  if (is.na(maintenance_landmark_weeks)) set_reason(d$maintenance, "maintenance design and no landmark given")
  v_usable <- is.na(v_reason)

  # Step 2: biomarker subgroups leave the all-comers update.
  set_reason(d$population_class == "subgroup",
             paste0("biomarker subgroup (", d$subgroup_kind, "); not used in the all-comers update"))
  set_reason(d$population_class == "other", "population is neither ITT nor pMMR/MSS")

  # Steps 3 and 4: one row per trial and endpoint.
  v_candidate <- which(is.na(v_reason))
  v_selected <- integer(0)
  v_selection <- character(0)
  v_key <- paste(d$trial, d$endpoint, sep = "|")
  for (key in unique(v_key[v_candidate])) {
    v_idx <- v_candidate[v_key[v_candidate] == key]
    v_pmmr <- v_idx[d$population_class[v_idx] == "pmmr"]
    v_itt <- integer(0)
    if (length(v_pmmr)) {
      v_itt <- setdiff(v_idx, v_pmmr)
      v_idx <- v_pmmr
    }
    v_cut <- as.numeric(d$data_cut[v_idx])
    v_full <- !is.na(d$read_from[v_idx]) & d$read_from[v_idx] == "full text"
    v_order <- v_idx[order(is.na(v_cut), -ifelse(is.na(v_cut), 0, v_cut), !(d$ci_level[v_idx] %in% 95),
                           !v_full, d$row[v_idx])]
    best <- v_order[1]
    v_reason[v_itt] <- paste0("ITT row not used: the pMMR/MSS row ", d$source_row[best],
                              " is used for this trial and endpoint")
    for (j in v_order[-1]) {
      v_reason[j] <- if (is.na(d$data_cut[j]) && !is.na(d$data_cut[best])) {
        paste0("data cut not stated; row ", d$source_row[best], " has the dated cut ", d$data_cut[best])
      } else if (!is.na(d$data_cut[j]) && !is.na(d$data_cut[best]) && d$data_cut[j] < d$data_cut[best]) {
        paste0("earlier data cut (", d$data_cut[j], "); row ", d$source_row[best], " has the latest (",
               d$data_cut[best], ")")
      } else if (!identical(d$ci_level[j], d$ci_level[best]) && identical(d$ci_level[best], 95)) {
        paste0(d$ci_level[j], "% interval; row ", d$source_row[best], " reports the same analysis with a 95% interval")
      } else {
        paste0("row ", d$source_row[best], " preferred (", if (v_full[match(best, v_idx)]) "full text" else
          "earlier table row", ")")
      }
    }
    pop <- if (d$population_class[best] == "pmmr") "pMMR/MSS population" else
      "ITT population (no usable pMMR/MSS row)"
    v_selected <- c(v_selected, best)
    v_selection <- c(v_selection, paste0(pop, "; ", if (is.na(d$data_cut[best])) "data cut not stated" else
      paste("data cut", d$data_cut[best]), "; ", d$ci_level[best], "% interval"))
  }

  # Output in the likelihood format.
  if (likelihood_input) {
    sel <- df[v_selected, , drop = FALSE]
  } else {
    sel <- suppressMessages(external_evidence_likelihood(df[v_selected, , drop = FALSE]))
    sel$source_row <- d$source_row[v_selected]
  }
  sel$population_class <- d$population_class[v_selected]
  sel$data_cut <- d$data_cut[v_selected]
  sel$ci_level <- d$ci_level[v_selected]
  v_lm <- if ("landmark_weeks" %in% names(sel)) sel$landmark_weeks else rep(NA_real_, nrow(sel))
  sel$landmark_weeks <- ifelse(sel$maintenance, ifelse(is.na(v_lm), maintenance_landmark_weeks, v_lm), 0)
  sel$dmmr_share_text <- ifelse(sel$population_class == "itt", d$dmmr_share_text[v_selected], NA_character_)
  sel$dmmr_share <- ifelse(sel$population_class == "itt", parse_dmmr_share(d$dmmr_share_text[v_selected]), NA_real_)
  sel$selection_reason <- v_selection
  sel <- sel[order(sel$trial, sel$endpoint), , drop = FALSE]
  rownames(sel) <- NULL
  if (anyDuplicated(paste(sel$trial, sel$endpoint))) stop("select_evidence_rows(): two rows for one trial and endpoint")

  v_excl <- setdiff(seq_len(nrow(d)), v_selected)
  excluded <- data.frame(source_row = d$source_row[v_excl], trial = d$trial[v_excl],
                         population = d$population[v_excl], endpoint = d$endpoint_text[v_excl],
                         exclusion_reason = v_reason[v_excl], stringsAsFactors = FALSE)
  if (anyNA(excluded$exclusion_reason)) stop("select_evidence_rows(): an excluded row has no reason")

  out_subgroups <- NULL
  if (subgroups) {
    v_sub <- which(d$population_class == "subgroup")
    match <- evidence_subgroup_match(d$subgroup_kind[v_sub], d$subgroup_definition[v_sub], d$population[v_sub],
                                     tmb_cutoff)
    out_subgroups <- data.frame(source_row = d$source_row[v_sub], trial = d$trial[v_sub],
                                population = d$population[v_sub], subgroup_definition = d$subgroup_definition[v_sub],
                                endpoint = d$endpoint[v_sub], log_hr = d$log_hr[v_sub], se = d$se[v_sub],
                                subgroup_kind = d$subgroup_kind[v_sub], usable = v_usable[v_sub],
                                unusable_reason = ifelse(v_usable[v_sub], NA_character_, v_reason[v_sub]),
                                match, stringsAsFactors = FALSE)
    rownames(out_subgroups) <- NULL
  }
  list(selected = sel, excluded = excluded, subgroups = out_subgroups,
       diagnostics = sel[sel$population_class == "itt", c("trial", "endpoint", "source_row", "dmmr_share",
                                                          "dmmr_share_text"), drop = FALSE])
}

#' Regenerate a sampling component from a posterior mean and covariance
#'
#' @param component Base joint sampling component (mvn_joint_v2).
#' @param posterior_mean,posterior_cov Posterior moments (joint names).
#' @param n_samples Number of draws (default: as the base component).
#' @param seed RNG seed (default: the base component's seed).
#' @param evidence_state List describing the state (stored as is).
#' @return Component in the sampling-cache layout: `method = "mvn_joint_v2"`,
#'   posterior-centred `original_os`/`original_pfs` (estimates and `cov`
#'   replaced), `draws`, `joint_covariance` (posterior covariance, prior
#'   covariance and source), `coefficient_summary`, `n_samples`, `n_failed`,
#'   distributions, `seed`, `evidence_state`, `creation_time`.
regenerate_evidence_component <- function(component, posterior_mean, posterior_cov,
                                          n_samples = component$n_samples, seed = component$seed,
                                          evidence_state = list()) {
  v_names <- joint_coefficient_names(component)
  if (!identical(names(posterior_mean), v_names) || !identical(rownames(posterior_cov), v_names)) {
    stop("Posterior moments must be named by joint_coefficient_names(component)")
  }
  if (min(eigen(posterior_cov, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("Posterior covariance is not positive definite")
  }
  centred <- lapply(c(os = "os", pfs = "pfs"), function(outcome) {
    fit <- component[[paste0("original_", outcome)]]
    v_full <- stats::setNames(fit$res.t[, "est"], rownames(fit$res))
    v_opt <- rownames(fit$res)[fit$optpars]
    v_full[v_opt] <- posterior_mean[paste0(outcome, "::", v_opt)]
    model <- build_sampled_flexsurv_model(fit, v_full)
    v_idx <- paste0(outcome, "::", v_opt)
    model$cov <- posterior_cov[v_idx, v_idx, drop = FALSE]
    dimnames(model$cov) <- dimnames(fit$cov)
    model$evidence_updated <- TRUE
    model
  })
  draws <- draw_joint_survival_coefficients(centred$os, centred$pfs, unname(posterior_cov), n_samples, seed)
  summarise_draws <- function(m, outcome) {
    data.frame(outcome = outcome, coefficient = colnames(m), mean = colMeans(m), sd = apply(m, 2, stats::sd),
               q2.5 = apply(m, 2, stats::quantile, 0.025), q97.5 = apply(m, 2, stats::quantile, 0.975),
               row.names = NULL, stringsAsFactors = FALSE)
  }
  list(method = "mvn_joint_v2",
       joint_covariance = list(covariance = posterior_cov,
                               prior_covariance = component$joint_covariance$covariance,
                               source = "conjugate evidence update (R/evidence_update.R)"),
       draws = draws, original_os = centred$os, original_pfs = centred$pfs,
       n_samples = as.integer(n_samples), n_failed = 0L,
       coefficient_summary = rbind(summarise_draws(draws$os, "os"), summarise_draws(draws$pfs, "pfs")),
       dist_os = component$dist_os, dist_pfs = component$dist_pfs, seed = seed,
       evidence_state = evidence_state, creation_time = Sys.time())
}

#' Build an evidence state: conjugate update, weight check, regenerated draws
#'
#' @param component Base joint sampling component (E0).
#' @param reference_population,weights Reference patients and fixed weights
#'   for the maps g_k.
#' @param evidence Data frame of external results: the raw evidence table
#'   (read_external_evidence()), external_evidence_likelihood() rows, or rows
#'   already selected by select_evidence_rows() (they carry
#'   `selection_reason`); or NULL.
#'   Optional columns override the defaults per row: `min_follow_up_weeks`,
#'   `accrual_weeks`, `landmark_weeks` (required for maintenance rows),
#'   `cells` (";"-separated codes).
#' @param select_rows Apply select_evidence_rows() to `evidence` unless it is
#'   already a selection (default TRUE; FALSE uses the rows as given).
#' @param maintenance_landmark_weeks Landmark of maintenance designs used by
#'   the selection.
#' @param sceptical List(coefs, sd) of a zero-centred sceptical prior, or NULL.
#' @param tau_shared,tau_trial,within_trial_correlation See evidence_covariance().
#' @param label Evidence-state label.
#' @param base_fingerprint Fingerprint of the base sampling cache.
#' @param evidence_md5 md5 of the evidence file (attr(read_external_evidence(), "md5")).
#' @param n_samples,seed Draws of the regenerated component.
#' @param tol Relinearisation tolerance.
#' @param check_args List of further importance_weight_check() arguments.
#' @param step Quadrature spacing of g_k (weeks).
#' @param biomarkers Biomarker columns.
#' @return List with `state` (sampling-cache layout: `joint`, `biomarkers`,
#'   `fingerprint`, `fingerprint_inputs`, `base_sampling_fingerprint`,
#'   `evidence_state`, `creation_time`), `update`, `check`, `observations`
#'   (data frame of the evidence rows with observed and predicted values) and
#'   `selection` (select_evidence_rows() output, NULL when not applied).
build_evidence_state <- function(component, reference_population, weights = NULL, evidence = NULL,
                                 select_rows = TRUE,
                                 maintenance_landmark_weeks = EVIDENCE_MAINTENANCE_LANDMARK_WEEKS,
                                 sceptical = NULL, tau_shared = 0, tau_trial = 0,
                                 within_trial_correlation = 0, label, base_fingerprint = NA_character_,
                                 evidence_md5 = NA_character_, n_samples = component$n_samples,
                                 seed = component$seed, tol = 1e-6, check_args = list(), step = 1,
                                 biomarkers = design_biomarkers()) {
  prior <- joint_prior_moments(component)
  v_names <- names(prior$mean)
  observations <- list()
  df_obs <- data.frame(label = character(), trial = character(), endpoint = character(), se = numeric(),
                       type = character(), stringsAsFactors = FALSE)
  specs <- list()
  selection <- NULL
  if (!is.null(evidence) && nrow(evidence) && select_rows && !"selection_reason" %in% names(evidence)) {
    selection <- select_evidence_rows(evidence, maintenance_landmark_weeks = maintenance_landmark_weeks)
    evidence <- selection$selected
  }
  if (!is.null(evidence) && nrow(evidence)) {
    for (i in seq_len(nrow(evidence))) {
      r <- evidence[i, , drop = FALSE]
      col_or <- function(col, default) if (col %in% names(r) && !is.na(r[[col]])) r[[col]] else default
      landmark <- col_or("landmark_weeks", 0)
      if (isTRUE(r$maintenance) && landmark == 0) {
        stop("Evidence row '", r$label, "' is a maintenance design: give landmark_weeks")
      }
      v_cells <- col_or("cells", NA)
      spec <- trial_summary_spec(r$label, r$endpoint, r$p_experimental,
                                 col_or("min_follow_up_weeks", r$follow_up_weeks), col_or("accrual_weeks", 0),
                                 landmark, if (is.na(v_cells)) NULL else strsplit(v_cells, ";", fixed = TRUE)[[1]],
                                 step)
      specs[[r$label]] <- spec
      observations[[length(observations) + 1]] <- list(
        label = r$label, y = r$log_hr,
        fn = make_trial_summary_map(spec, component, reference_population, weights, biomarkers))
      df_obs[nrow(df_obs) + 1, ] <- list(r$label, r$trial, r$endpoint, r$se, "external")
    }
  }
  if (!is.null(sceptical)) {
    terms <- sceptical_prior_terms(sceptical$coefs, v_names, sceptical$sd)
    for (j in seq_along(terms$y)) {
      observations[[length(observations) + 1]] <- list(
        label = terms$label[j], y = 0, a = stats::setNames(1, sceptical$coefs[j]))
      df_obs[nrow(df_obs) + 1, ] <- list(terms$label[j], "sceptical prior", sub("::.*$", "", sceptical$coefs[j]),
                                         sqrt(terms$V[j]), "prior")
    }
  }
  if (!length(observations)) stop("An evidence state needs external evidence or a sceptical prior")
  m_v <- evidence_covariance(df_obs, tau_shared, tau_trial, within_trial_correlation)
  update <- iterated_conjugate_update(prior$mean, prior$cov, observations, m_v, tol = tol)
  check <- do.call(importance_weight_check, c(list(
    draws = joint_draw_matrix(component), observations = observations, V = m_v,
    conjugate_mean = update$mean, prior_mean = prior$mean, conjugate_cov = update$cov), check_args))
  df_obs$observed <- update$observed
  df_obs$predicted_prior <- combine_observations(observations, v_names)$fn(prior$mean)
  df_obs$predicted_posterior <- update$predicted
  fingerprint <- evidence_state_fingerprint(
    base_fingerprint = base_fingerprint, label = label, evidence = evidence, evidence_md5 = evidence_md5,
    specs = specs, sceptical = sceptical, tau = list(shared = tau_shared, trial = tau_trial,
                                                     within_trial_correlation = within_trial_correlation),
    reference_population = reference_population, weights = weights, n_samples = n_samples, seed = seed,
    tol = tol, step = step)
  joint <- regenerate_evidence_component(
    component, update$mean, update$cov, n_samples, seed,
    evidence_state = list(label = label, prior_mean = prior$mean, posterior_mean = update$mean,
                          iterations = update$iterations, ess = check$ess, observations = df_obs))
  list(state = list(joint = joint, biomarkers = biomarkers, fingerprint = fingerprint$fingerprint,
                    fingerprint_inputs = fingerprint$inputs, base_sampling_fingerprint = base_fingerprint,
                    evidence_state = label, creation_time = Sys.time()),
       update = update, check = check, observations = df_obs, selection = selection)
}
