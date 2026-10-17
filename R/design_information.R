# ===============================================================================
# DESIGN INFORMATION OF THE JOINT GAMMA SURVIVAL MODELS (issue #202, paper 3 step 2)
# ===============================================================================
# flexsurv's gamma distribution ("gamma"; flexsurv.dists$gamma) has parameters
# shape k and rate lambda, with location parameter "rate": covariates act on
# log(rate), and both baseline parameters are optimised on the log scale. For
# patient i with covariate row x_i (no intercept) the optimisation-scale vector
# theta = (a, b, beta) gives
#
#   k = exp(a),  eta_i = b + x_i' beta,  lambda_i = exp(eta_i),
#   f(t) = lambda^k t^(k - 1) exp(-lambda t) / Gamma(k),  S(t) = Q(k, lambda t),
#
# with Q the regularised upper incomplete gamma function. Because the time
# scale is exp(-eta), exp(beta) is a rate ratio (accelerated failure time), not
# a hazard ratio. A patient's log-likelihood depends on theta only through
# (a, eta_i), so its negative Hessian is J_i' M_i J_i, where M_i is the 2 x 2
# negative Hessian in (a, eta) and J_i = d(a, eta_i) / d theta has rows
# e_shape and (0, 1, x_i'). The derivatives in eta are analytic for events and
# censored observations; the derivatives in a of a censored log Q(k, z) use
# central differences with Richardson extrapolation (gamma_censored_derivatives()).
# A left-truncated observation (entry at s > 0) subtracts the derivatives of
# log S(s).
#
# Three products:
#   1. stratum_observed_information(): the observed information of a fitted
#      gamma model at its estimates, split by arm x joint biomarker cell. The
#      strata sum to the inverse of the model's covariance (fit$cov, from
#      optim's numerical Hessian) within numerical tolerance; this is the
#      information decomposition paper 4 reports.
#   2. expected_gamma_information(): the expected observed information per
#      patient of a hypothetical design under a draw's own curves, integrated
#      over the event and censoring distribution (uniform accrual over A weeks,
#      administrative censoring at the analysis F weeks after the last
#      patient, optional randomisation landmark L with the clock of the model
#      kept and the likelihood left-truncated at L).
#   3. design_information(): expected information of n randomised patients
#      under a trial_design() (enrolled cells, cell proportions, allocation,
#      accrual, follow-up, randomisation week) for both endpoints, and V_d(n),
#      the covariance of the identifiable treatment-effect targets the design
#      delivers (design_treatment_targets()). An enriched design leaves
#      coefficient combinations unidentified (e.g. TMB/BRAF-positive only:
#      rate and tmb_braf, Rx and Rx:tmb_braf are aliased); the information is
#      singular, the null space is returned, and V_d(n) is computed for the
#      identifiable combinations only, which are checked to be orthogonal to
#      the null space.
#
# Randomisation at week L (the all-comers design S1 randomised at week 4 after
# CRP is known): only patients alive and progression-free at L are randomised;
# their cell mix is reweighted by P(PFS > L) under control (everyone receives
# FLOX before randomisation), and each endpoint's likelihood is left-truncated
# at L under the patient's arm. For OS this conditions on OS > L, not on
# PFS > L: the partitioned-survival model has no joint OS/PFS distribution.
#
# Cross-endpoint covariance: the OS and PFS estimates of one trial are
# correlated through its patients, which the marginal likelihoods do not show.
# design_information(cross_endpoint = "recolour") keeps the whitened
# cross-correlation of the targets under the current joint coefficient
# covariance (the #159 paired-bootstrap calibration) and recolours it with
# the design's marginal covariances, as calibrate_joint_covariance() does;
# "independent" sets the cross-block to zero.
#
# Requires R/pfs_os_violation_diagnostics.R (flexsurv_design_matrix()),
# R/joint_survival_sampling.R (covariance_root()) and
# R/voi_survival_parameters.R (survival_coefficient_map()).
# ===============================================================================

#' Step on the log-shape scale for the finite differences of censored terms
GAMMA_SHAPE_STEP <- 1e-3

#' Default joint-biomarker cell codes (biomarkers in get_biomarkers() order)
DESIGN_ALL_CELLS <- c("00", "01", "10", "11")

#' Default biomarker identifiers used when get_biomarkers() is unavailable
#'
#' @return Character vector of biomarker column names.
design_biomarkers <- function() {
  if (exists("get_biomarkers", mode = "function")) get_biomarkers() else c("crp", "tmb_braf")
}

#' Log-likelihood derivatives of observed gamma events in (log shape, eta)
#'
#' @param log_shape Log shape (scalar).
#' @param eta Linear predictor of the log rate (vector).
#' @param time Event times (vector, recycled with `eta`), strictly positive.
#' @return Matrix with columns `l`, `a`, `eta`, `aa`, `aeta`, `etaeta`: the
#'   log density and its first and second derivatives.
gamma_event_derivatives <- function(log_shape, eta, time) {
  k <- exp(log_shape)
  v_log_z <- eta + log(time)
  v_z <- exp(v_log_z)
  dig <- digamma(k)
  v_a <- k * (v_log_z - dig)
  cbind(l = k * eta + (k - 1) * log(time) - v_z - lgamma(k),
        a = v_a,
        eta = k - v_z,
        aa = v_a - k^2 * trigamma(k),
        aeta = rep(k, length(v_z)),
        etaeta = -v_z)
}

#' Analytic eta-derivatives of a censored gamma log-survival log Q(k, z)
#'
#' @param k Shape (scalar).
#' @param z Rate times time (vector, non-negative).
#' @return List with `l` (log Q), `eta` (first derivative in eta = log rate)
#'   and `etaeta` (second derivative).
gamma_censored_eta_terms <- function(k, z) {
  v_l <- stats::pgamma(z, shape = k, lower.tail = FALSE, log.p = TRUE)
  # h = z f(z) / Q(z): d log Q / d eta = -h, d^2 log Q / d eta^2 = -h (k - z) - h^2.
  v_h <- numeric(length(z))
  v_pos <- z > 0
  v_h[v_pos] <- exp(log(z[v_pos]) + stats::dgamma(z[v_pos], shape = k, log = TRUE) - v_l[v_pos])
  list(l = v_l, eta = -v_h, etaeta = -v_h * (k - z) - v_h^2)
}

#' Log-likelihood derivatives of censored gamma observations in (log shape, eta)
#'
#' Derivatives in the log shape use central differences at steps h and 2h
#' combined by Richardson extrapolation (truncation error O(h^4)).
#'
#' @param log_shape Log shape (scalar).
#' @param eta Linear predictor of the log rate (vector).
#' @param time Censoring times (vector, non-negative).
#' @param step Step h on the log-shape scale.
#' @return Matrix with columns `l`, `a`, `eta`, `aa`, `aeta`, `etaeta`.
gamma_censored_derivatives <- function(log_shape, eta, time, step = GAMMA_SHAPE_STEP) {
  v_z <- exp(eta) * time
  at <- function(offset) gamma_censored_eta_terms(exp(log_shape + offset), v_z)
  d0 <- at(0)
  dp1 <- at(step); dm1 <- at(-step); dp2 <- at(2 * step); dm2 <- at(-2 * step)
  richardson <- function(d1, d2) (4 * d1 - d2) / 3
  v_a <- richardson((dp1$l - dm1$l) / (2 * step), (dp2$l - dm2$l) / (4 * step))
  v_aa <- richardson((dp1$l - 2 * d0$l + dm1$l) / step^2, (dp2$l - 2 * d0$l + dm2$l) / (4 * step^2))
  v_aeta <- richardson((dp1$eta - dm1$eta) / (2 * step), (dp2$eta - dm2$eta) / (4 * step))
  cbind(l = d0$l, a = v_a, eta = d0$eta, aa = v_aa, aeta = v_aeta, etaeta = d0$etaeta)
}

#' Log-likelihood derivatives of gamma observations, events and censored
#'
#' @param log_shape Log shape (scalar).
#' @param eta Linear predictor of the log rate (vector).
#' @param time Observation times (vector).
#' @param event 1 for an event, 0 for a censored observation (vector).
#' @param step Finite-difference step for censored terms.
#' @return Matrix with columns `l`, `a`, `eta`, `aa`, `aeta`, `etaeta`.
gamma_loglik_derivatives <- function(log_shape, eta, time, event, step = GAMMA_SHAPE_STEP) {
  n <- max(length(eta), length(time), length(event))
  eta <- rep_len(eta, n); time <- rep_len(time, n); event <- rep_len(event, n)
  m_out <- matrix(NA_real_, n, 6, dimnames = list(NULL, c("l", "a", "eta", "aa", "aeta", "etaeta")))
  v_ev <- event == 1
  if (any(v_ev)) m_out[v_ev, ] <- gamma_event_derivatives(log_shape, eta[v_ev], time[v_ev])
  if (any(!v_ev)) m_out[!v_ev, ] <- gamma_censored_derivatives(log_shape, eta[!v_ev], time[!v_ev], step)
  m_out
}

#' Assemble a gamma information matrix from (a, eta) information contributions
#'
#' @param m_x1 Rate design matrix with the intercept column first (one row per
#'   patient or patient type); its first column is renamed "rate".
#' @param v_aa,v_aeta,v_etaeta Information contributions (negative second
#'   derivatives) per row, already multiplied by any row weights.
#' @param par_names Parameter order of the result (rownames of the fit's res).
#' @return Symmetric information matrix with dimnames `par_names`.
gamma_information_matrix <- function(m_x1, v_aa, v_aeta, v_etaeta, par_names) {
  colnames(m_x1)[1] <- "rate"
  v_ax <- colSums(v_aeta * m_x1)
  m_xx <- crossprod(m_x1, v_etaeta * m_x1)
  m_info <- rbind(c(sum(v_aa), v_ax), cbind(v_ax, m_xx))
  v_names <- c("shape", colnames(m_x1))
  dimnames(m_info) <- list(v_names, v_names)
  if (!setequal(v_names, par_names)) {
    stop("Information parameters (", paste(v_names, collapse = ", "),
         ") differ from the model parameters (", paste(par_names, collapse = ", "), ")")
  }
  m_info <- m_info[par_names, par_names, drop = FALSE]
  (m_info + t(m_info)) / 2
}

#' Check that a model is a flexsurvreg gamma fit without fixed parameters
#'
#' @param model flexsurvreg object.
#' @return The model, invisibly; stops otherwise.
assert_gamma_flexsurv <- function(model) {
  if (!inherits(model, "flexsurvreg") || !identical(model$dlist$name, "gamma")) {
    stop("A flexsurvreg fit with dist = \"gamma\" is required")
  }
  if (length(model$fixedpars)) stop("Fits with fixed parameters are not supported")
  if (!identical(model$dlist$location, "rate")) stop("Expected the gamma location parameter 'rate'")
  invisible(model)
}

#' Per-patient observed information contributions of a gamma fit
#'
#' @param model flexsurvreg gamma fit (right-censored or left-truncated data).
#' @param par Parameter vector on the optimisation scale named as
#'   rownames(model$res); defaults to the fitted estimates.
#' @param step Finite-difference step for censored terms.
#' @return List with `m_x1` (rate design matrix with intercept), `info`
#'   (matrix with columns aa, aeta, etaeta: negative second derivatives per
#'   patient), `event`, `par` and `frame` (the model frame).
gamma_patient_information <- function(model, par = NULL, step = GAMMA_SHAPE_STEP) {
  assert_gamma_flexsurv(model)
  if (is.null(par)) par <- stats::setNames(model$res.t[, "est"], rownames(model$res))
  m_y <- model$data$Y
  if (!attr(m_y, "type") %in% c("right", "counting")) {
    stop("Only right-censored or left-truncated data are supported; got ", attr(m_y, "type"))
  }
  m_x1 <- model$data$mml$rate
  v_cov <- colnames(m_x1)[-1]
  v_eta <- as.vector(par[["rate"]] + m_x1[, -1, drop = FALSE] %*% par[v_cov])
  v_event <- m_y[, "status"]
  m_d <- gamma_loglik_derivatives(par[["shape"]], v_eta, m_y[, "stop"], v_event, step)
  v_start <- m_y[, "start"]
  v_trunc <- v_start > 0
  if (any(v_trunc)) {
    m_s <- gamma_censored_derivatives(par[["shape"]], v_eta[v_trunc], v_start[v_trunc], step)
    m_d[v_trunc, c("aa", "aeta", "etaeta")] <- m_d[v_trunc, c("aa", "aeta", "etaeta")] -
      m_s[, c("aa", "aeta", "etaeta")]
  }
  list(m_x1 = m_x1, info = -m_d[, c("aa", "aeta", "etaeta"), drop = FALSE],
       event = v_event, par = par, frame = model$data$m)
}

#' Joint biomarker cell code of each row
#'
#' @param df Data frame with 0/1 (numeric or factor) biomarker columns.
#' @param biomarkers Biomarker columns, in cell-code order.
#' @return Character vector of cell codes such as "10".
design_cell_codes <- function(df, biomarkers = design_biomarkers()) {
  v_missing <- setdiff(biomarkers, names(df))
  if (length(v_missing)) stop("Biomarker columns missing: ", paste(v_missing, collapse = ", "))
  m_codes <- vapply(biomarkers, function(b) {
    v <- as.numeric(as.character(df[[b]]))
    if (anyNA(v) || !all(v %in% c(0, 1))) stop("Biomarker ", b, " must be 0/1 without missing values")
    as.character(v)
  }, character(nrow(df)))
  if (!is.matrix(m_codes)) m_codes <- matrix(m_codes, nrow = nrow(df))
  apply(m_codes, 1, paste, collapse = "")
}

#' Arm label ("control" or "experimental") from an Rx factor
#'
#' @param rx Factor whose first level is control and second experimental.
#' @return Character vector.
design_arm_labels <- function(rx) {
  if (!is.factor(rx) || nlevels(rx) != 2L) stop("Rx must be a two-level factor (control first)")
  c("control", "experimental")[as.integer(rx)]
}

#' Scale-free difference between two information matrices
#'
#' @param m_a,m_b Square matrices with the same dimnames.
#' @return Maximum absolute entry of D^(-1/2) (m_a - m_b) D^(-1/2), with D the
#'   diagonal of `m_b`, i.e. the largest difference relative to the geometric
#'   mean of the two diagonal entries.
scaled_matrix_difference <- function(m_a, m_b) {
  v_d <- sqrt(abs(diag(m_b)))
  v_d[v_d == 0] <- 1
  max(abs((m_a - m_b) / outer(v_d, v_d)))
}

#' Observed information of a gamma fit split by arm x joint biomarker cell
#'
#' @param model flexsurvreg gamma fit whose model frame contains `Rx` and the
#'   biomarker columns.
#' @param par Parameter vector (default: the estimates).
#' @param biomarkers Biomarker columns defining the cells.
#' @param step Finite-difference step for censored terms.
#' @return List with `strata` (named list of information matrices, names
#'   "<arm>|<cell>"), `total` (their sum), `counts` (data frame stratum, arm,
#'   cell, n, events), `inverse_vcov` (solve(model$cov)) and
#'   `scaled_difference` (scaled_matrix_difference(total, inverse_vcov)).
stratum_observed_information <- function(model, par = NULL, biomarkers = design_biomarkers(),
                                         step = GAMMA_SHAPE_STEP) {
  parts <- gamma_patient_information(model, par, step)
  df_frame <- parts$frame
  if (!"Rx" %in% names(df_frame)) stop("The model frame has no Rx column")
  v_arm <- design_arm_labels(df_frame$Rx)
  v_cell <- design_cell_codes(df_frame, biomarkers)
  v_stratum <- paste(v_arm, v_cell, sep = "|")
  v_par_names <- rownames(model$res)
  v_levels <- sort(unique(v_stratum), method = "radix")
  strata <- lapply(stats::setNames(v_levels, v_levels), function(s) {
    v_rows <- v_stratum == s
    gamma_information_matrix(parts$m_x1[v_rows, , drop = FALSE], parts$info[v_rows, "aa"],
                             parts$info[v_rows, "aeta"], parts$info[v_rows, "etaeta"], v_par_names)
  })
  m_total <- Reduce(`+`, strata)
  m_inv_vcov <- solve(model$cov)
  dimnames(m_inv_vcov) <- list(v_par_names, v_par_names)
  df_counts <- data.frame(stratum = v_levels,
                          arm = sub("\\|.*$", "", v_levels), cell = sub("^.*\\|", "", v_levels),
                          n = as.integer(table(factor(v_stratum, v_levels))),
                          events = as.numeric(tapply(parts$event, factor(v_stratum, v_levels), sum)),
                          stringsAsFactors = FALSE)
  list(strata = strata, total = m_total, counts = df_counts, inverse_vcov = m_inv_vcov,
       scaled_difference = scaled_matrix_difference(m_total, m_inv_vcov))
}

#' Survival of the uniform-accrual administrative censoring distribution
#'
#' Patients enter uniformly over `accrual_weeks`; the analysis is
#' `min_follow_up_weeks` after the last entry, so the potential follow-up C of
#' a patient (time from randomisation to analysis) is uniform on [F, F + A].
#'
#' @param s Time since randomisation (weeks).
#' @param min_follow_up_weeks F.
#' @param accrual_weeks A (0 gives fixed follow-up F).
#' @return P(C > s).
censoring_survival <- function(s, min_follow_up_weeks, accrual_weeks) {
  if (accrual_weeks > 0) {
    pmin(1, pmax(0, (min_follow_up_weeks + accrual_weeks - s) / accrual_weeks))
  } else {
    as.numeric(s < min_follow_up_weeks)
  }
}

#' Quadrature nodes on the clock since randomisation
#'
#' @param min_follow_up_weeks,accrual_weeks Follow-up F and accrual A.
#' @param step Node spacing (weeks); F and F + A are always nodes.
#' @param refine_origin Add geometrically spaced nodes in (0, step) when the
#'   clock starts at the model origin, where a gamma density with shape below
#'   one is unbounded and the log-time terms of the information are singular.
#' @return Sorted unique nodes from 0 to F + A.
design_time_nodes <- function(min_follow_up_weeks, accrual_weeks, step, refine_origin = FALSE) {
  if (!is.finite(min_follow_up_weeks) || min_follow_up_weeks <= 0) stop("min_follow_up_weeks must be positive")
  if (!is.finite(accrual_weeks) || accrual_weeks < 0) stop("accrual_weeks must be non-negative")
  end <- min_follow_up_weeks + accrual_weeks
  sort(unique(c(seq(0, min_follow_up_weeks, by = step), min_follow_up_weeks,
                if (refine_origin) step * 10^seq(-8, -0.1, length.out = 40),
                if (accrual_weeks > 0) seq(min_follow_up_weeks, end, by = step), end)),
       method = "radix")
}

#' Expected observed gamma information per patient of a hypothetical design
#'
#' Integrates the per-patient negative Hessian in (log shape, eta) over the
#' event and censoring distribution implied by the draw's own curves: events
#' observed in each quadrature interval (exact interval mass, derivatives at
#' the midpoint, censoring survival at the midpoint), administrative
#' censoring at C uniform on [F, F + A] (trapezoidal rule; a point mass at F
#' when A = 0), and left truncation at the randomisation landmark L.
#'
#' @param log_shape Log shape of the draw.
#' @param eta Linear predictors of the patient types (log rate).
#' @param min_follow_up_weeks,accrual_weeks Follow-up F and accrual A (weeks).
#' @param landmark_weeks Randomisation landmark L on the model clock.
#' @param step Quadrature node spacing (weeks).
#' @return Matrix, one row per patient type, columns `aa`, `aeta`, `etaeta`
#'   (expected information contributions) and `events` (expected probability
#'   of an observed event).
expected_gamma_information <- function(log_shape, eta, min_follow_up_weeks, accrual_weeks = 0,
                                       landmark_weeks = 0, step = 0.5) {
  k <- exp(log_shape)
  v_nodes <- design_time_nodes(min_follow_up_weeks, accrual_weeks, step,
                               refine_origin = landmark_weeks == 0)
  n_nodes <- length(v_nodes)
  n_types <- length(eta)
  v_rate <- exp(eta)
  v_s_land <- stats::pgamma(landmark_weeks, shape = k, rate = v_rate, lower.tail = FALSE)
  if (any(v_s_land <= 0)) stop("A patient type has zero survival at the landmark")
  m_surv <- outer(v_rate, landmark_weeks + v_nodes,
                  function(r, t) stats::pgamma(t, shape = k, rate = r, lower.tail = FALSE)) / v_s_land
  # Event part: P(event in interval j and observed) x information at the midpoint.
  v_mid <- (v_nodes[-1] + v_nodes[-n_nodes]) / 2
  m_mass <- (m_surv[, -n_nodes, drop = FALSE] - m_surv[, -1, drop = FALSE])
  m_obs <- sweep(m_mass, 2, censoring_survival(v_mid, min_follow_up_weeks, accrual_weeks), "*")
  m_ev <- gamma_loglik_derivatives(log_shape, rep(eta, times = n_nodes - 1),
                                   rep(landmark_weeks + v_mid, each = n_types), 1)
  info_sum <- function(m_w, m_d) {
    vapply(c("aa", "aeta", "etaeta"), function(j) rowSums(m_w * matrix(-m_d[, j], n_types)),
           numeric(n_types))
  }
  m_info <- info_sum(m_obs, m_ev)
  if (!is.matrix(m_info)) m_info <- matrix(m_info, nrow = n_types, dimnames = list(NULL, c("aa", "aeta", "etaeta")))
  # Censoring part: censored at c with density 1/A on [F, F + A] if still at risk.
  v_cens_nodes <- if (accrual_weeks > 0) v_nodes[v_nodes >= min_follow_up_weeks] else min_follow_up_weeks
  m_cens <- gamma_censored_derivatives(log_shape, rep(eta, times = length(v_cens_nodes)),
                                       rep(landmark_weeks + v_cens_nodes, each = n_types))
  m_s_cens <- m_surv[, match(v_cens_nodes, v_nodes), drop = FALSE]
  if (accrual_weeks > 0) {
    n_c <- length(v_cens_nodes)
    v_trap <- c(diff(v_cens_nodes), 0) / 2 + c(0, diff(v_cens_nodes)) / 2
    m_w_cens <- sweep(m_s_cens, 2, v_trap / accrual_weeks, "*")
    if (n_c < 2) stop("Too few censoring nodes; reduce step")
  } else {
    m_w_cens <- m_s_cens
  }
  m_info <- m_info + info_sum(m_w_cens, m_cens)
  if (landmark_weeks > 0) {
    m_trunc <- gamma_censored_derivatives(log_shape, eta, rep(landmark_weeks, n_types))
    m_info <- m_info - (-m_trunc[, c("aa", "aeta", "etaeta"), drop = FALSE])
  }
  cbind(m_info, events = rowSums(m_obs))
}

#' Specification of a candidate trial design
#'
#' @param name Design label.
#' @param enrolled_cells Cell codes enrolled (biomarkers in `biomarkers` order).
#' @param cell_proportions Named vector of cell proportions at registration
#'   (before any landmark selection); NULL uses the reference population's
#'   weighted mix. Restricted to `enrolled_cells` and renormalised.
#' @param p_experimental Allocation probability of the experimental arm.
#' @param min_follow_up_weeks Follow-up after the last patient enters.
#' @param accrual_weeks Accrual duration (weeks); give this or `accrual_rate`.
#' @param accrual_rate Randomised patients per week; the accrual duration is
#'   then n / accrual_rate.
#' @param randomisation_week Randomisation landmark L on the model clock (4 for
#'   randomisation after the week-4 CRP; 0 for randomisation at baseline).
#' @param biomarkers Biomarker columns defining the cells.
#' @return List of class "voi_trial_design".
trial_design <- function(name, enrolled_cells = DESIGN_ALL_CELLS, cell_proportions = NULL,
                         p_experimental = 0.5, min_follow_up_weeks, accrual_weeks = NULL,
                         accrual_rate = NULL, randomisation_week = 0,
                         biomarkers = design_biomarkers()) {
  if (is.null(accrual_weeks) == is.null(accrual_rate)) stop("Give exactly one of accrual_weeks and accrual_rate")
  if (!is.null(accrual_rate) && (!is.finite(accrual_rate) || accrual_rate <= 0)) stop("accrual_rate must be positive")
  if (!is.finite(p_experimental) || p_experimental <= 0 || p_experimental >= 1) stop("p_experimental must lie in (0, 1)")
  if (!length(enrolled_cells) || anyDuplicated(enrolled_cells) ||
      !all(nchar(enrolled_cells) == length(biomarkers)) || any(grepl("[^01]", enrolled_cells))) {
    stop("enrolled_cells must be distinct 0/1 codes of length ", length(biomarkers))
  }
  if (!is.null(cell_proportions) && (is.null(names(cell_proportions)) ||
      !all(enrolled_cells %in% names(cell_proportions)) || any(cell_proportions < 0))) {
    stop("cell_proportions must be non-negative and named by cell code, covering enrolled_cells")
  }
  if (!is.finite(randomisation_week) || randomisation_week < 0) stop("randomisation_week must be non-negative")
  structure(list(name = name, enrolled_cells = enrolled_cells, cell_proportions = cell_proportions,
                 p_experimental = p_experimental, min_follow_up_weeks = min_follow_up_weeks,
                 accrual_weeks = accrual_weeks, accrual_rate = accrual_rate,
                 randomisation_week = randomisation_week, biomarkers = biomarkers),
            class = "voi_trial_design")
}

#' Accrual duration of a design with n patients
#'
#' @param design trial_design().
#' @param n Number of randomised patients.
#' @return Weeks.
design_accrual_weeks <- function(design, n) {
  if (!is.null(design$accrual_weeks)) design$accrual_weeks else n / design$accrual_rate
}

#' Split a joint draw into endpoint coefficient vectors
#'
#' @param draw Either list(os = , pfs = ) of vectors named as rownames(fit$res),
#'   or one vector with names "os::<name>" and "pfs::<name>" (optimised
#'   parameters only; fixed ones are taken from the component's fits).
#' @param component Joint sampling component (original_os, original_pfs).
#' @return list(os = , pfs = ) full named vectors.
split_joint_draw <- function(draw, component) {
  originals <- list(os = component$original_os, pfs = component$original_pfs)
  lapply(stats::setNames(names(originals), names(originals)), function(outcome) {
    fit <- originals[[outcome]]
    v_names <- rownames(fit$res)
    if (is.list(draw)) {
      v <- draw[[outcome]]
      if (is.null(names(v))) names(v) <- v_names
      v <- v[v_names]
    } else {
      v <- stats::setNames(fit$res.t[, "est"], v_names)
      v_key <- paste0(outcome, "::", v_names)
      v_have <- v_key %in% names(draw)
      v[v_have] <- draw[v_key[v_have]]
    }
    if (anyNA(v)) stop("The ", outcome, " draw does not cover every model parameter")
    v
  })
}

#' Coefficient vectors of one cached draw
#'
#' @param component Joint sampling component with `draws$os`, `draws$pfs`.
#' @param idx Row of the draw matrices.
#' @return list(os = , pfs = ) named vectors.
component_draw <- function(component, idx) {
  list(os = component$draws$os[idx, ], pfs = component$draws$pfs[idx, ])
}

#' Fitted (mean) coefficient vectors of a sampling component
#'
#' @param component Joint sampling component.
#' @return list(os = , pfs = ) named vectors of the original estimates.
component_mean_draw <- function(component) {
  lapply(list(os = component$original_os, pfs = component$original_pfs), function(fit)
    stats::setNames(fit$res.t[, "est"], rownames(fit$res)))
}

#' Levels of the treatment factor of a flexsurvreg fit
#'
#' flexsurvreg stores factor levels in `covdata$xlev` (its `xlev` is NULL).
#'
#' @param model flexsurvreg fit with an `Rx` covariate.
#' @return Character vector: control level first, experimental second.
model_rx_levels <- function(model) {
  v_levels <- model$covdata$xlev$Rx
  if (is.null(v_levels) && is.factor(model$data$m$Rx)) v_levels <- levels(model$data$m$Rx)
  if (length(v_levels) != 2L) stop("The model has no two-level Rx factor")
  v_levels
}

#' Linear predictors of a gamma model for patients under a given arm
#'
#' @param model flexsurvreg gamma fit (formula and factor levels).
#' @param par Full named parameter vector.
#' @param newdata Patients; `Rx` is overwritten with `rx_level`.
#' @param rx_level Level of Rx to assign, or NULL to keep each patient's Rx.
#' @return Numeric vector of eta = log rate.
gamma_linear_predictor <- function(model, par, newdata, rx_level = NULL) {
  if (!is.null(rx_level)) newdata$Rx <- factor(rx_level, levels = model_rx_levels(model))
  m_x <- flexsurv_design_matrix(model, newdata)
  as.vector(par[["rate"]] + m_x %*% par[colnames(m_x)])
}

#' Patient mix of a design: reference patients by enrolled cell and arm
#'
#' @param design trial_design().
#' @param reference_population Data frame of reference patients (Age, sex,
#'   Rx factor, biomarkers); METIMMOX's prediction population in live runs.
#' @param weights Reference patient weights (default equal).
#' @param component Joint sampling component (for the landmark selection).
#' @param draw list(os, pfs) of the draw (for the landmark selection).
#' @return Data frame of patient types (reference rows x arm) with the
#'   reference covariates, `Rx`, `arm`, `cell`, `weight` (summing to one over
#'   randomised patients) and `selection` (P(PFS > L) under control).
design_patient_mix <- function(design, reference_population, weights = NULL, component, draw) {
  if (is.null(weights)) weights <- rep(1, nrow(reference_population))
  if (length(weights) != nrow(reference_population) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be non-negative and match the reference population")
  }
  v_cell <- design_cell_codes(reference_population, design$biomarkers)
  v_keep <- v_cell %in% design$enrolled_cells & weights > 0
  v_missing_cells <- setdiff(design$enrolled_cells, v_cell[v_keep])
  if (length(v_missing_cells)) stop("No reference patients in enrolled cell(s) ", paste(v_missing_cells, collapse = ", "))
  df_ref <- reference_population[v_keep, , drop = FALSE]
  v_w <- weights[v_keep]
  v_cell <- v_cell[v_keep]
  v_cell_mass <- tapply(v_w, v_cell, sum)
  v_target <- if (is.null(design$cell_proportions)) v_cell_mass else design$cell_proportions[names(v_cell_mass)]
  v_target <- v_target / sum(v_target)
  v_w <- v_w / v_cell_mass[v_cell] * v_target[v_cell]
  v_selection <- rep(1, length(v_w))
  if (design$randomisation_week > 0) {
    fit_pfs <- component$original_pfs
    v_eta <- gamma_linear_predictor(fit_pfs, draw$pfs, df_ref, model_rx_levels(fit_pfs)[1])
    v_selection <- stats::pgamma(design$randomisation_week, shape = exp(draw$pfs[["shape"]]),
                                 rate = exp(v_eta), lower.tail = FALSE)
  }
  v_w <- v_w * v_selection
  v_w <- v_w / sum(v_w)
  v_levels <- model_rx_levels(component$original_os)
  df_types <- rbind(df_ref, df_ref)
  df_types$Rx <- factor(rep(v_levels, each = nrow(df_ref)), levels = v_levels)
  df_types$arm <- rep(c("control", "experimental"), each = nrow(df_ref))
  df_types$cell <- rep(v_cell, 2)
  df_types$selection <- rep(v_selection, 2)
  df_types$weight <- c(v_w * (1 - design$p_experimental), v_w * design$p_experimental)
  rownames(df_types) <- NULL
  df_types
}

#' Treatment-effect targets a design can identify, per endpoint
#'
#' The treatment effect in cell (c, m) is Rx + c Rx:crp + m Rx:tmb_braf. The
#' targets span the effects of the enrolled cells; the basis prefers the
#' coefficients themselves (Rx, then the interactions) where they lie in that
#' span and completes it with cell effects. All four cells give
#' {Rx, Rx:crp, Rx:tmb_braf}; TMB/BRAF-positive cells only give
#' {Rx:crp, cell_effect_01}; CRP-positive cells only give
#' {Rx:tmb_braf, cell_effect_10}.
#'
#' @param coef_names Coefficient names of one endpoint (rownames(fit$res)).
#' @param outcome Endpoint label used in the row and column names.
#' @param enrolled_cells Enrolled cell codes.
#' @param biomarkers Biomarker identifiers (cell-code order).
#' @return Matrix, one row per target (names "<outcome>::<short>"), columns
#'   "<outcome>::<coefficient>".
design_treatment_targets <- function(coef_names, outcome, enrolled_cells, biomarkers = design_biomarkers()) {
  df_map <- survival_coefficient_map(coef_names, outcome, biomarkers)
  rx <- df_map$coefficient[df_map$role == "rx_main"]
  v_int <- vapply(biomarkers, function(b) {
    v <- df_map$coefficient[df_map$role == "interaction" & df_map$biomarker == b]
    if (length(v) != 1) stop("No unique ", outcome, " interaction for biomarker ", b)
    v
  }, character(1))
  v_cols <- c(rx, v_int)
  unit_vector <- function(j) { v <- numeric(length(v_cols)); v[j] <- 1; v }
  cell_vector <- function(code) c(1, as.numeric(strsplit(code, "")[[1]]))
  m_span <- t(vapply(enrolled_cells, cell_vector, numeric(length(v_cols))))
  rank_span <- qr(m_span)$rank
  in_span <- function(v) {
    fit <- qr.resid(qr(t(m_span)), v)
    sqrt(sum(fit^2)) < 1e-10
  }
  candidates <- c(stats::setNames(lapply(seq_along(v_cols), unit_vector),
                                  c("rx", paste0("rx_x_", biomarkers))),
                  stats::setNames(lapply(enrolled_cells, cell_vector),
                                  paste0("cell_effect_", enrolled_cells)))
  m_basis <- NULL
  for (nm in names(candidates)) {
    v <- candidates[[nm]]
    if (!in_span(v)) next
    m_try <- rbind(m_basis, stats::setNames(v, v_cols))
    if (qr(m_try)$rank > NROW(m_basis)) {
      rownames(m_try)[nrow(m_try)] <- nm
      m_basis <- m_try
    }
    if (NROW(m_basis) == rank_span) break
  }
  m_out <- matrix(0, nrow(m_basis), length(coef_names),
                  dimnames = list(paste0(outcome, "::", rownames(m_basis)), paste0(outcome, "::", coef_names)))
  m_out[, paste0(outcome, "::", v_cols)] <- m_basis
  m_out
}

#' Stack target matrices into one matrix over the joint coefficient columns
#'
#' Columns are matched by name (never by position), so endpoint blocks with
#' different column sets can be combined.
#'
#' @param l_targets List of target matrices (named rows and columns).
#' @param v_columns Joint column names ("os::<coef>", "pfs::<coef>").
#' @return Matrix with one row per target and columns `v_columns`; stops if a
#'   target loads on a column outside `v_columns` or row names repeat.
joint_target_matrix <- function(l_targets, v_columns) {
  m_out <- do.call(rbind, lapply(l_targets, function(m) {
    v_unknown <- setdiff(colnames(m)[colSums(abs(m)) > 0], v_columns)
    if (length(v_unknown)) stop("Target columns not in the model: ", paste(v_unknown, collapse = ", "))
    m_full <- matrix(0, nrow(m), length(v_columns), dimnames = list(rownames(m), v_columns))
    v_known <- intersect(colnames(m), v_columns)
    m_full[, v_known] <- m[, v_known, drop = FALSE]
    m_full
  }))
  if (is.null(rownames(m_out)) || anyDuplicated(rownames(m_out))) stop("Targets need unique row names")
  m_out
}

#' Moore-Penrose inverse, rank and null space of a symmetric information matrix
#'
#' @param m_info Symmetric positive semidefinite matrix.
#' @param tol Relative eigenvalue tolerance after diagonal scaling.
#' @return List with `inverse` (pseudo-inverse), `rank`, `null_space`
#'   (columns spanning the null space, on the original scale).
information_pseudo_inverse <- function(m_info, tol = 1e-9) {
  v_d <- sqrt(abs(diag(m_info)))
  v_d[v_d == 0] <- 1
  m_scaled <- m_info / outer(v_d, v_d)
  eig <- eigen((m_scaled + t(m_scaled)) / 2, symmetric = TRUE)
  v_keep <- eig$values > tol * max(eig$values)
  m_inv_scaled <- eig$vectors[, v_keep, drop = FALSE] %*%
    (t(eig$vectors[, v_keep, drop = FALSE]) / eig$values[v_keep])
  m_inverse <- m_inv_scaled / outer(v_d, v_d)
  dimnames(m_inverse) <- dimnames(m_info)
  # Null space of the scaled matrix mapped back: m_info v = 0 iff m_scaled (D v) = 0.
  m_null <- eig$vectors[, !v_keep, drop = FALSE] / v_d
  if (ncol(m_null)) m_null <- sweep(m_null, 2, sqrt(colSums(m_null^2)), "/")
  rownames(m_null) <- rownames(m_info)
  list(inverse = m_inverse, rank = sum(v_keep), null_space = m_null)
}

#' Expected information of n randomised patients under a design, both endpoints
#'
#' @param design trial_design().
#' @param n Number of randomised patients (scalar).
#' @param draw Draw whose curves define the expected information
#'   (split_joint_draw() formats; e.g. component_mean_draw(component)).
#' @param component Joint sampling component (fits and, for the recoloured
#'   cross-endpoint covariance, `joint_covariance$covariance`).
#' @param reference_population Reference patients (Age, sex, Rx, biomarkers).
#' @param weights Reference patient weights (default equal).
#' @param targets Joint target matrix (rows: targets; columns "os::<coef>",
#'   "pfs::<coef>"); NULL uses design_treatment_targets() for both endpoints.
#' @param cross_endpoint "recolour" or "independent" (see file header).
#' @param step Quadrature node spacing (weeks).
#' @return List of class "voi_design_information": `design`, `n`,
#'   `accrual_weeks`, `targets`, `covariance` (V_d(n) of the targets),
#'   `cross_endpoint`, `endpoints` (per endpoint: `information`, `rank`,
#'   `null_space`, `strata` (information by arm|cell), `expected_events`
#'   (data frame by arm and cell), `covariance` (pseudo-inverse)),
#'   `identifiable` (logical per target) and `patient_mix`.
design_information <- function(design, n, draw, component, reference_population, weights = NULL,
                               targets = NULL, cross_endpoint = c("recolour", "independent"),
                               step = 0.5) {
  if (!inherits(design, "voi_trial_design")) stop("design must come from trial_design()")
  cross_endpoint <- match.arg(cross_endpoint)
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) stop("n must be one positive number")
  draw <- split_joint_draw(draw, component)
  accrual <- design_accrual_weeks(design, n)
  df_mix <- design_patient_mix(design, reference_population, weights, component, draw)
  v_stratum <- paste(df_mix$arm, df_mix$cell, sep = "|")
  endpoints <- list()
  for (outcome in c("os", "pfs")) {
    fit <- component[[paste0("original_", outcome)]]
    assert_gamma_flexsurv(fit)
    v_par <- draw[[outcome]]
    v_names <- rownames(fit$res)
    df_new <- df_mix
    m_x <- flexsurv_design_matrix(fit, df_new)
    v_eta <- as.vector(v_par[["rate"]] + m_x %*% v_par[colnames(m_x)])
    m_e <- expected_gamma_information(v_par[["shape"]], v_eta, design$min_follow_up_weeks, accrual,
                                      design$randomisation_week, step)
    m_x1 <- cbind(rate = 1, m_x)
    v_w <- n * df_mix$weight
    v_levels <- sort(unique(v_stratum), method = "radix")
    strata <- lapply(stats::setNames(v_levels, v_levels), function(s) {
      v_rows <- v_stratum == s
      gamma_information_matrix(m_x1[v_rows, , drop = FALSE], v_w[v_rows] * m_e[v_rows, "aa"],
                               v_w[v_rows] * m_e[v_rows, "aeta"], v_w[v_rows] * m_e[v_rows, "etaeta"], v_names)
    })
    m_info <- Reduce(`+`, strata)
    pinv <- information_pseudo_inverse(m_info)
    df_events <- data.frame(stratum = v_levels, arm = sub("\\|.*$", "", v_levels),
                            cell = sub("^.*\\|", "", v_levels),
                            n = as.numeric(tapply(v_w, factor(v_stratum, v_levels), sum)),
                            events = as.numeric(tapply(v_w * m_e[, "events"], factor(v_stratum, v_levels), sum)),
                            stringsAsFactors = FALSE)
    endpoints[[outcome]] <- list(information = m_info, rank = pinv$rank, null_space = pinv$null_space,
                                 covariance = pinv$inverse, strata = strata, expected_events = df_events)
  }
  v_all_cols <- c(paste0("os::", rownames(component$original_os$res)),
                  paste0("pfs::", rownames(component$original_pfs$res)))
  if (is.null(targets)) {
    targets <- joint_target_matrix(list(
      design_treatment_targets(rownames(component$original_os$res), "os", design$enrolled_cells, design$biomarkers),
      design_treatment_targets(rownames(component$original_pfs$res), "pfs", design$enrolled_cells, design$biomarkers)),
      v_all_cols)
  }
  m_targets <- joint_target_matrix(list(targets), v_all_cols)
  # Block structure: each target must load on one endpoint.
  v_outcome_of_target <- vapply(rownames(m_targets), function(r) {
    v_out <- unique(sub("::.*$", "", v_all_cols[m_targets[r, ] != 0]))
    if (length(v_out) != 1) stop("Target ", r, " must load on exactly one endpoint")
    v_out
  }, character(1))
  v_identifiable <- logical(nrow(m_targets))
  names(v_identifiable) <- rownames(m_targets)
  v_blocks <- list()
  for (outcome in c("os", "pfs")) {
    v_rows <- which(v_outcome_of_target == outcome)
    if (!length(v_rows)) next
    v_cols <- paste0(outcome, "::", rownames(component[[paste0("original_", outcome)]]$res))
    m_a <- m_targets[v_rows, v_cols, drop = FALSE]
    m_null <- endpoints[[outcome]]$null_space
    if (ncol(m_null)) {
      v_proj <- apply(abs(m_a %*% m_null), 1, max) / pmax(1, sqrt(rowSums(m_a^2)))
      v_identifiable[v_rows] <- v_proj < 1e-6
    } else {
      v_identifiable[v_rows] <- TRUE
    }
    v_blocks[[outcome]] <- list(rows = v_rows, m_a = m_a,
                                covariance = m_a %*% endpoints[[outcome]]$covariance %*% t(m_a))
  }
  if (!all(v_identifiable)) {
    stop("Design ", design$name, " does not identify: ", paste(names(v_identifiable)[!v_identifiable], collapse = ", "))
  }
  m_cov <- matrix(0, nrow(m_targets), nrow(m_targets), dimnames = list(rownames(m_targets), rownames(m_targets)))
  for (blk in v_blocks) m_cov[blk$rows, blk$rows] <- blk$covariance
  if (cross_endpoint == "recolour" && length(v_blocks) == 2) {
    m_joint <- component$joint_covariance$covariance
    if (is.null(m_joint)) {
      message("No joint coefficient covariance in the component; cross-endpoint covariance set to zero")
      cross_endpoint <- "independent"
    } else {
      m_prior_t <- m_targets[, colnames(m_joint), drop = FALSE] %*% m_joint %*%
        t(m_targets[, colnames(m_joint), drop = FALSE])
      i_os <- v_blocks$os$rows; i_pfs <- v_blocks$pfs$rows
      m_r <- covariance_root(m_prior_t[i_os, i_os, drop = FALSE], inverse = TRUE) %*%
        m_prior_t[i_os, i_pfs, drop = FALSE] %*%
        covariance_root(m_prior_t[i_pfs, i_pfs, drop = FALSE], inverse = TRUE)
      m_cross <- covariance_root(m_cov[i_os, i_os, drop = FALSE]) %*% m_r %*%
        covariance_root(m_cov[i_pfs, i_pfs, drop = FALSE])
      m_cov[i_os, i_pfs] <- m_cross
      m_cov[i_pfs, i_os] <- t(m_cross)
    }
  }
  structure(list(design = design, n = n, accrual_weeks = accrual, targets = m_targets,
                 covariance = (m_cov + t(m_cov)) / 2, cross_endpoint = cross_endpoint,
                 endpoints = endpoints, identifiable = v_identifiable, patient_mix = df_mix),
            class = "voi_design_information")
}

#' Cox heuristic variance of a treatment-by-biomarker interaction
#'
#' For a log hazard-ratio interaction estimated from a 2 x 2 table of arm by
#' biomarker level, Var ~ sum over the four strata of 1 / expected events.
#'
#' @param info design_information() result.
#' @param outcome Endpoint.
#' @param biomarker Biomarker whose interaction is approximated.
#' @return Numeric scalar: sum(1 / d) over arm x biomarker level.
cox_interaction_heuristic <- function(info, outcome, biomarker) {
  df_ev <- info$endpoints[[outcome]]$expected_events
  j <- match(biomarker, info$design$biomarkers)
  if (is.na(j)) stop("Unknown biomarker ", biomarker)
  v_level <- substr(df_ev$cell, j, j)
  v_d <- tapply(df_ev$events, paste(df_ev$arm, v_level), sum)
  if (length(v_d) != 4) stop("The design does not have both levels of ", biomarker, " in both arms")
  sum(1 / v_d)
}

#' Simulate one endpoint of a trial under a design and a draw
#'
#' Patients are sampled from design_patient_mix(); each draws an entry time
#' (uniform over the accrual), a survival time from the draw's gamma model
#' conditional on surviving to the randomisation landmark, and is censored at
#' the analysis.
#'
#' @param design trial_design().
#' @param n Number of randomised patients.
#' @param draw Draw (split_joint_draw() formats).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param outcome Endpoint to simulate.
#' @return Data frame with the covariates, `start` (landmark), `stop`, `status`.
simulate_design_endpoint <- function(design, n, draw, component, reference_population,
                                     weights = NULL, outcome = "pfs") {
  draw <- split_joint_draw(draw, component)
  df_mix <- design_patient_mix(design, reference_population, weights, component, draw)
  fit <- component[[paste0("original_", outcome)]]
  v_par <- draw[[outcome]]
  v_rows <- sample.int(nrow(df_mix), n, replace = TRUE, prob = df_mix$weight)
  df_sim <- df_mix[v_rows, , drop = FALSE]
  v_eta <- gamma_linear_predictor(fit, v_par, df_sim)
  k <- exp(v_par[["shape"]])
  landmark <- design$randomisation_week
  v_s_land <- stats::pgamma(landmark, shape = k, rate = exp(v_eta), lower.tail = FALSE)
  v_time <- stats::qgamma(stats::runif(n) * v_s_land, shape = k, rate = exp(v_eta), lower.tail = FALSE)
  accrual <- design_accrual_weeks(design, n)
  v_follow <- design$min_follow_up_weeks + accrual - stats::runif(n, 0, accrual)
  df_sim$start <- landmark
  df_sim$stop <- pmin(v_time, landmark + v_follow)
  df_sim$status <- as.numeric(v_time <= landmark + v_follow)
  rownames(df_sim) <- NULL
  df_sim
}

#' Sampling covariance of design targets over simulated trials
#'
#' Simulates `n_trials` trials of one endpoint, refits the gamma model with
#' flexsurvreg (left-truncated at the landmark when it is positive) and returns
#' the empirical covariance of the target estimates.
#'
#' @param design,n,draw,component,reference_population,weights As in
#'   simulate_design_endpoint().
#' @param outcome Endpoint.
#' @param targets Target matrix restricted to this endpoint (columns
#'   "<outcome>::<coef>"); NULL uses design_treatment_targets().
#' @param n_trials Number of simulated trials.
#' @param seed RNG seed.
#' @return List with `covariance`, `mean`, `n_ok`, `n_failed` and `estimates`.
simulate_design_covariance <- function(design, n, draw, component, reference_population, weights = NULL,
                                       outcome = "pfs", targets = NULL, n_trials = 500, seed = 123L) {
  fit <- component[[paste0("original_", outcome)]]
  v_names <- rownames(fit$res)
  if (is.null(targets)) targets <- design_treatment_targets(v_names, outcome, design$enrolled_cells, design$biomarkers)
  m_a <- targets[, paste0(outcome, "::", v_names), drop = FALSE]
  v_rhs <- paste(deparse(stats::delete.response(stats::terms(fit$all.formulae$rate))[[2]]), collapse = " ")
  formula_sim <- stats::as.formula(paste(
    if (design$randomisation_week > 0) "survival::Surv(start, stop, status)" else "survival::Surv(stop, status)",
    "~", v_rhs))
  v_init <- split_joint_draw(draw, component)[[outcome]]
  set.seed(seed)
  m_est <- matrix(NA_real_, n_trials, nrow(m_a), dimnames = list(NULL, rownames(m_a)))
  for (b in seq_len(n_trials)) {
    df_sim <- simulate_design_endpoint(design, n, draw, component, reference_population, weights, outcome)
    refit <- tryCatch(suppressWarnings(flexsurv::flexsurvreg(
      formula_sim, data = df_sim, dist = "gamma",
      inits = c(exp(v_init[["shape"]]), exp(v_init[["rate"]]), v_init[-(1:2)]))),
      error = function(e) NULL)
    if (is.null(refit) || !identical(refit$opt$convergence, 0L)) next
    v_est <- stats::setNames(refit$res.t[, "est"], rownames(refit$res))
    if (!setequal(names(v_est), v_names)) next
    m_est[b, ] <- as.vector(m_a %*% v_est[v_names])
  }
  v_ok <- stats::complete.cases(m_est)
  list(covariance = stats::cov(m_est[v_ok, , drop = FALSE]), mean = colMeans(m_est[v_ok, , drop = FALSE]),
       n_ok = sum(v_ok), n_failed = sum(!v_ok), estimates = m_est[v_ok, , drop = FALSE])
}
