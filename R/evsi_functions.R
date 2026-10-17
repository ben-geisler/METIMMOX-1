# ===============================================================================
# EXPECTED VALUE OF SAMPLE INFORMATION BY DESIGN (issue #202, paper 3 step 4)
# ===============================================================================
# A design d with n patients delivers estimates T of its identifiable targets
# (design_treatment_targets(): the treatment-effect coefficient combinations
# the design can estimate) with T_i ~ N(A theta_i, V_d(n)) for draw i, V_d(n)
# from design_information() (R/design_information.R). Two estimators:
#
#   * "gam_regression" (Strong, Brennan, Oakley et al., Med Decis Making 2015):
#     regress each policy's incremental NMB (versus column 1 of the NMB matrix)
#     on T with a GAM and take E_T[max_d fitted] - max_d E[fitted], through
#     voi::evsi() with a custom data-generating function (voi 1.0.3 runs its
#     own GAM, the same as voi::evppi(), on the simulated T). The GAM formula
#     follows evppi_gam_formula() (tensor-product cubic splines up to four
#     targets, additive cubic splines beyond), and the Monte Carlo SE is the SD
#     over B draws of the regression coefficients (voi's `se = TRUE`), seeded
#     with `seed` (analysis_seed), as in calculate_evppi_regression().
#   * "gaussian_approximation" (Jalal & Alarid-Escudero, Med Decis Making 2018),
#     multivariate form: fit the EVPPI GAM of incremental NMB on the targets
#     A theta, then evaluate it at the preposterior means. With prior
#     covariance S of A theta, the posterior mean after observing T has
#     covariance S (S + V)^-1 S; draws with that covariance are
#     mu + M (A theta_i - mu) with M = R^(1/2) S^(-1/2), R = S (S + V)^-1 S.
#     Exact when NMB is linear in A theta and the targets are normal; its SE
#     comes from B draws of the GAM coefficients.
#
# The standard normal variates behind T are drawn once per call (seeded) and
# reused for every n (common random numbers), so EVSI(n) is smooth in n.
# Zero and missing results follow #168: EVPI exactly zero bounds every EVSI at
# zero (method "zero_evpi_bound", SE 0); a small positive EVPI (<= min_evpi)
# gives NA with method "not_estimated"; n = 0 gives EVSI 0 (method
# "zero_information").
#
# Requires R/evppi_functions.R (calculate_evpi_from_nmb(), evppi_gam_formula()),
# R/joint_survival_sampling.R (covariance_root()) and, for
# evsi_design_from_information(), R/design_information.R.
# ===============================================================================

#' EVSI methods implemented here
EVSI_METHODS <- c("gam_regression", "gaussian_approximation")

#' Specification of a design for the EVSI engine
#'
#' @param name Design label.
#' @param targets Matrix of identifiable targets (rows) over the joint
#'   coefficient columns of `theta_draws` ("os::<coef>", "pfs::<coef>").
#' @param covariances Named list over sample sizes (names as.character(n)) of
#'   V_d(n): a q x q matrix common to every draw, or an array
#'   c(n_draws, q, q) with one covariance per draw.
#' @return List of class "voi_evsi_design".
evsi_design <- function(name, targets, covariances) {
  if (!is.matrix(targets) || is.null(rownames(targets)) || is.null(colnames(targets))) {
    stop("targets must be a matrix with row and column names")
  }
  q <- nrow(targets)
  for (nm in names(covariances)) {
    v_dim <- dim(covariances[[nm]])
    if (!(identical(v_dim, c(q, q)) || (length(v_dim) == 3 && identical(v_dim[2:3], c(q, q))))) {
      stop("Covariance for n = ", nm, " must be ", q, " x ", q, " or n_draws x ", q, " x ", q)
    }
  }
  structure(list(name = name, targets = targets, covariances = covariances), class = "voi_evsi_design")
}

#' Build an EVSI design from design_information() at a grid of sample sizes
#'
#' @param design trial_design().
#' @param n_grid Positive sample sizes.
#' @param draw Draw defining the expected information (typically
#'   component_mean_draw(component)).
#' @param component,reference_population,weights As in design_information().
#' @param ... Further design_information() arguments.
#' @return evsi_design() with the design's targets and V_d(n) at each n;
#'   attribute `information` holds the design_information() results.
evsi_design_from_information <- function(design, n_grid, draw, component, reference_population,
                                         weights = NULL, ...) {
  n_grid <- n_grid[n_grid > 0]
  infos <- lapply(n_grid, function(n) design_information(design, n, draw, component, reference_population,
                                                          weights, ...))
  names(infos) <- as.character(n_grid)
  out <- evsi_design(design$name, infos[[1]]$targets, lapply(infos, function(x) x$covariance))
  attr(out, "information") <- infos
  out
}

#' Target values A theta_i of every draw
#'
#' @param theta_draws Matrix of joint coefficient draws (row-aligned with the
#'   NMB matrix; columns named as the targets' columns).
#' @param targets Target matrix.
#' @return Matrix n_draws x q with the targets' row names as column names.
evsi_target_values <- function(theta_draws, targets) {
  v_used <- colnames(targets)[colSums(abs(targets)) > 0]
  v_missing <- setdiff(v_used, colnames(theta_draws))
  if (length(v_missing)) stop("theta_draws lacks target column(s): ", paste(v_missing, collapse = ", "))
  m <- theta_draws[, v_used, drop = FALSE] %*% t(targets[, v_used, drop = FALSE])
  colnames(m) <- rownames(targets)
  m
}

#' Simulated summary statistics T = A theta + L z for one sample size
#'
#' @param m_theta Target values (n_draws x q).
#' @param covariance V_d(n): q x q matrix or n_draws x q x q array.
#' @param m_z Standard normal variates (n_draws x q).
#' @return Data frame with columns T1..Tq.
evsi_summary_draws <- function(m_theta, covariance, m_z) {
  q <- ncol(m_theta)
  if (length(dim(covariance)) == 2) {
    m_l <- chol((covariance + t(covariance)) / 2)  # upper: t(m_l) %*% m_l = V
    m_t <- m_theta + m_z %*% m_l
  } else {
    if (dim(covariance)[1] != nrow(m_theta)) stop("Per-draw covariances must have one slice per draw")
    m_t <- m_theta
    for (i in seq_len(nrow(m_theta))) {
      m_t[i, ] <- m_theta[i, ] + as.vector(m_z[i, , drop = FALSE] %*% chol(covariance[i, , ]))
    }
  }
  df <- as.data.frame(m_t)
  names(df) <- paste0("T", seq_len(q))
  df
}

#' GAM right-hand side for q summary statistics (the EVPPI rule)
#'
#' @param v_names Variable names.
#' @return NULL (voi's tensor-product default, up to four variables) or the
#'   additive formula of evppi_gam_formula().
evsi_gam_formula <- function(v_names) evppi_gam_formula(v_names)

#' Regression EVSI for one design and sample size via voi::evsi()
#'
#' @param nmb_matrix NMB per draw (rows) and policy (columns; column 1 the
#'   reference).
#' @param m_theta Target values (n_draws x q).
#' @param covariance V_d(n).
#' @param m_z Standard normal variates (n_draws x q).
#' @param n Sample size (passed to the data-generating function).
#' @param seed Seed for the SE draws.
#' @param B Number of regression-coefficient draws for the SE.
#' @return List with `evsi`, `se`.
evsi_regression_one <- function(nmb_matrix, m_theta, covariance, m_z, n, seed, B) {
  df_inputs <- as.data.frame(m_theta)
  names(df_inputs) <- paste0("theta", seq_len(ncol(m_theta)))
  df_t <- evsi_summary_draws(m_theta, covariance, m_z)
  datagen <- function(inputs, n = 1) {
    if (nrow(inputs) != nrow(df_t)) stop("datagen: unexpected number of input rows")
    df_t
  }
  v_formula <- evsi_gam_formula(names(df_t))
  args <- list(outputs = nmb_matrix, inputs = df_inputs, datagen_fn = datagen, pars = names(df_inputs),
               n = n, method = "gam", se = TRUE, B = B)
  if (!is.null(v_formula)) args$gam_formula <- v_formula
  set.seed(seed)
  res <- do.call(voi::evsi, args)
  list(evsi = as.numeric(res$evsi[1]), se = as.numeric(res$se[1]))
}

#' GAM formula string used by the Gaussian approximation
#'
#' @param v_names Variable names.
#' @return Character right-hand side (voi's tensor default or additive splines).
evsi_ga_rhs <- function(v_names) {
  v_formula <- evsi_gam_formula(v_names)
  if (!is.null(v_formula)) return(v_formula)
  sprintf("te(%s, bs='cr'%s)", paste(v_names, collapse = ", "), if (length(v_names) >= 4) ", k=4" else "")
}

#' Fit the EVPPI GAMs of incremental NMB on the targets (Gaussian approximation)
#'
#' @param nmb_matrix NMB matrix (column 1 the reference).
#' @param m_theta Target values.
#' @return List of mgcv::gam fits, one per non-reference policy, and the
#'   variable names.
evsi_ga_models <- function(nmb_matrix, m_theta) {
  df <- as.data.frame(m_theta)
  v_names <- paste0("theta", seq_len(ncol(m_theta)))
  names(df) <- v_names
  f <- stats::as.formula(paste("y ~", evsi_ga_rhs(v_names)))
  models <- lapply(seq_len(ncol(nmb_matrix))[-1], function(d) {
    df$y <- nmb_matrix[, d] - nmb_matrix[, 1]
    mgcv::gam(f, data = df)
  })
  list(models = models, names = v_names)
}

#' Gaussian-approximation EVSI for one sample size
#'
#' @param ga evsi_ga_models() result.
#' @param m_theta Target values.
#' @param covariance V_d(n) (per-draw arrays are averaged over draws).
#' @param seed Seed for the SE draws.
#' @param B Number of GAM-coefficient draws for the SE.
#' @return List with `evsi`, `se`.
evsi_gaussian_one <- function(ga, m_theta, covariance, seed, B) {
  if (length(dim(covariance)) == 3) covariance <- apply(covariance, c(2, 3), mean)
  v_mu <- colMeans(m_theta)
  m_s <- stats::cov(m_theta)
  m_r <- m_s %*% solve(m_s + covariance) %*% m_s
  m_r <- (m_r + t(m_r)) / 2
  m_m <- covariance_root(m_r) %*% covariance_root(m_s, inverse = TRUE)
  m_pre <- sweep(sweep(m_theta, 2, v_mu) %*% t(m_m), 2, v_mu, "+")
  df_pre <- as.data.frame(m_pre)
  names(df_pre) <- ga$names
  evsi_from_fitted <- function(m_fit) {
    m_all <- cbind(0, m_fit)
    mean(apply(m_all, 1, max)) - max(colMeans(m_all))
  }
  l_x <- lapply(ga$models, function(model) stats::predict(model, newdata = df_pre, type = "lpmatrix"))
  m_fit <- vapply(seq_along(ga$models), function(d) as.vector(l_x[[d]] %*% stats::coef(ga$models[[d]])),
                  numeric(nrow(m_pre)))
  if (!is.matrix(m_fit)) m_fit <- matrix(m_fit, ncol = length(ga$models))
  set.seed(seed)
  l_beta <- lapply(ga$models, function(model) mvtnorm::rmvnorm(B, stats::coef(model), model$Vp))
  v_rep <- vapply(seq_len(B), function(b) {
    evsi_from_fitted(vapply(seq_along(ga$models), function(d) as.vector(l_x[[d]] %*% l_beta[[d]][b, ]),
                            numeric(nrow(m_pre))))
  }, numeric(1))
  list(evsi = evsi_from_fitted(m_fit), se = stats::sd(v_rep))
}

#' EVSI of a design over a grid of sample sizes
#'
#' @param nmb_matrix NMB per draw (rows) and policy (columns, column 1 the
#'   reference, e.g. standard of care).
#' @param theta_draws Joint coefficient draws row-aligned with `nmb_matrix`
#'   (bound by `model_idx`, never by draw number).
#' @param design evsi_design().
#' @param n_grid Sample sizes (0 allowed; others need a covariance in the design).
#' @param methods Subset of EVSI_METHODS.
#' @param seed Seed of the normal variates and the SE draws (analysis_seed).
#' @param B Coefficient draws for the SEs (1,000 as for EVPPI).
#' @param min_evpi EVPI at or below which (but above zero) EVSI is not estimated.
#' @return Data frame with `design`, `n`, `evsi`, `mcse`, `method`, `evpi`,
#'   `error`.
evsi_by_design <- function(nmb_matrix, theta_draws, design, n_grid, methods = EVSI_METHODS,
                           seed = get0("analysis_seed", ifnotfound = 123L), B = 1000, min_evpi = 0.01) {
  if (!inherits(design, "voi_evsi_design")) stop("design must come from evsi_design()")
  methods <- match.arg(methods, EVSI_METHODS, several.ok = TRUE)
  nmb_matrix <- as.matrix(nmb_matrix)
  if (nrow(nmb_matrix) != nrow(theta_draws)) stop("nmb_matrix and theta_draws must have the same rows")
  if (ncol(nmb_matrix) < 2) stop("At least two policies are needed")
  if (any(!is.finite(nmb_matrix))) stop("nmb_matrix must be finite")
  if (any(!is.finite(n_grid)) || any(n_grid < 0)) stop("n_grid must be non-negative")
  evpi <- calculate_evpi_from_nmb(nmb_matrix)
  make_rows <- function(n, evsi, mcse, method, error = "") {
    data.frame(design = design$name, n = n, evsi = evsi, mcse = mcse, method = method, evpi = evpi,
               error = error, stringsAsFactors = FALSE)
  }
  grid_rows <- function(method, evsi, mcse, error) {
    do.call(rbind, lapply(n_grid, function(n) make_rows(n, evsi, mcse, method, error)))
  }
  if (isTRUE(evpi == 0)) return(grid_rows("zero_evpi_bound", 0, 0, ""))
  if (!is.finite(evpi) || evpi <= min_evpi) {
    return(grid_rows("not_estimated", NA_real_, NA_real_, "EVPI below the estimator's threshold; not estimated"))
  }
  m_theta <- evsi_target_values(theta_draws, design$targets)
  set.seed(seed)
  m_z <- matrix(stats::rnorm(nrow(m_theta) * ncol(m_theta)), nrow(m_theta))
  ga <- if ("gaussian_approximation" %in% methods) evsi_ga_models(nmb_matrix, m_theta) else NULL
  rows <- list()
  for (n in n_grid) {
    if (n == 0) {
      rows[[length(rows) + 1]] <- make_rows(0, 0, 0, "zero_information")
      next
    }
    covariance <- design$covariances[[as.character(n)]]
    if (is.null(covariance)) stop("Design ", design$name, " has no covariance for n = ", n)
    for (method in methods) {
      result <- tryCatch(
        if (method == "gam_regression") evsi_regression_one(nmb_matrix, m_theta, covariance, m_z, n, seed, B)
        else evsi_gaussian_one(ga, m_theta, covariance, seed, B),
        error = function(e) e)
      rows[[length(rows) + 1]] <- if (inherits(result, "error")) {
        make_rows(n, NA_real_, NA_real_, method, conditionMessage(result))
      } else {
        make_rows(n, result$evsi, result$se, method)
      }
    }
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}
