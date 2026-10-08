# ===============================================================================
# TEST: regression-based EVPPI estimator (issue #152)
# ===============================================================================
# The kNN estimator this replaces reported neighbourhood re-weighting bias as
# value of information: a pure-noise parameter returned a non-zero EVPPI, a
# group could come out below its own members, and every value was then capped
# at EVPI. These checks pin the replacement (voi::evppi, GAM regression) to
# known answers on a synthetic decision problem. No trial data or caches are
# needed; the voi package must be installed.
#
#   1. Analytic case. Two strategies whose incremental NMB is theta + noise,
#      with theta ~ N(mu, sigma^2) independent of the noise. The true
#      EVPPI(theta) is E[max(theta, 0)] - max(E[theta], 0)
#      = mu * Phi(mu/sigma) + sigma * phi(mu/sigma) - max(mu, 0),
#      recovered within Monte Carlo error. The noise inflates EVPI above
#      EVPPI, so the estimate must also stay below EVPI.
#   2. Pure noise. A parameter unrelated to the model has EVPPI within Monte
#      Carlo error of zero and is reported as estimated, never floored.
#   3. Groups. {theta, noise} carries the same information as theta alone, so
#      the joint EVPPI matches within Monte Carlo error and the
#      group-consistency check passes; a constructed violation is rejected.
#      Groups with more than four members use the additive GAM formula.
#   4. Seeds. The point estimate is deterministic; the standard error is
#      reproducible under the shared seed and changes with it.
#   5. Failures are NA with an error message, not zero.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        tests/test_evppi_estimator.R
# ===============================================================================

if (!requireNamespace("voi", quietly = TRUE)) {
  stop("The voi package is required for the EVPPI estimator tests.")
}
source(here::here("R", "evppi_functions.R"))

#' Evaluate `expr` with its printed output suppressed and return its value.
quiet <- function(expr) {
  out <- NULL
  utils::capture.output(out <- expr)
  out
}

# ---- Synthetic decision problem -------------------------------------------
set.seed(20260904L)
n <- 4000L
mu <- 300
sigma <- 1000
noise_sd <- 300

theta <- rnorm(n, mu, sigma)
v_inc_nmb <- theta + rnorm(n, 0, noise_sd)  # incremental NMB, strategy two vs one
df_psa_params <- data.frame(
  theta = theta,
  noise = rnorm(n),
  n2 = rnorm(n), n3 = rnorm(n), n4 = rnorm(n)
)
# With wtp = 1 and zero costs, effect equals NMB.
psa_obj <- list(
  cost = cbind(control = rep(0, n), experimental = rep(0, n)),
  effect = cbind(control = rep(0, n), experimental = v_inc_nmb)
)

analytic_evppi <- mu * pnorm(mu / sigma) + sigma * dnorm(mu / sigma) - max(mu, 0)
# Same quantity on the realised theta draws (removes the sampling error of
# theta itself, which the estimator cannot see).
sample_evppi <- mean(pmax(theta, 0)) - max(mean(theta), 0)

# ---- 1. Analytic case ------------------------------------------------------
res_theta <- quiet(calculate_evppi_regression(
  psa_obj, df_psa_params, "theta", wtp = 1, B = 500, seed = 123L
))
stopifnot(
  is.null(res_theta$error),
  is.finite(res_theta$evppi),
  is.finite(res_theta$evppi_se), res_theta$evppi_se > 0,
  identical(res_theta$method, "gam"),
  res_theta$n_sim == n
)
tol_theta <- 4 * res_theta$evppi_se + 0.02 * sample_evppi
stopifnot(
  abs(res_theta$evppi - sample_evppi) < tol_theta,
  abs(res_theta$evppi - analytic_evppi) < 0.10 * analytic_evppi,
  res_theta$evppi < res_theta$evpi
)

# ---- 2. Pure noise ---------------------------------------------------------
res_noise <- quiet(calculate_evppi_regression(
  psa_obj, df_psa_params, "noise", wtp = 1, B = 500, seed = 123L
))
stopifnot(
  is.null(res_noise$error),
  is.finite(res_noise$evppi),
  abs(res_noise$evppi) < max(4 * res_noise$evppi_se, 0.02 * res_noise$evpi)
)

# ---- 3. Groups -------------------------------------------------------------
res_group <- quiet(calculate_evppi_regression(
  psa_obj, df_psa_params, c("theta", "noise"), wtp = 1, B = 500, seed = 123L
))
stopifnot(
  is.null(res_group$error),
  abs(res_group$evppi - res_theta$evppi) <
    4 * sqrt(res_group$evppi_se^2 + res_theta$evppi_se^2) + 0.02 * res_theta$evpi
)

# More than four parameters switch to the additive formula and still recover
# the information carried by theta.
stopifnot(
  is.null(evppi_gam_formula(c("a", "b", "c", "d"))),
  identical(
    evppi_gam_formula(c("a", "b", "c", "d", "e")),
    "s(a, bs = 'cr') + s(b, bs = 'cr') + s(c, bs = 'cr') + s(d, bs = 'cr') + s(e, bs = 'cr')"
  )
)
res_five <- quiet(calculate_evppi_regression(
  psa_obj, df_psa_params, names(df_psa_params), wtp = 1, B = 500, seed = 123L
))
stopifnot(
  is.null(res_five$error),
  grepl("s(theta, bs = 'cr')", res_five$gam_formula, fixed = TRUE),
  abs(res_five$evppi - res_theta$evppi) <
    4 * sqrt(res_five$evppi_se^2 + res_theta$evppi_se^2) + 0.02 * res_theta$evpi
)

groups <- list(both = c("theta", "noise"), all_five = names(df_psa_params))
df_results <- quiet(run_evppi_analysis(
  psa_obj, df_psa_params, wtp = 1,
  evppi_params = names(df_psa_params), param_groups = groups,
  seed = 123L, B = 200
))
v_expected_cols <- c("parameter", "evppi", "evppi_se", "evpi",
                   "evppi_percent_of_evpi", "method", "n_params", "n_sim", "error")
stopifnot(
  nrow(df_results) == length(df_psa_params) + length(groups),
  all(v_expected_cols %in% names(df_results)),
  all(is.finite(df_results$evppi)),
  all(is.finite(df_results$evppi_se)),
  all(!nzchar(df_results$error)),
  all(paste0("[GROUP] ", names(groups)) %in% df_results$parameter),
  df_results$parameter[1] == "theta"  # ranked first
)
df_consistency <- attr(df_results, "group_consistency")
stopifnot(
  is.data.frame(df_consistency),
  nrow(df_consistency) == length(groups),
  all(df_consistency$max_member == "theta"),
  !any(df_consistency$violation)
)

# A group reported below its largest member beyond tolerance must be rejected.
df_bad <- df_results
df_bad$evppi[df_bad$parameter == "[GROUP] both"] <-
  df_bad$evppi[df_bad$parameter == "theta"] - 0.5 * df_bad$evpi[1]
df_violations <- check_evppi_group_consistency(df_bad, groups)
stopifnot(
  df_violations$violation[df_violations$group == "both"],
  !df_violations$violation[df_violations$group == "all_five"]
)
assertion <- tryCatch(assert_evppi_group_consistency(df_bad, groups),
                      error = function(e) e)
stopifnot(
  inherits(assertion, "error"),
  grepl("below its largest member", conditionMessage(assertion), fixed = TRUE)
)
warned <- tryCatch(
  quiet(run_evppi_analysis(
    psa_obj, df_psa_params, wtp = 1, evppi_params = "theta",
    param_groups = list(both = c("theta", "noise")), seed = 123L, B = 50,
    assert_groups = FALSE
  )),
  warning = function(w) w
)
stopifnot(!inherits(warned, "warning"))  # consistent groups do not warn

# ---- 4. Seeds --------------------------------------------------------------
#' EVPPI of theta (B = 100) with the given seed after unrelated RNG use; returns
#' the calculate_evppi_regression() result list.
seed_run <- function(seed) {
  invisible(runif(20))  # unrelated RNG use must not matter
  quiet(calculate_evppi_regression(psa_obj, df_psa_params, "theta", wtp = 1,
                                   B = 100, seed = seed))
}
run_a <- seed_run(7L)
run_b <- seed_run(7L)
run_c <- seed_run(8L)
stopifnot(
  identical(run_a$evppi, run_b$evppi),
  identical(run_a$evppi_se, run_b$evppi_se),
  identical(run_a$evppi, run_c$evppi),
  !identical(run_a$evppi_se, run_c$evppi_se)
)

# ---- 5. Failures are NA, never zero ----------------------------------------
res_missing <- quiet(calculate_evppi_regression(
  psa_obj, df_psa_params, "not_a_parameter", wtp = 1
))
stopifnot(!is.null(res_missing$error), is.na(res_missing$evppi))

df_psa_const <- df_psa_params
df_psa_const$const <- 1
res_const <- quiet(calculate_evppi_regression(psa_obj, df_psa_const, "const", wtp = 1))
stopifnot(!is.null(res_const$error), is.na(res_const$evppi))

# A dampack-style object stores effects as `effectiveness`.
psa_dampack <- list(cost = psa_obj$cost, effectiveness = psa_obj$effect)
res_dampack <- quiet(calculate_evppi_regression(
  psa_dampack, df_psa_params, "theta", wtp = 1, B = 50, seed = 123L
))
stopifnot(identical(res_dampack$evppi, res_theta$evppi))

cat(sprintf(paste0(
  "EVPPI estimator tests passed.\n",
  "  analytic EVPPI(theta) = %.2f, sample = %.2f, estimate = %.2f (SE %.2f), EVPI = %.2f\n",
  "  noise EVPPI = %.3f (SE %.3f); group {theta, noise} = %.2f (SE %.2f); ",
  "five-parameter additive = %.2f (SE %.2f)\n"),
  analytic_evppi, sample_evppi, res_theta$evppi, res_theta$evppi_se, res_theta$evpi,
  res_noise$evppi, res_noise$evppi_se, res_group$evppi, res_group$evppi_se,
  res_five$evppi, res_five$evppi_se
))
