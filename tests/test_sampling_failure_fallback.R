# Regression test for the survival coefficient sampler's failure handling.
#
# Before issue #156 the nonparametric bootstrap substituted the original fit
# when a resampled fit failed (flagged `failed = TRUE`). The multivariate-normal
# sampler cannot fail draw by draw; its only failure mode is a fit without a
# usable covariance matrix, which must stop loudly rather than fall back.

# Evaluate only the sampler so this focused test does not run the
# cache-generation side effects in 06_sampling.R.
sampling_expressions <- parse("analysis/06_sampling.R")
#' The single top-level `name <- ...` expression of analysis/06_sampling.R, so a
#' production helper can be evaluated without running the script.
extract_assignment <- function(name) {
  v_matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(v_matches) == 1L)
  sampling_expressions[[which(v_matches)]]
}

test_env <- new.env(parent = globalenv())
eval(extract_assignment("SAMPLING_METHOD"), envir = test_env)
eval(extract_assignment("sample_survival_coefficients"), envir = test_env)

#' Stand-in for flexsurvreg(): a three-parameter fit (shape, rate, predictor)
#' carrying the covariance matrix `cov`.
mock_fit <- function(cov) {
  v_par_names <- c("shape", "rate", "predictor")
  v_est <- c(0.1, -3, 0.2)
  m_res <- matrix(v_est, ncol = 4, nrow = 3,
                dimnames = list(v_par_names, c("est", "L95%", "U95%", "se")))
  list(
    res = m_res, res.t = m_res, coefficients = setNames(v_est, v_par_names),
    opt = list(par = v_est), cov = cov,
    optpars = 1:3, fixedpars = integer(0),
    dlist = list(pars = c("shape", "rate"), inv.transforms = list(exp, exp))
  )
}

# 1. Converged fits: every draw is produced, none fails, nothing is substituted.
test_env$flexsurvreg <- function(...) mock_fit(diag(c(0.01, 0.04, 0.02)))
result <- test_env$sample_survival_coefficients(
  formula_os = response ~ predictor,
  formula_pfs = response ~ predictor,
  data = data.frame(predictor = seq_len(3)),
  n_samples = 4
)
stopifnot(
  result$n_failed == 0L,
  is.null(result$samples),
  identical(dim(result$draws$os), c(4L, 3L)),
  identical(dim(result$draws$pfs), c(4L, 3L)),
  identical(colnames(result$draws$os), c("shape", "rate", "predictor")),
  all(is.finite(result$draws$os)), all(is.finite(result$draws$pfs)),
  # OS and PFS draws are independent, not copies of one another.
  !identical(result$draws$os, result$draws$pfs),
  nrow(result$coefficient_summary) == 6L
)

# 2. A fit without a covariance matrix (non-converged) must stop, naming the
#    outcome, instead of silently substituting the point estimate.
test_env$fit_calls <- 0L
test_env$flexsurvreg <- with(test_env, function(...) {
  fit_calls <<- fit_calls + 1L
  mock_fit(if (fit_calls == 2L) matrix(NA_real_, 3, 3) else diag(3) * 0.01)
})
failure <- tryCatch(
  test_env$sample_survival_coefficients(
    formula_os = response ~ predictor,
    formula_pfs = response ~ predictor,
    data = data.frame(predictor = seq_len(3)),
    n_samples = 2
  ),
  error = identity
)
stopifnot(
  inherits(failure, "error"),
  grepl("PFS model", conditionMessage(failure)),
  grepl("non-converged", conditionMessage(failure))
)

cat("Sampling failure handling regression test passed.\n")
