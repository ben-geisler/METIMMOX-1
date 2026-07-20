# Regression test for the correlated-resampling failure fallback.

# Evaluate only sample_correlated_survival() so this focused test does not run
# the cache-generation side effects in 06_sampling.R.
sampling_expressions <- parse("scripts/R/analysis/06_sampling.R")
is_sampling_function <- vapply(sampling_expressions, function(expr) {
  is.call(expr) &&
    identical(expr[[1]], as.name("<-")) &&
    identical(expr[[2]], as.name("sample_correlated_survival"))
}, logical(1))

stopifnot(sum(is_sampling_function) == 1)

test_env <- new.env(parent = baseenv())
eval(sampling_expressions[[which(is_sampling_function)]], envir = test_env)

# Let the two original fits succeed, then force the first resampled fit to fail.
test_env$fit_calls <- 0L
test_env$flexsurvreg <- with(test_env, function(...) {
  fit_calls <<- fit_calls + 1L
  if (fit_calls > 2L) {
    stop("forced resample failure")
  }
  list(
    coefficients = c(intercept = fit_calls),
    fit_id = fit_calls
  )
})

result <- test_env$sample_correlated_survival(
  formula_os = response ~ predictor,
  formula_pfs = response ~ predictor,
  data = data.frame(predictor = seq_len(3)),
  n_samples = 1
)

stopifnot(result$n_failed == 1)
stopifnot(length(result$samples) == 1)
stopifnot(!is.null(result$samples[[1]]))
stopifnot(isTRUE(result$samples[[1]]$failed))
stopifnot(identical(result$samples[[1]]$os$model, result$original_os))
stopifnot(identical(result$samples[[1]]$pfs$model, result$original_pfs))

cat("Sampling failure fallback regression test passed.\n")
