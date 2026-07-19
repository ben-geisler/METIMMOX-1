# White-box regression tests for WB-011 / CS-008.

source("scripts/R/functions/psa_functions.R")

make_params <- function() {
  list(
    c_test_CRP = 1,
    c_test_NGS = 1,
    c_test_biomarker = list(crp = 1, tmb_braf = 1)
  )
}

psa_params <- data.frame(sim = seq_len(4))

# Successful iterations report and return a zero fallback rate.
model_fun <- function(params, ...) {
  data.frame(Strategy = "control", Cost = 1, Effect = 1)
}
output <- capture.output({
  result <- run_psa_analysis(
    psa_params = psa_params,
    l_params_base = make_params(),
    param_distributions = list(),
    strategies = "control",
    time_horizon = 1,
    cl = 1,
    n_sim = 4
  )
})
stopifnot(result$fallback_count == 0)
stopifnot(result$fallback_rate == 0)
stopifnot(result$fallback_threshold == 0.02)
stopifnot(any(grepl("PSA initial-pass fallback rate: 0/4 \\(0.00%", output)))

# One fallback in four iterations exceeds the default 2% threshold and stops
# before the replacement phase can conceal the high initial failure rate.
model_fun <- function(params, sim_idx, ...) {
  result <- data.frame(Strategy = "control", Cost = 1, Effect = 1)
  attr(result, "fallback_used") <- sim_idx == 1
  result
}
output <- capture.output({
  fallback_error <- tryCatch(
    run_psa_analysis(
      psa_params = psa_params,
      l_params_base = make_params(),
      param_distributions = list(),
      strategies = "control",
      time_horizon = 1,
      cl = 1,
      n_sim = 4
    ),
    error = identity
  )
})
stopifnot(inherits(fallback_error, "error"))
stopifnot(grepl("25.00% \\(1/4\\) exceeds", conditionMessage(fallback_error)))
stopifnot(any(grepl("PSA initial-pass fallback rate: 1/4 \\(25.00%", output)))

cat("PSA fallback reporting and threshold tests passed.\n")
