# White-box regression tests for WB-011 / CS-008.

source("R/model_configs.R")
source("R/psa_functions.R")

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
stopifnot(result$dropped_count == 0)
stopifnot(result$n_sim == 4)
stopifnot(identical(result$retained_iterations, 1:4))
stopifnot(identical(result$model_idx, 1:4))
stopifnot(identical(result$failed_draw_policy, psa_failed_draw_policy()))
stopifnot(identical(
  result$replacement_model_policy,
  "cyclic_next_10_other_cached_models_v1"
))
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

# Replacement candidates explicitly cover every other cached survival model
# once. There is no out-of-range "unused" model bank.
stopifnot(
  identical(psa_replacement_candidates(1L, 4L), 2:4),
  identical(psa_replacement_candidates(3L, 4L), c(4L, 1L, 2L)),
  identical(psa_replacement_candidates(1L, 20L), 2:11),
  identical(psa_replacement_candidates(1L, 1L), integer(0))
)

# Every failed draw gets an independent candidate scan. The first failed draw
# below exhausts all three alternatives; the second must still try model 4 and
# succeed. This also verifies that replacement starts after each failed draw's
# original model rather than wrapping n_sim + 1 back to model 1.
psa_params_with_two_failures <- data.frame(
  sim = seq_len(4),
  draw_id = seq_len(4)
)
total_calls <- 0L
replacement_calls <- integer(0)
model_fun <- function(params, sim_idx, ...) {
  total_calls <<- total_calls + 1L
  initial_pass <- total_calls <= 4L
  if (!initial_pass) {
    replacement_calls <<- c(replacement_calls, sim_idx)
  }

  result <- data.frame(
    Strategy = "control",
    Cost = 10 * params$draw_id + sim_idx,
    Effect = params$draw_id
  )
  attr(result, "fallback_used") <- if (initial_pass) {
    params$draw_id %in% c(1L, 3L)
  } else {
    params$draw_id == 1L
  }
  result
}
replacement_warnings <- character(0)
output <- capture.output({
  result <- withCallingHandlers(
    run_psa_analysis(
      psa_params = psa_params_with_two_failures,
      l_params_base = c(make_params(), list(draw_id = 0L)),
      param_distributions = list(draw_id = list()),
      strategies = "control",
      time_horizon = 1,
      cl = 1,
      n_sim = 4,
      fallback_threshold = 1
    ),
    warning = function(w) {
      replacement_warnings <<- c(
        replacement_warnings,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )
})
stopifnot(
  identical(replacement_calls, c(2L, 3L, 4L, 4L)),
  identical(result$dropped_iterations, 1L),
  identical(result$retained_iterations, 2:4),
  # Draw 3 was re-run with cached model 4; model_idx records that (issue #156).
  identical(result$model_idx, c(2L, 4L, 4L)),
  identical(as.numeric(result$cost[, "control"]), c(22, 34, 44)),
  any(grepl("iteration 1 after 3 attempts", replacement_warnings)),
  any(grepl("Replacement model attempts: 4", output))
)

# An iteration that remains invalid after replacement is dropped as a whole;
# it is never replaced with column means. Use a permissive threshold here so
# the test reaches the replacement/drop phase.
psa_params_with_failure <- data.frame(
  sim = seq_len(4),
  invalid = c(TRUE, FALSE, FALSE, FALSE)
)
model_fun <- function(params, ...) {
  if (isTRUE(params$invalid)) {
    stop("deliberate unrecoverable draw")
  }
  data.frame(Strategy = "control", Cost = 1, Effect = 1)
}
drop_warnings <- character(0)
output <- capture.output({
  result <- withCallingHandlers(
    run_psa_analysis(
      psa_params = psa_params_with_failure,
      l_params_base = c(make_params(), list(invalid = FALSE)),
      param_distributions = list(invalid = list()),
      strategies = "control",
      time_horizon = 1,
      cl = 1,
      n_sim = 4,
      fallback_threshold = 1
    ),
    warning = function(w) {
      drop_warnings <<- c(drop_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
})
stopifnot(result$fallback_count == 1)
stopifnot(result$dropped_count == 1)
stopifnot(identical(result$dropped_iterations, 1L))
stopifnot(identical(result$retained_iterations, 2:4))
stopifnot(identical(result$model_idx, 2:4))
stopifnot(result$n_sim == 3)
stopifnot(nrow(result$cost) == 3, nrow(result$effect) == 3)
stopifnot(!anyNA(result$cost), !anyNA(result$effect))
stopifnot(any(grepl("Dropping 1 PSA iteration", drop_warnings)))
stopifnot(!any(grepl("column mean", c(output, drop_warnings), ignore.case = TRUE)))

cat("PSA fallback reporting and threshold tests passed.\n")
