# White-box regression tests for WB-011 / CS-008.

source("R/model_configs.R")
source("R/psa_functions.R")

#' Minimal base parameter list (diagnostic prices only) for run_psa_analysis().
make_params <- function() {
  list(
    c_test_CRP = 1,
    c_test_NGS = 1,
    c_test_biomarker = list(crp = 1, tmb_braf = 1)
  )
}

df_psa_params <- data.frame(sim = seq_len(4))

# Successful iterations report and return a zero fallback rate.
#' Stub model: one control row with unit cost and effect, never a fallback.
model_fun <- function(params, ...) {
  data.frame(Strategy = "control", Cost = 1, Effect = 1)
}
v_output <- capture.output({
  result <- run_psa_analysis(
    psa_params = df_psa_params,
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
stopifnot(any(grepl("PSA initial-pass fallback rate: 0/4 \\(0.00%", v_output)))

# One fallback in four iterations exceeds the default 2% threshold and stops
# before the replacement phase can conceal the high initial failure rate.
#' Stub model: unit cost and effect, flagged as a fallback for draw 1 only.
model_fun <- function(params, sim_idx, ...) {
  df_result <- data.frame(Strategy = "control", Cost = 1, Effect = 1)
  attr(df_result, "fallback_used") <- sim_idx == 1
  df_result
}
v_output <- capture.output({
  fallback_error <- tryCatch(
    run_psa_analysis(
      psa_params = df_psa_params,
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
stopifnot(any(grepl("PSA initial-pass fallback rate: 1/4 \\(25.00%", v_output)))

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
df_psa_params_with_two_failures <- data.frame(
  sim = seq_len(4),
  draw_id = seq_len(4)
)
total_calls <- 0L
v_replacement_calls <- integer(0)
#' Stub model that records replacement calls; economic draws 1 and 3 fail
#' initially and draw 1 also fails with every replacement model.
model_fun <- function(params, sim_idx, ...) {
  total_calls <<- total_calls + 1L
  initial_pass <- total_calls <= 4L
  if (!initial_pass) {
    v_replacement_calls <<- c(v_replacement_calls, sim_idx)
  }

  df_result <- data.frame(
    Strategy = "control",
    Cost = 10 * params$draw_id + sim_idx,
    Effect = params$draw_id
  )
  attr(df_result, "fallback_used") <- if (initial_pass) {
    params$draw_id %in% c(1L, 3L)
  } else {
    params$draw_id == 1L
  }
  df_result
}
v_replacement_warnings <- character(0)
v_output <- capture.output({
  result <- withCallingHandlers(
    run_psa_analysis(
      psa_params = df_psa_params_with_two_failures,
      l_params_base = c(make_params(), list(draw_id = 0L)),
      param_distributions = list(draw_id = list()),
      strategies = "control",
      time_horizon = 1,
      cl = 1,
      n_sim = 4,
      fallback_threshold = 1
    ),
    warning = function(w) {
      v_replacement_warnings <<- c(
        v_replacement_warnings,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )
})
stopifnot(
  identical(v_replacement_calls, c(2L, 3L, 4L, 4L)),
  identical(result$dropped_iterations, 1L),
  identical(result$retained_iterations, 2:4),
  # Draw 3 was re-run with cached model 4; model_idx records that (issue #156).
  identical(result$model_idx, c(2L, 4L, 4L)),
  identical(as.numeric(result$cost[, "control"]), c(22, 34, 44)),
  any(grepl("iteration 1 after 3 attempts", v_replacement_warnings)),
  any(grepl("Replacement model attempts: 4", v_output))
)

# An iteration that remains invalid after replacement is dropped as a whole;
# it is never replaced with column means. Use a permissive threshold here so
# the test reaches the replacement/drop phase.
df_psa_params_with_failure <- data.frame(
  sim = seq_len(4),
  invalid = c(TRUE, FALSE, FALSE, FALSE)
)
#' Stub model that errors for the draw flagged `invalid` and otherwise returns
#' unit cost and effect.
model_fun <- function(params, ...) {
  if (isTRUE(params$invalid)) {
    stop("deliberate unrecoverable draw")
  }
  data.frame(Strategy = "control", Cost = 1, Effect = 1)
}
v_drop_warnings <- character(0)
v_output <- capture.output({
  result <- withCallingHandlers(
    run_psa_analysis(
      psa_params = df_psa_params_with_failure,
      l_params_base = c(make_params(), list(invalid = FALSE)),
      param_distributions = list(invalid = list()),
      strategies = "control",
      time_horizon = 1,
      cl = 1,
      n_sim = 4,
      fallback_threshold = 1
    ),
    warning = function(w) {
      v_drop_warnings <<- c(v_drop_warnings, conditionMessage(w))
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
stopifnot(any(grepl("Dropping 1 PSA iteration", v_drop_warnings)))
stopifnot(!any(grepl("column mean", c(v_output, v_drop_warnings), ignore.case = TRUE)))

# Non-finite values, duplicate/missing rows and missing columns all use the
# same threshold and replacement path on initial and replacement calls (#171).
df_valid <- data.frame(Strategy = c("control", "crp"), Cost = c(1, 2), Effect = c(1, 1.5))
invalid_results <- list()
for (column in c("Cost", "Effect")) for (value in c(NA_real_, NaN, Inf, -Inf)) {
  df_bad <- df_valid
  df_bad[[column]][2] <- value
  invalid_results[[length(invalid_results) + 1L]] <- df_bad
}
invalid_results <- c(invalid_results, list(df_valid[1, ], df_valid[c(1, 1), ],
  df_valid[c(1, 2, 2), ], df_valid[, c("Strategy", "Effect")]))
for (df_bad in invalid_results) {
  calls <- 0L
  model_fun <- function(params, sim_idx, ...) {
    calls <<- calls + 1L
    if (sim_idx <= 3) df_bad else df_valid
  }
  invisible(capture.output(err <- tryCatch(run_psa_analysis(
    data.frame(sim = 1:100), make_params(), list(), c("control", "crp"), 1, 1, 100),
    error = identity)))
  stopifnot(inherits(err, "error"), grepl("3.00% \\(3/100\\) exceeds", conditionMessage(err)),
            calls == 100L)

  # Exactly 2% is allowed, but both failures must be replaced; model 2 is
  # still invalid as a replacement and model 3 is used for both economic rows.
  calls <- 0L
  model_fun <- function(params, sim_idx, ...) {
    calls <<- calls + 1L
    if (sim_idx <= 2) df_bad else df_valid[2:1, ]
  }
  invisible(capture.output(result <- run_psa_analysis(
    data.frame(sim = 1:100), make_params(), list(), c("control", "crp"), 1, 1, 100)))
  stopifnot(result$fallback_count == 2, result$fallback_rate == 0.02,
            result$dropped_count == 0, result$n_sim == 100, calls == 103L,
            identical(result$model_idx[1:3], c(3L, 3L, 3L)),
            all(result$cost[, "control"] == 1), all(result$cost[, "crp"] == 2))
}

cat("PSA fallback reporting and threshold tests passed.\n")
