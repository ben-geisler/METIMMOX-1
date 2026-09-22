# Regression test for WB-001 / BB-PSM-ORD-crp / BB-PSM-ORD-tmb_braf.
# Run from the repository root with Rscript.

library(here)

run_survival_ordering_test <- function() {
  caller_owned_value <- "preserved"
  source("analysis/02_setup_and_global_variables.R", local = environment())
  source("analysis/03_biomarker_strategies.R", local = environment())
  source("analysis/04_parametric_survival_analysis.R", local = environment())

  stopifnot(
    identical(caller_owned_value, "preserved"),
    length(time_points) == 521L,
    isTRUE(basecase_ordering_check$ordered),
    identical(basecase_ordering_check$n_violations, 0L),
    nrow(basecase_ordering_check$details) == 5L,
    all(basecase_ordering_check$details$ordered)
  )

  source("analysis/05_basecase_input_parameters.R", local = environment())
  source("R/model_fun.R", local = environment())
  source("R/calculate_outcomes.R", local = environment())

  basecase_result <- model_fun(l_params_base, determpsa = "det")
  stopifnot(
    nrow(basecase_result) == 3L,
    all(is.finite(basecase_result$Cost)),
    all(is.finite(basecase_result$Effect))
  )

  # The clamp is forbidden in deterministic mode; supplied curves opt in explicitly.
  crossed_params <- l_params_base
  crossed_params$p_pfs$control_PFS[2] <-
    crossed_params$p_os$control_OS[2] + 0.01

  deterministic_error <- tryCatch(
    {
      model_fun(crossed_params, determpsa = "det")
      NULL
    },
    error = identity
  )
  stopifnot(
    inherits(deterministic_error, "error"),
    grepl("Base-case survival ordering invariant failed",
          conditionMessage(deterministic_error), fixed = TRUE)
  )

  psa_warning <- NULL
  psa_result <- withCallingHandlers(
    model_fun(crossed_params, determpsa = "curves"),
    warning = function(w) {
      if (grepl("PFS > OS constraint enforced", conditionMessage(w),
                fixed = TRUE)) {
        psa_warning <<- w
        invokeRestart("muffleWarning")
      }
    }
  )
  stopifnot(
    !is.null(psa_warning),
    inherits(psa_warning, "survival_ordering_warning"),
    identical(psa_warning$strategy, get_control_strategy()),
    identical(psa_warning$subgroup, "control"),
    grepl("max excess=", conditionMessage(psa_warning), fixed = TRUE),
    all(is.finite(psa_result$Cost)),
    all(is.finite(psa_result$Effect))
  )

  cat("PASS: OS >= PFS at all 521 points for control and four subgroups; ",
      "the clamp requires PSA or explicit supplied-curve mode.\n", sep = "")
}

run_survival_ordering_test()
