# Regression test for WB-001 / BB-PSM-ORD-crp / BB-PSM-ORD-tmb_braf.
# Run from the repository root with Rscript.

library(here)

#' Source scripts 02-05 into a local environment and check base-case OS >= PFS
#' ordering, the deterministic error, the PSA clamp warning and the OS > 1 guard.
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

  df_basecase_result <- model_fun(l_params_base, determpsa = "det")
  stopifnot(
    nrow(df_basecase_result) == 3L,
    all(is.finite(df_basecase_result$Cost)),
    all(is.finite(df_basecase_result$Effect))
  )

  # The clamp is forbidden in deterministic mode; supplied curves opt in explicitly.
  # The crossed PFS curve is itself a valid survival curve (starts at 1,
  # nonincreasing, within [0, 1]), so the ordering check, not the survival-curve
  # contract of issue #60, is what rejects it.
  crossed_params <- l_params_base
  crossed_params$p_pfs$control_PFS <-
    pmin(1, crossed_params$p_os$control_OS + 0.01)

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
  df_psa_result <- withCallingHandlers(
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
    all(is.finite(df_psa_result$Cost)),
    all(is.finite(df_psa_result$Effect))
  )

  # An OS value above 1 used to pass (dead occupancy -0.1, finite QALYs); the
  # survival-curve contract now rejects it in every mode (issue #60).
  above_one <- l_params_base
  above_one$p_os$control_OS[2] <- 1.1
  for (mode in c("det", "curves")) {
    above_one_error <- tryCatch({ model_fun(above_one, determpsa = mode); NULL },
                                error = identity)
    stopifnot(
      inherits(above_one_error, "error"),
      grepl("Survival curve 'control_OS' leaves [0, 1]",
            conditionMessage(above_one_error), fixed = TRUE)
    )
  }

  cat("PASS: OS >= PFS at all 521 points for control and four subgroups; ",
      "the clamp requires PSA or explicit supplied-curve mode; ",
      "OS above 1 is rejected.\n", sep = "")
}

run_survival_ordering_test()
