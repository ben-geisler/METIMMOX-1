# ===============================================================================
# TEST: PSA strategy means are centred on the base case (issue #151)
# ===============================================================================
# The base case and the PSA must use ONE survival model. Before issue #151 the
# PSA control arm came from a separate age/sex-only bootstrap fitted on the 36
# control-arm patients while the base case predicted the joint model with
# Rx = control; the PSA control QALY mean sat 16 standard errors below the base
# case and every incremental quantity, CEAC, EVPI and EVPPI inherited the shift.
#
# Pass criterion (per strategy, for both cost and QALYs):
#   |PSA mean - base-case value| <= 5 * SE(PSA mean)
# where SE is the Monte Carlo standard error of the PSA mean. Five standard
# errors leaves room for the mild nonlinearity bias of a correctly centred
# partitioned survival model while catching a mis-specified control arm.
#
# Also checks, without any cache, that the PSA control-curve helper forces
# Rx = control for every patient.
#
# SKIPS (exit 0 with a notice) when the PSA cache is absent so the test is
# usable on a fresh clone. Requires the trial data and scripts 02-05.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" \
#        scripts/R/tests/test_psa_basecase_alignment.R
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here, survival, flexsurv, dplyr, tidyr, dampack)

SE_TOLERANCE <- 5

# -----------------------------------------------------------------------------
# 1. Structural check: the PSA control curve is the joint model at Rx = control
# -----------------------------------------------------------------------------
sampling_expressions <- parse(here::here("scripts/R/analysis/06_sampling.R"))
extract_assignment <- function(name) {
  matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(matches) == 1L)
  sampling_expressions[[which(matches)]]
}

mock_env <- new.env(parent = baseenv())
mock_env$extract_all_survival_probabilities <- function(prediction) {
  do.call(cbind, lapply(prediction$.pred, function(x) x$.pred_survival))
}
# A model whose prediction depends only on the treatment level, so any leak of
# a patient's actual arm into the control curve is visible.
mock_env$predict <- function(object, newdata, type, times) {
  value <- ifelse(newdata$Rx == "experimental", 0.9, 0.5)
  data.frame(.pred = I(lapply(value, function(v) {
    data.frame(.pred_survival = rep(v, length(times)))
  })))
}
eval(extract_assignment("generate_psa_population_averaged_predictions"),
     envir = mock_env)

patient_data <- data.frame(
  Age = c(50, 60, 70, 80),
  Rx = factor(c("control", "experimental", "control", "experimental"),
              levels = c("control", "experimental"))
)
control_curve <- mock_env$generate_psa_population_averaged_predictions(
  list(samples = list(list(os = list(model = list())))),
  outcome = "os", data_original = patient_data, time_points = 0:3
)
stopifnot(identical(control_curve, rep(0.5, 4)))
cat("PASS: PSA control curve is predicted with Rx = control for every patient\n")

# -----------------------------------------------------------------------------
# 2. Numerical check against the cached PSA
# -----------------------------------------------------------------------------
suppressMessages({
  source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
})
psa_file <- psa_obj_path(utility_source_label)
if (!file.exists(psa_file)) {
  cat("SKIP: PSA cache absent (", basename(psa_file), "); numerical ",
      "alignment check not run.\n", sep = "")
  quit(status = 0)
}

suppressMessages({
  source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
  source(here::here("scripts/R/analysis/04_parametric_survival_analysis.R"))
  source(here::here("scripts/R/analysis/05_basecase_input_parameters.R"))
  source(here::here("scripts/R/functions/model_fun.R"))
  source(here::here("scripts/R/functions/calculate_outcomes.R"))
  source(here::here("scripts/R/functions/cea_helpers.R"))
})

base_results <- run_basecase(verbose = FALSE)$base_results
psa_obj <- load_psa_cache(util_label = utility_source_label, verbose = FALSE)
stopifnot(!is.null(psa_obj))

comparison <- create_psa_basecase_comparison(psa_obj, base_results)
stopifnot(
  nrow(comparison) == 2L * length(get_strategies()),
  setequal(comparison$Strategy, get_strategies()),
  all(is.finite(comparison$Diff_SE))
)

cat("\nPSA mean versus base case (difference in Monte Carlo SEs):\n")
print(
  transform(
    comparison,
    Base_Case = signif(Base_Case, 6),
    PSA_Mean = signif(PSA_Mean, 6),
    PSA_SE = signif(PSA_SE, 4),
    Diff_SE = round(Diff_SE, 2)
  ),
  row.names = FALSE
)

# Informational: increments versus control. A shift common to every strategy
# (bootstrap nonlinearity bias of the shared joint model) cancels here; a
# control-specific shift, such as the separate control model removed in
# issue #151, does not.
incremental <- create_psa_incremental_comparison(psa_obj, base_results)
cat("\nPSA incremental means versus base-case increments (vs control):\n")
print(
  transform(
    incremental,
    Base_Case = signif(Base_Case, 6),
    PSA_Mean = signif(PSA_Mean, 6),
    PSA_SE = signif(PSA_SE, 4),
    Diff_SE = round(Diff_SE, 2)
  ),
  row.names = FALSE
)

violations <- comparison[abs(comparison$Diff_SE) > SE_TOLERANCE, , drop = FALSE]
if (nrow(violations) > 0L) {
  cat("\nFAIL: PSA means further than", SE_TOLERANCE,
      "standard errors from the base case:\n")
  print(violations, row.names = FALSE)
  qaly_shift <- comparison$Diff[comparison$Outcome == "QALYs"]
  if (all(qaly_shift > 0) || all(qaly_shift < 0)) {
    cat("\nNote: the QALY shift has the same sign for every strategy ",
        "(range ", sprintf("%+.4f", min(qaly_shift)), " to ",
        sprintf("%+.4f", max(qaly_shift)), "), i.e. a shift common to the ",
        "shared joint model rather than a control-specific one; see the ",
        "incremental table above.\n", sep = "")
  }
  stop("PSA strategy means are not within ", SE_TOLERANCE,
       " standard errors of the base case (issue #151 criterion).")
}

cat("\nPASS: all PSA strategy means (cost and QALYs) are within ",
    SE_TOLERANCE, " standard errors of the base case (n_sim = ",
    psa_obj$n_sim, ").\n", sep = "")
