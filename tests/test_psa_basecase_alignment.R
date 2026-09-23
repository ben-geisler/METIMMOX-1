# Structural PSA/base-case contract (issues #151, #180).
# The historical |PSA mean - deterministic value| <= 5 MCSE criterion is retired:
# E[f(theta)] need not equal f(E[theta]). Zero-uncertainty equality is checked
# by test_psa_zero_uncertainty.R; uncertain mean differences remain report diagnostics.
# No trial data or production caches are needed for this structural check.

# -----------------------------------------------------------------------------
# 1. Structural check: the PSA control curve is the joint model at Rx = control
# -----------------------------------------------------------------------------
sampling_expressions <- parse(here::here("analysis/06_sampling.R"))
extract_assignment <- function(name) {
  matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(matches) == 1L)
  sampling_expressions[[which(matches)]]
}

mock_env <- new.env(parent = baseenv())
# Since issue #156 the prediction helper resolves the draw through the cache
# accessor; a legacy-shaped sample list is returned as-is here.
mock_env$sampled_survival_models <- function(component, idx) component$samples[[idx]]
mock_env$weighted_survival_average <- function(curves, weights) as.vector(curves %*% (weights / sum(weights)))
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
