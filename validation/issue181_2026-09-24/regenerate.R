# Run from the repository root after committing the endpoint source changes.
pacman::p_load(devtools, readxl, dplyr, tableone, ggplot2, flexsurv,
               survival, survminer, gems, mstate, tidyverse, xtable,
               darthtools, dampack, mvtnorm, Matrix, here, voi)
for (script in c("02_setup_and_global_variables", "03_biomarker_strategies",
                 "04_parametric_survival_analysis", "05_basecase_input_parameters")) {
  source(paste0("analysis/", script, ".R"))
}
stopifnot(models$best_fit$os_distribution == "gamma",
          models$best_fit$pfs_distribution == "gamma")
for (helper in c("model_fun", "calculate_outcomes", "prediction_functions")) {
  source(paste0("R/", helper, ".R"))
}
for (script in c("06_sampling", "07_traces", "08_basecase_analysis",
                 "08b_enriched_population_analysis", "09_DSA", "10_PSA",
                 "11_EVPPIs", "12_scenario_EVPPIs")) {
  cat("\nISSUE181 START", script, format(Sys.time()), "\n")
  source(paste0("analysis/", script, ".R"))
  cat("\nISSUE181 COMPLETE", script, format(Sys.time()), "\n")
}
