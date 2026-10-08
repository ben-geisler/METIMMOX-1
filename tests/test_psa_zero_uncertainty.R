# Degenerate-covariance check required by #159: identical coefficients, economic
# parameters and population weights must reproduce every deterministic result.
# Trial data only; no production sampling/PSA cache is read or written.
if (!file.exists("data/tidy/METIMMOX.rds")) {
  cat("SKIP: confidential trial data absent.\n")
  quit(status = 0)
}
for (f in c("02_setup_and_global_variables.R", "03_biomarker_strategies.R",
            "04_parametric_survival_analysis.R", "05_basecase_input_parameters.R"))
  source(file.path("analysis", f))
for (f in c("joint_survival_sampling", "model_fun", "calculate_outcomes"))
  source(file.path("R", paste0(f, ".R")))
expressions <- parse("analysis/06_sampling.R")
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("generate_psa_population_averaged_predictions"))) eval(expr)
}
os <- models$best_fit$os
pfs <- models$best_fit$pfs
p <- length(os$opt$par) + length(pfs$opt$par)
draws <- draw_joint_survival_coefficients(os, pfs, matrix(0, p, p), 3L, 123L)
sampling_models <- list(joint = list(method = "mvn_joint_v2", draws = draws,
  original_os = os, original_pfs = pfs, n_samples = 3L,
  dist_os = os$dlist$name, dist_pfs = pfs$dlist$name))
df_base <- model_fun(l_params_base, determpsa = "det")
for (i in 1:3) {
  df_result <- model_fun(l_params_base, determpsa = "psa", sim_idx = i)
  stopifnot(!isTRUE(attr(df_result, "fallback_used")),
            identical(df_base$Strategy, df_result$Strategy),
            max(abs(df_base$Cost - df_result$Cost)) < 1e-8,
            max(abs(df_base$Effect - df_result$Effect)) < 1e-12)
}
cat("PASS: zero covariance and fixed parameters/weights reproduce all deterministic costs and QALYs.\n")
