# Trial-data confirmation of #166 and draw-aligned ordering diagnostics.
# Sampling uses a temporary cache; no production cache is read or written.
if (!file.exists("data/tidy/METIMMOX.rds")) {
  cat("SKIP: confidential trial data absent.\n")
  quit(status = 0)
}
#' Source the pipeline (scripts 02-06) into a local environment with a temporary
#' cache directory and run both trial-data checks; stops on failure.
run_test <- function() {
  directory <- tempfile("population-contract-")
  dir.create(directory)
  old <- options(metimmox.cache_dir = directory)
  on.exit(options(old), add = TRUE)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  for (file in c("02_setup_and_global_variables.R", "03_biomarker_strategies.R",
                 "04_parametric_survival_analysis.R", "05_basecase_input_parameters.R", "06_sampling.R"))
    source(file.path("analysis", file), local = environment())
  for (file in c("model_fun", "calculate_outcomes", "psa_functions", "pfs_os_violation_diagnostics"))
    source(file.path("R", paste0(file, ".R")), local = environment())

  # Every subgroup receives the joint model's control treatment, equal drug
  # schedules, and zero screening cost: the review's identical-treatment probe.
  curves <- l_params_base$population_curves
  for (bm in get_biomarkers()) for (part in c("positive", "negative")) {
    v_rows <- curves[[bm]][[part]]$rows
    curves[[bm]][[part]]$os <- curves$control$os[, v_rows, drop = FALSE]
    curves[[bm]][[part]]$pfs <- curves$control$pfs[, v_rows, drop = FALSE]
  }
  p <- set_population_predictions(l_params_base, standardize_population_curves(curves,
    l_params_base$prediction_population, l_params_base$population_weights))
  p$l_nivo[] <- 0
  p$l_FLOX_exp <- p$l_FLOX_control
  p$c_test_CRP <- p$c_test_NGS <- 0
  for (bm in get_biomarkers()) for (mult in c(0.8, 1, 1.2)) {
    changed <- p; key <- paste0("p_", bm); changed[[key]] <- p[[key]] * mult
    df_result <- model_fun(changed)
    stopifnot(max(abs(df_result$Cost - df_result$Cost[1])) < 1e-10,
              max(abs(df_result$Effect - df_result$Effect[1])) < 1e-10)
  }
  cat("PASS: trial-data identical-treatment increments are zero at all prevalence DSA endpoints.\n")

  df_draws <- generate_psa_samples(param_distributions, 3L, analysis_seed)
  # Repeat a survival model under different population/utility draws, as a
  # replacement would: diagnostics must follow draw IDs, not model IDs alone.
  df_draws$model_idx <- c(1L, 3L, 3L)
  diagnostic <- run_pfs_os_violation_diagnostics(sampling_models, data_complete,
    data_complete, l_params_base, time_horizon, cache_path = file.path(directory, "diagnostics.rds"),
    n_draws = 3L, verbose = FALSE, psa_params = df_draws)
  for (i in 1:3) {
    p <- l_params_base
    for (nm in names(param_distributions)) p[[nm]] <- df_draws[[nm]][i]
    df_result <- suppressWarnings(model_fun(p, determpsa = "psa", sim_idx = df_draws$model_idx[i]))
    df_rows <- diagnostic$per_curve[diagnostic$per_curve$draw == i, ]
    v_qaly <- setNames(df_rows$qaly_clamped, df_rows$curve)
    v_expected <- c(control = v_qaly[["control"]], vapply(get_biomarkers(), function(bm)
      p[[paste0("p_", bm)]] * v_qaly[[paste0(bm, "_pos")]] +
      (1 - p[[paste0("p_", bm)]]) * v_qaly[[paste0(bm, "_neg")]], numeric(1)))
    stopifnot(max(abs(df_result$Effect - v_expected[df_result$Strategy])) < 1e-10)
  }
  cat("PASS: weighted diagnostic curves match predict() and reproduce PSA QALYs with repeated model indices.\n")
}
run_test()
