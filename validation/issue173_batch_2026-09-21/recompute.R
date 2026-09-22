# Recompute economic outcomes for #173/#171/#163 using the unchanged saved
# coefficient/parameter draws. Uses the existing validated closed-form gamma
# predictor and the explicit supplied-curve mode; never substitutes old outcomes.
# Run from the repository root BEFORE taking the fixed #173 snapshot.
options(metimmox.sampling_allow_regenerate = FALSE)
for (f in c("02_setup_and_global_variables.R", "03_biomarker_strategies.R",
            "04_parametric_survival_analysis.R", "05_basecase_input_parameters.R",
            "06_sampling.R")) source(file.path("analysis", f))
for (f in c("model_fun", "calculate_outcomes", "psa_functions", "cea_helpers",
            "evppi_functions", "scenario_analysis", "pfs_os_violation_diagnostics"))
  source(file.path("R", paste0(f, ".R")))

sampling_before <- tools::md5sum(sampling_cache_path())
old_psa <- readRDS("data/output/snapshots/psa_173_baseline_1b50e8d.rds")
draws <- readRDS(psa_params_path())
# Only the parameter values must match the baseline: a rerun can read the
# already refreshed pair. Rebind metadata in memory for the baseline check.
stopifnot(identical(cache_fingerprint(draws), old_psa$pair_metadata$params_hash))
attr(draws, "pair_metadata") <- old_psa$pair_metadata
validate_psa_pair(old_psa, draws)
stopifnot(old_psa$n_sim == n_sim, old_psa$fallback_count == 0,
          old_psa$dropped_count == 0, identical(draws$sim, seq_len(n_sim)),
          identical(draws$model_idx, seq_len(n_sim)),
          identical(old_psa$fingerprint_inputs$sampling_fingerprint, sampling_models$fingerprint))
# Confirm that saved economic/population draws equal today's seeded draws.
fresh_draws <- generate_psa_samples(param_distributions, n_sim, analysis_seed)
stopifnot(isTRUE(all.equal(draws[, names(fresh_draws)], fresh_draws,
                           check.attributes = FALSE, tolerance = 0)))

# Independent before-fix calculation verifies the fast curves against ALL
# 5,000 saved PSA rows, not just a few prediction examples.
before <- new.env(parent = globalenv())
for (file in c("R/model_fun.R", "R/calculate_outcomes.R")) {
  code <- system2("git", c("show", paste0("1b50e8d:", file)), stdout = TRUE)
  eval(parse(text = code), envir = before)
}
component <- get_joint_sampling_models(sampling_models)
population <- l_params_base$prediction_population
checked_draws <- unique(as.integer(round(seq(1, n_sim, length.out = 12))))
curve_check <- validate_direct_curves(component, psa_violation_subgroups(population),
  population, population, time_points, draws = checked_draws, tol = 1e-12,
  psa_params = draws, params = l_params_base)
cat("Maximum curve difference versus standard prediction:", max(curve_check$max_abs_diff), "\n")

scenarios <- define_scenarios(base_wtp = WTP,
  base_c_drug_nivo = l_params_base$c_drug_nivo, decreased_c_drug_nivo = 4641)
prices <- unique(scenarios$c_drug_nivo)
costs <- lapply(prices, function(x) matrix(NA_real_, n_sim, length(strategies),
  dimnames = list(NULL, strategies)))
effects <- costs[[1]]
max_old_cost_error <- max_old_effect_error <- max_standard_cost_error <- max_standard_effect_error <- 0
for (i in seq_len(n_sim)) {
  p <- l_params_base
  for (nm in names(param_distributions)) p[[nm]] <- draws[[nm]][i]
  weights <- model_population_weights(p, population)
  curves <- psa_draw_subgroup_curves(component, draws$model_idx[i],
    psa_violation_subgroups(population, weights = weights), time_points)
  for (nm in names(curves)) {
    p$p_os[[paste0(nm, "_OS")]] <- curves[[nm]]$os
    p$p_pfs[[paste0(nm, "_PFS")]] <- curves[[nm]]$pfs
  }
  p$population_weights <- weights
  p$population_curves <- NULL
  old_p <- sync_biomarker_test_costs(p) # legacy calculation requires its lookup
  original <- suppressWarnings(before$model_fun(old_p, determpsa = "psa"))
  old_cost_error <- max(abs(original$Cost - as.numeric(old_psa$cost[i, ])))
  old_effect_error <- max(abs(original$Effect - as.numeric(old_psa$effect[i, ])))
  stopifnot(old_cost_error < 1e-7, old_effect_error < 1e-11)
  max_old_cost_error <- max(max_old_cost_error, old_cost_error)
  max_old_effect_error <- max(max_old_effect_error, old_effect_error)
  for (j in seq_along(prices)) {
    p$c_drug_nivo <- prices[j]
    result <- suppressWarnings(model_fun(p, determpsa = "curves"))
    stopifnot(psa_outcomes_complete(result, strategies), !isTRUE(attr(result, "fallback_used")))
    costs[[j]][i, ] <- result$Cost
    if (j == 1L) effects[i, ] <- result$Effect
    else stopifnot(identical(effects[i, ], setNames(result$Effect, strategies)))
    if (i %in% checked_draws) {
      standard <- suppressWarnings(model_fun(p, determpsa = "psa", sim_idx = draws$model_idx[i]))
      cost_error <- max(abs(result$Cost - standard$Cost))
      effect_error <- max(abs(result$Effect - standard$Effect))
      stopifnot(!isTRUE(attr(standard, "fallback_used")), cost_error < 1e-7, effect_error < 1e-11)
      max_standard_cost_error <- max(max_standard_cost_error, cost_error)
      max_standard_effect_error <- max(max_standard_effect_error, effect_error)
    }
  }
  if (i %% 500L == 0L) cat("Recomputed", i, "of", n_sim, "draws for both prices\n")
}
provenance <- list(method = "validated_direct_gamma_supplied_curves",
  issues = c(173L, 171L, 163L), baseline_commit = "1b50e8d",
  sampling_md5 = unname(sampling_before), checked_standard_draws = checked_draws,
  max_curve_error = max(curve_check$max_abs_diff),
  max_old_cost_error = max_old_cost_error, max_old_effect_error = max_old_effect_error,
  max_standard_cost_error = max_standard_cost_error, max_standard_effect_error = max_standard_effect_error)
print(provenance)
jsonlite::write_json(provenance,
  "validation/issue173_batch_2026-09-21/recomputation_checks.json",
  auto_unbox = TRUE, pretty = TRUE, digits = 17)

make_pair <- function(j, additional = NULL) {
  params <- l_params_base
  params$c_drug_nivo <- prices[j]
  obj <- dampack::make_psa_obj(as.data.frame(costs[[j]]), as.data.frame(effects),
    strategies = strategies, currency = "EUR")
  obj$requested_n_sim <- n_sim
  obj$fallback_count <- obj$dropped_count <- 0L
  obj$dropped_iterations <- integer(0)
  obj$model_idx <- draws$model_idx
  obj$failed_draw_policy <- psa_failed_draw_policy()
  obj$replacement_model_policy <- "cyclic_next_10_other_cached_models_v1"
  obj$sampling_method <- component$method
  identity <- psa_cache_fingerprint(sampling_models$fingerprint, params,
    param_distributions, strategies, n_sim, analysis_seed, time_horizon, cl)
  obj$fingerprint <- identity$fingerprint
  obj$fingerprint_inputs <- identity$inputs
  obj$recomputation <- provenance
  pp <- draws
  attr(pp, "pair_metadata") <- NULL
  if (!is.null(additional)) pp <- cbind(pp, additional[pp$model_idx, , drop = FALSE])
  attr(pp, "seed") <- analysis_seed
  bind_psa_pair(obj, pp)
}
pair <- make_pair(match(l_params_base$c_drug_nivo, prices))
saveRDS(pair$psa_obj, psa_obj_path())
saveRDS(pair$psa_params, psa_params_path())
source("analysis/11_EVPPIs.R")

# Use the canonical scenario EVPPI functions and script-12 save/summary code.
config <- configure_evppi_analysis(param_distributions, param_groups)
interactions <- extract_interaction_coefficients(sampling_models, n_sim)
all_scenario_results <- list()
for (j in seq_along(prices)) {
  scenario_pair <- make_pair(j, interactions)
  for (k in which(scenarios$c_drug_nivo == prices[j])) {
    all_scenario_results[[scenarios$scenario_id[k]]] <- run_scenario_evppi(
      scenarios[k, ], scenario_pair$psa_obj, scenario_pair$psa_params,
      config$params, config$groups, seed = analysis_seed)
  }
}
script <- readLines("analysis/12_scenario_EVPPIs.R", warn = FALSE)
start <- grep("^evppi_all_scenarios <- compile_evppi_results", script)
stopifnot(length(start) == 1L)
eval(parse(text = script[start:length(script)]))
diagnostics <- run_pfs_os_violation_diagnostics(sampling_models, population, data, l_params_base,
  time_horizon, psa_params = pair$psa_params)
stopifnot(identical(sampling_before, tools::md5sum(sampling_cache_path())))
cat("All economic caches recomputed; survival sampling file unchanged.\n")
