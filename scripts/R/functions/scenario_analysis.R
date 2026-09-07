# Scenario analysis functions for running multiple WTP and cost scenarios
# Optimized to avoid redundant PSA runs: scenarios sharing the same cost
# parameters reuse PSA results and only recalculate EVPPI at different WTPs.

#' Define scenarios for analysis
#'
#' @param base_wtp Base case WTP threshold
#' @param base_c_drug_nivo Base case nivolumab cost
#' @param decreased_c_drug_nivo Decreased nivolumab cost for scenarios
#' @return Data frame with scenario definitions
define_scenarios <- function(base_wtp = 51000,
                             base_c_drug_nivo = NULL,
                             decreased_c_drug_nivo = 4641) {

  scenarios <- data.frame(
    scenario_id = c("base", "s1_decreased_nivo",
                    "s2_wtp100k", "s3_wtp100k_decreased",
                    "s4_wtp150k", "s5_wtp150k_decreased"),
    scenario_name = c("Base Case",
                      "Decreased Nivolumab Cost",
                      "WTP €100,000",
                      "WTP €100,000 + Decreased Nivolumab",
                      "WTP €150,000",
                      "WTP €150,000 + Decreased Nivolumab"),
    wtp = c(base_wtp, base_wtp,
            100000, 100000,
            150000, 150000),
    c_drug_nivo = c(base_c_drug_nivo, decreased_c_drug_nivo,
                    base_c_drug_nivo, decreased_c_drug_nivo,
                    base_c_drug_nivo, decreased_c_drug_nivo),
    stringsAsFactors = FALSE
  )

  return(scenarios)
}

#' Run PSA for a unique cost configuration
#'
#' Generates PSA samples and runs the full PSA simulation for a given
#' nivolumab cost level. Returns PSA object and parameters that can be
#' reused across scenarios sharing the same cost configuration.
#'
#' @param c_drug_nivo Nivolumab cost for this configuration
#' @param l_params_base Base parameter list
#' @param param_distributions Parameter distributions
#' @param strategies Vector of strategy names
#' @param time_horizon Time horizon
#' @param cl Cluster for parallel processing
#' @param n_sim Number of simulations
#' @param seed RNG seed used for PSA and EVPPI stochastic blocks
#' @return List with psa_obj and psa_params
run_scenario_psa <- function(c_drug_nivo, l_params_base, param_distributions,
                             strategies, time_horizon, cl, n_sim,
                             seed = 123L) {

  cat("\n  Running PSA for nivolumab cost:", c_drug_nivo, "\n")

  # Update base parameters for this cost configuration
  scenario_params <- l_params_base
  scenario_params$c_drug_nivo <- c_drug_nivo

  # Update parameter distributions for this cost configuration. Since issue #154
  # the nivolumab price is a fixed input rather than a sampled one, so there is
  # normally no distribution to rescale; only re-centre it if the parameter set
  # in force does sample it.
  scenario_param_dist <- param_distributions
  if (!is.null(scenario_param_dist$c_drug_nivo)) {
    cv_costs <- 0.2
    scenario_param_dist$c_drug_nivo <- list(
      dist = "gamma",
      shape = 1 / cv_costs^2,
      rate = (1 / cv_costs^2) / c_drug_nivo
    )
  }

  # Add draw-aligned interaction coefficients when sampling models are available.
  interaction_coefs <- NULL
  if (exists("extract_interaction_coefficients") &&
      exists("sampling_models") && !is.null(sampling_models)) {
    interaction_coefs <- extract_interaction_coefficients(
      sampling_models = sampling_models,
      n_sim = n_sim
    )
    if (!is.null(interaction_coefs) && nrow(interaction_coefs) != n_sim) {
      interaction_coefs <- NULL
    }
  }

  psa_build <- build_psa_obj(
    l_params_base = scenario_params,
    param_distributions = scenario_param_dist,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim,
    seed = seed,
    additional_params = interaction_coefs
  )
  psa_obj <- psa_build$psa_obj
  scenario_psa_params <- psa_build$psa_params

  cat("  PSA Summary:\n")
  print(summary(psa_obj))

  return(list(
    psa_obj = psa_obj,
    psa_params = scenario_psa_params
  ))
}

#' Run EVPPI analysis for a single scenario using pre-computed PSA results
#'
#' @param scenario_row Single row from scenarios data frame
#' @param psa_obj Pre-computed PSA object
#' @param psa_params Pre-computed PSA parameter samples
#' @param evppi_params Parameters to analyze for EVPPI (should be fully
#'   assembled by caller, including prevalence and interaction params)
#' @param param_groups Optional named list of parameter groups for joint EVPPI
#' @param seed RNG seed used for grouped-parameter subsampling
#' @return List with PSA object, PSA params, EVPPI results, and scenario info
run_scenario_evppi <- function(scenario_row, psa_obj, psa_params,
                               evppi_params, param_groups = NULL,
                               seed = 123L) {

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Running EVPPI for Scenario:", scenario_row$scenario_name, "\n")
  cat("  WTP:", scenario_row$wtp, "\n")
  cat("  Nivolumab cost:", scenario_row$c_drug_nivo, "\n")
  cat("  (PSA results reused from cost group)\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Run EVPPI analysis with this scenario's WTP
  evppi_results <- run_evppi_analysis(
    psa_obj = psa_obj,
    psa_params = psa_params,
    wtp = scenario_row$wtp,
    evppi_params = evppi_params,
    param_groups = param_groups,
    seed = seed
  )

  # Add scenario information to results
  if (nrow(evppi_results) > 0) {
    evppi_results$scenario_id <- scenario_row$scenario_id
    evppi_results$scenario_name <- scenario_row$scenario_name
    evppi_results$wtp <- scenario_row$wtp
  }

  return(list(
    psa_obj = psa_obj,
    psa_params = psa_params,
    evppi_results = evppi_results,
    scenario_info = scenario_row
  ))
}

#' Run PSA and EVPPI analysis for a single scenario (legacy interface)
#'
#' Retained for backward compatibility. Runs a full PSA + EVPPI for one
#' scenario. For batch processing, prefer run_all_scenarios() which
#' automatically deduplicates PSA runs across scenarios.
#'
#' @param scenario_row Single row from scenarios data frame
#' @param psa_params PSA parameter samples (unused; generated internally)
#' @param l_params_base Base parameter list
#' @param param_distributions Parameter distributions
#' @param strategies Vector of strategy names
#' @param time_horizon Time horizon
#' @param cl Cluster for parallel processing
#' @param n_sim Number of simulations
#' @param evppi_params Parameters to analyze for EVPPI
#' @param param_groups Optional named list of parameter groups for joint EVPPI
#' @param seed RNG seed used for PSA and EVPPI stochastic blocks
#' @return List with PSA object and EVPPI results
run_scenario_analysis <- function(scenario_row, psa_params, l_params_base,
                                  param_distributions, strategies,
                                  time_horizon, cl, n_sim, evppi_params,
                                  param_groups = NULL, seed = 123L) {

  # Run PSA for this scenario's cost configuration
  psa_result <- run_scenario_psa(
    c_drug_nivo = scenario_row$c_drug_nivo,
    l_params_base = l_params_base,
    param_distributions = param_distributions,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim,
    seed = seed
  )

  # Run EVPPI using the PSA results
  return(run_scenario_evppi(
    scenario_row = scenario_row,
    psa_obj = psa_result$psa_obj,
    psa_params = psa_result$psa_params,
    evppi_params = evppi_params,
    param_groups = param_groups,
    seed = seed
  ))
}

#' Run all scenarios with optimized PSA reuse
#'
#' Groups scenarios by unique cost configurations (c_drug_nivo), runs PSA
#' once per group, then computes EVPPI for each WTP level within the group.
#' This avoids redundant PSA runs for scenarios that differ only in WTP.
#'
#' @param scenarios Data frame with scenario definitions
#' @param psa_params Ignored (retained for backward compatibility)
#' @param l_params_base Base parameter list
#' @param param_distributions Parameter distributions
#' @param strategies Vector of strategy names
#' @param time_horizon Time horizon
#' @param cl Cluster for parallel processing
#' @param n_sim Number of simulations
#' @param evppi_params Parameters to analyze for EVPPI
#' @param param_groups Optional named list of parameter groups for joint EVPPI
#' @param seed RNG seed used for PSA and EVPPI stochastic blocks
#' @return List with all scenario results
run_all_scenarios <- function(scenarios, psa_params = NULL, l_params_base,
                              param_distributions, strategies, time_horizon,
                              cl, n_sim, evppi_params, param_groups = NULL,
                              seed = 123L) {

  all_results <- list()

  # Identify unique cost configurations
  unique_nivo_costs <- unique(scenarios$c_drug_nivo)
  n_unique <- length(unique_nivo_costs)
  n_total <- nrow(scenarios)

  cat("Scenario optimization: ", n_total, " scenarios with ", n_unique,
      " unique cost configuration(s)\n", sep = "")
  cat("PSA runs needed: ", n_unique, " (saved ",
      n_total - n_unique, " redundant PSA runs)\n\n", sep = "")

  # Run PSA once per unique cost configuration, then EVPPI per scenario
  for (nivo_cost in unique_nivo_costs) {

    # Get all scenarios sharing this cost configuration
    cost_group <- scenarios[scenarios$c_drug_nivo == nivo_cost, ]

    cat("--- Cost group: c_drug_nivo =", nivo_cost, "---\n")
    cat("Scenarios in this group:", paste(cost_group$scenario_id, collapse = ", "), "\n")

    # Run PSA once for this cost group
    psa_result <- run_scenario_psa(
      c_drug_nivo = nivo_cost,
      l_params_base = l_params_base,
      param_distributions = param_distributions,
      strategies = strategies,
      time_horizon = time_horizon,
      cl = cl,
      n_sim = n_sim,
      seed = seed
    )

    # Run EVPPI for each scenario in this cost group (varying WTP)
    for (j in seq_len(nrow(cost_group))) {
      scenario_results <- run_scenario_evppi(
        scenario_row = cost_group[j, ],
        psa_obj = psa_result$psa_obj,
        psa_params = psa_result$psa_params,
        evppi_params = evppi_params,
        param_groups = param_groups,
        seed = seed
      )

      all_results[[cost_group$scenario_id[j]]] <- scenario_results
    }
  }

  return(all_results)
}

#' Compile EVPPI results across all scenarios
#'
#' @param all_results List of scenario results from run_all_scenarios
#' @return Combined data frame with all EVPPI results
compile_evppi_results <- function(all_results) {

  evppi_combined <- data.frame()

  for (scenario_id in names(all_results)) {
    scenario_evppi <- all_results[[scenario_id]]$evppi_results
    if (nrow(scenario_evppi) > 0) {
      evppi_combined <- rbind(evppi_combined, scenario_evppi)
    }
  }

  return(evppi_combined)
}

#' Deterministic alternative-utility scenarios (Table S8)
#'
#' Run deterministically on the base-case survival curves by table_s8.qmd;
#' they are not part of the PSA/EVPPI scenario set in define_scenarios().
#' The literature utility pair has NO citation yet: `source` is a placeholder
#' that must be replaced with the reference before publication (issue #157).
#' `u_decrement` is derived so the pair stays consistent with the derived-
#' utility convention of issue #154.
#'
#' @return Data frame with scenario_id, scenario_name, u_np, u_p, source.
define_utility_scenarios <- function() {
  data.frame(
    scenario_id = "u_literature",
    scenario_name = "Literature utilities (u_np 0.80, u_p 0.65)",
    u_np = 0.80,
    u_p = 0.65,
    source = "[CITATION PLACEHOLDER: reference for u_np 0.80 / u_p 0.65 to be added]",
    stringsAsFactors = FALSE
  )
}
