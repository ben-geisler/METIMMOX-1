# Scenario analysis functions for running multiple WTP and cost scenarios

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

#' Run PSA and EVPPI analysis for a single scenario
#' 
#' @param scenario_row Single row from scenarios data frame
#' @param psa_params PSA parameter samples
#' @param l_params_base Base parameter list
#' @param param_distributions Parameter distributions
#' @param strategies Vector of strategy names
#' @param time_horizon Time horizon
#' @param cl Cluster for parallel processing
#' @param n_sim Number of simulations
#' @param evppi_params Parameters to analyze for EVPPI
#' @return List with PSA object and EVPPI results
run_scenario_analysis <- function(scenario_row, psa_params, l_params_base, 
                                  param_distributions, strategies, 
                                  time_horizon, cl, n_sim, evppi_params) {
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Running Scenario:", scenario_row$scenario_name, "\n")
  cat("  WTP:", scenario_row$wtp, "\n")
  cat("  Nivolumab cost:", scenario_row$c_drug_nivo, "\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Update base parameters for this scenario
  scenario_params <- l_params_base
  scenario_params$c_drug_nivo <- scenario_row$c_drug_nivo
  
  # Update parameter distributions for this scenario
  scenario_param_dist <- param_distributions
  # Recalculate nivolumab cost distribution
  cv_costs <- 0.2
  nivo_gamma <- list(
    shape = 1/cv_costs^2,
    rate = (1/cv_costs^2) / scenario_row$c_drug_nivo
  )
  scenario_param_dist$c_drug_nivo <- c(list(dist = "gamma"), nivo_gamma)
  
  # Generate new PSA samples with updated distributions
  scenario_psa_params <- generate_psa_samples(scenario_param_dist, n_sim)
  
  # Run PSA
  psa_results <- run_psa_analysis(
    psa_params = scenario_psa_params,
    l_params_base = scenario_params,
    param_distributions = scenario_param_dist,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim
  )
  
  # Create PSA object
  psa_obj <- dampack::make_psa_obj(
    cost = as.data.frame(psa_results$cost),
    effect = as.data.frame(psa_results$effect),
    strategies = strategies,
    currency = "€"
  )
  
  cat("\nPSA Summary:\n")
  print(summary(psa_obj))
  
  # Run EVPPI analysis
  evppi_results <- run_evppi_analysis(
    psa_obj = psa_obj,
    psa_params = scenario_psa_params,
    wtp = scenario_row$wtp,
    evppi_params = evppi_params
  )
  
  # Add scenario information to results
  if (nrow(evppi_results) > 0) {
    evppi_results$scenario_id <- scenario_row$scenario_id
    evppi_results$scenario_name <- scenario_row$scenario_name
    evppi_results$wtp <- scenario_row$wtp
  }
  
  return(list(
    psa_obj = psa_obj,
    psa_params = scenario_psa_params,
    evppi_results = evppi_results,
    scenario_info = scenario_row
  ))
}

#' Run all scenarios and compile results
#' 
#' @param scenarios Data frame with scenario definitions
#' @param ... Additional parameters passed to run_scenario_analysis
#' @return List with all scenario results
run_all_scenarios <- function(scenarios, ...) {
  
  all_results <- list()
  
  for (i in 1:nrow(scenarios)) {
    scenario_results <- run_scenario_analysis(
      scenario_row = scenarios[i, ],
      ...
    )
    
    all_results[[scenarios$scenario_id[i]]] <- scenario_results
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