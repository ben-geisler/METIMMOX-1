# Multi-scenario EVPPI analysis
# This script runs PSA and EVPPI analysis across multiple scenarios
# varying WTP thresholds and nivolumab costs

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel, ggplot2, tidyr)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/psa_functions.R"))
source(here::here("scripts/R/functions/evppi_functions.R"))
source(here::here("scripts/R/functions/scenario_analysis.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}

# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# ===============================================================================
# DEFINE SCENARIOS
# ===============================================================================

# Get base case nivolumab cost from parameters
base_c_drug_nivo <- l_params_base$c_drug_nivo

# Define all scenarios
scenarios <- define_scenarios(
  base_wtp = WTP,  # Assuming WTP = 51000 from environment
  base_c_drug_nivo = base_c_drug_nivo,
  decreased_c_drug_nivo = 4641
)

cat("=== SCENARIO DEFINITIONS ===\n")
print(scenarios)
cat("\n")

# ===============================================================================
# DEFINE PARAMETERS FOR EVPPI ANALYSIS
# ===============================================================================

evppi_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS", 
                  "c_test_CT", "u_np", "u_p", "c_other_last",
                  "c_test_blood", "c_other_visit", "c_other_baseline", "c_other_follow")

# ===============================================================================
# RUN ALL SCENARIOS
# ===============================================================================

cat("\n=== RUNNING ALL SCENARIOS ===\n\n")

all_scenario_results <- run_all_scenarios(
  scenarios = scenarios,
  psa_params = NULL,  # Will be generated fresh for each scenario
  l_params_base = l_params_base,
  param_distributions = param_distributions,
  strategies = strategies,
  time_horizon = time_horizon,
  cl = cl,
  n_sim = n_sim,
  evppi_params = evppi_params
)

# ===============================================================================
# COMPILE AND ANALYZE RESULTS
# ===============================================================================

cat("\n=== COMPILING RESULTS ACROSS SCENARIOS ===\n")

# Compile all EVPPI results
evppi_all_scenarios <- compile_evppi_results(all_scenario_results)

# Summary statistics
cat("\nEVPPI Results Summary:\n")
if (nrow(evppi_all_scenarios) > 0) {
  summary_table <- evppi_all_scenarios %>%
    group_by(scenario_name, wtp) %>%
    summarise(
      n_parameters = n(),
      total_evpi = first(evpi),
      mean_evppi = mean(evppi),
      max_evppi = max(evppi),
      top_parameter = parameter[which.max(evppi)],
      .groups = "drop"
    )
  
  print(summary_table)
  
  # Detailed results by scenario
  cat("\n=== DETAILED EVPPI RESULTS BY SCENARIO ===\n")
  for (scenario_id in scenarios$scenario_id) {
    scenario_data <- evppi_all_scenarios %>%
      filter(scenario_id == !!scenario_id) %>%
      arrange(desc(evppi))
    
    if (nrow(scenario_data) > 0) {
      cat("\n", scenario_data$scenario_name[1], "\n")
      cat("WTP: €", format(scenario_data$wtp[1], big.mark = ","), "\n", sep = "")
      cat("EVPI: €", round(scenario_data$evpi[1], 2), "\n", sep = "")
      cat("\nTop 5 Parameters:\n")
      print(scenario_data[1:min(5, nrow(scenario_data)), 
                          c("parameter", "evppi", "evppi_percent_of_evpi")])
    }
  }
} else {
  cat("No EVPPI results generated (likely due to zero EVPI in all scenarios)\n")
}

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

# Save all results
save(all_scenario_results, evppi_all_scenarios, scenarios,
     file = here::here("data/tidy/scenario_evppi_results.RData"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to: data/tidy/scenario_evppi_results.RData\n")