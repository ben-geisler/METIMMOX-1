# Multi-scenario EVPPI analysis
# This script runs PSA and EVPPI analysis across multiple scenarios
# varying WTP thresholds and nivolumab costs for the single joint economic model.

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel, ggplot2, tidyr, voi)

# Load functions
source(here::here("R/model_fun.R"))
source(here::here("R/calculate_outcomes.R"))
source(here::here("R/psa_functions.R"))
source(here::here("R/evppi_functions.R"))
source(here::here("R/scenario_analysis.R"))

# Validate the base PSA against the current sampling, code and economic inputs.
source(here::here("R/cea_helpers.R"))
psa_obj <- load_current_psa_cache()

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
# DEFINE BASE PARAMETERS FOR EVPPI ANALYSIS
# ===============================================================================
# Cost, utility, prevalence, and interaction parameters for the single economic model.
if (!exists("param_distributions") || !exists("param_groups")) {
  stop("Parameter distributions not found. Run 06_sampling.R first.")
}

# ===============================================================================
# RUN SCENARIOS
# ===============================================================================

# Assemble full EVPPI parameters and groups using the same configuration as script 11.
evppi_config <- configure_evppi_analysis(param_distributions, param_groups)
v_evppi_params <- evppi_config$params
scenario_param_groups <- evppi_config$groups
rm(evppi_config)
cat("Scenario EVPPI parameters:", length(v_evppi_params),
    "| parameter groups:", length(scenario_param_groups), "\n")

cat("\n=== RUNNING ALL SCENARIOS FOR THE JOINT ECONOMIC MODEL ===\n\n")

# run_all_scenarios groups scenarios by unique cost configurations,
# runs PSA once per group, then computes EVPPI at each WTP level.
# For 6 scenarios (3 WTP x 2 nivo costs), this runs 2 PSA instead of 6.
all_scenario_results <- run_all_scenarios(
  scenarios = scenarios,
  l_params_base = l_params_base,
  param_distributions = param_distributions,
  strategies = strategies,
  time_horizon = time_horizon,
  cl = cl,
  n_sim = n_sim,
  evppi_params = v_evppi_params,
  param_groups = scenario_param_groups,
  seed = analysis_seed
)

evppi_all_scenarios <- compile_evppi_results(all_scenario_results)

# ===============================================================================
# COMPILE AND ANALYZE RESULTS
# ===============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("COMPILING RESULTS ACROSS SCENARIOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# ===============================================================================
# ADD POPULATION-LEVEL EVPPI
# ===============================================================================

# Add population-level EVPPI to all scenario results
if (nrow(evppi_all_scenarios) > 0) {
  evppi_all_scenarios <- evppi_all_scenarios %>%
    mutate(
      evppi_population_millions = calculate_population_evppi(
        evppi,
        annual_incidence = annual_incidence_norway,
        research_horizon = research_horizon_years,
        discount_rate = discount_rate_research
      ),
      evpi_population_millions = calculate_population_evppi(
        evpi,
        annual_incidence = annual_incidence_norway,
        research_horizon = research_horizon_years,
        discount_rate = discount_rate_research
      )
    )

  # Print scenario-specific population EVPPI summary
  cat("\n=== Population-Level EVPPI by Scenario ===\n")
  cat(sprintf("Annual incidence: %d patients/year\n", annual_incidence_norway))
  cat(sprintf("Research horizon: %d years\n", research_horizon_years))
  cat(sprintf("Discount rate: %.1f%%\n\n", discount_rate_research * 100))

  for (scen in unique(evppi_all_scenarios$scenario_name)) {
    cat(sprintf("%s:\n", scen))
    df_scen_data <- evppi_all_scenarios %>%
      filter(scenario_name == scen) %>%
      arrange(desc(evppi_population_millions)) %>%
      head(3)

    if (nrow(df_scen_data) > 0) {
      for (i in seq_len(nrow(df_scen_data))) {
        cat(sprintf("  %s: €%.2fM\n",
                    df_scen_data$parameter[i],
                    df_scen_data$evppi_population_millions[i]))
      }
    }
    cat("\n")
  }
}

# Summary statistics
cat("\nEVPPI Results Summary:\n")
if (nrow(evppi_all_scenarios) > 0) {
  df_summary_table <- evppi_all_scenarios %>%
    group_by(scenario_name, wtp) %>%
    summarise(
      n_parameters = n(),
      total_evpi = first(evpi),
      mean_evppi = mean(evppi),
      max_evppi = max(evppi),
      top_parameter = parameter[which.max(evppi)],
      .groups = "drop"
    )

  print(df_summary_table)

  # Detailed results by scenario
  cat("\n=== DETAILED EVPPI RESULTS BY SCENARIO ===\n")

  for (scenario_id in scenarios$scenario_id) {
    df_scenario_data <- evppi_all_scenarios %>%
      filter(scenario_id == !!scenario_id) %>%
      arrange(desc(evppi))

    if (nrow(df_scenario_data) > 0) {
      cat("\n  ", df_scenario_data$scenario_name[1], "\n")
      cat("  WTP: EUR", format(df_scenario_data$wtp[1], big.mark = ","), "\n")
      cat("  EVPI: EUR", round(df_scenario_data$evpi[1], 2), "\n")
      cat("\n  Top 3 Parameters:\n")
      df_top3 <- df_scenario_data[1:min(3, nrow(df_scenario_data)),
                            c("parameter", "evppi", "evppi_percent_of_evpi")]
      for (i in seq_len(nrow(df_top3))) {
        cat("    ", i, ". ", df_top3$parameter[i], ": EUR",
            round(df_top3$evppi[i], 2), " (", round(df_top3$evppi_percent_of_evpi[i], 1), "%)\n", sep = "")
      }
    }
  }
} else {
  cat("No EVPPI results generated (likely due to zero EVPI in all scenarios)\n")
}

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

# Save all results as a list in .rds format (consistent with other cache files)
scenario_identity <- current_scenario_fingerprint()
evppi_cache <- list(
  fingerprint = scenario_identity$fingerprint,
  fingerprint_inputs = scenario_identity$inputs,
  all_scenario_results = all_scenario_results,
  evppi_all_scenarios = evppi_all_scenarios,
  scenarios = scenarios,
  failed_draw_policy = psa_failed_draw_policy(),
  # Provenance (issue #156): the sampling cache these scenario PSAs drew from.
  sampling_fingerprint = sampling_models$fingerprint,
  sampling_method = get_joint_sampling_models(sampling_models)$method
)

scenario_evppi_cache_file <- scenario_evppi_path()
evppi_cache$result_hash <- scenario_result_hash(evppi_cache)
saveRDS(evppi_cache, file = scenario_evppi_cache_file)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", scenario_evppi_cache_file, "\n")
cat("\nSaved objects (in list):\n")
cat("  - all_scenario_results: Scenario PSA and EVPPI results\n")
cat("  - evppi_all_scenarios: Combined EVPPI data across scenarios\n")
cat("  - scenarios: Scenario definitions\n")
