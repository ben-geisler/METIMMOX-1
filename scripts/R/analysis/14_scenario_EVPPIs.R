# Multi-scenario EVPPI analysis
# This script runs PSA and EVPPI analysis across multiple scenarios
# varying WTP thresholds and nivolumab costs for the single joint economic model.

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
# DEFINE BASE PARAMETERS FOR EVPPI ANALYSIS
# ===============================================================================
# Cost, utility, prevalence, and interaction parameters for the single economic model.

evppi_params_base <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS",
                       "c_test_CT", "u_np", "u_p", "c_other_last",
                       "c_test_blood", "c_other_visit", "c_other_baseline",
                       "c_other_follow")

# ===============================================================================
# RUN SCENARIOS
# ===============================================================================

# Assemble full evppi_params (matches 13_EVPPIs.R)
evppi_params <- c(evppi_params_base, "p_crp", "p_tmb_braf")

# Add interaction parameters if available
if (exists("get_interaction_evppi_params")) {
  interaction_evppi_params <- get_interaction_evppi_params()
  if (length(interaction_evppi_params) > 0) {
    evppi_params <- c(evppi_params, interaction_evppi_params)
    cat("Scenario EVPPI includes", length(interaction_evppi_params),
        "interaction parameters\n")
  }
}

# Build parameter groups for joint EVPPI analysis
# (matches the definitions in 08_sampling.R + 13_EVPPIs.R)
scenario_param_groups <- list(
  drug_costs = c("c_drug_nivo", "c_drug_FLOX"),
  test_costs = c("c_test_CT", "c_test_blood", "c_test_NGS"),
  other_costs = c("c_other_visit", "c_other_baseline",
                  "c_other_follow", "c_other_last"),
  all_costs = c("c_drug_nivo", "c_drug_FLOX", "c_test_CT",
                "c_test_blood", "c_test_NGS", "c_other_visit",
                "c_other_baseline", "c_other_follow", "c_other_last"),
  utilities = c("u_np", "u_p"),
  prevalence = c("p_crp", "p_tmb_braf")
)

# Add interaction parameter groups if available
if (exists("add_interaction_param_groups")) {
  scenario_param_groups <- add_interaction_param_groups(scenario_param_groups)
  cat("Total param_groups:", length(scenario_param_groups), "\n")
}

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
  evppi_params = evppi_params,
  param_groups = scenario_param_groups
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
    scen_data <- evppi_all_scenarios %>%
      filter(scenario_name == scen) %>%
      arrange(desc(evppi_population_millions)) %>%
      head(3)

    if (nrow(scen_data) > 0) {
      for (i in 1:nrow(scen_data)) {
        cat(sprintf("  %s: €%.2fM\n",
                    scen_data$parameter[i],
                    scen_data$evppi_population_millions[i]))
      }
    }
    cat("\n")
  }
}

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
      cat("\n  ", scenario_data$scenario_name[1], "\n")
      cat("  WTP: EUR", format(scenario_data$wtp[1], big.mark = ","), "\n")
      cat("  EVPI: EUR", round(scenario_data$evpi[1], 2), "\n")
      cat("\n  Top 3 Parameters:\n")
      top3 <- scenario_data[1:min(3, nrow(scenario_data)),
                            c("parameter", "evppi", "evppi_percent_of_evpi")]
      for (i in 1:nrow(top3)) {
        cat("    ", i, ". ", top3$parameter[i], ": EUR",
            round(top3$evppi[i], 2), " (", round(top3$evppi_percent_of_evpi[i], 1), "%)\n", sep = "")
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
evppi_cache <- list(
  all_scenario_results = all_scenario_results,
  evppi_all_scenarios = evppi_all_scenarios,
  scenarios = scenarios
)

scenario_evppi_cache_file <- here::here("data/tidy", paste0("scenario_evppi_results_", utility_source_label, ".rds"))
saveRDS(evppi_cache, file = scenario_evppi_cache_file)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", scenario_evppi_cache_file, "\n")
cat("\nSaved objects (in list):\n")
cat("  - all_scenario_results: Scenario PSA and EVPPI results\n")
cat("  - evppi_all_scenarios: Combined EVPPI data across scenarios\n")
cat("  - scenarios: Scenario definitions\n")
