# Multi-scenario EVPPI analysis with multi-model structure support
# This script runs PSA and EVPPI analysis across multiple scenarios
# varying WTP thresholds and nivolumab costs, for all three model structures.
#
# Model structures:
#   MODEL_STRUCTURE = 0: Joint (Model A) - all biomarkers + all interactions
#   MODEL_STRUCTURE = 1: Focused (Model B) - all biomarkers + one interaction
#   MODEL_STRUCTURE = 2: Separate (Model C) - one biomarker + one interaction
#
# Set RUN_ALL_MODELS <- TRUE to run for all model structures (generates 18 result sets)
# Set RUN_ALL_MODELS <- FALSE to run for current MODEL_STRUCTURE only (default)

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
source(here::here("scripts/R/functions/multi_model_cea.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}

# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# ===============================================================================
# MULTI-MODEL CONFIGURATION
# ===============================================================================

# Set to TRUE to run for all model structures, FALSE for current MODEL_STRUCTURE only
if (!exists("RUN_ALL_MODELS")) {
  RUN_ALL_MODELS <- FALSE
}

# Define model structures to process
if (RUN_ALL_MODELS) {
  model_structures_to_run <- c(0, 1, 2)
  model_labels <- c("joint", "focused", "separate")
  model_names <- c("Model A: Joint", "Model B: Focused", "Model C: Separate")
  cat("=== MULTI-MODEL MODE: Running for all 3 model structures ===\n\n")
} else {
  model_structures_to_run <- MODEL_STRUCTURE
  model_labels <- c("joint", "focused", "separate")[MODEL_STRUCTURE + 1]
  model_names <- c("Model A: Joint", "Model B: Focused", "Model C: Separate")[MODEL_STRUCTURE + 1]
  cat("=== SINGLE MODEL MODE: Running for MODEL_STRUCTURE =", MODEL_STRUCTURE, "===\n\n")
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
# Cost and utility parameters are common across all model structures.
# Prevalence, interaction params, and param_groups are model-structure-dependent
# and are assembled inside the model loop below.

evppi_params_base <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS",
                       "c_test_CT", "u_np", "u_p", "c_other_last",
                       "c_test_blood", "c_other_visit", "c_other_baseline",
                       "c_other_follow")

# ===============================================================================
# RUN SCENARIOS FOR EACH MODEL STRUCTURE
# ===============================================================================

# Initialize storage for all results across models
all_model_scenario_results <- list()
all_model_evppi_results <- data.frame()

for (m_idx in seq_along(model_structures_to_run)) {
  current_model_structure <- model_structures_to_run[m_idx]
  current_model_label <- model_labels[m_idx]
  current_model_name <- model_names[m_idx]

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("PROCESSING MODEL STRUCTURE:", current_model_structure, "(", current_model_label, ")\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Set MODEL_STRUCTURE globally
  MODEL_STRUCTURE <<- current_model_structure

  # Re-source scripts for this MODEL_STRUCTURE
  cat("Re-sourcing analysis scripts for model structure", current_model_structure, "...\n")
  suppressMessages(suppressWarnings({
    source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
    source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))
  }))

  # Check for existing PSA cache
  psa_cache_file <- here::here("data", "tidy", paste0("psa_obj_", current_model_label, ".rds"))

  if (file.exists(psa_cache_file)) {
    cat("Loading PSA cache from:", psa_cache_file, "\n")
    psa_obj_cached <- readRDS(psa_cache_file)
  } else {
    cat("WARNING: PSA cache not found for", current_model_label, "model.\n")
    cat("Please run 12_PSA.R with RUN_ALL_MODELS = TRUE first.\n")
    next
  }

  # Update base nivolumab cost from current parameters
  base_c_drug_nivo <- l_params_base$c_drug_nivo

  # Assemble full evppi_params for this model structure
  # (matches the logic in 13_EVPPIs.R)
  evppi_params <- c(evppi_params_base, "p_crp", "p_tmb_braf")
  if (current_model_structure != 1) {
    # Model B (focused) excludes TLR; Models A and C include it
    evppi_params <- c(evppi_params, "p_tlr")
  }

  # Add interaction parameters if available
  if (exists("get_interaction_evppi_params")) {
    interaction_evppi_params <- get_interaction_evppi_params(
      current_model_structure
    )
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
    prevalence = c("p_crp", "p_tlr", "p_tmb_braf")
  )

  # Add interaction parameter groups if available
  if (exists("add_interaction_param_groups")) {
    scenario_param_groups <- add_interaction_param_groups(
      scenario_param_groups, current_model_structure
    )
    cat("Total param_groups:", length(scenario_param_groups), "\n")
  }

  cat("\n=== RUNNING ALL SCENARIOS FOR", current_model_name, "===\n\n")

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

  # Store results for this model
  all_model_scenario_results[[current_model_label]] <- all_scenario_results

  # Compile EVPPI results and add model structure column
  evppi_results <- compile_evppi_results(all_scenario_results)

  if (nrow(evppi_results) > 0) {
    evppi_results$model_structure <- current_model_structure
    evppi_results$model_label <- current_model_label
    evppi_results$model_name <- current_model_name

    all_model_evppi_results <- rbind(all_model_evppi_results, evppi_results)
  }

  cat("\nCompleted scenarios for", current_model_name, "\n")
}

# ===============================================================================
# COMPILE AND ANALYZE RESULTS
# ===============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("COMPILING RESULTS ACROSS ALL MODELS AND SCENARIOS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Rename for consistency with EVPPIs.qmd expectations
evppi_all_scenarios <- all_model_evppi_results

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
      filter(scenario_name == scen, model_label == "joint") %>%
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
    group_by(model_name, scenario_name, wtp) %>%
    summarise(
      n_parameters = n(),
      total_evpi = first(evpi),
      mean_evppi = mean(evppi),
      max_evppi = max(evppi),
      top_parameter = parameter[which.max(evppi)],
      .groups = "drop"
    )

  print(summary_table)

  # Detailed results by model and scenario
  cat("\n=== DETAILED EVPPI RESULTS BY MODEL AND SCENARIO ===\n")

  for (model_label in unique(evppi_all_scenarios$model_label)) {
    model_data <- evppi_all_scenarios %>% filter(model_label == !!model_label)
    model_name <- unique(model_data$model_name)

    cat("\n", rep("-", 60), "\n", sep = "")
    cat(model_name, "\n")
    cat(rep("-", 60), "\n", sep = "")

    for (scenario_id in scenarios$scenario_id) {
      scenario_data <- model_data %>%
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
  }
} else {
  cat("No EVPPI results generated (likely due to zero EVPI in all scenarios)\n")
}

# ===============================================================================
# SAVE RESULTS
# ===============================================================================

# For backward compatibility, also create all_scenario_results from first/primary model
if (length(all_model_scenario_results) > 0) {
  all_scenario_results <- all_model_scenario_results[[1]]
}

# Save all results as a list in .rds format (consistent with other cache files)
evppi_cache <- list(
  all_scenario_results = all_scenario_results,
  all_model_scenario_results = all_model_scenario_results,
  evppi_all_scenarios = evppi_all_scenarios,
  scenarios = scenarios
)

saveRDS(evppi_cache, file = here::here("data/tidy/scenario_evppi_results.rds"))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to: data/tidy/scenario_evppi_results.rds\n")
cat("\nSaved objects (in list):\n")
cat("  - all_scenario_results: Results for primary model (backward compatibility)\n")
cat("  - all_model_scenario_results: Results for all model structures\n")
cat("  - evppi_all_scenarios: Combined EVPPI data with model_structure column\n")
cat("  - scenarios: Scenario definitions\n")

if (RUN_ALL_MODELS) {
  cat("\nProcessed", length(model_structures_to_run), "model structures x",
      nrow(scenarios), "scenarios =",
      length(model_structures_to_run) * nrow(scenarios), "total result sets\n")
}
