# Quick preview version of scenario analysis for testing
# Runs with reduced PSA iterations (500 instead of 5000)
# This provides faster results for testing the biosimilar scenario report

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel, ggplot2, tidyr)

# Load setup scripts first
cat("Loading analysis scripts...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))
source(here::here("scripts/R/analysis/08_sampling.R"))

# Then load functions (after setup scripts, which may define dependencies)
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/psa_functions.R"))
source(here::here("scripts/R/functions/evppi_functions.R"))
source(here::here("scripts/R/functions/scenario_analysis.R"))
source(here::here("scripts/R/functions/multi_model_cea.R"))

# ===============================================================================
# PREVIEW CONFIGURATION
# ===============================================================================

# Reduce PSA iterations for quick preview
n_sim_preview <- 500  # Instead of 5000
cat("\n=== PREVIEW MODE ===\n")
cat("PSA iterations reduced to", n_sim_preview, "for faster testing\n\n")

# Run only Model A (joint) for quick preview
model_structures_to_run <- 0
model_labels <- "joint"
model_names <- "Model A: Joint"

# Focus on just base case and biosimilar scenarios
scenarios_preview <- data.frame(
  scenario_id = c("base", "s1_decreased_nivo"),
  scenario_name = c("Base Case", "Decreased Nivolumab Cost"),
  wtp = c(51000, 51000),
  c_drug_nivo = c(l_params_base$c_drug_nivo, 4641),
  stringsAsFactors = FALSE
)

cat("=== SCENARIO DEFINITIONS (PREVIEW) ===\n")
print(scenarios_preview)
cat("\n")

# ===============================================================================
# DEFINE PARAMETERS FOR EVPPI ANALYSIS
# ===============================================================================

evppi_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS",
                  "c_test_CT", "u_np", "u_p", "c_other_last",
                  "c_test_blood", "c_other_visit", "c_other_baseline", "c_other_follow")

# ===============================================================================
# RUN PREVIEW SCENARIOS
# ===============================================================================

cat("\n=== RUNNING PREVIEW SCENARIOS FOR Model A: Joint ===\n\n")

all_scenario_results_preview <- list()

for (i in 1:nrow(scenarios_preview)) {
  scenario_row <- scenarios_preview[i, ]

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
  cv_costs <- 0.2
  nivo_gamma <- list(
    shape = 1/cv_costs^2,
    rate = (1/cv_costs^2) / scenario_row$c_drug_nivo
  )
  scenario_param_dist$c_drug_nivo <- c(list(dist = "gamma"), nivo_gamma)

  # Generate new PSA samples with updated distributions (reduced iterations)
  scenario_psa_params <- generate_psa_samples(scenario_param_dist, n_sim_preview)

  # Run PSA
  psa_results <- run_psa_analysis(
    psa_params = scenario_psa_params,
    l_params_base = scenario_params,
    param_distributions = scenario_param_dist,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim_preview
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

  # Store results
  all_scenario_results_preview[[scenario_row$scenario_id]] <- list(
    psa_obj = psa_obj,
    psa_params = scenario_psa_params,
    evppi_results = evppi_results,
    scenario_info = scenario_row
  )
}

# ===============================================================================
# SAVE PREVIEW RESULTS
# ===============================================================================

# Create preview cache structure (compatible with biosimilar_scenario.qmd)
evppi_all_scenarios_preview <- do.call(rbind, lapply(all_scenario_results_preview, function(x) x$evppi_results))

if (nrow(evppi_all_scenarios_preview) > 0) {
  evppi_all_scenarios_preview$model_structure <- 0
  evppi_all_scenarios_preview$model_label <- "joint"
  evppi_all_scenarios_preview$model_name <- "Model A: Joint"

  # Add population-level EVPPI (matches 14_scenario_EVPPIs.R)
  evppi_all_scenarios_preview <- evppi_all_scenarios_preview %>%
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
}

# Create multi-model structure for compatibility
all_model_scenario_results_preview <- list(
  joint = all_scenario_results_preview
)

evppi_cache_preview <- list(
  all_scenario_results = all_scenario_results_preview,
  all_model_scenario_results = all_model_scenario_results_preview,
  evppi_all_scenarios = evppi_all_scenarios_preview,
  scenarios = scenarios_preview
)

# Save to preview cache file (utility-source-aware name, matching 14_scenario_EVPPIs.R)
preview_cache_file <- here::here("data/tidy",
                                 paste0("scenario_evppi_results_", utility_source_label,
                                        "_PREVIEW.rds"))
saveRDS(evppi_cache_preview, file = preview_cache_file)

cat("\n=== PREVIEW COMPLETE ===\n")
cat("Results saved to:", preview_cache_file, "\n")

cat("Summary:\n")
cat("- PSA iterations:", n_sim_preview, "(10% of full analysis)\n")
cat("- Model structures:", "1 (Model A only)\n")
cat("- Scenarios:", nrow(scenarios_preview), "(Base case + Biosimilar)\n")
cat("- Expected runtime:", "~5-10 minutes total\n")
