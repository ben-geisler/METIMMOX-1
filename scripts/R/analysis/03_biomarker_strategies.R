#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(dplyr)

# Create binary biomarker variables with simplified names
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))

# Convert sex to factor BEFORE creating subsets
# This ensures all subset datasets (data_control, data_crp, etc.) have sex as factor
data$sex <- as.factor(data$sex)

# Create limited dataset with only the variables we need
data <- data %>%
  select(ID, PFSwk, Progression, OSwk, Death, Rx, crp, tlr, tmb_braf, Age, sex)

# Extract control and experimental groups
data_control <- subset(data, Rx == levels(data$Rx)[1])
data_experimental <- subset(data, Rx == levels(data$Rx)[2])

# Calculate prevalence rates
p_crp <- mean(data$crp, na.rm = TRUE)
p_tlr <- mean(data$tlr, na.rm = TRUE)
p_tmb_braf <- mean(data$tmb_braf, na.rm = TRUE)

# Create biomarker-guided datasets
# CRP combined dataset: CRP+ from experimental, CRP- from control
data_crp <- rbind(
  transform(subset(data_experimental, crp == 1), biomarker_status = "positive"),
  transform(subset(data_control, crp == 0), biomarker_status = "negative")
)

# TLR combined dataset: TLR+ from experimental, TLR- from control
data_tlr <- rbind(
  transform(subset(data_experimental, tlr == 1), biomarker_status = "positive"),
  transform(subset(data_control, tlr == 0), biomarker_status = "negative")
)

# TMB/BRAF combined dataset: TMB/BRAF+ from experimental, TMB/BRAF- from control
data_tmb_braf <- rbind(
  transform(subset(data_experimental, tmb_braf == 1), biomarker_status = "positive"),
  transform(subset(data_control, tmb_braf == 0), biomarker_status = "negative")
)

# Add strategy labels for model identification
data_control$strategy <- "control"
data_crp$strategy <- "crp"
data_tlr$strategy <- "tlr"
data_tmb_braf$strategy <- "tmb_braf"

# Define key vectors from the central economic-model configuration.
strategies <- get_strategies()
biomarkers <- get_biomarkers()

# Join keyed metadata and calculated values by strategy ID. This remains correct
# if the configured strategy order changes.
strategy_metadata <- get_strategy_metadata(strategies)
prevalence_by_strategy <- c(
  setNames(1, get_control_strategy()),
  setNames(vapply(biomarkers, function(x) get(paste0("p_", x)), numeric(1)),
           biomarkers)
)
n_by_strategy <- setNames(
  vapply(strategies, function(x) nrow(get(paste0("data_", x))), integer(1)),
  strategies
)
strategies_df <- transform(
  strategy_metadata,
  prevalence = unname(prevalence_by_strategy[id]),
  n_patients = unname(n_by_strategy[id])
)

# Print the final dataframe
print(strategies_df)

# Clean up intermediate variables
rm(list = setdiff(ls(), c(
  # Main datasets
  "data", "data_control", "data_crp", "data_tlr", "data_tmb_braf",
  # Model parameters and structure
  "strategies", "biomarkers", "strategies_df",
  # Biomarker prevalence
  "p_crp", "p_tlr", "p_tmb_braf",
  # Other essential variables from 02_setup_and_global_variables.R
  "time_horizon", "cl", "WTP", "DSA_mult", "n_samples", "n_sim",
  "analysis_seed", "dr",
  # Population parameters for EVPPI scaling
  "annual_incidence_norway", "research_horizon_years", "discount_rate_research",
  # Model configuration switches
  "UTILITY_SOURCE", "utility_source_label",
  # Model config helper functions
  "get_strategies", "get_control_strategy", "get_biomarkers",
  "get_model_configs", "get_current_model_config",
  "get_strategy_metadata", "strategy_display_name", "biomarker_cost_key",
  "biomarker_prevalence_key", "sync_biomarker_test_costs",
  "get_strategy_formula", "get_control_formula", "get_model_formulas",
  "print_model_config", "CONTROL_STRATEGY", "ALL_STRATEGIES", "ALL_BIOMARKERS",
  "STRATEGY_METADATA", "BIOMARKER_METADATA",
  # Cache-path helper functions
  "resolve_util_label", "sampling_cache_path", "psa_obj_path", "psa_params_path",
  "evppi_path", "scenario_evppi_path"
)))
