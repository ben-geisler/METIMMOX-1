# ===============================================================================
# PARAMETRIC SURVIVAL ANALYSIS WITH MODEL SELECTION
# ===============================================================================
# This script performs comprehensive parametric survival analysis by:
# 1. Fitting multiple parametric distributions across two model structures
# 2. Automatically selecting best-fitting models based on AIC
# 3. Generating population-level predictions for all treatment strategies
# ===============================================================================

# Load parametric model fitting functions
source("scripts/R/functions/para_model_fit.R")

# Load prediction functions
source("scripts/R/functions/prediction_functions.R")

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(survival, flexsurv, dplyr)

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

# Create time points sequence for prediction (from setup file)
time_points <- seq(0, time_horizon, by = 1)

# Ensure categorical variables are properly coded as factors
data$sex <- as.factor(data$sex)
data$crp <- as.factor(data$crp)
data$tlr <- as.factor(data$tlr)
data$tmb_braf <- as.factor(data$tmb_braf)

# Remove rows with missing values for modeling
data_complete <- data[complete.cases(data[, c("Age", "sex", "Rx", "crp", "tlr", "tmb_braf", "OSwk", "Death", "PFSwk", "Progression")]), ]

cat("Data preparation complete. Using", nrow(data_complete), "of", nrow(data), "observations.\n")

# ===============================================================================
# MODEL STRUCTURE DEFINITIONS
# ===============================================================================

# Define the two model structures for comparison:
# 1. Full model: includes age, sex, treatment, and biomarker interactions
# 2. No age/sex model: includes only treatment and biomarker interactions
model_formulas <- list(
  full = list(
    os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
    pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
  ),
  no_age_sex = list(
    os = Surv(OSwk, Death) ~ Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
    pfs = Surv(PFSwk, Progression) ~ Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
  )
)

# Define parametric distributions to test
distributions_to_test <- c("exponential", "weibull", "weibullph", "llogis", 
                           "lognormal", "gamma", "gompertz", "gengamma", "genf")

# ===============================================================================
# MODEL FITTING
# ===============================================================================

cat("Starting parametric model fitting process...\n")

# Initialize storage structure:
# models$[structure]$[outcome]$[distribution] = fitted flexsurvreg object
# models$[structure]$[outcome]_ic = information criteria table
# models$best_fit$[outcome] = best model overall
models <- list()
models$full <- list()
models$no_age_sex <- list()

# Fit models for each structure
for (structure_name in names(model_formulas)) {
  cat("Fitting", structure_name, "models...\n")
  
  # Fit OS models using fit_all_direct function
  cat("  - OS models\n")
  os_results <- fit_all_direct(
    fit_data = data_complete,
    fit_formula = model_formulas[[structure_name]]$os,
    fit_dists = distributions_to_test
  )
  
  # Fit PFS models using fit_all_direct function
  cat("  - PFS models\n")
  pfs_results <- fit_all_direct(
    fit_data = data_complete,
    fit_formula = model_formulas[[structure_name]]$pfs,
    fit_dists = distributions_to_test
  )
  
  # Store fitted models in organized structure
  models[[structure_name]]$os <- os_results$fitted_models
  models[[structure_name]]$pfs <- pfs_results$fitted_models
  
  # Store Kaplan-Meier fits for reference
  models[[structure_name]]$km_os <- os_results$km_all
  models[[structure_name]]$km_pfs <- pfs_results$km_all
  
  # Extract and store information criteria (AIC/BIC)
  os_ic <- extract_ic_single(os_results)
  pfs_ic <- extract_ic_single(pfs_results)
  
  models[[structure_name]]$os_ic <- os_ic
  models[[structure_name]]$pfs_ic <- pfs_ic
  
  # Identify best models within this structure
  best_os <- find_best_model(os_ic, criterion = "AIC")
  best_pfs <- find_best_model(pfs_ic, criterion = "AIC")
  
  if (!is.null(best_os)) {
    models[[structure_name]]$best_os_dist <- best_os$distribution
    models[[structure_name]]$best_os_aic <- best_os$criterion_value
    cat("    Best OS distribution:", best_os$distribution, "(AIC =", round(best_os$criterion_value, 2), ")\n")
  }
  
  if (!is.null(best_pfs)) {
    models[[structure_name]]$best_pfs_dist <- best_pfs$distribution
    models[[structure_name]]$best_pfs_aic <- best_pfs$criterion_value
    cat("    Best PFS distribution:", best_pfs$distribution, "(AIC =", round(best_pfs$criterion_value, 2), ")\n")
  }
}

# ===============================================================================
# GLOBAL BEST MODEL SELECTION
# ===============================================================================

cat("\nDetermining overall best models across all structures...\n")

# Combine AIC values from both structures to find global optimum
all_os_aic <- rbind(
  data.frame(structure = "full", models$full$os_ic),
  data.frame(structure = "no_age_sex", models$no_age_sex$os_ic)
)

all_pfs_aic <- rbind(
  data.frame(structure = "full", models$full$pfs_ic),
  data.frame(structure = "no_age_sex", models$no_age_sex$pfs_ic)
)

# Select globally best OS model
best_os_overall <- all_os_aic[!is.na(all_os_aic$AIC), ]
if (nrow(best_os_overall) > 0) {
  best_os_idx <- which.min(best_os_overall$AIC)
  best_os_structure <- best_os_overall$structure[best_os_idx]
  best_os_dist <- best_os_overall$Distribution[best_os_idx]
  
  # Store best model information
  models$best_fit$os <- models[[best_os_structure]]$os[[best_os_dist]]
  models$best_fit$os_structure <- best_os_structure
  models$best_fit$os_distribution <- best_os_dist
  models$best_fit$os_aic <- best_os_overall$AIC[best_os_idx]
  
  cat("Overall best OS model:", best_os_structure, "-", best_os_dist, "(AIC =", round(models$best_fit$os_aic, 2), ")\n")
}

# Select globally best PFS model
best_pfs_overall <- all_pfs_aic[!is.na(all_pfs_aic$AIC), ]
if (nrow(best_pfs_overall) > 0) {
  best_pfs_idx <- which.min(best_pfs_overall$AIC)
  best_pfs_structure <- best_pfs_overall$structure[best_pfs_idx]
  best_pfs_dist <- best_pfs_overall$Distribution[best_pfs_idx]
  
  # Store best model information
  models$best_fit$pfs <- models[[best_pfs_structure]]$pfs[[best_pfs_dist]]
  models$best_fit$pfs_structure <- best_pfs_structure
  models$best_fit$pfs_distribution <- best_pfs_dist
  models$best_fit$pfs_aic <- best_pfs_overall$AIC[best_pfs_idx]
  
  cat("Overall best PFS model:", best_pfs_structure, "-", best_pfs_dist, "(AIC =", round(models$best_fit$pfs_aic, 2), ")\n")
}

# ===============================================================================
# POPULATION-LEVEL PREDICTIONS
# ===============================================================================

cat("\nGenerating population-level predictions using best-fitting models...\n")

# Proceed only if both best models were successfully identified
if (!is.null(models$best_fit$os) && !is.null(models$best_fit$pfs)) {
  
  # Create model list for prediction functions
  best_models <- list(
    os = models$best_fit$os,
    pfs = models$best_fit$pfs
  )
  
  # Extract treatment level references
  exp_rx <- levels(data_complete$Rx)[2]  # Experimental treatment
  ctrl_rx <- levels(data_complete$Rx)[1]  # Control treatment
  
  # Calculate reference values for population-representative predictions
  ref_age <- mean(data_complete$Age, na.rm = TRUE)
  ref_sex <- get_modal_category(data_complete$sex)
  ref_crp <- get_modal_category(data_complete$crp)
  ref_tlr <- get_modal_category(data_complete$tlr)
  ref_tmb_braf <- get_modal_category(data_complete$tmb_braf)
  
  # Generate predictions for all strategies using prediction functions
  predictions <- generate_all_strategy_predictions(
    models = best_models,
    strategies_df = strategies_df,
    ref_age = ref_age,
    ref_sex = ref_sex,
    ref_crp = ref_crp,
    ref_tlr = ref_tlr,
    ref_tmb_braf = ref_tmb_braf,
    exp_rx = exp_rx,
    ctrl_rx = ctrl_rx,
    time_points = time_points,
    data_complete = data_complete
  )
  
  cat("Population-level predictions generated for all", length(strategies_df$id), "strategies.\n")
  
} else {
  cat("ERROR: Could not determine best-fitting models. Check model fitting results.\n")
  predictions <- list()
}

# ===============================================================================
# GENERATE OUTPUT FILES
# ===============================================================================

# Create reference values output
source("scripts/R/functions/ref_values_emm.R")

# Create parametric model fit comparison table
source("scripts/R/functions/para_model_fit_table.R")

# Create survival plots
source("scripts/R/functions/survival_plots.R")

# ===============================================================================
# SUMMARY
# ===============================================================================

cat("\n=== PARAMETRIC SURVIVAL ANALYSIS SUMMARY ===\n")
cat("Models fitted:\n")
cat("- Full model (Age + sex + biomarker interactions):", length(distributions_to_test), "distributions\n")
cat("- No age/sex model (biomarker interactions only):", length(distributions_to_test), "distributions\n")
cat("- Total models fitted:", 2 * 2 * length(distributions_to_test), "(2 structures × 2 outcomes × distributions)\n")

if (length(predictions) > 0) {
  cat("\nBest-fitting parametric models selected and predictions generated for all strategies.\n")
  cat("Output files generated:\n")
  cat("- scripts/R/functions/ref_values_emm.R: Reference values for population predictions\n")
  cat("- scripts/R/functions/para_model_fit_table.R: Parametric model comparison table with AIC/BIC\n") 
  cat("- scripts/R/functions/survival_plots.R: Biomarker-specific survival plots\n")
} else {
  cat("\nWarning: Best-fitting parametric models could not be determined.\n")
}

cat("\nFinal storage structure:\n")
cat("- models$full: Full model results for all distributions\n")
cat("- models$no_age_sex: No age/sex model results for all distributions\n")
cat("- models$best_fit: Best-fitting parametric models with metadata\n")
cat("- predictions: Strategy-level predictions using best parametric models\n")