# ===============================================================================
# PARAMETRIC SURVIVAL ANALYSIS WITH MODEL SELECTION
# ===============================================================================
# This script performs comprehensive parametric survival analysis by:
# 1. Fitting multiple parametric distributions across model structures
# 2. Automatically selecting best-fitting models based on AIC
# 3. Generating population-level predictions for all treatment strategies
# 
# Testing and validation outputs are handled in scripts/R/tests/para_models.R
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman, here)
p_load(survival, flexsurv, dplyr)

# ===============================================================================
# MODEL SELECTION CONFIGURATION
# ===============================================================================

# CONFIGURATION SWITCH: Controls which model structures to consider
# 0 = Full model (age- and sex-adjusted) only
# 1 = Both full model and model without age- and sex-adjustments
USE_BOTH_MODELS <- 0

# Load parametric model fitting functions
para_model_fit_path <- here::here("scripts", "R", "functions", "para_model_fit.R")
if (file.exists(para_model_fit_path)) {
  source(para_model_fit_path)
} else {
  message("File not found: ", para_model_fit_path)
}
rm(para_model_fit_path)

# Load prediction functions
prediction_functions_path <- here::here("scripts", "R", "functions", "prediction_functions.R")
if (file.exists(prediction_functions_path)) {
  source(prediction_functions_path)
} else {
  message("File not found: ", prediction_functions_path)
}
rm(prediction_functions_path)

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

# Create time points sequence for prediction (from setup file)
time_points <- seq(0, time_horizon, by = 1)

# Ensure categorical variables are properly coded as factors
if ("sex" %in% names(data) && !is.null(data$sex)) {
  data$sex <- as.factor(data$sex)
}
if ("crp" %in% names(data) && !is.null(data$crp)) {
  data$crp <- as.factor(data$crp)
}
if ("tlr" %in% names(data) && !is.null(data$tlr)) {
  data$tlr <- as.factor(data$tlr)
}
if ("tmb_braf" %in% names(data) && !is.null(data$tmb_braf)) {
  data$tmb_braf <- as.factor(data$tmb_braf)
}

# Remove rows with missing values for modeling
data_complete <- data[complete.cases(data[, c("Age", "sex", "Rx", "crp", "tlr", "tmb_braf", "OSwk", "Death", "PFSwk", "Progression")]), ]

# ===============================================================================
# MODEL STRUCTURE DEFINITIONS
# ===============================================================================

# Define model structures based on configuration switch
if (USE_BOTH_MODELS == 0) {
  # Only use full model (age- and sex-adjusted)
  model_formulas <- list(
    full = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
    )
  )
  cat("Configuration: Using full model (age- and sex-adjusted) only\n")
  
} else if (USE_BOTH_MODELS == 1) {
  # Use both full model and model without age/sex adjustments
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
  cat("Configuration: Using both full model and model without age/sex adjustments\n")
  
} else {
  # Default to full model if invalid switch value
  model_formulas <- list(
    full = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
    )
  )
  cat("Warning: Invalid USE_BOTH_MODELS value. Defaulting to full model only\n")
}

# Print which models will be fitted
cat("Model structures to be fitted: ", paste(names(model_formulas), collapse = ", "), "\n")

# Define parametric distributions to test
distributions_to_test <- c("exponential", "weibull", "weibullph", "llogis", 
                           "lognormal", "gamma", "gompertz", "gengamma", "genf")

# ===============================================================================
# MODEL FITTING
# ===============================================================================

# Initialize storage structure:
# models$[structure]$[outcome]$[distribution] = fitted flexsurvreg object
# models$[structure]$[outcome]_ic = information criteria table
# models$best_fit$[outcome] = best model overall
models <- list()

# Initialize model structure lists based on configuration
for (structure_name in names(model_formulas)) {
  models[[structure_name]] <- list()
}

# Fit models for each selected structure
for (structure_name in names(model_formulas)) {
  
  cat("Fitting models for structure:", structure_name, "\n")
  
  # Fit OS models using fit_all_direct function
  os_results <- fit_all_direct(
    fit_data = data_complete,
    fit_formula = model_formulas[[structure_name]]$os,
    fit_dists = distributions_to_test
  )
  
  # Fit PFS models using fit_all_direct function
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
  }
  
  if (!is.null(best_pfs)) {
    models[[structure_name]]$best_pfs_dist <- best_pfs$distribution
    models[[structure_name]]$best_pfs_aic <- best_pfs$criterion_value
  }
}

# ===============================================================================
# GLOBAL BEST MODEL SELECTION
# ===============================================================================

# Combine AIC values from all fitted structures to find global optimum
all_os_aic <- do.call(rbind, lapply(names(model_formulas), function(structure_name) {
  if (!is.null(models[[structure_name]]$os_ic)) {
    data.frame(structure = structure_name, models[[structure_name]]$os_ic)
  }
}))

all_pfs_aic <- do.call(rbind, lapply(names(model_formulas), function(structure_name) {
  if (!is.null(models[[structure_name]]$pfs_ic)) {
    data.frame(structure = structure_name, models[[structure_name]]$pfs_ic)
  }
}))

# Select globally best OS model
if (!is.null(all_os_aic) && nrow(all_os_aic) > 0) {
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
    
    cat("Best OS model: ", best_os_structure, " structure, ", best_os_dist, " distribution (AIC: ", 
        round(best_os_overall$AIC[best_os_idx], 2), ")\n")
  }
}

# Select globally best PFS model
if (!is.null(all_pfs_aic) && nrow(all_pfs_aic) > 0) {
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
    
    cat("Best PFS model: ", best_pfs_structure, " structure, ", best_pfs_dist, " distribution (AIC: ", 
        round(best_pfs_overall$AIC[best_pfs_idx], 2), ")\n")
  }
}

# ===============================================================================
# POPULATION-LEVEL PREDICTIONS
# ===============================================================================

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
  
} else {
  # If model selection failed, create empty predictions object
  predictions <- list()
  warning("Model fitting failed - predictions object is empty")
}

# ===============================================================================
# FINAL STORAGE STRUCTURE DOCUMENTATION
# ===============================================================================
# 
# CONFIGURATION:
# - USE_BOTH_MODELS: Numeric switch controlling model selection
#   * 0 = Only fit full model with age and sex adjustments (recommended for decision analysis)
#   * 1 = Fit both full model and model without age/sex adjustments
#
# MODELS OBJECT STRUCTURE:
# ├── models$[selected_structures]      # Based on USE_BOTH_MODELS setting
# │   ├── $os                          # OS models by distribution
# │   │   ├── $weibull                # flexsurvreg object
# │   │   ├── $exponential            # flexsurvreg object
# │   │   └── $[other_distributions]
# │   ├── $pfs                         # PFS models by distribution  
# │   ├── $os_ic                       # AIC/BIC table for OS models
# │   ├── $pfs_ic                      # AIC/BIC table for PFS models
# │   ├── $km_os                       # Kaplan-Meier fit for OS
# │   └── $km_pfs                      # Kaplan-Meier fit for PFS
# └── models$best_fit                  # Globally best models from selected structures
#     ├── $os                          # Best OS model object
#     ├── $os_structure               # Structure name
#     ├── $os_distribution            # e.g., "weibull"
#     ├── $os_aic                     # AIC value
#     └── [same for pfs]
#
# PREDICTIONS OBJECT STRUCTURE:
# ├── predictions$control              # Standard of care strategy
# │   ├── $os                         # Population OS curve
# │   └── $pfs                        # Population PFS curve
# ├── predictions$crp                 # CRP-guided strategy
# │   ├── $os                         # Population-weighted OS curve
# │   ├── $pfs                        # Population-weighted PFS curve
# │   ├── $biomarker_positive         # CRP+ subgroup predictions
# │   ├── $biomarker_negative         # CRP- subgroup predictions
# │   └── $prevalence                 # CRP prevalence used for weighting
# └── predictions$[tlr/tmb_braf]      # Same structure for other biomarkers
# ===============================================================================