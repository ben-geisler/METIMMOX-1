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
#we will set this now with the global variables, so it can be overridden later more easily

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
# Note: sex is converted to factor in script 03 before subsets are created
# This prevents type mismatch in resampled models (issue #72)
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
# MODEL_STRUCTURE controls which model formula structure is used:
#   0 = "joint"    - All biomarkers + all treatment interactions in ONE model (Model A)
#   1 = "focused"  - All biomarkers as main effects + ONE interaction per model (Model B)
#   2 = "separate" - Only ONE biomarker + its interaction per model (Model C)
#
# USE_BOTH_MODELS controls age/sex adjustment:
#   0 = Full model with age and sex adjustments only
#   1 = Both full model and model without age/sex adjustments
# ===============================================================================

# Validate MODEL_STRUCTURE
if (!exists("MODEL_STRUCTURE")) {
  MODEL_STRUCTURE <- 0
  cat("MODEL_STRUCTURE not defined, defaulting to 0 (joint)\n")
}

if (MODEL_STRUCTURE == 0) {
  # =========================================================================
  # JOINT (Model A): All biomarkers + all treatment interactions in ONE model
  # =========================================================================
  model_type_label <- "joint"
  cat("Configuration: Using JOINT model (Model A - all biomarkers + all interactions)\n")

  # Define model formulas based on USE_BOTH_MODELS
  if (USE_BOTH_MODELS == 0) {
    model_formulas <- list(
      full = list(
        os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
        pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
      )
    )
  } else if (USE_BOTH_MODELS == 1) {
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
  } else {
    model_formulas <- list(
      full = list(
        os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx,
        pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx
      )
    )
    cat("Warning: Invalid USE_BOTH_MODELS value. Defaulting to full model only\n")
  }

} else if (MODEL_STRUCTURE == 1) {
  # =========================================================================
  # FOCUSED (Model B): All biomarkers as main effects + ONE interaction per model
  # =========================================================================
  # Each biomarker strategy uses a separate model with:
  # - All biomarkers as main effects (adjusting for their prognostic value)
  # - Only the relevant biomarker's treatment interaction
  # This allows each biomarker's predictive effect to be estimated independently
  # while still adjusting for other biomarkers' prognostic effects.
  # =========================================================================
  model_type_label <- "focused"
  cat("Configuration: Using FOCUSED models (Model B - all biomarkers + one interaction each)\n")
  cat("Fitting 3 separate models (one per biomarker strategy).\n")

  model_formulas <- list(
    crp = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tlr + tmb_braf + crp:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tlr + tmb_braf + crp:Rx
    ),
    tlr = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tlr + tmb_braf + tlr:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tlr + tmb_braf + tlr:Rx
    ),
    tmb_braf = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tlr + tmb_braf + tmb_braf:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tlr + tmb_braf + tmb_braf:Rx
    )
  )

} else if (MODEL_STRUCTURE == 2) {
  # =========================================================================
  # SEPARATE (Model C): Only ONE biomarker + its interaction per model
  # =========================================================================
  # Each biomarker strategy uses a separate model with:
  # - Only the relevant biomarker's treatment interaction
  # - No adjustment for other biomarkers
  # This is the most parsimonious approach, estimating each biomarker's
  # predictive effect in isolation.
  # =========================================================================
  model_type_label <- "separate"
  cat("Configuration: Using SEPARATE models (Model C - one biomarker + its interaction only)\n")
  cat("Fitting 3 separate models (one per biomarker strategy).\n")

  model_formulas <- list(
    crp = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx
    ),
    tlr = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + tlr:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tlr:Rx
    ),
    tmb_braf = list(
      os = Surv(OSwk, Death) ~ Age + sex + Rx + tmb_braf:Rx,
      pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tmb_braf:Rx
    )
  )

} else {
  stop("Invalid MODEL_STRUCTURE value: ", MODEL_STRUCTURE,
       ". Must be 0 (joint), 1 (focused), or 2 (separate)")
}

cat("Model type label:", model_type_label, "\n")

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
# For Model A (joint): Select globally best OS and PFS models across structures
# For Models B/C (focused/separate): Select best model for EACH biomarker
# ===============================================================================

# Check if this is a biomarker-specific model structure (Models B/C)
is_biomarker_specific <- MODEL_STRUCTURE %in% c(1, 2)
biomarker_names <- c("crp", "tlr", "tmb_braf")

if (is_biomarker_specific) {
  # =========================================================================
  # BIOMARKER-SPECIFIC BEST MODEL SELECTION (Models B/C)
  # =========================================================================
  # Each biomarker strategy has its own model; select best distribution for each

  cat("Selecting best models for each biomarker strategy...\n")

  for (biomarker in biomarker_names) {
    if (!is.null(models[[biomarker]])) {
      # Select best OS model for this biomarker
      if (!is.null(models[[biomarker]]$os_ic)) {
        os_ic <- models[[biomarker]]$os_ic
        os_ic_valid <- os_ic[!is.na(os_ic$AIC), ]
        if (nrow(os_ic_valid) > 0) {
          best_os_idx <- which.min(os_ic_valid$AIC)
          best_os_dist <- os_ic_valid$Distribution[best_os_idx]

          models$best_fit[[biomarker]]$os <- models[[biomarker]]$os[[best_os_dist]]
          models$best_fit[[biomarker]]$os_distribution <- best_os_dist
          models$best_fit[[biomarker]]$os_aic <- os_ic_valid$AIC[best_os_idx]

          cat("  Best OS model for ", biomarker, ": ", best_os_dist,
              " (AIC: ", round(os_ic_valid$AIC[best_os_idx], 2), ")\n", sep = "")
        }
      }

      # Select best PFS model for this biomarker
      if (!is.null(models[[biomarker]]$pfs_ic)) {
        pfs_ic <- models[[biomarker]]$pfs_ic
        pfs_ic_valid <- pfs_ic[!is.na(pfs_ic$AIC), ]
        if (nrow(pfs_ic_valid) > 0) {
          best_pfs_idx <- which.min(pfs_ic_valid$AIC)
          best_pfs_dist <- pfs_ic_valid$Distribution[best_pfs_idx]

          models$best_fit[[biomarker]]$pfs <- models[[biomarker]]$pfs[[best_pfs_dist]]
          models$best_fit[[biomarker]]$pfs_distribution <- best_pfs_dist
          models$best_fit[[biomarker]]$pfs_aic <- pfs_ic_valid$AIC[best_pfs_idx]

          cat("  Best PFS model for ", biomarker, ": ", best_pfs_dist,
              " (AIC: ", round(pfs_ic_valid$AIC[best_pfs_idx], 2), ")\n", sep = "")
        }
      }
    }
  }

} else {
  # =========================================================================
  # SINGLE MODEL BEST SELECTION (Model A - Joint)
  # =========================================================================
  # Select globally best OS and PFS models across all structures

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
}

# ===============================================================================
# POPULATION-LEVEL PREDICTIONS
# ===============================================================================
# Uses population averaging approach to account for heterogeneity in patient
# characteristics (age, sex, biomarkers). Predictions are generated for all
# patients using their actual covariate values, then averaged within relevant
# subgroups. This is the methodologically rigorous approach for CEA.
# See GitHub Issue #69 for methodology discussion.
#
# For Model A (joint): Use single best model for all strategies
# For Models B/C (focused/separate): Use biomarker-specific models
# ===============================================================================

if (is_biomarker_specific) {
  # =========================================================================
  # BIOMARKER-SPECIFIC PREDICTIONS (Models B/C)
  # =========================================================================
  # Each biomarker strategy uses its own fitted model

  # Check that all biomarker models were successfully fitted
  all_biomarkers_fitted <- all(sapply(biomarker_names, function(bm) {
    !is.null(models$best_fit[[bm]]$os) && !is.null(models$best_fit[[bm]]$pfs)
  }))

  if (all_biomarkers_fitted) {
    # Generate predictions using biomarker-specific models
    predictions <- generate_population_averaged_predictions_multimodel(
      biomarker_models = models$best_fit,
      strategies_df = strategies_df,
      data_complete = data_complete,
      time_points = time_points
    )
  } else {
    predictions <- list()
    warning("Some biomarker models failed to fit - predictions object is empty")
  }

} else {
  # =========================================================================
  # SINGLE MODEL PREDICTIONS (Model A - Joint)
  # =========================================================================
  # Use same model for all strategies

  if (!is.null(models$best_fit$os) && !is.null(models$best_fit$pfs)) {
    # Create model list for prediction functions
    best_models <- list(
      os = models$best_fit$os,
      pfs = models$best_fit$pfs
    )

    # Generate population-averaged predictions for all strategies
    # This function:
    # 1. Predicts survival for EACH patient using their actual age, sex, biomarkers
    # 2. For each strategy, assigns treatment based on biomarker status
    # 3. Averages predictions within relevant subgroups (biomarker+, biomarker-)
    # 4. Computes population-weighted averages using observed prevalence
    predictions <- generate_population_averaged_predictions(
      models = best_models,
      strategies_df = strategies_df,
      data_complete = data_complete,
      time_points = time_points
    )
  } else {
    # If model selection failed, create empty predictions object
    predictions <- list()
    warning("Model fitting failed - predictions object is empty")
  }
}

# ===============================================================================
# FINAL STORAGE STRUCTURE DOCUMENTATION
# ===============================================================================
#
# CONFIGURATION:
# - MODEL_STRUCTURE: Controls which model formula structure is used
#   * 0 = "joint" (Model A): All biomarkers + all interactions in one model
#   * 1 = "focused" (Model B): All biomarkers as main effects + one interaction each
#   * 2 = "separate" (Model C): Only one biomarker + its interaction per model
#
# - USE_BOTH_MODELS: Controls age/sex adjustment (only for Model A)
#   * 0 = Only fit full model with age and sex adjustments
#   * 1 = Fit both full model and model without age/sex adjustments
#
# MODELS OBJECT STRUCTURE:
#
# For MODEL_STRUCTURE = 0 (Model A - Joint):
# ├── models$full                       # Single model with all biomarker interactions
# │   ├── $os                          # OS models by distribution
# │   │   ├── $weibull                # flexsurvreg object
# │   │   ├── $exponential            # flexsurvreg object
# │   │   └── $[other_distributions]
# │   ├── $pfs                         # PFS models by distribution
# │   ├── $os_ic                       # AIC/BIC table for OS models
# │   ├── $pfs_ic                      # AIC/BIC table for PFS models
# │   ├── $km_os                       # Kaplan-Meier fit for OS
# │   └── $km_pfs                      # Kaplan-Meier fit for PFS
# └── models$best_fit                  # Globally best models
#     ├── $os                          # Best OS model object
#     ├── $os_structure               # Structure name ("full")
#     ├── $os_distribution            # e.g., "gamma"
#     ├── $os_aic                     # AIC value
#     └── [same for pfs]
#
# For MODEL_STRUCTURE = 1 or 2 (Models B/C - Biomarker-specific):
# ├── models$crp                        # CRP-specific model
# │   ├── $os, $pfs, $os_ic, $pfs_ic, $km_os, $km_pfs  # Same structure as above
# ├── models$tlr                        # TLR-specific model
# │   ├── ... (same structure)
# ├── models$tmb_braf                   # TMB/BRAF-specific model
# │   ├── ... (same structure)
# └── models$best_fit                   # Best models PER BIOMARKER
#     ├── $crp$os, $crp$pfs            # Best CRP models
#     ├── $crp$os_distribution, $crp$os_aic
#     ├── $tlr$os, $tlr$pfs            # Best TLR models
#     ├── ... (same structure)
#     └── $tmb_braf$os, $tmb_braf$pfs  # Best TMB/BRAF models
#
# PREDICTIONS OBJECT STRUCTURE (same for all MODEL_STRUCTURE values):
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