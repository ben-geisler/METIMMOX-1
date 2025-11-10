#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr, tidyr)

# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# ===============================================================================
# BOOTSTRAP CACHE CONFIGURATION
# ===============================================================================

# Define cache directory and file paths
cache_dir <- here("data", "bootstrap_cache")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
  cat("Created bootstrap cache directory:", cache_dir, "\n")
}

# Create cache file path based on n_samples and model configuration
cache_file <- here(cache_dir, paste0("boot_models_n", n_samples, "_", 
                                     ifelse(USE_BOTH_MODELS == 0, "full", "both"), 
                                     ".rds"))

# ===============================================================================
# MODEL FORMULA SELECTION BASED ON PARAMETRIC SURVIVAL ANALYSIS SWITCH
# ===============================================================================

# Use the same switch value from 06_parametric_survival analysis.R
# USE_BOTH_MODELS: 0 = full model only, 1 = both models
# Note: We assume the value is already set in the environment

# Define model formulas based on the switch (matching 06_parametric_survival analysis.R)
if (USE_BOTH_MODELS == 0) {
  # Only use full model (age- and sex-adjusted) - matching the parametric analysis
  model_type <- "full"
  cat("Bootstrap sampling: Using full model (age- and sex-adjusted) to match parametric analysis\n")
  
  # Control model formula (age and sex adjusted)
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex
  
  # Function to create biomarker model formulas (age and sex adjusted)
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    }
  }
  
} else if (USE_BOTH_MODELS == 1) {
  # For compatibility, we'll use the full model structure when both are available
  # This ensures consistency with whichever model was selected as best in parametric analysis
  model_type <- "full"  # Default to full for bootstrap sampling
  cat("Bootstrap sampling: Using full model structure (age- and sex-adjusted) for consistency\n")
  
  # Control model formula (age and sex adjusted)
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex
  
  # Function to create biomarker model formulas (age and sex adjusted)
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    }
  }
  
} else {
  # Default to full model if invalid switch value
  model_type <- "full"
  cat("Warning: Invalid USE_BOTH_MODELS value in bootstrap sampling. Defaulting to full model\n")
  
  # Control model formula (age and sex adjusted)
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex
  
  # Function to create biomarker model formulas (age and sex adjusted)
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    }
  }
}

# ===============================================================================
# CORRELATED BOOTSTRAP FUNCTION
# ===============================================================================
# This function implements non-parametric bootstrap that maintains correlation
# between PFS and OS by resampling patients and fitting both models to the same
# resampled dataset.
# ===============================================================================

bootstrap_correlated_survival <- function(formula_os, formula_pfs, data, 
                                          dist = "weibull", n_boot = n_samples) {
  
  # Get the best-fitting distribution from parametric analysis if available
  if (exists("models") && !is.null(models$best_fit$os_distribution)) {
    dist <- models$best_fit$os_distribution
    cat("Using distribution from parametric analysis:", dist, "\n")
  }
  
  n_patients <- nrow(data)
  boot_samples <- vector("list", n_boot)
  n_failed <- 0
  
  cat("Generating", n_boot, "correlated bootstrap samples for", n_patients, "patients...\n")
  cat("This preserves PFS-OS correlation by fitting both models to the same resampled data\n")
  
  # Fit original models for reference
  original_os <- flexsurvreg(formula_os, data = data, dist = dist)
  original_pfs <- flexsurvreg(formula_pfs, data = data, dist = dist)
  
  # Set up progress reporting
  start_time <- Sys.time()
  
  for (i in 1:n_boot) {
    # Resample patients with replacement (non-parametric bootstrap)
    boot_idx <- sample(1:n_patients, size = n_patients, replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    # Fit BOTH models to the SAME resampled dataset
    # This is critical: PFS and OS share the same patient cohort
    tryCatch({
      os_model <- flexsurvreg(formula_os, data = boot_data, dist = dist)
      pfs_model <- flexsurvreg(formula_pfs, data = boot_data, dist = dist)
      
      boot_samples[[i]] <- list(
        os = list(
          model = os_model,
          coefficients = os_model$coefficients,
          dist = dist
        ),
        pfs = list(
          model = pfs_model,
          coefficients = pfs_model$coefficients,
          dist = dist
        ),
        boot_idx = boot_idx  # Store for reproducibility/debugging
      )
    }, error = function(e) {
      # If bootstrap sample fails, use original model
      cat("Warning: Bootstrap sample", i, "failed:", conditionMessage(e), "\n")
      boot_samples[[i]] <- list(
        os = list(
          model = original_os,
          coefficients = original_os$coefficients,
          dist = dist
        ),
        pfs = list(
          model = original_pfs,
          coefficients = original_pfs$coefficients,
          dist = dist
        ),
        failed = TRUE
      )
      n_failed <<- n_failed + 1
    })
    
    # Progress reporting
    if (i %% 500 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- i / elapsed
      remaining <- (n_boot - i) / rate
      cat(sprintf("Completed %d/%d bootstrap samples (%.1f%%) - %d failures - ETA: %.1f min\n", 
                  i, n_boot, 100*i/n_boot, n_failed, remaining))
    }
  }
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("Bootstrap complete: %d successful, %d failed in %.1f minutes\n", 
              n_boot - n_failed, n_failed, total_time))
  
  return(list(
    samples = boot_samples,
    original_os = original_os,
    original_pfs = original_pfs,
    n_boot = n_boot,
    n_failed = n_failed,
    dist = dist,
    creation_time = Sys.time()
  ))
}

# ===============================================================================
# CHECK CACHE OR GENERATE BOOTSTRAP SAMPLES
# ===============================================================================

if (file.exists(cache_file)) {
  cat("\n=== Loading bootstrap samples from cache ===\n")
  cat("Cache file:", cache_file, "\n")
  
  boot_models <- readRDS(cache_file)
  
  # Validate cache
  if (is.list(boot_models) && 
      "control" %in% names(boot_models) &&
      all(biomarkers %in% names(boot_models))) {
    
    cat("Cache loaded successfully!\n")
    cat("- Control samples:", boot_models$control$n_boot, "\n")
    cat("- Biomarkers:", paste(biomarkers, collapse = ", "), "\n")
    cat("- Created:", format(boot_models$control$creation_time), "\n")
    cat("- Distribution:", boot_models$control$dist, "\n")
    
    # Check if n_samples matches
    if (boot_models$control$n_boot != n_samples) {
      cat("WARNING: Cached samples (", boot_models$control$n_boot, 
          ") != requested (", n_samples, ")\n")
      cat("Will use cached samples. Delete cache file to regenerate.\n")
    }
    
  } else {
    cat("Cache file corrupted or incomplete. Regenerating...\n")
    boot_models <- NULL
  }
  
} else {
  cat("\n=== No cache found - will generate bootstrap samples ===\n")
  boot_models <- NULL
}

# ===============================================================================
# BOOTSTRAP MODEL FITTING WITH CORRELATION (if not cached)
# ===============================================================================

if (is.null(boot_models)) {
  
  # Initialize list to store bootstrapped models
  boot_models <- list()
  
  cat("\n=== Starting correlated bootstrap sampling ===\n")
  cat("This will take some time but only needs to run once.\n")
  cat("Results will be cached to:", cache_file, "\n\n")
  
  # Bootstrap control models with age and sex adjustments
  cat("\nBootstrapping control strategy...\n")
  boot_models$control <- bootstrap_correlated_survival(
    formula_os = control_os_formula,
    formula_pfs = control_pfs_formula,
    data = data_control,
    n_boot = n_samples
  )
  
  # Bootstrap biomarker models with age and sex adjustments
  for (biomarker in biomarkers) {
    cat("\nBootstrapping", biomarker, "strategy...\n")
    
    # OS model with age, sex, biomarker, treatment, and interaction
    os_formula <- create_biomarker_formula("os", biomarker)
    pfs_formula <- create_biomarker_formula("pfs", biomarker)
    
    boot_models[[biomarker]] <- bootstrap_correlated_survival(
      formula_os = os_formula,
      formula_pfs = pfs_formula,
      data = data,
      n_boot = n_samples
    )
  }
  
  # Save to cache
  cat("\n=== Saving bootstrap samples to cache ===\n")
  tryCatch({
    saveRDS(boot_models, cache_file)
    cat("Bootstrap samples saved successfully to:", cache_file, "\n")
    cat("File size:", format(file.info(cache_file)$size / 1024^2, digits = 2), "MB\n")
  }, error = function(e) {
    cat("WARNING: Failed to save cache:", conditionMessage(e), "\n")
  })
}

# Print confirmation of model structures
cat("\n=== Bootstrap model structures ready ===\n")
cat("- Control models: Age + sex adjusted\n")
cat("- Biomarker models: Age + sex + biomarker + treatment + interaction\n")
cat("- Model type used:", model_type, "\n")
cat("- Number of bootstrap samples:", boot_models$control$n_boot, "\n")
cat("- PFS and OS are CORRELATED within each bootstrap sample\n")

# ===============================================================================
# PREDICTION FUNCTION FOR BOOTSTRAP SAMPLES - CORRECTED VERSION
# ===============================================================================

generate_bootstrap_predictions <- function(boot_model_list, outcome = "os", 
                                           sample_idx = 1, newdata = NULL, 
                                           time_points) {
  # Make sure time_points passed in has correct length
  if (length(time_points) != time_points_length) {
    warning("time_points length mismatch in bootstrap predictions")
    time_points <- time_points[1:min(length(time_points), time_points_length)]
  }
  
  # Extract the specific bootstrap sample
  boot_sample <- boot_model_list$samples[[sample_idx]]
  
  if (is.null(boot_sample)) {
    warning("Bootstrap sample ", sample_idx, " is NULL, using original model")
    model_obj <- if (outcome == "os") boot_model_list$original_os else boot_model_list$original_pfs
  } else {
    # Get the model for the requested outcome (os or pfs)
    model_obj <- boot_sample[[outcome]]$model
  }
  
  # Generate predictions using the bootstrap-fitted model
  tryCatch({
    if (is.null(newdata)) {
      # No newdata - predict at mean covariate values
      pred <- predict(model_obj, type = "survival", times = time_points)
    } else {
      # With newdata - predict for specific covariate combination
      pred <- predict(model_obj, newdata = newdata, type = "survival", times = time_points)
    }
    
    # ========================================================================
    # CRITICAL FIX: Extract survival predictions from nested tibble structure
    # ========================================================================
    # predict.flexsurvreg returns a tibble with a list column called .pred
    # Each element of .pred is itself a tibble with columns:
    #   .eval_time (the time points) and .pred_survival (the survival probabilities)
    
    if (is.data.frame(pred) && ".pred" %in% names(pred)) {
      # This is the nested tibble structure from flexsurv
      
      if (is.null(newdata)) {
        # Without newdata: multiple rows (one per patient in original data)
        # We want to take the MEAN across all patients
        # Extract all survival predictions and average them
        all_surv_probs <- lapply(pred$.pred, function(x) x$.pred_survival)
        
        # Convert list to matrix (each column is a patient, each row is a time point)
        surv_matrix <- do.call(cbind, all_surv_probs)
        
        # Take row means (average across patients at each time point)
        mean_survival <- rowMeans(surv_matrix, na.rm = TRUE)
        
        return(mean_survival)
        
      } else {
        # With newdata: typically just one row
        # Extract the survival predictions for that row
        if (nrow(pred) == 1) {
          # Single prediction
          survival_pred <- pred$.pred[[1]]$.pred_survival
          return(survival_pred)
        } else {
          # Multiple rows in newdata - take the first one
          warning("Multiple rows in newdata, using first row")
          survival_pred <- pred$.pred[[1]]$.pred_survival
          return(survival_pred)
        }
      }
    } else {
      warning("Unexpected prediction structure from flexsurvreg")
      return(rep(NA, length(time_points)))
    }
    
  }, error = function(e) {
    warning("Prediction failed for bootstrap sample ", sample_idx, ": ", conditionMessage(e))
    # Return NA vector of correct length
    return(rep(NA, length(time_points)))
  })
}

# ===============================================================================
# PARAMETER DISTRIBUTIONS FOR UNCERTAINTY ANALYSIS (NON-SURVIVAL PARAMETERS)
# ===============================================================================

# Create parameter distributions for other model parameters
# Cost parameters (assume 20% coefficient of variation for costs)
cv_costs <- 0.2

# Function to calculate gamma parameters from mean and CV
get_gamma_params <- function(mean_val, cv) {
  shape <- 1/cv^2
  rate <- shape/mean_val
  return(list(shape = shape, rate = rate))
}

# Drug costs - gamma distributions
dist_c_drug_nivo <- c(list(dist = "gamma"), 
                      get_gamma_params(l_params_base$c_drug_nivo, cv_costs))

dist_c_drug_FLOX <- c(list(dist = "gamma"), 
                      get_gamma_params(l_params_base$c_drug_FLOX, cv_costs))

# Test costs - gamma distributions
dist_c_test_CT <- c(list(dist = "gamma"), 
                    get_gamma_params(l_params_base$c_test_CT, cv_costs))

dist_c_test_blood <- c(list(dist = "gamma"), 
                       get_gamma_params(l_params_base$c_test_blood, cv_costs))

dist_c_test_NGS <- c(list(dist = "gamma"), 
                     get_gamma_params(l_params_base$c_test_NGS, cv_costs))

# Other costs - gamma distributions
dist_c_other_visit <- c(list(dist = "gamma"), 
                        get_gamma_params(l_params_base$c_other_visit, cv_costs))

dist_c_other_baseline <- c(list(dist = "gamma"), 
                           get_gamma_params(l_params_base$c_other_baseline, cv_costs))

dist_c_other_follow <- c(list(dist = "gamma"), 
                         get_gamma_params(l_params_base$c_other_follow, cv_costs))

dist_c_other_last <- c(list(dist = "gamma"), 
                       get_gamma_params(l_params_base$c_other_last, cv_costs))


# Utility parameters - beta distributions (bounded between 0 and 1)
# Using method of moments to get alpha and beta parameters
get_beta_params <- function(mean_val, cv) {
  var_val <- (mean_val * cv)^2
  alpha <- mean_val * (mean_val * (1 - mean_val) / var_val - 1)
  beta <- (1 - mean_val) * (mean_val * (1 - mean_val) / var_val - 1)
  return(list(shape1 = alpha, shape2 = beta))
}

cv_utilities <- 0.15  # 15% CV for utilities
dist_u_np <- c(list(dist = "beta"), get_beta_params(l_params_base$u_np, cv_utilities))
dist_u_p <- c(list(dist = "beta"), get_beta_params(l_params_base$u_p, cv_utilities))

# Create comprehensive parameter distributions list
# Note: Survival curves are NOT in this list - they come from bootstrap samples
param_distributions <- list(
  # Cost parameters
  c_drug_nivo = dist_c_drug_nivo,
  c_drug_FLOX = dist_c_drug_FLOX,
  c_test_CT = dist_c_test_CT,
  c_test_blood = dist_c_test_blood,
  c_test_NGS = dist_c_test_NGS,
  c_other_visit = dist_c_other_visit,
  c_other_baseline = dist_c_other_baseline,
  c_other_follow = dist_c_other_follow,
  c_other_last = dist_c_other_last,
  
  # Utility parameters
  u_np = dist_u_np,
  u_p = dist_u_p
)

cat("\n=== Parameter distributions configured ===\n")
cat("- Cost parameters: Gamma distributions (CV =", cv_costs, ")\n")
cat("- Utility parameters: Beta distributions (CV =", cv_utilities, ")\n")
cat("- Survival parameters: Correlated bootstrap samples (n =", boot_models$control$n_boot, ")\n")
cat("\nReady for PSA analysis\n")