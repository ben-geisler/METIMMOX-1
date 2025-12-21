#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr, tidyr)

# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# ===============================================================================
# SAMPLING CACHE CONFIGURATION
# ===============================================================================

# Define cache directory and file paths
cache_dir <- here("data", "tidy")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
  cat("Created sampling cache directory:", cache_dir, "\n")
}

# Create cache file path based on n_samples and model configuration
# Include MODEL_STRUCTURE in filename to ensure separate caches per structure
model_structure_label <- c("joint", "focused", "separate")[MODEL_STRUCTURE + 1]
cache_file <- here(cache_dir, paste0("sampling_models_n", n_samples, "_",
                                     ifelse(USE_BOTH_MODELS == 0, "full", "both"), "_",
                                     model_structure_label,
                                     ".rds"))

# Check for old cache files in deprecated location
old_cache_dir <- here("data", "bootstrap_cache")
old_cache_file <- here(old_cache_dir, paste0("boot_models_n", n_samples, "_",
                                              ifelse(USE_BOTH_MODELS == 0, "full", "both"),
                                              ".rds"))
if (file.exists(old_cache_file) && !file.exists(cache_file)) {
  cat("\n*** NOTICE: Old cache file detected ***\n")
  cat("Old location:", old_cache_file, "\n")
  cat("New location:", cache_file, "\n")
  cat("Please regenerate sampling models. Old cache will not be used.\n\n")
}

# ===============================================================================
# MODEL FORMULA SELECTION BASED ON GLOBAL MODEL_STRUCTURE SWITCH
# ===============================================================================
# Matches clinical_effectiveness.qmd Model A/B/C structure
# This ensures PSA uses the SAME model structure as base case
#
# MODEL_STRUCTURE options:
#   0 = "joint"    - All biomarkers + all treatment interactions in ONE model (Model A)
#   1 = "focused"  - All biomarkers as main effects + ONE interaction per model (Model B)
#   2 = "separate" - Only ONE biomarker + its interaction per model (Model C)
# ===============================================================================

if (MODEL_STRUCTURE == 0) {
  # =========================================================================
  # JOINT (Model A): All biomarkers + all treatment interactions
  # =========================================================================
  model_type <- "joint"
  cat("Resampling: Using JOINT model (Model A - all biomarkers + all interactions)\n")

  # Control model - full model with all terms, fitted to ALL patients
  # Note: For control strategy, we still use simpler model since control patients
  # don't have treatment variation
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex

  # Biomarker formulas - JOINT model with ALL biomarker interactions
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx)
    } else {
      return(Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx)
    }
  }

} else if (MODEL_STRUCTURE == 1) {
  # =========================================================================
  # FOCUSED (Model B): All biomarkers as main effects + ONE interaction per model
  # =========================================================================
  model_type <- "focused"
  cat("Resampling: Using FOCUSED models (Model B - all biomarkers + one interaction)\n")

  # Control model - all biomarkers as main effects, no interactions
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex + crp + tlr + tmb_braf
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex + crp + tlr + tmb_braf

  # Biomarker formulas - all biomarkers + only THIS biomarker's interaction
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + crp + tlr + tmb_braf + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tlr + tmb_braf + ", biomarker, ":Rx")))
    }
  }

} else if (MODEL_STRUCTURE == 2) {
  # =========================================================================
  # SEPARATE (Model C): Only ONE biomarker + its interaction per model
  # =========================================================================
  model_type <- "separate"
  cat("Resampling: Using SEPARATE models (Model C - one biomarker only)\n")

  # Control model - simple age + sex (fitted to control patients only)
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex

  # Biomarker formulas - only THIS biomarker + its interaction
  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    }
  }

} else {
  stop("Invalid MODEL_STRUCTURE value: ", MODEL_STRUCTURE,
       ". Must be 0 (joint), 1 (focused), or 2 (separate)")
}

# ===============================================================================
# CORRELATED SAMPLING FUNCTION
# ===============================================================================
# This function implements non-parametric resampling that maintains correlation
# between PFS and OS by resampling patients and fitting both models to the same
# resampled dataset.
# ===============================================================================

sample_correlated_survival <- function(formula_os, formula_pfs, data,
                                       dist = "weibull", n_samples = n_samples) {
  
  # Get the best-fitting distribution from parametric analysis if available
  if (exists("models") && !is.null(models$best_fit$os_distribution)) {
    dist <- models$best_fit$os_distribution
    cat("Using distribution from parametric analysis:", dist, "\n")
  }
  
  n_patients <- nrow(data)
  sampled_models <- vector("list", n_samples)
  n_failed <- 0

  cat("Generating", n_samples, "correlated resampled models for", n_patients, "patients...\n")
  cat("This preserves PFS-OS correlation by fitting both models to the same resampled data\n")
  
  # Fit original models for reference
  original_os <- flexsurvreg(formula_os, data = data, dist = dist)
  original_pfs <- flexsurvreg(formula_pfs, data = data, dist = dist)
  
  # Set up progress reporting
  start_time <- Sys.time()
  
  for (i in 1:n_samples) {
    # Resample patients with replacement (non-parametric resampling)
    resample_idx <- sample(1:n_patients, size = n_patients, replace = TRUE)
    resampled_data <- data[resample_idx, ]

    # Fit BOTH models to the SAME resampled dataset
    # This is critical: PFS and OS share the same patient cohort
    tryCatch({
      os_model <- flexsurvreg(formula_os, data = resampled_data, dist = dist)
      pfs_model <- flexsurvreg(formula_pfs, data = resampled_data, dist = dist)

      sampled_models[[i]] <- list(
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
        resample_idx = resample_idx  # Store for reproducibility/debugging
      )
    }, error = function(e) {
      # If resampled model fails, use original model
      cat("Warning: Resampled model", i, "failed:", conditionMessage(e), "\n")
      sampled_models[[i]] <- list(
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
      remaining <- (n_samples - i) / rate
      cat(sprintf("Completed %d/%d resampled models (%.1f%%) - %d failures - ETA: %.1f min\n",
                  i, n_samples, 100*i/n_samples, n_failed, remaining))
    }
  }

  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("Resampling complete: %d successful, %d failed in %.1f minutes\n",
              n_samples - n_failed, n_failed, total_time))

  return(list(
    samples = sampled_models,
    original_os = original_os,
    original_pfs = original_pfs,
    n_samples = n_samples,
    n_failed = n_failed,
    dist = dist,
    creation_time = Sys.time()
  ))
}

# ===============================================================================
# CHECK CACHE OR GENERATE SAMPLING MODELS
# ===============================================================================

if (file.exists(cache_file)) {
  cat("\n=== Loading sampling models from cache ===\n")
  cat("Cache file:", cache_file, "\n")

  sampling_models <- readRDS(cache_file)

  # Validate cache
  if (is.list(sampling_models) &&
      "control" %in% names(sampling_models) &&
      all(biomarkers %in% names(sampling_models))) {

    cat("Cache loaded successfully!\n")
    cat("- Control samples:", sampling_models$control$n_samples, "\n")
    cat("- Biomarkers:", paste(biomarkers, collapse = ", "), "\n")
    cat("- Created:", format(sampling_models$control$creation_time), "\n")
    cat("- Distribution:", sampling_models$control$dist, "\n")

    # Check if n_samples matches
    if (sampling_models$control$n_samples != n_samples) {
      cat("WARNING: Cached models (", sampling_models$control$n_samples,
          ") != requested (", n_samples, ")\n")
      cat("Will use cached models. Delete cache file to regenerate.\n")
    }

  } else {
    cat("Cache file corrupted or incomplete. Regenerating...\n")
    sampling_models <- NULL
  }

} else {
  cat("\n=== No cache found - will generate sampling models ===\n")
  sampling_models <- NULL
}

# ===============================================================================
# SURVIVAL MODEL RESAMPLING WITH CORRELATION (if not cached)
# ===============================================================================

if (is.null(sampling_models)) {

  # Initialize list to store resampled models
  sampling_models <- list()

  cat("\n=== Starting correlated survival resampling ===\n")
  cat("This will take some time but only needs to run once.\n")
  cat("Results will be cached to:", cache_file, "\n\n")

  # Sample control models with age and sex adjustments
  cat("\nSampling control strategy...\n")
  sampling_models$control <- sample_correlated_survival(
    formula_os = control_os_formula,
    formula_pfs = control_pfs_formula,
    data = data_control,
    n_samples = n_samples
  )

  # Sample biomarker models with age and sex adjustments
  for (biomarker in biomarkers) {
    cat("\nSampling", biomarker, "strategy...\n")

    # OS model with age, sex, biomarker, treatment, and interaction
    os_formula <- create_biomarker_formula("os", biomarker)
    pfs_formula <- create_biomarker_formula("pfs", biomarker)

    sampling_models[[biomarker]] <- sample_correlated_survival(
      formula_os = os_formula,
      formula_pfs = pfs_formula,
      data = data,
      n_samples = n_samples
    )
  }

  # Save to cache
  cat("\n=== Saving sampling models to cache ===\n")
  tryCatch({
    saveRDS(sampling_models, cache_file)
    cat("Sampling models saved successfully to:", cache_file, "\n")
    cat("File size:", format(file.info(cache_file)$size / 1024^2, digits = 2), "MB\n")
  }, error = function(e) {
    cat("WARNING: Failed to save cache:", conditionMessage(e), "\n")
  })
}

# Print confirmation of model structures
cat("\n=== Sampling model structures ready ===\n")
cat("- MODEL_STRUCTURE:", MODEL_STRUCTURE, "(", model_type, ")\n")
if (MODEL_STRUCTURE == 0) {
  cat("- Control models: Age + sex adjusted\n")
  cat("- Biomarker models: Age + sex + Rx + crp:Rx + tlr:Rx + tmb_braf:Rx (JOINT)\n")
} else if (MODEL_STRUCTURE == 1) {
  cat("- Control models: Age + sex + crp + tlr + tmb_braf\n")
  cat("- Biomarker models: Age + sex + Rx + all biomarkers + [biomarker]:Rx (FOCUSED)\n")
} else if (MODEL_STRUCTURE == 2) {
  cat("- Control models: Age + sex adjusted\n")
  cat("- Biomarker models: Age + sex + Rx + [biomarker]:Rx (SEPARATE)\n")
}
cat("- Number of resampled models:", sampling_models$control$n_samples, "\n")
cat("- PFS and OS are CORRELATED within each resampled model\n")

# ===============================================================================
# PREDICTION FUNCTION FOR SAMPLED MODELS - CORRECTED VERSION
# ===============================================================================

generate_sampling_predictions <- function(sampling_model_list, outcome = "os",
                                         sample_idx = 1, newdata = NULL,
                                         time_points) {
  # Make sure time_points passed in has correct length
  if (length(time_points) != time_points_length) {
    warning("time_points length mismatch in sampling predictions")
    time_points <- time_points[1:min(length(time_points), time_points_length)]
  }

  # Extract the specific sampled model
  sampled_model <- sampling_model_list$samples[[sample_idx]]

  if (is.null(sampled_model)) {
    warning("Sampled model ", sample_idx, " is NULL, using original model")
    model_obj <- if (outcome == "os") sampling_model_list$original_os else sampling_model_list$original_pfs
  } else {
    # Get the model for the requested outcome (os or pfs)
    model_obj <- sampled_model[[outcome]]$model
  }

  # Generate predictions using the resampled model
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
    warning("Prediction failed for resampled model ", sample_idx, ": ", conditionMessage(e))
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

# Prevalence parameters - beta distributions (bounded between 0 and 1)
cv_prevalence <- 0.15  # 15% CV for prevalence (same as utilities)
dist_p_crp <- c(list(dist = "beta"), get_beta_params(l_params_base$p_crp, cv_prevalence))
dist_p_tlr <- c(list(dist = "beta"), get_beta_params(l_params_base$p_tlr, cv_prevalence))
dist_p_tmb_braf <- c(list(dist = "beta"), get_beta_params(l_params_base$p_tmb_braf, cv_prevalence))

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
  u_p = dist_u_p,

  # Prevalence parameters
  p_crp = dist_p_crp,
  p_tlr = dist_p_tlr,
  p_tmb_braf = dist_p_tmb_braf
)

# Parameter groups for sensitivity analysis (EVPPI)
param_groups <- list(
  drug_costs = c("c_drug_nivo", "c_drug_FLOX"),
  test_costs = c("c_test_CT", "c_test_blood", "c_test_NGS"),
  other_costs = c("c_other_visit", "c_other_baseline",
                  "c_other_follow", "c_other_last"),
  all_costs = c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                "c_test_NGS", "c_other_visit", "c_other_baseline",
                "c_other_follow", "c_other_last"),
  utilities = c("u_np", "u_p"),
  prevalence = c("p_crp", "p_tlr", "p_tmb_braf")
)

# ===============================================================================
# POPULATION AVERAGING FUNCTION FOR PSA (BIOMARKER STRATEGIES)
# ===============================================================================
# This function implements population averaging for biomarker strategies in PSA.
# Uses the ORIGINAL population (not resampled cohort) for predictions.
# This ensures:
#   1. Consistent population across all PSA iterations
#   2. True biomarker prevalence is maintained
#   3. External validity - predicting for actual patient population
# The resampled model captures parameter uncertainty; population averaging
# integrates out patient heterogeneity.
# See GitHub Issues #73, #74, #87 for methodology discussion.
# ===============================================================================

generate_psa_population_averaged_predictions <- function(sampling_model_list,
                                                         biomarker_name,
                                                         outcome = "os",
                                                         sample_idx = 1,
                                                         data_original,
                                                         time_points) {

  # Extract the resampled model
  sampled_model <- sampling_model_list$samples[[sample_idx]]

  if (is.null(sampled_model)) {
    warning("Sampled model ", sample_idx, " is NULL - returning NULL")
    return(NULL)
  }

  # Get the model object for this outcome
  model_obj <- sampled_model[[outcome]]$model

  if (is.null(model_obj)) {
    warning("Model object for ", outcome, " is NULL in sample ", sample_idx)
    return(NULL)
  }

  # Get treatment levels
  exp_rx <- levels(data_original$Rx)[2]
  ctrl_rx <- levels(data_original$Rx)[1]

  # Initialize outputs
  survival_pos <- NULL
  survival_neg <- NULL

  # -----------------------------------------------------------------------
  # BIOMARKER-POSITIVE SUBGROUP: All biomarker+ patients with experimental Rx
  # -----------------------------------------------------------------------
  # Use ORIGINAL population, not resampled cohort
  # Use explicit type conversion for robustness (works with factor, character, or numeric)
  biomarker_pos_data <- data_original[
    as.numeric(as.character(data_original[[biomarker_name]])) == 1, ]

  if (nrow(biomarker_pos_data) > 0) {
    # Assign experimental treatment to all biomarker+ patients
    biomarker_pos_data$Rx <- factor(exp_rx, levels = levels(data_original$Rx))

    # Predict for ALL biomarker+ patients in original population
    tryCatch({
      pred_pos <- predict(model_obj, newdata = biomarker_pos_data,
                         type = "survival", times = time_points)

      # Extract and average survival probabilities across all patients
      all_surv_probs_pos <- lapply(pred_pos$.pred, function(x) x$.pred_survival)
      surv_matrix_pos <- do.call(cbind, all_surv_probs_pos)
      survival_pos <- rowMeans(surv_matrix_pos, na.rm = TRUE)

    }, error = function(e) {
      warning("Prediction failed for ", biomarker_name, "+ in sim ", sample_idx,
              ": ", conditionMessage(e))
      survival_pos <<- rep(NA, length(time_points))
    })
  } else {
    warning("No biomarker+ patients in original data for ", biomarker_name)
    survival_pos <- rep(NA, length(time_points))
  }

  # -----------------------------------------------------------------------
  # BIOMARKER-NEGATIVE SUBGROUP: All biomarker- patients with control Rx
  # -----------------------------------------------------------------------
  # Use ORIGINAL population, not resampled cohort
  # Use explicit type conversion for robustness (works with factor, character, or numeric)
  biomarker_neg_data <- data_original[
    as.numeric(as.character(data_original[[biomarker_name]])) == 0, ]

  if (nrow(biomarker_neg_data) > 0) {
    # Assign control treatment to all biomarker- patients
    biomarker_neg_data$Rx <- factor(ctrl_rx, levels = levels(data_original$Rx))

    # Predict for ALL biomarker- patients in original population
    tryCatch({
      pred_neg <- predict(model_obj, newdata = biomarker_neg_data,
                         type = "survival", times = time_points)

      # Extract and average survival probabilities across all patients
      all_surv_probs_neg <- lapply(pred_neg$.pred, function(x) x$.pred_survival)
      surv_matrix_neg <- do.call(cbind, all_surv_probs_neg)
      survival_neg <- rowMeans(surv_matrix_neg, na.rm = TRUE)

    }, error = function(e) {
      warning("Prediction failed for ", biomarker_name, "- in sim ", sample_idx,
              ": ", conditionMessage(e))
      survival_neg <<- rep(NA, length(time_points))
    })
  } else {
    warning("No biomarker- patients in original data for ", biomarker_name)
    survival_neg <- rep(NA, length(time_points))
  }

  # Return both subgroup predictions
  return(list(
    positive = survival_pos,
    negative = survival_neg
  ))
}

# ===============================================================================
# CONTROL STRATEGY PSA PREDICTION FUNCTION
# ===============================================================================
# This function generates population-averaged predictions for the CONTROL strategy
# using the ORIGINAL control population (not the resampled cohort).
# This matches the methodology used for biomarker strategies (Issue #97).
# ===============================================================================

generate_psa_control_predictions <- function(sampling_model_list,
                                              outcome = "os",
                                              sample_idx = 1,
                                              data_control_original,
                                              time_points) {
  # Extract the resampled model
  sampled_model <- sampling_model_list$samples[[sample_idx]]

  if (is.null(sampled_model)) {
    warning("Sampled model ", sample_idx, " is NULL - returning NULL")
    return(NULL)
  }

  model_obj <- sampled_model[[outcome]]$model

  if (is.null(model_obj)) {
    warning("Model object for ", outcome, " is NULL in sample ", sample_idx)
    return(NULL)
  }

  # Predict for ALL patients in original control population
  tryCatch({
    pred <- predict(model_obj, newdata = data_control_original,
                    type = "survival", times = time_points)

    # Extract and average survival probabilities across all patients
    all_surv_probs <- lapply(pred$.pred, function(x) x$.pred_survival)
    surv_matrix <- do.call(cbind, all_surv_probs)
    survival_avg <- rowMeans(surv_matrix, na.rm = TRUE)

    return(survival_avg)

  }, error = function(e) {
    warning("Prediction failed for control in sim ", sample_idx,
            ": ", conditionMessage(e))
    return(rep(NA, length(time_points)))
  })
}

cat("\n=== Parameter distributions configured ===\n")
cat("- Cost parameters: Gamma distributions (CV =", cv_costs, ")\n")
cat("- Utility parameters: Beta distributions (CV =", cv_utilities, ")\n")
cat("- Survival parameters: Correlated resampled models (n =", sampling_models$control$n_samples, ")\n")
cat("- PSA population averaging function ready\n")
cat("\nReady for PSA analysis\n")