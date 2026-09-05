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
sampling_cache_dir <- cache_dir()
if (!dir.exists(sampling_cache_dir)) {
  dir.create(sampling_cache_dir, recursive = TRUE)
  cat("Created sampling cache directory:", sampling_cache_dir, "\n")
}

# Create cache file path based on n_samples
cache_file <- sampling_cache_path(n_samples)
sampling_seed <- analysis_seed

# ===============================================================================
# MODEL FORMULA SELECTION
# ===============================================================================
# Model formulas are defined in scripts/R/functions/model_configs.R
# This ensures PSA uses the same single joint model as base case.
#
# Economic biomarker strategies are CRP and TMB/BRAF. Both share one joint
# formula, so one bootstrap of the joint model serves every strategy. The
# control arm is NOT a separate age/sex model: it is the same joint fit
# predicted with Rx forced to the control level, exactly as in the base case
# (issue #151). The cache therefore holds one component per biomarker, all
# pointing at the same joint bootstrap; a legacy "control" component, if
# present in an older cache, is ignored.
# ===============================================================================

# Source model configurations if not already loaded
if (!exists("get_model_configs")) {
  source(here::here("scripts/R/functions/model_configs.R"))
}
if (!exists("extract_all_survival_probabilities", mode = "function")) {
  source(here::here("scripts/R/functions/prediction_functions.R"))
}

# Get model configuration
model_config <- get_current_model_config()
cat("Resampling: Using", model_config$label, "\n")
cat("Description:", model_config$description, "\n")

# Create biomarker formula function using central config
create_biomarker_formula <- function(outcome, biomarker) {
  return(get_strategy_formula(biomarker, outcome))
}

# ===============================================================================
# CORRELATED SAMPLING FUNCTION
# ===============================================================================
# This function implements non-parametric resampling that maintains correlation
# between PFS and OS by resampling patients and fitting both models to the same
# resampled dataset.
# ===============================================================================

sample_correlated_survival <- function(formula_os, formula_pfs, data,
                                       dist_os = "weibull", dist_pfs = "weibull",
                                       n_samples = n_samples, seed = 123L) {

  n_patients <- nrow(data)
  sampled_models <- vector("list", n_samples)
  n_failed <- 0

  cat("Generating", n_samples, "correlated resampled models for", n_patients, "patients...\n")
  cat("Using OS distribution:", dist_os, "| PFS distribution:", dist_pfs, "\n")
  cat("This preserves PFS-OS correlation by fitting both models to the same resampled data\n")

  # Fit original models for reference
  original_os <- flexsurvreg(formula_os, data = data, dist = dist_os)
  original_pfs <- flexsurvreg(formula_pfs, data = data, dist = dist_pfs)

  # Make bootstrap indices independent of prior RNG use and cache branches.
  set.seed(seed)

  # Set up progress reporting
  start_time <- Sys.time()

  for (i in seq_len(n_samples)) {
    # Resample patients with replacement (non-parametric resampling)
    resample_idx <- sample(1:n_patients, size = n_patients, replace = TRUE)
    resampled_data <- data[resample_idx, ]

    # Fit BOTH models to the SAME resampled dataset
    # This is critical: PFS and OS share the same patient cohort
    sampled_models[[i]] <- tryCatch({
      os_model <- flexsurvreg(formula_os, data = resampled_data, dist = dist_os)
      pfs_model <- flexsurvreg(formula_pfs, data = resampled_data, dist = dist_pfs)

      list(
        os = list(
          model = os_model,
          coefficients = os_model$coefficients,
          dist = dist_os
        ),
        pfs = list(
          model = pfs_model,
          coefficients = pfs_model$coefficients,
          dist = dist_pfs
        ),
        resample_idx = resample_idx  # Store for reproducibility/debugging
      )
    }, error = function(e) {
      # If resampled model fails, use original model
      cat("Warning: Resampled model", i, "failed:", conditionMessage(e), "\n")
      list(
        os = list(
          model = original_os,
          coefficients = original_os$coefficients,
          dist = dist_os
        ),
        pfs = list(
          model = original_pfs,
          coefficients = original_pfs$coefficients,
          dist = dist_pfs
        ),
        failed = TRUE
      )
    })

    if (isTRUE(sampled_models[[i]]$failed)) {
      n_failed <- n_failed + 1
    }
    
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
    dist_os = dist_os,
    dist_pfs = dist_pfs,
    seed = seed,
    creation_time = Sys.time()
  ))
}

# Return TRUE only when every requested cache component records the expected seed.
sampling_cache_seed_matches <- function(sampling_models, components,
                                        seed = 123L) {
  all(vapply(
    components,
    function(component) identical(sampling_models[[component]]$seed, seed),
    logical(1)
  ))
}

# ===============================================================================
# DISTRIBUTION RESOLUTION HELPER
# ===============================================================================
# Resolves the AIC-best parametric distributions from the models object.
# Returns a list with $os and $pfs distribution names.
# ===============================================================================

resolve_best_distributions <- function(models) {
  default_dist <- "weibull"

  if (is.null(models)) {
    cat("WARNING: 'models' is NULL. Using default distribution:", default_dist, "\n")
    return(list(os = default_dist, pfs = default_dist))
  }

  os_dist <- models$best_fit$os_distribution
  pfs_dist <- models$best_fit$pfs_distribution

  if (is.null(os_dist)) {
    cat("WARNING: No selected OS distribution found. Using default:", default_dist, "\n")
    os_dist <- default_dist
  }
  if (is.null(pfs_dist)) {
    cat("WARNING: No selected PFS distribution found. Using default:", default_dist, "\n")
    pfs_dist <- default_dist
  }

  cat("Resolved distributions - OS:", os_dist, "| PFS:", pfs_dist, "\n")
  return(list(os = os_dist, pfs = pfs_dist))
}

# ===============================================================================
# CHECK CACHE OR GENERATE SAMPLING MODELS
# ===============================================================================

if (file.exists(cache_file)) {
  cat("\n=== Loading sampling models from cache ===\n")
  cat("Cache file:", cache_file, "\n")

  sampling_models <- readRDS(cache_file)

  # Validate cache contains the single economic model biomarkers. Each
  # biomarker component holds the same joint bootstrap, which also supplies the
  # control arm (issue #151); a legacy "control" component is not required.
  biomarkers_to_validate <- get_biomarkers()
  if (is.list(sampling_models) &&
      all(biomarkers_to_validate %in% names(sampling_models))) {
    expected_dists <- resolve_best_distributions(models)
    cache_components <- biomarkers_to_validate
    distributions_match <- all(vapply(
      cache_components,
      function(component) {
        cached <- sampling_models[[component]]
        identical(cached$dist_os, expected_dists$os) &&
          identical(cached$dist_pfs, expected_dists$pfs)
      },
      logical(1)
    ))
    seed_matches <- sampling_cache_seed_matches(
      sampling_models,
      cache_components,
      sampling_seed
    )

    if (!distributions_match) {
      cat("Cache distributions do not match the ordering-constrained base case ",
          "(OS: ", expected_dists$os, ", PFS: ", expected_dists$pfs,
          "). Regenerating...\n", sep = "")
      sampling_models <- NULL
    } else if (!seed_matches) {
      cat("Cache RNG seed is missing or does not match the requested seed (",
          sampling_seed, "). Regenerating...\n", sep = "")
      sampling_models <- NULL
    } else {
      cat("Cache loaded successfully!\n")
      joint_models <- get_joint_sampling_models(sampling_models)
      cat("- Joint-model samples:", joint_models$n_samples, "\n")
      cat("- Biomarkers:", paste(biomarkers_to_validate, collapse = ", "), "\n")
      cat("- Created:", format(joint_models$creation_time), "\n")
      cat("- Distributions: OS", joint_models$dist_os,
          "| PFS", joint_models$dist_pfs, "\n")
      cat("- RNG seed:", joint_models$seed, "\n")
      if (get_control_strategy() %in% names(sampling_models)) {
        cat("- Legacy control component present and ignored: the control arm ",
            "is predicted from the joint model with Rx = control (issue #151)\n",
            sep = "")
      }

      # Check if n_samples matches
      if (joint_models$n_samples != n_samples) {
        cat("WARNING: Cached models (", joint_models$n_samples,
            ") != requested (", n_samples, ")\n")
        cat("Will use cached models. Delete cache file to regenerate.\n")
      }
      rm(joint_models)
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

  # Bootstrap the joint economic model (age, sex, treatment, and the
  # biomarker-by-treatment interactions). No separate control bootstrap is
  # fitted: the control arm is this same joint fit predicted with Rx = control,
  # matching the base case (issue #151). Biomarker strategies that share a
  # formula share one bootstrap, stored under each biomarker key.
  selected_dists <- resolve_best_distributions(models)
  biomarkers_to_sample <- get_biomarkers()
  fitted_formula_sets <- list()
  for (biomarker in biomarkers_to_sample) {
    os_formula <- create_biomarker_formula("os", biomarker)
    pfs_formula <- create_biomarker_formula("pfs", biomarker)
    formula_key <- paste(deparse(os_formula), deparse(pfs_formula), sep = " || ")

    if (!is.null(fitted_formula_sets[[formula_key]])) {
      cat("\n", biomarker, " strategy shares the joint model bootstrap already ",
          "fitted for ", fitted_formula_sets[[formula_key]], "\n", sep = "")
      sampling_models[[biomarker]] <-
        sampling_models[[fitted_formula_sets[[formula_key]]]]
      next
    }

    cat("\nSampling joint model bootstrap for", biomarker, "strategy...\n")
    sampling_models[[biomarker]] <- sample_correlated_survival(
      formula_os = os_formula,
      formula_pfs = pfs_formula,
      data = data,
      dist_os = selected_dists$os,
      dist_pfs = selected_dists$pfs,
      n_samples = n_samples,
      seed = sampling_seed
    )
    fitted_formula_sets[[formula_key]] <- biomarker
  }
  rm(fitted_formula_sets)

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

# Print confirmation of model structure
cat("\n=== Sampling model structure ready ===\n")
cat("- Model:", model_config$label, "\n")
cat("- Biomarker strategies: CRP, TMB/BRAF\n")
cat("- Formulas defined in: scripts/R/functions/model_configs.R\n")
cat("- Number of resampled models:",
    get_joint_sampling_models(sampling_models)$n_samples, "\n")
cat("- PFS and OS are CORRELATED within each resampled model\n")
cat("- Control arm: joint model predicted with Rx = control (issue #151)\n")

# ===============================================================================
# PARAMETER DISTRIBUTIONS FOR UNCERTAINTY ANALYSIS (NON-SURVIVAL PARAMETERS)
# ===============================================================================

# Survival curves are excluded because they come from correlated resampling.
parameter_config <- configure_parameter_distributions(l_params_base)
param_distributions <- parameter_config$distributions
param_groups <- parameter_config$groups
parameter_spec <- parameter_config$spec
rm(parameter_config)

# ===============================================================================
# POPULATION AVERAGING FUNCTION FOR PSA
# ===============================================================================
# With biomarker_name = NULL, predicts the control curve over all patients with
# Rx forced to the control level, i.e. the joint model's standard-of-care
# prediction, mirroring generate_population_averaged_predictions() in the base
# case (issue #151). Otherwise predicts the biomarker-positive/experimental and
# biomarker-negative/control subgroup curves using one shared code path.

generate_psa_population_averaged_predictions <- function(sampling_model_list,
                                                         biomarker_name = NULL,
                                                         outcome = "os",
                                                         sample_idx = 1,
                                                         data_original,
                                                         time_points) {
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

  average_prediction <- function(subgroup_data, rx_level = NULL, label) {
    if (nrow(subgroup_data) == 0) {
      warning("No patients in ", label, " subgroup")
      return(rep(NA_real_, length(time_points)))
    }
    if (!is.null(rx_level)) {
      subgroup_data$Rx <- factor(rx_level, levels = levels(data_original$Rx))
    }
    tryCatch({
      prediction <- predict(
        model_obj, newdata = subgroup_data,
        type = "survival", times = time_points
      )
      rowMeans(extract_all_survival_probabilities(prediction), na.rm = TRUE)
    }, error = function(e) {
      warning("Prediction failed for ", label, " in sim ", sample_idx,
              ": ", conditionMessage(e))
      rep(NA_real_, length(time_points))
    })
  }

  if (is.null(biomarker_name)) {
    return(average_prediction(
      data_original,
      rx_level = levels(data_original$Rx)[1],
      label = "control"
    ))
  }

  status <- as.numeric(as.character(data_original[[biomarker_name]]))
  subgroup_spec <- list(
    positive = list(status = 1, rx = levels(data_original$Rx)[2], suffix = "+"),
    negative = list(status = 0, rx = levels(data_original$Rx)[1], suffix = "-")
  )
  lapply(subgroup_spec, function(subgroup) {
    average_prediction(
      data_original[status == subgroup$status, , drop = FALSE],
      rx_level = subgroup$rx,
      label = paste0(biomarker_name, subgroup$suffix)
    )
  })
}

# Report what the PSA actually samples. Since issue #154 the specification also
# holds parameters that are fixed in the PSA (unit prices) and one that is
# derived rather than drawn (u_p), so the summary is built per group instead of
# assuming a single cost CV and a single utility CV.
sampled_spec <- parameter_spec[parameter_spec$psa & !parameter_spec$derived, ,
                               drop = FALSE]
fixed_spec <- parameter_spec[!parameter_spec$psa, , drop = FALSE]
derived_spec <- parameter_spec[parameter_spec$derived, , drop = FALSE]

cat("\n=== Parameter distributions configured ===\n")
for (group_name in unique(sampled_spec$group)) {
  rows <- sampled_spec[sampled_spec$group == group_name, , drop = FALSE]
  cat(sprintf("- %s: %s (CV = %s): %s\n",
              group_name,
              paste(unique(rows$distribution), collapse = "/"),
              paste(unique(rows$cv), collapse = "/"),
              paste(rows$parameter, collapse = ", ")))
}
if (nrow(derived_spec) > 0) {
  cat("- Derived in the PSA (not drawn):",
      paste(derived_spec$parameter, collapse = ", "), "\n")
}
if (nrow(fixed_spec) > 0) {
  cat("- Fixed in the PSA, varied in the DSA only (issue #154):",
      paste(fixed_spec$parameter, collapse = ", "), "\n")
}
cat("- Survival parameters: Correlated resampled models (n =",
    sampling_models[[get_control_strategy()]]$n_samples, ")\n")
cat("- PSA population averaging function ready\n")
cat("\nReady for PSA analysis\n")
