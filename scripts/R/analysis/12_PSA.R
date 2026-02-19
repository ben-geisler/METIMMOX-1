# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/psa_functions.R"))

# ===============================================================================
# RUN_ALL_MODELS SWITCH
# ===============================================================================
# When RUN_ALL_MODELS = TRUE, this script will loop through all 3 MODEL_STRUCTURE
# values (0, 1, 2) and generate PSA caches for each. This is useful for pre-generating
# all caches needed for multi-model CEA comparison.
#
# When RUN_ALL_MODELS = FALSE (default), the script uses the current MODEL_STRUCTURE
# value and generates a single PSA cache file.
# ===============================================================================

if (!exists("RUN_ALL_MODELS")) {
  RUN_ALL_MODELS <- FALSE
}

if (RUN_ALL_MODELS) {
  cat("\n=== RUN_ALL_MODELS mode: Generating PSA caches for all model structures ===\n")
  model_structures_to_run <- c(0, 1, 2)
  original_model_structure <- MODEL_STRUCTURE
} else {
  model_structures_to_run <- MODEL_STRUCTURE
}

# Loop through model structures (single iteration if RUN_ALL_MODELS = FALSE)
for (current_model_structure in model_structures_to_run) {

  # Set MODEL_STRUCTURE for this iteration
  MODEL_STRUCTURE <- current_model_structure
  model_structure_label <- c("joint", "focused", "separate")[MODEL_STRUCTURE + 1]

  if (RUN_ALL_MODELS) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("Processing MODEL_STRUCTURE =", MODEL_STRUCTURE, "(", model_structure_label, ")\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")

    # Re-source scripts 06, 07, 08 to update model and sampling for this MODEL_STRUCTURE
    cat("Re-sourcing analysis scripts for new MODEL_STRUCTURE...\n")
    source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
    source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))
    # Don't re-source 08 - use existing sampling_models for each structure
    # The sampling cache is loaded based on MODEL_STRUCTURE in 08_sampling.R
    # Instead, load the appropriate sampling cache directly
    sampling_cache_file <- here("data", "tidy", paste0("sampling_models_n", n_samples, "_",
                                                       ifelse(USE_BOTH_MODELS == 0, "full", "both"), "_",
                                                       model_structure_label, ".rds"))
    if (file.exists(sampling_cache_file)) {
      sampling_models <- readRDS(sampling_cache_file)
      cat("Loaded sampling cache:", sampling_cache_file, "\n")
    } else {
      cat("ERROR: Sampling cache not found:", sampling_cache_file, "\n")
      cat("Please run 08_sampling.R with RUN_ALL_MODELS = TRUE first.\n")
      next
    }
  }

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}

# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Define cache file paths with model structure label
cache_file_obj <- here("data", "tidy", paste0("psa_obj_", model_structure_label, ".rds"))
cache_file_params <- here("data", "tidy", paste0("psa_params_", model_structure_label, ".rds"))

# Check if PSA cache exists and is valid
psa_cached <- FALSE
if (file.exists(cache_file_obj) && file.exists(cache_file_params)) {
  cat("PSA cache files detected. Attempting to load...\n")

  # Try to load cached PSA results
  tryCatch({
    psa_obj <- readRDS(cache_file_obj)
    psa_params <- readRDS(cache_file_params)

    # Validate cache
    cache_valid <- TRUE

    # Check if n_sim matches
    if (psa_obj$n_sim != n_sim) {
      cat("  Cache validation failed: n_sim mismatch (cached:", psa_obj$n_sim,
          ", expected:", n_sim, ")\n")
      cache_valid <- FALSE
    }

    # Check if strategies match
    if (!all(psa_obj$strategies == strategies)) {
      cat("  Cache validation failed: strategy mismatch\n")
      cache_valid <- FALSE
    }

    # Check if psa_params has correct number of rows
    if (nrow(psa_params) != n_sim) {
      cat("  Cache validation failed: psa_params row count mismatch\n")
      cache_valid <- FALSE
    }

    if (cache_valid) {
      cat("  Cache validation successful!\n")
      cat("  Loaded PSA results from cache:\n")
      cat("    - File:", cache_file_obj, "\n")
      cat("    - Modified:", format(file.info(cache_file_obj)$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
      cat("    - Simulations:", psa_obj$n_sim, "\n")
      cat("    - Strategies:", psa_obj$n_strategies, "\n")
      psa_cached <- TRUE
    } else {
      cat("  Cache invalid. PSA will be regenerated.\n")
      psa_cached <- FALSE
    }

  }, error = function(e) {
    cat("  Error loading cache:", e$message, "\n")
    cat("  PSA will be regenerated.\n")
    psa_cached <- FALSE
  })
}

# Only run PSA if cache was not loaded
if (!psa_cached) {
  cat("\n=== Generating new PSA results ===\n")

  # Generate PSA samples
  cat("Generating PSA samples for", n_sim, "simulations\n")
  psa_params <- generate_psa_samples(param_distributions, n_sim)

  # Run the PSA
  cat("Starting PSA with", n_sim, "simulations\n")
  psa_results <- run_psa_analysis(
    psa_params = psa_params,
    l_params_base = l_params_base,
    param_distributions = param_distributions,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim
  )

  # Create the PSA object using dampack's make_psa_obj function
  psa_obj <- dampack::make_psa_obj(
    cost = as.data.frame(psa_results$cost),
    effect = as.data.frame(psa_results$effect),
    strategies = strategies,
    currency = "€"
  )

  cat("\n=== PSA generation complete ===\n")
}

# Verify the PSA object structure
cat("PSA object created with:\n")
cat(" - Simulations:", psa_obj$n_sim, "\n")
cat(" - Strategies:", psa_obj$n_strategies, "\n")
cat(" - Cost matrix dimensions:", dim(psa_obj$cost), "\n")
cat(" - Effect matrix dimensions:", dim(psa_obj$effect), "\n")

psa_summary <- summary(psa_obj)
print(psa_summary)

psa_summary_c95ci <- summary(psa_obj, calc_sds = TRUE)
print(psa_summary_c95ci)

# Calculate ICERs from mean PSA costs and effects
mean_costs <- colMeans(psa_obj$cost, na.rm = TRUE)
mean_effects <- colMeans(psa_obj$effect, na.rm = TRUE)
cea_psa <- calculate_icers(
  cost = mean_costs,
  effect = mean_effects,
  strategies = psa_obj$strategies
)
print(cea_psa)

# ICE scatter plot
plot(psa_obj, frontier = TRUE, points = TRUE)

# Create cost-effectiveness acceptability curves
c_wtp <- seq(from = 0, to = 2e5, by = 5e3)
ceac_obj <- ceac(c_wtp, psa_obj)
ceac_sum <- summary(ceac_obj)
print(ceac_sum)
plot(ceac_obj, frontier = TRUE, points = TRUE, currency = "€")

# Save results (only if newly generated)
if (!psa_cached) {
  cat("\nSaving PSA results to cache...\n")
  tryCatch({
    saveRDS(psa_params, cache_file_params)
    saveRDS(psa_obj, cache_file_obj)
    cat("  PSA results saved successfully to:\n")
    cat("    - ", cache_file_params, "\n")
    cat("    - ", cache_file_obj, "\n")
  }, error = function(e) {
    warning("Failed to save PSA cache: ", e$message)
  })
} else {
  cat("\nUsing cached PSA results (not saving).\n")
}

} # End of RUN_ALL_MODELS loop

# Restore original MODEL_STRUCTURE if we were in RUN_ALL_MODELS mode
if (RUN_ALL_MODELS) {
  MODEL_STRUCTURE <- original_model_structure
  cat("\n=== RUN_ALL_MODELS complete ===\n")
  cat("Generated PSA caches for MODEL_STRUCTURE: 0 (joint), 1 (focused), 2 (separate)\n")
  cat("Restored MODEL_STRUCTURE to:", MODEL_STRUCTURE, "\n")
}