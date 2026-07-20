# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/psa_functions.R"))

# Define cache file paths with utility source label
cache_file_obj <- psa_obj_path()
cache_file_params <- psa_params_path()
psa_seed <- analysis_seed

# Check if PSA cache exists and is valid
psa_cached <- FALSE
if (file.exists(cache_file_obj) && file.exists(cache_file_params)) {
  cat("PSA cache files detected. Attempting to load...\n")

  # Try to load cached PSA results
  psa_cached <- tryCatch({
    psa_obj <- readRDS(cache_file_obj)
    psa_params <- readRDS(cache_file_params)

    # Validate cache
    cache_valid <- TRUE

    # The effective n_sim may be smaller when failed draws were dropped.
    cached_requested_n_sim <- if (!is.null(psa_obj$requested_n_sim)) {
      psa_obj$requested_n_sim
    } else {
      psa_obj$n_sim
    }
    if (cached_requested_n_sim != n_sim) {
      cat("  Cache validation failed: requested n_sim mismatch (cached:", cached_requested_n_sim,
          ", expected:", n_sim, ")\n")
      cache_valid <- FALSE
    }

    # Check if strategies match
    if (!all(psa_obj$strategies == strategies)) {
      cat("  Cache validation failed: strategy mismatch\n")
      cache_valid <- FALSE
    }

    # Cost/effect outcomes and parameter samples must remain row-aligned.
    if (nrow(psa_params) != psa_obj$n_sim) {
      cat("  Cache validation failed: psa_params/PSA row count mismatch\n")
      cache_valid <- FALSE
    }

    # Legacy caches may contain mean-imputed draws or the old shared-budget
    # replacement loop and must be rebuilt.
    if (!identical(psa_obj$failed_draw_policy, psa_failed_draw_policy())) {
      cat("  Cache validation failed: legacy failed-draw handling policy\n")
      cache_valid <- FALSE
    }

    # Legacy caches without seed metadata, and caches made with another seed,
    # cannot guarantee reproducibility against a fresh run.
    if (!psa_samples_seed_matches(psa_params, psa_seed)) {
      cat("  Cache validation failed: RNG seed missing or mismatched (expected:",
          psa_seed, ")\n")
      cache_valid <- FALSE
    }

    # Parameter-set changes (including new diagnostic costs) invalidate old caches.
    missing_param_columns <- setdiff(names(param_distributions), names(psa_params))
    if (length(missing_param_columns) > 0) {
      cat("  Cache validation failed: missing PSA parameters:",
          paste(missing_param_columns, collapse = ", "), "\n")
      cache_valid <- FALSE
    }

    if (cache_valid) {
      cat("  Cache validation successful!\n")
      cat("  Loaded PSA results from cache:\n")
      cat("    - File:", cache_file_obj, "\n")
      cat("    - Modified:", format(file.info(cache_file_obj)$mtime, "%Y-%m-%d %H:%M:%S"), "\n")
      cat("    - Simulations:", psa_obj$n_sim, "\n")
      cat("    - Strategies:", psa_obj$n_strategies, "\n")
    } else {
      cat("  Cache invalid. PSA will be regenerated.\n")
    }

    cache_valid

  }, error = function(e) {
    cat("  Error loading cache:", e$message, "\n")
    cat("  PSA will be regenerated.\n")
    FALSE
  })
}

# Only run PSA if cache was not loaded
if (!psa_cached) {
  cat("\n=== Generating new PSA results ===\n")

  cat("Starting PSA with", n_sim, "simulations\n")
  psa_build <- build_psa_obj(
    l_params_base = l_params_base,
    param_distributions = param_distributions,
    strategies = strategies,
    time_horizon = time_horizon,
    cl = cl,
    n_sim = n_sim,
    seed = psa_seed
  )
  psa_obj <- psa_build$psa_obj
  psa_params <- psa_build$psa_params
  psa_results <- psa_build$psa_results
  rm(psa_build)

  cat(sprintf(
    "PSA fallback diagnostic: %d/%d initial iterations (%.2f%%; maximum permitted %.2f%%)\n",
    psa_results$fallback_count, n_sim, 100 * psa_results$fallback_rate,
    100 * psa_results$fallback_threshold
  ))

  cat(sprintf(
    "PSA dropped-draw diagnostic: %d unrecoverable iteration(s); effective n_sim = %d\n",
    psa_results$dropped_count, psa_results$n_sim
  ))

  cat("\n=== PSA generation complete ===\n")
} else {
  cat(paste0(
    "PSA fallback diagnostic: unavailable for the existing cache ",
    "(the PSA loop was not run).\n"
  ))
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
