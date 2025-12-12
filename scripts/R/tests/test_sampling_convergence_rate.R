# ==============================================================================
# Sampling Convergence Rate Analysis Script
# ==============================================================================
# This script quantifies how many bootstrap resamples fail to converge during
# the sample_correlated_survival() function in 08_sampling.R.
#
# Two modes:
#   1. Inspect existing cache - Check stored failure counts (fast)
#   2. Fresh sampling test - Re-run sampling with configurable n (slower)
#
# USAGE: Adjust settings below (after dependencies), then run the entire script.
# ==============================================================================

# =============================================================================
# SETUP - Load dependencies first (these clear the environment)
# =============================================================================
cat("=== Sampling Convergence Rate Analysis ===\n\n")

# Set working directory
setwd("c:/Users/benjampg/git/METIMMOX-1")

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr)

# Source required setup scripts (NOTE: 02_setup clears the environment!)
cat("Loading dependencies...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))

cat("Dependencies loaded.\n\n")

# =============================================================================
# CONFIGURABLE PARAMETERS - Set AFTER dependencies (which clear environment)
# =============================================================================
n_test_samples <- 100   # Number of fresh bootstrap samples to test
run_fresh_test <- FALSE # Set to TRUE to run fresh sampling test (slower)

# =============================================================================
# PART 1: INSPECT EXISTING CACHE
# =============================================================================
cat("===============================================================================\n")
cat("                    PART 1: EXISTING CACHE INSPECTION                          \n")
cat("===============================================================================\n\n")

# Define cache file path (matching 08_sampling.R)
cache_file <- here("data", "tidy", paste0("sampling_models_n", n_samples, "_",
                                          ifelse(USE_BOTH_MODELS == 0, "full", "both"),
                                          ".rds"))

if (file.exists(cache_file)) {
  cat("Cache file found:", cache_file, "\n")
  cat("File size:", format(file.info(cache_file)$size / 1024^2, digits = 2), "MB\n")
  cat("Modified:", format(file.info(cache_file)$mtime, "%Y-%m-%d %H:%M:%S"), "\n\n")

  # Load cache
  sampling_models <- readRDS(cache_file)

  # Strategy names
  strategy_names <- c("control", biomarkers)

  cat("FAILURE COUNTS FROM CACHE:\n")
  cat("--------------------------\n")

  total_samples <- 0
  total_failed <- 0

  for (strategy in strategy_names) {
    if (!is.null(sampling_models[[strategy]])) {
      n_samp <- sampling_models[[strategy]]$n_samples
      n_fail <- sampling_models[[strategy]]$n_failed

      # Also count failed = TRUE markers in individual samples
      n_fail_markers <- sum(sapply(sampling_models[[strategy]]$samples, function(x) {
        isTRUE(x$failed)
      }))

      # Use the larger of the two counts (in case n_failed wasn't updated correctly)
      n_fail_actual <- max(n_fail, n_fail_markers)

      cat(sprintf("  %-15s %5d / %5d failed (%5.2f%%)\n",
                  paste0(strategy, ":"), n_fail_actual, n_samp,
                  100 * n_fail_actual / n_samp))

      total_samples <- total_samples + n_samp
      total_failed <- total_failed + n_fail_actual
    } else {
      cat(sprintf("  %-15s [NOT FOUND IN CACHE]\n", paste0(strategy, ":")))
    }
  }

  cat("--------------------------\n")
  cat(sprintf("  %-15s %5d / %5d failed (%5.2f%%)\n",
              "TOTAL:", total_failed, total_samples,
              100 * total_failed / total_samples))

  cat("\nCache created:", format(sampling_models$control$creation_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Distribution used:", sampling_models$control$dist, "\n")

} else {
  cat("No cache file found at:", cache_file, "\n")
  cat("Run 08_sampling.R first to generate the sampling cache.\n")
}

# =============================================================================
# PART 2: FRESH SAMPLING TEST (OPTIONAL)
# =============================================================================

# Safeguard: define default if not set
if (!exists("run_fresh_test")) run_fresh_test <- FALSE
if (!exists("n_test_samples")) n_test_samples <- 100

if (run_fresh_test) {
  cat("\n")
  cat("===============================================================================\n")
  cat("                    PART 2: FRESH SAMPLING TEST                                \n")
  cat("===============================================================================\n\n")

  cat("Running fresh sampling test with n =", n_test_samples, "samples per strategy\n")
  cat("This tests current convergence rates with smaller sample size.\n\n")

  # Get distribution from parametric analysis
  if (exists("models") && !is.null(models$best_fit$os_distribution)) {
    test_dist <- models$best_fit$os_distribution
  } else {
    test_dist <- "gengamma"  # Default
  }
  cat("Distribution:", test_dist, "\n\n")

  # Define formulas (matching 08_sampling.R)
  control_os_formula <- Surv(OSwk, Death) ~ Age + sex
  control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex

  create_biomarker_formula <- function(outcome, biomarker) {
    if (outcome == "os") {
      return(as.formula(paste0("Surv(OSwk, Death) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    } else {
      return(as.formula(paste0("Surv(PFSwk, Progression) ~ Age + sex + Rx + ", biomarker, ":Rx")))
    }
  }

  # Function to test sampling convergence (simplified version)
  test_sampling_convergence <- function(formula_os, formula_pfs, data,
                                        dist, n_samples, strategy_name) {
    n_patients <- nrow(data)
    n_failed <- 0
    n_os_failed <- 0
    n_pfs_failed <- 0
    error_messages <- character(0)

    cat("Testing", strategy_name, "(n =", n_patients, "patients)...\n")

    # Fit original models for fallback
    original_os <- tryCatch({
      flexsurvreg(formula_os, data = data, dist = dist)
    }, error = function(e) {
      cat("  ERROR: Original OS model failed to fit!\n")
      return(NULL)
    })

    original_pfs <- tryCatch({
      flexsurvreg(formula_pfs, data = data, dist = dist)
    }, error = function(e) {
      cat("  ERROR: Original PFS model failed to fit!\n")
      return(NULL)
    })

    if (is.null(original_os) || is.null(original_pfs)) {
      return(list(
        n_samples = n_samples,
        n_failed = n_samples,
        n_os_failed = n_samples,
        n_pfs_failed = n_samples,
        error_messages = "Original model fitting failed"
      ))
    }

    start_time <- Sys.time()

    for (i in 1:n_samples) {
      # Resample patients
      resample_idx <- sample(1:n_patients, size = n_patients, replace = TRUE)
      resampled_data <- data[resample_idx, ]

      os_ok <- TRUE
      pfs_ok <- TRUE

      # Try OS model
      tryCatch({
        os_model <- flexsurvreg(formula_os, data = resampled_data, dist = dist)
      }, error = function(e) {
        os_ok <<- FALSE
        n_os_failed <<- n_os_failed + 1
        error_messages <<- c(error_messages, paste0("OS[", i, "]: ", conditionMessage(e)))
      })

      # Try PFS model
      tryCatch({
        pfs_model <- flexsurvreg(formula_pfs, data = resampled_data, dist = dist)
      }, error = function(e) {
        pfs_ok <<- FALSE
        n_pfs_failed <<- n_pfs_failed + 1
        error_messages <<- c(error_messages, paste0("PFS[", i, "]: ", conditionMessage(e)))
      })

      # Count as failed if either model failed
      if (!os_ok || !pfs_ok) {
        n_failed <- n_failed + 1
      }

      # Progress
      if (i %% 25 == 0) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        cat(sprintf("  Progress: %d/%d - %d failures so far (%.1fs)\n",
                    i, n_samples, n_failed, elapsed))
      }
    }

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("  Complete: %d/%d failed (%.2f%%) in %.1fs\n",
                n_failed, n_samples, 100 * n_failed / n_samples, elapsed))

    return(list(
      n_samples = n_samples,
      n_failed = n_failed,
      n_os_failed = n_os_failed,
      n_pfs_failed = n_pfs_failed,
      error_messages = error_messages
    ))
  }

  # Test each strategy
  fresh_results <- list()

  cat("\n--- Testing Control Strategy ---\n")
  fresh_results$control <- test_sampling_convergence(
    formula_os = control_os_formula,
    formula_pfs = control_pfs_formula,
    data = data_control,
    dist = test_dist,
    n_samples = n_test_samples,
    strategy_name = "Control"
  )

  for (biomarker in biomarkers) {
    cat("\n--- Testing", toupper(biomarker), "Strategy ---\n")
    fresh_results[[biomarker]] <- test_sampling_convergence(
      formula_os = create_biomarker_formula("os", biomarker),
      formula_pfs = create_biomarker_formula("pfs", biomarker),
      data = data,
      dist = test_dist,
      n_samples = n_test_samples,
      strategy_name = biomarker
    )
  }

  # Summary
  cat("\n")
  cat("FRESH SAMPLING TEST RESULTS:\n")
  cat("----------------------------\n")

  total_test_samples <- 0
  total_test_failed <- 0

  for (strategy in c("control", biomarkers)) {
    res <- fresh_results[[strategy]]
    cat(sprintf("  %-15s %3d / %3d failed (%5.2f%%)  [OS: %d, PFS: %d]\n",
                paste0(strategy, ":"), res$n_failed, res$n_samples,
                100 * res$n_failed / res$n_samples,
                res$n_os_failed, res$n_pfs_failed))

    total_test_samples <- total_test_samples + res$n_samples
    total_test_failed <- total_test_failed + res$n_failed
  }

  cat("----------------------------\n")
  cat(sprintf("  %-15s %3d / %3d failed (%5.2f%%)\n",
              "TOTAL:", total_test_failed, total_test_samples,
              100 * total_test_failed / total_test_samples))

  # Show sample error messages if any
  all_errors <- unlist(lapply(fresh_results, function(x) x$error_messages))
  if (length(all_errors) > 0) {
    cat("\nSAMPLE ERROR MESSAGES (first 10):\n")
    for (msg in head(all_errors, 10)) {
      cat("  ", msg, "\n")
    }
    if (length(all_errors) > 10) {
      cat("  ... and", length(all_errors) - 10, "more\n")
    }
  }
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n")
cat("===============================================================================\n")
cat("                           ANALYSIS COMPLETE                                   \n")
cat("===============================================================================\n")

if (exists("total_failed") && total_samples > 0) {
  cat("\nCACHE SUMMARY: ", total_failed, "/", total_samples, " samples failed (",
      sprintf("%.2f%%", 100 * total_failed / total_samples), ")\n", sep = "")
}

if (run_fresh_test && exists("total_test_failed")) {
  cat("FRESH TEST:   ", total_test_failed, "/", total_test_samples, " samples failed (",
      sprintf("%.2f%%", 100 * total_test_failed / total_test_samples), ")\n", sep = "")
}

cat("\nNote: Failed models fall back to the original (non-resampled) model.\n")
cat("This means parameter uncertainty is underestimated for those iterations.\n")
