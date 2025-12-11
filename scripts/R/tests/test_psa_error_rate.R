# ==============================================================================
# PSA Error Rate Analysis Script
# ==============================================================================
# This script quantifies how many PSA iterations lead to errors (prediction
# failures, model convergence issues, etc.) and outputs error counts with
# percentages.
#
# USAGE: Adjust n_test_sim below, then run the entire script.
# ==============================================================================

# =============================================================================
# CONFIGURABLE PARAMETER - Adjust this value as needed
# =============================================================================
n_test_sim <- 100  # Number of PSA iterations to test

# =============================================================================
# SETUP
# =============================================================================
cat("=== PSA Error Rate Analysis ===\n")
cat("Testing", n_test_sim, "PSA iterations for errors\n\n")

# Set working directory
setwd("c:/Users/benjampg/git/METIMMOX-1")

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr)

# Source required analysis scripts
cat("Loading dependencies...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))

# Override n_samples for testing (to ensure sampling cache is available)
# We need the cache but will only test n_test_sim iterations
source(here::here("scripts/R/analysis/08_sampling.R"))

# Source model functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/prediction_functions.R"))

cat("Dependencies loaded successfully.\n\n")

# =============================================================================
# ERROR TRACKING INITIALIZATION
# =============================================================================
error_tracking <- list(
  # Prediction failures (model fell back to base case curves)
  control_failures = integer(0),      # Control strategy prediction failures
  crp_failures = integer(0),          # CRP biomarker prediction failures
  tlr_failures = integer(0),          # TLR biomarker prediction failures
  tmb_braf_failures = integer(0),     # TMB/BRAF biomarker prediction failures

  # Complete iteration failures (model_fun crashed)
  complete_failures = integer(0),

  # PFS > OS constraint violations (not errors, but tracked for context)
  pfs_os_violations = list(
    control = integer(0),
    crp_pos = integer(0),
    crp_neg = integer(0),
    tlr_pos = integer(0),
    tlr_neg = integer(0),
    tmb_braf_pos = integer(0),
    tmb_braf_neg = integer(0)
  )
)

# =============================================================================
# PSA ERROR TESTING LOOP
# =============================================================================
cat("Running PSA iterations with error tracking...\n")
start_time <- Sys.time()

for (i in 1:n_test_sim) {
  # Create parameter set for this simulation (use base params)
  sim_params <- l_params_base

  # Run model_fun with error and warning capture
  tryCatch({
    withCallingHandlers({
      sim_results <- model_fun(
        params = sim_params,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE,
        sim_idx = i
      )
    }, warning = function(w) {
      msg <- conditionMessage(w)

      # Track prediction failures (fallback to base case)
      if (grepl("Failed to generate control curves", msg) ||
          grepl("NA values in control survival", msg)) {
        error_tracking$control_failures <<- c(error_tracking$control_failures, i)
      }

      if (grepl("Failed to generate crp curves", msg) ||
          grepl("NA/NULL values in crp survival", msg) ||
          grepl("NA/NULL values in crp", msg)) {
        error_tracking$crp_failures <<- c(error_tracking$crp_failures, i)
      }

      if (grepl("Failed to generate tlr curves", msg) ||
          grepl("NA/NULL values in tlr survival", msg) ||
          grepl("NA/NULL values in tlr", msg)) {
        error_tracking$tlr_failures <<- c(error_tracking$tlr_failures, i)
      }

      if (grepl("Failed to generate tmb_braf curves", msg) ||
          grepl("NA/NULL values in tmb_braf survival", msg) ||
          grepl("NA/NULL values in tmb_braf", msg)) {
        error_tracking$tmb_braf_failures <<- c(error_tracking$tmb_braf_failures, i)
      }

      # Track PFS > OS violations
      if (grepl("PFS > OS constraint enforced", msg)) {
        if (grepl("\\(control\\)", msg)) {
          error_tracking$pfs_os_violations$control <<-
            c(error_tracking$pfs_os_violations$control, i)
        } else if (grepl("\\(crp\\+\\)", msg)) {
          error_tracking$pfs_os_violations$crp_pos <<-
            c(error_tracking$pfs_os_violations$crp_pos, i)
        } else if (grepl("\\(crp-\\)", msg)) {
          error_tracking$pfs_os_violations$crp_neg <<-
            c(error_tracking$pfs_os_violations$crp_neg, i)
        } else if (grepl("\\(tlr\\+\\)", msg)) {
          error_tracking$pfs_os_violations$tlr_pos <<-
            c(error_tracking$pfs_os_violations$tlr_pos, i)
        } else if (grepl("\\(tlr-\\)", msg)) {
          error_tracking$pfs_os_violations$tlr_neg <<-
            c(error_tracking$pfs_os_violations$tlr_neg, i)
        } else if (grepl("\\(tmb_braf\\+\\)", msg)) {
          error_tracking$pfs_os_violations$tmb_braf_pos <<-
            c(error_tracking$pfs_os_violations$tmb_braf_pos, i)
        } else if (grepl("\\(tmb_braf-\\)", msg)) {
          error_tracking$pfs_os_violations$tmb_braf_neg <<-
            c(error_tracking$pfs_os_violations$tmb_braf_neg, i)
        }
        # Suppress PFS > OS warnings (expected behavior)
        invokeRestart("muffleWarning")
      }

      # Let other warnings through but don't duplicate tracking
    })

  }, error = function(e) {
    # Complete iteration failure
    error_tracking$complete_failures <<- c(error_tracking$complete_failures, i)
    cat("  ERROR in iteration", i, ":", conditionMessage(e), "\n")
  })

  # Progress reporting
  if (i %% 10 == 0 || i == n_test_sim) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate <- i / elapsed
    remaining <- (n_test_sim - i) / rate
    cat(sprintf("  Progress: %d/%d (%.0f%%) - Elapsed: %.1fs - ETA: %.1fs\n",
                i, n_test_sim, 100 * i / n_test_sim, elapsed, remaining))
  }
}

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

# =============================================================================
# RESULTS SUMMARY
# =============================================================================
cat("\n")
cat("===============================================================================\n")
cat("                        PSA ERROR RATE ANALYSIS RESULTS                        \n")
cat("===============================================================================\n")
cat(sprintf("Simulations tested: %d\n", n_test_sim))
cat(sprintf("Total runtime: %.1f seconds\n", total_time))
cat("\n")

# Deduplicate iteration lists (in case same iteration had multiple warnings)
control_failures <- unique(error_tracking$control_failures)
crp_failures <- unique(error_tracking$crp_failures)
tlr_failures <- unique(error_tracking$tlr_failures)
tmb_braf_failures <- unique(error_tracking$tmb_braf_failures)
complete_failures <- unique(error_tracking$complete_failures)

# Calculate unique iterations with any prediction failure
all_prediction_failures <- unique(c(control_failures, crp_failures,
                                     tlr_failures, tmb_braf_failures))

cat("PREDICTION FAILURES (model fell back to base case curves):\n")
cat(sprintf("  Control strategy:     %3d iterations (%5.1f%%)\n",
            length(control_failures), 100 * length(control_failures) / n_test_sim))
cat(sprintf("  CRP strategy:         %3d iterations (%5.1f%%)\n",
            length(crp_failures), 100 * length(crp_failures) / n_test_sim))
cat(sprintf("  TLR strategy:         %3d iterations (%5.1f%%)\n",
            length(tlr_failures), 100 * length(tlr_failures) / n_test_sim))
cat(sprintf("  TMB/BRAF strategy:    %3d iterations (%5.1f%%)\n",
            length(tmb_braf_failures), 100 * length(tmb_braf_failures) / n_test_sim))
cat(sprintf("  Any prediction fail:  %3d iterations (%5.1f%%)\n",
            length(all_prediction_failures), 100 * length(all_prediction_failures) / n_test_sim))
cat("\n")

cat("COMPLETE FAILURES (model_fun crashed):\n")
cat(sprintf("  Total crashes:        %3d iterations (%5.1f%%)\n",
            length(complete_failures), 100 * length(complete_failures) / n_test_sim))
if (length(complete_failures) > 0) {
  cat(sprintf("  Failed iterations: %s\n", paste(complete_failures, collapse = ", ")))
}
cat("\n")

# Overall error summary
all_failures <- unique(c(all_prediction_failures, complete_failures))
cat("OVERALL ERROR SUMMARY:\n")
cat(sprintf("  Iterations with any error: %3d of %d (%5.1f%%)\n",
            length(all_failures), n_test_sim, 100 * length(all_failures) / n_test_sim))
cat(sprintf("  Clean iterations:          %3d of %d (%5.1f%%)\n",
            n_test_sim - length(all_failures), n_test_sim,
            100 * (n_test_sim - length(all_failures)) / n_test_sim))
cat("\n")

# PFS > OS violations (not errors, but useful context)
all_violations <- unique(unlist(error_tracking$pfs_os_violations))
if (length(all_violations) > 0) {
  cat("PFS > OS CONSTRAINT VIOLATIONS (not errors, corrected automatically):\n")
  cat(sprintf("  Iterations affected:  %3d of %d (%5.1f%%)\n",
              length(all_violations), n_test_sim, 100 * length(all_violations) / n_test_sim))
  cat("  Breakdown by subgroup:\n")

  subgroup_names <- c(
    control = "    Control",
    crp_pos = "    CRP+",
    crp_neg = "    CRP-",
    tlr_pos = "    TLR+",
    tlr_neg = "    TLR-",
    tmb_braf_pos = "    TMB/BRAF+",
    tmb_braf_neg = "    TMB/BRAF-"
  )

  for (subgroup in names(error_tracking$pfs_os_violations)) {
    n_violations <- length(unique(error_tracking$pfs_os_violations[[subgroup]]))
    if (n_violations > 0) {
      cat(sprintf("%s: %d iterations (%.1f%%)\n",
                  subgroup_names[subgroup], n_violations,
                  100 * n_violations / n_test_sim))
    }
  }
  cat("  (PFS capped at OS to ensure valid state occupancy - expected behavior)\n")
}

cat("\n===============================================================================\n")
cat("                               ANALYSIS COMPLETE                               \n")
cat("===============================================================================\n")

# Return results object for programmatic use
psa_error_results <- list(
  n_sim = n_test_sim,
  runtime_seconds = total_time,
  prediction_failures = list(
    control = control_failures,
    crp = crp_failures,
    tlr = tlr_failures,
    tmb_braf = tmb_braf_failures,
    any = all_prediction_failures
  ),
  complete_failures = complete_failures,
  all_failures = all_failures,
  pfs_os_violations = error_tracking$pfs_os_violations,
  error_rate = 100 * length(all_failures) / n_test_sim
)

cat("\nResults object 'psa_error_results' available for further analysis.\n")
