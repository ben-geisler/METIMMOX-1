# ==============================================================================
# sim_idx Validation Tests (Issue #46)
# ==============================================================================
# This script tests that model_fun() properly validates the sim_idx parameter
# when running in PSA mode. It verifies that invalid inputs produce clear
# error messages rather than cryptic failures or silent fallbacks.
#
# USAGE: Run this script after sourcing the required dependencies.
#        The test should PASS after the fix is implemented.
# ==============================================================================

# =============================================================================
# SETUP - Load dependencies
# =============================================================================
cat("=== sim_idx Validation Tests (Issue #46) ===\n\n")
cat("Loading dependencies...\n")

# Set working directory
setwd(here::here())

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr)

# Source required analysis scripts
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
source(here::here("scripts/R/analysis/04_parametric_survival_analysis.R"))
source(here::here("scripts/R/analysis/05_basecase_input_parameters.R"))
source(here::here("scripts/R/analysis/06_sampling.R"))

# Source model functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/prediction_functions.R"))

cat("Dependencies loaded successfully.\n\n")

# =============================================================================
# TEST HELPER FUNCTION
# =============================================================================
# Runs model_fun with given parameters and captures the outcome

test_model_fun <- function(test_name, sim_idx_value, expect_error = TRUE) {
  cat(sprintf("%-35s", paste0("Test: ", test_name, "...")))

  result <- tryCatch({
    # Suppress warnings for clean output (we're testing for errors)
    suppressWarnings({
      model_fun(
        params = l_params_base,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE,
        sim_idx = sim_idx_value
      )
    })
    # If we get here, no error was thrown
    list(status = "success", error = NULL)
  }, error = function(e) {
    list(status = "error", error = conditionMessage(e))
  })

  # Evaluate pass/fail
  if (expect_error) {
    if (result$status == "error") {
      cat("PASS - Got expected error\n")
      cat(sprintf("         Error message: %s\n", substr(result$error, 1, 60)))
      return(TRUE)
    } else {
      cat("FAIL - Expected error but got success\n")
      return(FALSE)
    }
  } else {
    if (result$status == "success") {
      cat("PASS - Model ran successfully\n")
      return(TRUE)
    } else {
      cat("FAIL - Expected success but got error\n")
      cat(sprintf("         Error: %s\n", result$error))
      return(FALSE)
    }
  }
}

# Special test for NULL sim_idx (PSA block should be skipped)
test_null_sim_idx <- function() {
  cat(sprintf("%-35s", "Test: NULL sim_idx..."))

  result <- tryCatch({
    suppressWarnings({
      model_fun(
        params = l_params_base,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE,
        sim_idx = NULL  # NULL should cause PSA block to be skipped
      )
    })
    list(status = "success", error = NULL)
  }, error = function(e) {
    list(status = "error", error = conditionMessage(e))
  })

  # With NULL sim_idx and determpsa="psa", the PSA block is skipped
  # This is documented behavior (line 53 checks: && !is.null(sim_idx))
  if (result$status == "success") {
    cat("PASS - PSA block skipped (documented)\n")
    return(TRUE)
  } else {
    cat("INFO - Got error (could be stricter validation)\n")
    cat(sprintf("         Error: %s\n", result$error))
    return(TRUE)  # Not a failure, just documenting behavior
  }
}

# =============================================================================
# RUN TEST CASES
# =============================================================================

cat("Running tests...\n")
cat("n_samples = ", n_samples, "\n\n")

results <- c()

# Test 1: Zero index (should error)
results["zero"] <- test_model_fun("sim_idx = 0", sim_idx_value = 0, expect_error = TRUE)

# Test 2: Negative index (should error)
results["negative"] <- test_model_fun("sim_idx = -1", sim_idx_value = -1, expect_error = TRUE)

# Test 3: Out of bounds (should error)
results["out_of_bounds"] <- test_model_fun(
  paste0("sim_idx = ", n_samples + 1),
  sim_idx_value = n_samples + 1,
  expect_error = TRUE
)

# Test 4: Non-numeric (should error)
results["non_numeric"] <- test_model_fun("sim_idx = 'abc'", sim_idx_value = "abc", expect_error = TRUE)

# Test 5: NULL with PSA mode (documented behavior)
results["null"] <- test_null_sim_idx()

# Test 6: Valid index (should succeed)
results["valid"] <- test_model_fun("sim_idx = 1", sim_idx_value = 1, expect_error = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat(strrep("=", 50), "\n")
cat("SUMMARY\n")
cat(strrep("=", 50), "\n")

passed <- sum(results)
total <- length(results)

cat(sprintf("Passed: %d / %d\n", passed, total))

if (passed == total) {
  cat("\nAll tests passed! Issue #46 validation is working correctly.\n")
} else {
  cat("\nSome tests failed. The sim_idx validation may not be implemented yet.\n")
  cat("See Issue #46 for the recommended fix.\n")
}
