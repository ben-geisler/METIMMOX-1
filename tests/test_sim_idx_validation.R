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

source("R/model_configs.R")
source("R/model_fun.R")
source("R/calculate_outcomes.R")
source("R/prediction_functions.R")

n_samples <- 2L
time_horizon <- 52L
cl <- 1 / 52
curve <- exp(-seq(0, time_horizon) / 100)
l_params_base <- list(
  dr_costs = 0.04, dr_effects = 0.04, u_np = 0.73, u_p = 0.59,
  c_drug_nivo = 100, c_drug_FLOX = 10, c_test_CT = 1, c_test_blood = 1,
  c_test_CRP = 1, c_test_NGS = 1, c_test_biomarker = list(crp = 1, tmb_braf = 1),
  c_other_visit = 1, c_other_baseline = 1, c_other_follow = 1, c_other_last = 1,
  p_crp = 0.5, p_tmb_braf = 0.5)
for (name in c("l_nivo", "l_FLOX_exp", "l_FLOX_control", "l_CT", "l_blood", "l_visit"))
  l_params_base[[name]] <- rep(0, time_horizon + 1L)
for (outcome in c("OS", "PFS")) {
  keys <- paste0(c("control", "crp_pos", "crp_neg", "tmb_braf_pos", "tmb_braf_neg"), "_", outcome)
  l_params_base[[paste0("p_", tolower(outcome))]] <- setNames(rep(list(curve), length(keys)), keys)
}
data <- data_complete <- data.frame(ID = 1:2, crp = 0:1, tmb_braf = 0:1)
l_params_base$prediction_population <- data_complete
l_params_base$population_weights <- c(0.5, 0.5)
sampling_models <- list(joint = list(samples = rep(list(list(failed = FALSE)), n_samples)))
# Index validation is the unit under test; deterministic prediction fixture
# exercises the remaining calculation without fitting or accessing any cache.
generate_psa_population_averaged_predictions <- function(sampling_model_list,
    biomarker_name = NULL, outcome, sample_idx, data_original, time_points, weights = NULL) {
  stopifnot(sample_idx %in% seq_len(n_samples))
  if (is.null(biomarker_name)) curve else list(positive = curve, negative = curve)
}

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
      output <- model_fun(
        params = l_params_base,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE,
        sim_idx = sim_idx_value
      )
      stopifnot(!isTRUE(attr(output, "fallback_used")),
                all(is.finite(output$Cost)), all(is.finite(output$Effect)))
    })
    # If we get here, no error was thrown
    list(status = "success", error = NULL)
  }, error = function(e) {
    list(status = "error", error = conditionMessage(e))
  })

  # Evaluate pass/fail
  if (expect_error) {
    if (result$status == "error" && grepl("sim_idx", result$error, fixed = TRUE)) {
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
    cat("FAIL - Expected the documented NULL behaviour\n")
    cat(sprintf("         Error: %s\n", result$error))
    return(FALSE)  # The documented NULL behaviour is an asserted contract.
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

stopifnot(all(results))
