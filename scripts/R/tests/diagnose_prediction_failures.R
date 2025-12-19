# ==============================================================================
# Diagnostic Script: PSA Prediction Failure Analysis
# ==============================================================================
# This script investigates WHY 75% of PSA iterations have prediction failures
# despite only 0.01% of bootstrap models failing to converge.
#
# Key Question: What causes predictions to return NA values when models converge?
# ==============================================================================

cat("=== PSA Prediction Failure Diagnosis ===\n\n")

# Set working directory
setwd("c:/Users/benjampg/git/METIMMOX-1")

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr)

# Source required setup scripts
cat("Loading dependencies...\n")
source(here::here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here::here("scripts/R/analysis/03_biomarker_strategies.R"))
source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))
source(here::here("scripts/R/analysis/08_sampling.R"))

cat("Dependencies loaded.\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
n_test <- 100  # Number of models to test in detail

# =============================================================================
# PART 1: INSPECT SAMPLED MODEL STRUCTURE
# =============================================================================
cat("===============================================================================\n")
cat("                    PART 1: MODEL STRUCTURE INSPECTION                         \n")
cat("===============================================================================\n\n")

# Check first sampled model
cat("Checking structure of sampled model #1...\n")
model1 <- sampling_models$control$samples[[1]]
cat("  Model 1 is NULL?", is.null(model1), "\n")
cat("  Model 1 failed?", isTRUE(model1$failed), "\n")
cat("  OS model class:", class(model1$os$model), "\n")
cat("  PFS model class:", class(model1$pfs$model), "\n")
cat("  Distribution:", model1$os$dist, "\n")

# Check coefficient structure
cat("\n  OS coefficients:\n")
print(model1$os$coefficients)

cat("\n  OS model summary:\n")
cat("    Data nobs:", model1$os$model$N, "\n")
cat("    Model distribution:", model1$os$model$dlist$name, "\n")

# =============================================================================
# PART 2: TEST PREDICTION DIRECTLY
# =============================================================================
cat("\n===============================================================================\n")
cat("                    PART 2: DIRECT PREDICTION TEST                             \n")
cat("===============================================================================\n\n")

# Test prediction with model #1 (should work)
test_model <- model1$os$model

cat("Testing predict() directly on model #1...\n")

# Test 1: Predict without newdata (uses original model data)
cat("\n  Test 1: predict() without newdata...\n")
tryCatch({
  pred1 <- predict(test_model, type = "survival", times = c(0, 52, 104, 156, 208, 260))
  cat("    Result structure: ", class(pred1), " with names:", names(pred1), "\n")
  cat("    Number of rows:", nrow(pred1), "\n")
  cat("    .pred structure:", class(pred1$.pred[[1]]), "\n")
  cat("    First survival values:", head(pred1$.pred[[1]]$.pred_survival), "\n")
  cat("    Any NA in first prediction?", any(is.na(pred1$.pred[[1]]$.pred_survival)), "\n")
}, error = function(e) {
  cat("    ERROR:", conditionMessage(e), "\n")
})

# Test 2: Predict with newdata from original data
cat("\n  Test 2: predict() with newdata from data_control...\n")
tryCatch({
  # Use first patient from original data
  newdata <- data_control[1, , drop = FALSE]
  cat("    newdata columns:", paste(names(newdata), collapse = ", "), "\n")
  cat("    newdata$sex class:", class(newdata$sex), "\n")
  cat("    newdata$Rx class:", class(newdata$Rx), "\n")

  pred2 <- predict(test_model, newdata = newdata, type = "survival", times = c(0, 52, 104))
  cat("    Prediction successful!\n")
  cat("    Survival values:", pred2$.pred[[1]]$.pred_survival, "\n")
}, error = function(e) {
  cat("    ERROR:", conditionMessage(e), "\n")
})

# Test 3: Predict with newdata from full data (used in PSA for biomarker strategies)
cat("\n  Test 3: predict() with newdata from full 'data' object...\n")
biomarker_model <- sampling_models$crp$samples[[1]]$os$model
tryCatch({
  # Subset to biomarker+ patients
  crp_pos_data <- data[data$crp == 1, ]
  cat("    Number of CRP+ patients:", nrow(crp_pos_data), "\n")
  cat("    crp_pos_data$sex class:", class(crp_pos_data$sex), "\n")
  cat("    crp_pos_data$Rx class:", class(crp_pos_data$Rx), "\n")

  # Assign experimental treatment
  exp_rx <- levels(data$Rx)[2]
  crp_pos_data$Rx <- factor(exp_rx, levels = levels(data$Rx))

  pred3 <- predict(biomarker_model, newdata = crp_pos_data, type = "survival", times = c(0, 52, 104))
  cat("    Prediction successful!\n")
  cat("    Number of predictions:", nrow(pred3), "\n")

  # Check for NA values in predictions
  all_surv <- lapply(pred3$.pred, function(x) x$.pred_survival)
  any_na <- any(sapply(all_surv, function(x) any(is.na(x))))
  cat("    Any NA in predictions?", any_na, "\n")

  if (any_na) {
    # Find which patients have NA
    na_patients <- which(sapply(all_surv, function(x) any(is.na(x))))
    cat("    Patients with NA:", head(na_patients), "\n")
  }
}, error = function(e) {
  cat("    ERROR:", conditionMessage(e), "\n")
})

# =============================================================================
# PART 3: TEST MULTIPLE MODELS FOR NA PATTERNS
# =============================================================================
cat("\n===============================================================================\n")
cat("                    PART 3: MULTI-MODEL NA ANALYSIS                            \n")
cat("===============================================================================\n\n")

cat("Testing first", n_test, "sampled models for NA patterns...\n\n")

na_analysis <- data.frame(
  model_idx = integer(),
  failed_flag = logical(),
  os_pred_error = logical(),
  os_has_na = logical(),
  pfs_pred_error = logical(),
  pfs_has_na = logical()
)

for (i in 1:n_test) {
  model <- sampling_models$control$samples[[i]]

  failed_flag <- isTRUE(model$failed)
  os_pred_error <- FALSE
  os_has_na <- FALSE
  pfs_pred_error <- FALSE
  pfs_has_na <- FALSE

  # Test OS prediction
  tryCatch({
    pred_os <- predict(model$os$model, type = "survival", times = seq(0, time_horizon))
    all_surv <- lapply(pred_os$.pred, function(x) x$.pred_survival)
    surv_matrix <- do.call(cbind, all_surv)
    mean_surv <- rowMeans(surv_matrix, na.rm = TRUE)
    os_has_na <- any(is.na(mean_surv))
  }, error = function(e) {
    os_pred_error <<- TRUE
  })

  # Test PFS prediction
  tryCatch({
    pred_pfs <- predict(model$pfs$model, type = "survival", times = seq(0, time_horizon))
    all_surv <- lapply(pred_pfs$.pred, function(x) x$.pred_survival)
    surv_matrix <- do.call(cbind, all_surv)
    mean_surv <- rowMeans(surv_matrix, na.rm = TRUE)
    pfs_has_na <- any(is.na(mean_surv))
  }, error = function(e) {
    pfs_pred_error <<- TRUE
  })

  na_analysis <- rbind(na_analysis, data.frame(
    model_idx = i,
    failed_flag = failed_flag,
    os_pred_error = os_pred_error,
    os_has_na = os_has_na,
    pfs_pred_error = pfs_pred_error,
    pfs_has_na = pfs_has_na
  ))
}

cat("Results for control strategy:\n")
cat("  Models with failed flag:", sum(na_analysis$failed_flag), "/", n_test, "\n")
cat("  OS prediction errors:", sum(na_analysis$os_pred_error), "/", n_test, "\n")
cat("  OS predictions with NA:", sum(na_analysis$os_has_na), "/", n_test, "\n")
cat("  PFS prediction errors:", sum(na_analysis$pfs_pred_error), "/", n_test, "\n")
cat("  PFS predictions with NA:", sum(na_analysis$pfs_has_na), "/", n_test, "\n")

# =============================================================================
# PART 4: TEST BIOMARKER STRATEGY PREDICTIONS
# =============================================================================
cat("\n===============================================================================\n")
cat("                    PART 4: BIOMARKER PREDICTION TEST                          \n")
cat("===============================================================================\n\n")

cat("Testing generate_psa_population_averaged_predictions()...\n\n")

# Test with first few models
biomarker_na_counts <- list()

for (biomarker in c("crp", "tlr", "tmb_braf")) {
  cat("Testing", toupper(biomarker), "strategy...\n")

  na_count <- 0
  error_count <- 0

  for (i in 1:n_test) {
    tryCatch({
      preds <- generate_psa_population_averaged_predictions(
        sampling_model_list = sampling_models[[biomarker]],
        biomarker_name = biomarker,
        outcome = "os",
        sample_idx = i,
        data_original = data,
        time_points = seq(0, time_horizon)
      )

      if (is.null(preds) ||
          any(is.na(preds$positive)) ||
          any(is.na(preds$negative))) {
        na_count <- na_count + 1
      }
    }, error = function(e) {
      error_count <<- error_count + 1
    })
  }

  cat("  Predictions with NULL/NA:", na_count, "/", n_test,
      "(", sprintf("%.1f%%", 100 * na_count / n_test), ")\n")
  cat("  Prediction errors:", error_count, "/", n_test,
      "(", sprintf("%.1f%%", 100 * error_count / n_test), ")\n")

  biomarker_na_counts[[biomarker]] <- list(na = na_count, error = error_count)
}

# =============================================================================
# PART 5: DETAILED ERROR MESSAGE CAPTURE
# =============================================================================
cat("\n===============================================================================\n")
cat("                    PART 5: DETAILED ERROR CAPTURE                             \n")
cat("===============================================================================\n\n")

cat("Capturing detailed error messages for first 5 failures...\n\n")

error_messages <- character()
idx <- 1

while (length(error_messages) < 5 && idx <= min(n_test * 5, n_samples)) {
  tryCatch({
    suppressWarnings({
      preds <- generate_psa_population_averaged_predictions(
        sampling_model_list = sampling_models$crp,
        biomarker_name = "crp",
        outcome = "os",
        sample_idx = idx,
        data_original = data,
        time_points = seq(0, time_horizon)
      )

      if (is.null(preds) || any(is.na(preds$positive)) || any(is.na(preds$negative))) {
        # Try to get more detail
        model <- sampling_models$crp$samples[[idx]]$os$model
        crp_pos <- data[data$crp == 1, ]
        crp_pos$Rx <- factor(levels(data$Rx)[2], levels = levels(data$Rx))

        pred_detail <- predict(model, newdata = crp_pos, type = "survival", times = c(0, 52))
        all_surv <- lapply(pred_detail$.pred, function(x) x$.pred_survival)

        na_per_patient <- sapply(all_surv, function(x) sum(is.na(x)))
        if (any(na_per_patient > 0)) {
          error_messages <- c(error_messages,
                              paste0("Model ", idx, ": NA in ", sum(na_per_patient > 0), " patients"))
        }
      }
    })
  }, error = function(e) {
    error_messages <<- c(error_messages, paste0("Model ", idx, ": ", conditionMessage(e)))
  })

  idx <- idx + 1
}

if (length(error_messages) > 0) {
  cat("Error/NA messages:\n")
  for (msg in error_messages) {
    cat("  ", msg, "\n")
  }
} else {
  cat("No errors captured in first", idx - 1, "models\n")
}

# =============================================================================
# PART 6: TEST model_fun() DIRECTLY (matching PSA error rate test)
# =============================================================================
cat("\n===============================================================================\n")
cat("                    PART 6: MODEL_FUN() DIRECT TEST                            \n")
cat("===============================================================================\n\n")

# Source model functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))

cat("Testing model_fun() with PSA mode for", n_test, "iterations...\n\n")

fallback_count <- 0
error_count <- 0
warning_details <- list()

for (i in 1:n_test) {
  tryCatch({
    withCallingHandlers({
      sim_results <- model_fun(
        params = l_params_base,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE,
        sim_idx = i
      )

      # Check fallback attribute
      if (isTRUE(attr(sim_results, "fallback_used"))) {
        fallback_count <- fallback_count + 1
        warning_details[[length(warning_details) + 1]] <- list(
          idx = i,
          type = "fallback"
        )
      }

    }, warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("NA values|Failed to generate|NA/NULL values", msg)) {
        if (length(warning_details) < 10) {
          warning_details[[length(warning_details) + 1]] <<- list(
            idx = i,
            message = msg
          )
        }
      }
      invokeRestart("muffleWarning")
    })
  }, error = function(e) {
    error_count <<- error_count + 1
  })

  if (i %% 25 == 0) {
    cat(sprintf("  Progress: %d/%d - Fallbacks: %d, Errors: %d\n",
                i, n_test, fallback_count, error_count))
  }
}

cat("\nmodel_fun() RESULTS:\n")
cat("  Fallback used (predictions failed):", fallback_count, "/", n_test,
    "(", sprintf("%.1f%%", 100 * fallback_count / n_test), ")\n")
cat("  Complete errors:", error_count, "/", n_test,
    "(", sprintf("%.1f%%", 100 * error_count / n_test), ")\n")

if (length(warning_details) > 0) {
  cat("\nSample warning messages (first 5):\n")
  for (j in 1:min(5, length(warning_details))) {
    detail <- warning_details[[j]]
    if (!is.null(detail$message)) {
      cat("  Iteration", detail$idx, ":", substr(detail$message, 1, 80), "\n")
    }
  }
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n===============================================================================\n")
cat("                           DIAGNOSIS COMPLETE                                  \n")
cat("===============================================================================\n\n")

cat("SUMMARY OF FINDINGS:\n")
cat("-------------------\n")
cat("1. Control strategy NA (generate_sampling_predictions):",
    sum(na_analysis$os_has_na), "/", n_test, "\n")
cat("2. Biomarker CRP NA (generate_psa_population_averaged_predictions):",
    biomarker_na_counts$crp$na, "/", n_test, "\n")
cat("3. model_fun() fallbacks:", fallback_count, "/", n_test, "\n")

if (fallback_count > 0 && biomarker_na_counts$crp$na == 0) {
  cat("\nKEY INSIGHT: Fallbacks in model_fun() but not in individual prediction tests\n")
  cat("This suggests the issue is in how model_fun checks for NA values.\n")
}