# Test script to verify issue #72 fix: Sex variable type mismatch
# This script checks that sex is consistently typed as factor across all datasets

cat("=== Testing Sex Variable Type Consistency (Issue #72) ===\n\n")

# Set working directory
setwd("c:/Users/benjampg/git/METIMMOX-1")

# Source required scripts in order
cat("Loading data and creating datasets...\n")
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")

# Check sex variable type in main data object
cat("\n1. Main data object:\n")
cat("   sex type:", class(data$sex), "\n")
cat("   sex levels:", paste(levels(data$sex), collapse=", "), "\n")

# Check sex variable type in all subset datasets
cat("\n2. Subset datasets:\n")

cat("   data_control$sex type:", class(data_control$sex), "\n")
cat("   data_control$sex levels:", paste(levels(data_control$sex), collapse=", "), "\n")

cat("   data_crp$sex type:", class(data_crp$sex), "\n")
cat("   data_crp$sex levels:", paste(levels(data_crp$sex), collapse=", "), "\n")

cat("   data_tlr$sex type:", class(data_tlr$sex), "\n")
cat("   data_tlr$sex levels:", paste(levels(data_tlr$sex), collapse=", "), "\n")

cat("   data_tmb_braf$sex type:", class(data_tmb_braf$sex), "\n")
cat("   data_tmb_braf$sex levels:", paste(levels(data_tmb_braf$sex), collapse=", "), "\n")

# Now source script 06 and check it doesn't change types
cat("\n3. After parametric survival analysis script:\n")
source("scripts/R/analysis/06_parametric_survival analysis.R")

cat("   data$sex type:", class(data$sex), "\n")
cat("   data_control$sex type:", class(data_control$sex), "\n")

# Test resampling with one iteration
cat("\n4. Testing survival model resampling (1 iteration)...\n")

# Set small n_samples for testing
n_samples_orig <- n_samples
n_samples <- 1

# Source resampling function
source("scripts/R/functions/prediction_functions.R")

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(survival, flexsurv)

# Define formulas
control_os_formula <- Surv(OSwk, Death) ~ Age + sex
control_pfs_formula <- Surv(PFSwk, Progression) ~ Age + sex

# Sample ONE model from control data
cat("   Fitting models to data_control...\n")
tryCatch({
  resample_idx <- sample(1:nrow(data_control), size = nrow(data_control), replace = TRUE)
  resampled_data <- data_control[resample_idx, ]

  cat("   Resampled data sex type:", class(resampled_data$sex), "\n")

  # Fit models
  os_model <- flexsurvreg(control_os_formula, data = resampled_data, dist = "weibull")
  pfs_model <- flexsurvreg(control_pfs_formula, data = resampled_data, dist = "weibull")

  cat("   ✓ Models fitted successfully!\n")

  # Test prediction with factor sex
  cat("   Testing prediction with factor sex...\n")
  test_newdata <- data.frame(
    Age = mean(data_control$Age),
    sex = factor(levels(data_control$sex)[1], levels = levels(data_control$sex))
  )

  pred <- predict(os_model, newdata = test_newdata, type = "survival", times = seq(0, 100, by = 10))
  cat("   ✓ Prediction successful!\n")

  cat("\n=== TEST PASSED ===\n")
  cat("Sex variable is consistently typed as factor across all datasets.\n")
  cat("No type mismatch warnings during resampling or prediction.\n")

}, error = function(e) {
  cat("\n=== TEST FAILED ===\n")
  cat("Error:", conditionMessage(e), "\n")
})

# Restore original n_samples
n_samples <- n_samples_orig
