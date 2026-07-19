# ===============================================================================
# SENSITIVITY ANALYSIS: SEPARATE BIOMARKER-SPECIFIC SURVIVAL MODELS
# ===============================================================================
# This script fits separate survival models for each biomarker strategy,
# where each model only includes the interaction term for that biomarker.
#
# Purpose: Compare clinical outcomes between:
#   - All-in-One model: ~ Age + sex + Rx + crp*Rx + tlr*Rx + tmb_braf*Rx
#   - Separate models: ~ Age + sex + Rx + [biomarker]:Rx (one per biomarker)
#
# All models use the gamma distribution (selected by AIC in the primary analysis).
#
# Issue #90: https://github.com/ben-geisler/METIMMOX-1/issues/90
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr, tidyr)

# ===============================================================================
# SETUP
# ===============================================================================

# Source required analysis scripts if not already loaded
if (!exists("data_complete") || !exists("strategies_df")) {
  cat("Loading required analysis scripts...\n")
  source(here("scripts/R/analysis/02_setup_and_global_variables.R"))
  source(here("scripts/R/analysis/03_biomarker_strategies.R"))
}

# Source prediction functions
source(here("scripts/R/functions/prediction_functions.R"))

# Source script 06 to get the full models and best-fit distribution
if (!exists("models") || is.null(models$best_fit)) {
  cat("Loading parametric survival analysis...\n")
  source(here("scripts/R/analysis/06_parametric_survival_analysis.R"))
}

# Get best-fit distribution from main analysis
best_os_dist <- models$best_fit$os_distribution
best_pfs_dist <- models$best_fit$pfs_distribution
cat("Using distributions:", best_os_dist, "(OS),", best_pfs_dist, "(PFS)\n")

# ===============================================================================
# DEFINE SEPARATE MODEL FORMULAS
# ===============================================================================

# Each model includes only ONE biomarker interaction term
separate_formulas <- list(
  crp = list(
    os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx,
    pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx
  ),
  tlr = list(
    os = Surv(OSwk, Death) ~ Age + sex + Rx + tlr:Rx,
    pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tlr:Rx
  ),
  tmb_braf = list(
    os = Surv(OSwk, Death) ~ Age + sex + Rx + tmb_braf:Rx,
    pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tmb_braf:Rx
  )
)

# ===============================================================================
# FIT SEPARATE MODELS
# ===============================================================================

cat("\n========================================\n")
cat("Fitting separate biomarker models\n")
cat("========================================\n\n")

separate_models <- list()

for (biomarker in names(separate_formulas)) {
  cat("Fitting models for", toupper(biomarker), "strategy...\n")

  # Fit OS model
  os_model <- flexsurvreg(
    formula = separate_formulas[[biomarker]]$os,
    data = data_complete,
    dist = best_os_dist
  )

  # Fit PFS model
  pfs_model <- flexsurvreg(
    formula = separate_formulas[[biomarker]]$pfs,
    data = data_complete,
    dist = best_pfs_dist
  )

  separate_models[[biomarker]] <- list(
    os = os_model,
    pfs = pfs_model
  )

  cat("  - OS model AIC:", round(AIC(os_model), 2), "\n")
  cat("  - PFS model AIC:", round(AIC(pfs_model), 2), "\n")
}

cat("\nSeparate model fitting complete.\n")

# ===============================================================================
# GENERATE PREDICTIONS FROM SEPARATE MODELS
# ===============================================================================

cat("\n========================================\n")
cat("Generating predictions from separate models\n")
cat("========================================\n\n")

# Extract treatment level references
exp_rx <- levels(data_complete$Rx)[2]  # Experimental treatment
ctrl_rx <- levels(data_complete$Rx)[1]  # Control treatment

separate_predictions <- list()

for (biomarker in names(separate_models)) {
  cat("Generating predictions for", toupper(biomarker), "strategy...\n")

  # Get biomarker-positive patients
  # Use explicit type conversion for robustness (works with factor, character, or numeric)
  biomarker_pos_data <- data_complete[
    as.numeric(as.character(data_complete[[biomarker]])) == 1, ]
  biomarker_pos_data$Rx <- factor(exp_rx, levels = levels(data_complete$Rx))

  # Get biomarker-negative patients
  # Use explicit type conversion for robustness (works with factor, character, or numeric)
  biomarker_neg_data <- data_complete[
    as.numeric(as.character(data_complete[[biomarker]])) == 0, ]
  biomarker_neg_data$Rx <- factor(ctrl_rx, levels = levels(data_complete$Rx))

  # Predict OS for biomarker-positive (experimental treatment)
  os_pred_pos <- predict(separate_models[[biomarker]]$os,
                         newdata = biomarker_pos_data,
                         type = "survival", times = time_points)
  os_matrix_pos <- extract_all_survival_probabilities(os_pred_pos)
  biomarker_pos_os <- rowMeans(os_matrix_pos, na.rm = TRUE)

  # Predict PFS for biomarker-positive (experimental treatment)
  pfs_pred_pos <- predict(separate_models[[biomarker]]$pfs,
                          newdata = biomarker_pos_data,
                          type = "survival", times = time_points)
  pfs_matrix_pos <- extract_all_survival_probabilities(pfs_pred_pos)
  biomarker_pos_pfs <- rowMeans(pfs_matrix_pos, na.rm = TRUE)

  # Predict OS for biomarker-negative (control treatment)
  os_pred_neg <- predict(separate_models[[biomarker]]$os,
                         newdata = biomarker_neg_data,
                         type = "survival", times = time_points)
  os_matrix_neg <- extract_all_survival_probabilities(os_pred_neg)
  biomarker_neg_os <- rowMeans(os_matrix_neg, na.rm = TRUE)

  # Predict PFS for biomarker-negative (control treatment)
  pfs_pred_neg <- predict(separate_models[[biomarker]]$pfs,
                          newdata = biomarker_neg_data,
                          type = "survival", times = time_points)
  pfs_matrix_neg <- extract_all_survival_probabilities(pfs_pred_neg)
  biomarker_neg_pfs <- rowMeans(pfs_matrix_neg, na.rm = TRUE)

  # Calculate prevalence and weighted average
  prevalence <- mean(as.numeric(as.character(data_complete[[biomarker]])), na.rm = TRUE)
  weighted_os <- prevalence * biomarker_pos_os + (1 - prevalence) * biomarker_neg_os
  weighted_pfs <- prevalence * biomarker_pos_pfs + (1 - prevalence) * biomarker_neg_pfs

  separate_predictions[[biomarker]] <- list(
    strategy = biomarker,
    os = weighted_os,
    pfs = weighted_pfs,
    biomarker_positive = list(
      os = biomarker_pos_os,
      pfs = biomarker_pos_pfs
    ),
    biomarker_negative = list(
      os = biomarker_neg_os,
      pfs = biomarker_neg_pfs
    ),
    prevalence = prevalence
  )

  cat("  - Prevalence:", round(prevalence * 100, 1), "%\n")
  cat("  - Biomarker+ patients:", nrow(biomarker_pos_data), "\n")
  cat("  - Biomarker- patients:", nrow(biomarker_neg_data), "\n")
}

cat("\nPrediction generation complete.\n")

# ===============================================================================
# GENERATE CONTROL PREDICTIONS (for reference)
# ===============================================================================

cat("\n========================================\n")
cat("Generating control strategy predictions\n")
cat("========================================\n\n")

# Control: all patients receive control treatment
# Use predictions from script 06 if available, otherwise generate here
if (exists("predictions") && !is.null(predictions$control)) {
  control_os <- predictions$control$os
  control_pfs <- predictions$control$pfs
  cat("Using control predictions from main analysis.\n")
} else {
  # Generate control predictions using the all-in-one model
  control_data <- data_complete
  control_data$Rx <- factor(ctrl_rx, levels = levels(data_complete$Rx))

  os_pred_ctrl <- predict(models$best_fit$os, newdata = control_data,
                          type = "survival", times = time_points)
  pfs_pred_ctrl <- predict(models$best_fit$pfs, newdata = control_data,
                           type = "survival", times = time_points)

  os_matrix_ctrl <- extract_all_survival_probabilities(os_pred_ctrl)
  pfs_matrix_ctrl <- extract_all_survival_probabilities(pfs_pred_ctrl)

  control_os <- rowMeans(os_matrix_ctrl, na.rm = TRUE)
  control_pfs <- rowMeans(pfs_matrix_ctrl, na.rm = TRUE)
  cat("Generated control predictions from all-in-one model.\n")
}

cat("Control predictions available.\n")

# ===============================================================================
# CREATE COEFFICIENT COMPARISON TABLE
# ===============================================================================

cat("\n========================================\n")
cat("Extracting model coefficients\n")
cat("========================================\n\n")

# Extract coefficients from all-in-one model (from script 06)
allinone_os_model <- models$best_fit$os
allinone_pfs_model <- models$best_fit$pfs

# Function to extract interaction coefficient from model
# Uses flexible pattern matching to find biomarker:Rx interaction terms
extract_interaction_coef <- function(model, biomarker_name) {
  coefs <- as.data.frame(model$res.t)
  coef_names <- rownames(coefs)

  # Debug: print available coefficient names (comment out after testing)
  # cat("  Available coefficients:", paste(coef_names, collapse = ", "), "\n")

  # Flexible pattern: look for biomarker name followed by :Rx (case-insensitive)
  # Pattern matches: "crp1:RxExp", "crp:RxExp", "RxExp:crp1", etc.
  pattern1 <- paste0(biomarker_name, ".*:Rx")  # biomarker....:Rx...
  pattern2 <- paste0("Rx.*:", biomarker_name)  # Rx....:biomarker...

  matching_rows <- grep(pattern1, coef_names, value = TRUE, ignore.case = TRUE)
  if (length(matching_rows) == 0) {
    matching_rows <- grep(pattern2, coef_names, value = TRUE, ignore.case = TRUE)
  }

  if (length(matching_rows) > 0) {
    row_name <- matching_rows[1]
    return(list(
      estimate = coefs[row_name, "est"],
      se = coefs[row_name, "se"],
      hr = exp(coefs[row_name, "est"]),
      hr_lower = exp(coefs[row_name, "est"] - 1.96 * coefs[row_name, "se"]),
      hr_upper = exp(coefs[row_name, "est"] + 1.96 * coefs[row_name, "se"])
    ))
  } else {
    warning(paste("No interaction coefficient found for", biomarker_name))
    return(list(estimate = NA, se = NA, hr = NA, hr_lower = NA, hr_upper = NA))
  }
}

# Build comparison table
coef_comparison <- data.frame(
  Biomarker = character(),
  Outcome = character(),
  Model_Type = character(),
  HR = numeric(),
  HR_Lower = numeric(),
  HR_Upper = numeric(),
  stringsAsFactors = FALSE
)

for (biomarker in c("crp", "tlr", "tmb_braf")) {
  # All-in-One model coefficients
  allinone_os_coef <- extract_interaction_coef(allinone_os_model, biomarker)
  allinone_pfs_coef <- extract_interaction_coef(allinone_pfs_model, biomarker)

  # Separate model coefficients
  separate_os_coef <- extract_interaction_coef(separate_models[[biomarker]]$os, biomarker)
  separate_pfs_coef <- extract_interaction_coef(separate_models[[biomarker]]$pfs, biomarker)

  # Add rows to comparison table
  coef_comparison <- rbind(coef_comparison, data.frame(
    Biomarker = toupper(biomarker),
    Outcome = "OS",
    Model_Type = "All-in-One",
    HR = allinone_os_coef$hr,
    HR_Lower = allinone_os_coef$hr_lower,
    HR_Upper = allinone_os_coef$hr_upper
  ))

  coef_comparison <- rbind(coef_comparison, data.frame(
    Biomarker = toupper(biomarker),
    Outcome = "OS",
    Model_Type = "Separate",
    HR = separate_os_coef$hr,
    HR_Lower = separate_os_coef$hr_lower,
    HR_Upper = separate_os_coef$hr_upper
  ))

  coef_comparison <- rbind(coef_comparison, data.frame(
    Biomarker = toupper(biomarker),
    Outcome = "PFS",
    Model_Type = "All-in-One",
    HR = allinone_pfs_coef$hr,
    HR_Lower = allinone_pfs_coef$hr_lower,
    HR_Upper = allinone_pfs_coef$hr_upper
  ))

  coef_comparison <- rbind(coef_comparison, data.frame(
    Biomarker = toupper(biomarker),
    Outcome = "PFS",
    Model_Type = "Separate",
    HR = separate_pfs_coef$hr,
    HR_Lower = separate_pfs_coef$hr_lower,
    HR_Upper = separate_pfs_coef$hr_upper
  ))
}

# Print coefficient comparison
cat("Hazard Ratios for Treatment × Biomarker Interactions:\n\n")
print(coef_comparison, digits = 3, row.names = FALSE)

# ===============================================================================
# CALCULATE UNADJUSTED LIFE-YEARS COMPARISON
# ===============================================================================

cat("\n========================================\n")
cat("Unadjusted Life-Years Comparison (Area Under Survival Curve)\n")
cat("========================================\n\n")

cat("Note: These are UNADJUSTED life-years (no utility weighting applied).\n")
cat("      Life-years are calculated as the area under the survival curve.\n\n")

# Function to calculate unadjusted life-years (restricted mean survival time)
# Uses trapezoidal rule to integrate survival curve
calculate_life_years <- function(survival_probs, time_points) {
  # Convert time from weeks to years
  time_years <- time_points / 52

  # Trapezoidal integration
  n <- length(time_years)
  life_years <- sum((survival_probs[-n] + survival_probs[-1]) / 2 * diff(time_years))

  return(life_years)
}

# Calculate control strategy life-years first (reference for all comparisons)
control_os_ly <- calculate_life_years(control_os, time_points)
control_pfs_ly <- calculate_life_years(control_pfs, time_points)

# Build THREE life-years comparison tables:
# 1. Biomarker-positive subgroups
# 2. Biomarker-negative subgroups
# 3. Overall strategy (population-weighted)

# ---- Table 1: Biomarker-Positive Subgroups ----
ly_biomarker_positive <- data.frame(
  Biomarker = character(),
  Outcome = character(),
  Model_Type = character(),
  Life_Years = numeric(),
  stringsAsFactors = FALSE
)

# ---- Table 2: Biomarker-Negative Subgroups ----
ly_biomarker_negative <- data.frame(
  Biomarker = character(),
  Outcome = character(),
  Model_Type = character(),
  Life_Years = numeric(),
  stringsAsFactors = FALSE
)

# ---- Table 3: Overall Strategy (Population-Weighted) ----
ly_overall_strategy <- data.frame(
  Biomarker = character(),
  Outcome = character(),
  Model_Type = character(),
  Life_Years = numeric(),
  stringsAsFactors = FALSE
)

for (biomarker in c("crp", "tlr", "tmb_braf")) {
  prevalence <- separate_predictions[[biomarker]]$prevalence

  # === BIOMARKER-POSITIVE ===
  # All-in-One model
  allinone_pos_os_ly <- calculate_life_years(predictions[[biomarker]]$biomarker_positive$os, time_points)
  allinone_pos_pfs_ly <- calculate_life_years(predictions[[biomarker]]$biomarker_positive$pfs, time_points)
  # Separate model
  separate_pos_os_ly <- calculate_life_years(separate_predictions[[biomarker]]$biomarker_positive$os, time_points)
  separate_pos_pfs_ly <- calculate_life_years(separate_predictions[[biomarker]]$biomarker_positive$pfs, time_points)

  ly_biomarker_positive <- rbind(ly_biomarker_positive,
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "All-in-One", Life_Years = round(allinone_pos_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "Separate", Life_Years = round(separate_pos_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "All-in-One", Life_Years = round(allinone_pos_pfs_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "Separate", Life_Years = round(separate_pos_pfs_ly, 3))
  )

  # === BIOMARKER-NEGATIVE ===
  # All-in-One model
  allinone_neg_os_ly <- calculate_life_years(predictions[[biomarker]]$biomarker_negative$os, time_points)
  allinone_neg_pfs_ly <- calculate_life_years(predictions[[biomarker]]$biomarker_negative$pfs, time_points)
  # Separate model
  separate_neg_os_ly <- calculate_life_years(separate_predictions[[biomarker]]$biomarker_negative$os, time_points)
  separate_neg_pfs_ly <- calculate_life_years(separate_predictions[[biomarker]]$biomarker_negative$pfs, time_points)

  ly_biomarker_negative <- rbind(ly_biomarker_negative,
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "All-in-One", Life_Years = round(allinone_neg_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "Separate", Life_Years = round(separate_neg_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "All-in-One", Life_Years = round(allinone_neg_pfs_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "Separate", Life_Years = round(separate_neg_pfs_ly, 3))
  )

  # === OVERALL STRATEGY (Population-Weighted) ===
  # All-in-One model (use existing predictions which are already population-weighted)
  allinone_overall_os_ly <- calculate_life_years(predictions[[biomarker]]$os, time_points)
  allinone_overall_pfs_ly <- calculate_life_years(predictions[[biomarker]]$pfs, time_points)
  # Separate model (calculate population-weighted from subgroups)
  separate_overall_os <- prevalence * separate_predictions[[biomarker]]$biomarker_positive$os +
                         (1 - prevalence) * separate_predictions[[biomarker]]$biomarker_negative$os
  separate_overall_pfs <- prevalence * separate_predictions[[biomarker]]$biomarker_positive$pfs +
                          (1 - prevalence) * separate_predictions[[biomarker]]$biomarker_negative$pfs
  separate_overall_os_ly <- calculate_life_years(separate_overall_os, time_points)
  separate_overall_pfs_ly <- calculate_life_years(separate_overall_pfs, time_points)

  ly_overall_strategy <- rbind(ly_overall_strategy,
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "All-in-One", Life_Years = round(allinone_overall_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "OS", Model_Type = "Separate", Life_Years = round(separate_overall_os_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "All-in-One", Life_Years = round(allinone_overall_pfs_ly, 3)),
    data.frame(Biomarker = toupper(biomarker), Outcome = "PFS", Model_Type = "Separate", Life_Years = round(separate_overall_pfs_ly, 3))
  )
}

cat("Unadjusted Life-Years - Biomarker-Positive Subgroups:\n\n")
print(ly_biomarker_positive, row.names = FALSE)

cat("\nUnadjusted Life-Years - Biomarker-Negative Subgroups:\n\n")
print(ly_biomarker_negative, row.names = FALSE)

cat("\nUnadjusted Life-Years - Overall Strategy (Population-Weighted):\n\n")
print(ly_overall_strategy, row.names = FALSE)

# Store control life-years for use in reporting
control_life_years <- list(
  os = control_os_ly,
  pfs = control_pfs_ly
)

# Keep legacy ly_comparison for backward compatibility
ly_comparison <- ly_biomarker_positive

# ===============================================================================
# SUMMARY OUTPUT
# ===============================================================================

cat("\n========================================\n")
cat("SUMMARY: Separate vs All-in-One Model Comparison\n")
cat("========================================\n\n")

cat("Parametric distribution used:", best_os_dist, "(OS),", best_pfs_dist, "(PFS)\n\n")

cat("Objects created:\n")
cat("  - separate_models: List of fitted flexsurvreg objects (one per biomarker)\n")
cat("  - separate_predictions: List of survival predictions from separate models\n")
cat("  - control_os, control_pfs: Control strategy survival curves\n")
cat("  - coef_comparison: Hazard ratio comparison table\n")
cat("  - ly_comparison: Unadjusted life-years comparison table\n\n")

cat("Key finding: Compare HR values between 'All-in-One' and 'Separate' models.\n")
cat("If HRs differ substantially, biomarker correlation is affecting treatment effect estimates.\n")

# For backward compatibility, also create aliases
independent_models <- separate_models
independent_predictions <- separate_predictions
