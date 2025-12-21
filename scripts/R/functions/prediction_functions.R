# ===============================================================================
# PREDICTION FUNCTIONS FOR PARAMETRIC SURVIVAL MODELS
# ===============================================================================
# This file contains functions for generating individual and population-level
# survival predictions from fitted parametric models
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(survival, flexsurv, dplyr, tidyr)

# ===============================================================================
# POPULATION AVERAGING FUNCTIONS
# ===============================================================================
# These functions implement methodologically rigorous population-averaged
# predictions that account for the full distribution of prognostic factors
# (age, sex, biomarkers) in the patient population.
# ===============================================================================

# Helper function to extract survival probabilities from flexsurv predict objects
# Handles the nested structure: pred$.pred[[i]]$.pred_survival
extract_all_survival_probabilities <- function(pred_object) {
  # pred_object is the result of predict(flexsurvreg_model, newdata = data, times = time_points)
  # Structure: tibble with $.pred column containing list of predictions
  # Each element in $.pred has $.pred_survival (vector of survival probs at each time point)

  # Extract all survival probability vectors (one per patient)
  all_surv_probs <- lapply(pred_object$.pred, function(x) x$.pred_survival)

  # Convert to matrix: rows = time points, columns = patients
  surv_matrix <- do.call(cbind, all_surv_probs)

  return(surv_matrix)
}

# Main function for population-averaged predictions across all strategies
# Uses actual patient-level covariate distributions instead of reference patient
generate_population_averaged_predictions <- function(models, strategies_df,
                                                     data_complete, time_points) {

  # Extract treatment level references
  exp_rx <- levels(data_complete$Rx)[2]  # Experimental treatment
  ctrl_rx <- levels(data_complete$Rx)[1]  # Control treatment

  # Initialize predictions list
  predictions <- list()

  # -------------------------------------------------------------------------
  # CONTROL STRATEGY: All patients receive standard of care
  # -------------------------------------------------------------------------

  cat("Generating population-averaged predictions for control strategy...\n")

  # Create dataset with all patients assigned to control treatment
  control_data <- data_complete
  control_data$Rx <- factor(ctrl_rx, levels = levels(data_complete$Rx))

  # Predict for all patients
  os_pred_control <- predict(models$os, newdata = control_data,
                              type = "survival", times = time_points)
  pfs_pred_control <- predict(models$pfs, newdata = control_data,
                               type = "survival", times = time_points)

  # Extract and average across population
  os_matrix <- extract_all_survival_probabilities(os_pred_control)
  pfs_matrix <- extract_all_survival_probabilities(pfs_pred_control)

  control_os <- rowMeans(os_matrix, na.rm = TRUE)
  control_pfs <- rowMeans(pfs_matrix, na.rm = TRUE)

  predictions$control <- list(
    strategy = "control",
    os = control_os,
    pfs = control_pfs
  )

  # -------------------------------------------------------------------------
  # BIOMARKER-GUIDED STRATEGIES: Treatment assigned by biomarker status
  # -------------------------------------------------------------------------

  # Loop through each biomarker strategy
  for (biomarker_name in c("crp", "tlr", "tmb_braf")) {

    cat("Generating population-averaged predictions for", biomarker_name, "strategy...\n")

    # Get strategy info
    strategy_info <- strategies_df[strategies_df$id == biomarker_name, ]

    # -----------------------------------------------------------------------
    # BIOMARKER-POSITIVE SUBGROUP: Receive experimental treatment
    # -----------------------------------------------------------------------

    # Subset to biomarker-positive patients
    # Use explicit type conversion for robustness (works with factor, character, or numeric)
    biomarker_pos_data <- data_complete[as.numeric(as.character(data_complete[[biomarker_name]])) == 1, ]

    if (nrow(biomarker_pos_data) > 0) {
      # Assign experimental treatment
      biomarker_pos_data$Rx <- factor(exp_rx, levels = levels(data_complete$Rx))

      # Predict for all biomarker-positive patients
      os_pred_pos <- predict(models$os, newdata = biomarker_pos_data,
                             type = "survival", times = time_points)
      pfs_pred_pos <- predict(models$pfs, newdata = biomarker_pos_data,
                              type = "survival", times = time_points)

      # Extract and average
      os_matrix_pos <- extract_all_survival_probabilities(os_pred_pos)
      pfs_matrix_pos <- extract_all_survival_probabilities(pfs_pred_pos)

      biomarker_pos_os <- rowMeans(os_matrix_pos, na.rm = TRUE)
      biomarker_pos_pfs <- rowMeans(pfs_matrix_pos, na.rm = TRUE)

    } else {
      warning(paste("No biomarker-positive patients for", biomarker_name))
      biomarker_pos_os <- rep(NA, length(time_points))
      biomarker_pos_pfs <- rep(NA, length(time_points))
    }

    # -----------------------------------------------------------------------
    # BIOMARKER-NEGATIVE SUBGROUP: Receive control treatment
    # -----------------------------------------------------------------------

    # Subset to biomarker-negative patients
    # Use explicit type conversion for robustness (works with factor, character, or numeric)
    biomarker_neg_data <- data_complete[as.numeric(as.character(data_complete[[biomarker_name]])) == 0, ]

    if (nrow(biomarker_neg_data) > 0) {
      # Assign control treatment
      biomarker_neg_data$Rx <- factor(ctrl_rx, levels = levels(data_complete$Rx))

      # Predict for all biomarker-negative patients
      os_pred_neg <- predict(models$os, newdata = biomarker_neg_data,
                             type = "survival", times = time_points)
      pfs_pred_neg <- predict(models$pfs, newdata = biomarker_neg_data,
                              type = "survival", times = time_points)

      # Extract and average
      os_matrix_neg <- extract_all_survival_probabilities(os_pred_neg)
      pfs_matrix_neg <- extract_all_survival_probabilities(pfs_pred_neg)

      biomarker_neg_os <- rowMeans(os_matrix_neg, na.rm = TRUE)
      biomarker_neg_pfs <- rowMeans(pfs_matrix_neg, na.rm = TRUE)

    } else {
      warning(paste("No biomarker-negative patients for", biomarker_name))
      biomarker_neg_os <- rep(NA, length(time_points))
      biomarker_neg_pfs <- rep(NA, length(time_points))
    }

    # -----------------------------------------------------------------------
    # POPULATION-WEIGHTED AVERAGE: Combine using prevalence
    # -----------------------------------------------------------------------

    # Calculate prevalence (proportion biomarker-positive)
    prevalence <- mean(as.numeric(as.character(data_complete[[biomarker_name]])), na.rm = TRUE)

    # Weighted average
    weighted_os <- prevalence * biomarker_pos_os + (1 - prevalence) * biomarker_neg_os
    weighted_pfs <- prevalence * biomarker_pos_pfs + (1 - prevalence) * biomarker_neg_pfs

    # Store results in same format as current predict_strategy() output
    predictions[[biomarker_name]] <- list(
      strategy = biomarker_name,
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
  }

  cat("Population-averaged predictions complete for all strategies.\n")

  return(predictions)
}