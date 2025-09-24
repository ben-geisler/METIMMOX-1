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

# Function for individual patient-level predictions (PSA-ready)
# Used for probabilistic sensitivity analysis and individual patient predictions
predict_patient <- function(models, age, sex, rx_treatment, crp_val, tlr_val, tmb_braf_val, 
                            time_points, data_complete) {
  
  # Create patient data with proper factor levels
  patient_data <- data.frame(
    Age = age,
    sex = factor(sex, levels = levels(data_complete$sex)),
    Rx = factor(rx_treatment, levels = levels(data_complete$Rx)),
    crp = factor(crp_val, levels = levels(data_complete$crp)),
    tlr = factor(tlr_val, levels = levels(data_complete$tlr)),
    tmb_braf = factor(tmb_braf_val, levels = levels(data_complete$tmb_braf))
  )
  
  # Generate predictions
  os_pred <- predict(models$os, newdata = patient_data, type = "survival", times = time_points)
  pfs_pred <- predict(models$pfs, newdata = patient_data, type = "survival", times = time_points)
  
  # Return survival probabilities
  return(list(
    os = unnest(os_pred, .pred)$.pred_survival,
    pfs = unnest(pfs_pred, .pred)$.pred_survival
  ))
}

# Function for strategy-level predictions with population marginalization
# Handles both control strategy and biomarker-guided strategies
predict_strategy <- function(models, strategy_name, strategies_df, 
                             ref_age, ref_sex, ref_crp, ref_tlr, ref_tmb_braf,
                             exp_rx, ctrl_rx, time_points, data_complete) {
  
  # Get strategy info from strategies dataframe
  strategy_info <- strategies_df[strategies_df$id == strategy_name, ]
  
  if (strategy_name == "control") {
    # Standard of care: all patients get control treatment
    # Uses modal categories for representative population
    pred <- predict_patient(
      models = models,
      age = ref_age,
      sex = ref_sex,
      rx_treatment = ctrl_rx,
      crp_val = ref_crp,
      tlr_val = ref_tlr,
      tmb_braf_val = ref_tmb_braf,
      time_points = time_points,
      data_complete = data_complete
    )
    
    return(list(
      strategy = strategy_name,
      os = pred$os,
      pfs = pred$pfs
    ))
    
  } else {
    # Biomarker-guided strategies: marginalize over biomarker status
    biomarker_name <- strategy_name  # crp, tlr, or tmb_braf
    prevalence <- strategy_info$prevalence
    
    # Biomarker-positive patients: get experimental treatment
    pred_positive <- predict_patient(
      models = models,
      age = ref_age,
      sex = ref_sex,
      rx_treatment = exp_rx,
      crp_val = ifelse(biomarker_name == "crp", "1", ref_crp),
      tlr_val = ifelse(biomarker_name == "tlr", "1", ref_tlr),
      tmb_braf_val = ifelse(biomarker_name == "tmb_braf", "1", ref_tmb_braf),
      time_points = time_points,
      data_complete = data_complete
    )
    
    # Biomarker-negative patients: get control treatment
    pred_negative <- predict_patient(
      models = models,
      age = ref_age,
      sex = ref_sex,
      rx_treatment = ctrl_rx,
      crp_val = ifelse(biomarker_name == "crp", "0", ref_crp),
      tlr_val = ifelse(biomarker_name == "tlr", "0", ref_tlr),
      tmb_braf_val = ifelse(biomarker_name == "tmb_braf", "0", ref_tmb_braf),
      time_points = time_points,
      data_complete = data_complete
    )
    
    # Population-weighted marginalization using observed prevalence
    os_marginalized <- prevalence * pred_positive$os + (1 - prevalence) * pred_negative$os
    pfs_marginalized <- prevalence * pred_positive$pfs + (1 - prevalence) * pred_negative$pfs
    
    return(list(
      strategy = strategy_name,
      os = os_marginalized,
      pfs = pfs_marginalized,
      biomarker_positive = pred_positive,
      biomarker_negative = pred_negative,
      prevalence = prevalence
    ))
  }
}

# Helper function for batch predictions (useful for PSA)
# Generates predictions for multiple parameter sets at once
batch_predict_patients <- function(models, param_matrix, time_points, data_complete) {
  
  # param_matrix should have columns: age, sex, rx_treatment, crp_val, tlr_val, tmb_braf_val
  # Each row represents one parameter set
  
  results <- list()
  
  for (i in 1:nrow(param_matrix)) {
    params <- param_matrix[i, ]
    
    result <- predict_patient(
      models = models,
      age = params$age,
      sex = params$sex,
      rx_treatment = params$rx_treatment,
      crp_val = params$crp_val,
      tlr_val = params$tlr_val,
      tmb_braf_val = params$tmb_braf_val,
      time_points = time_points,
      data_complete = data_complete
    )
    
    results[[i]] <- result
  }
  
  return(results)
}

# Function to generate predictions for all strategies at once
# Useful for ensuring consistency across strategies
generate_all_strategy_predictions <- function(models, strategies_df, 
                                              ref_age, ref_sex, ref_crp, ref_tlr, ref_tmb_braf,
                                              exp_rx, ctrl_rx, time_points, data_complete) {
  
  predictions <- list()
  
  for (strategy in strategies_df$id) {
    predictions[[strategy]] <- predict_strategy(
      models = models, 
      strategy_name = strategy, 
      strategies_df = strategies_df,
      ref_age = ref_age,
      ref_sex = ref_sex, 
      ref_crp = ref_crp,
      ref_tlr = ref_tlr,
      ref_tmb_braf = ref_tmb_braf,
      exp_rx = exp_rx,
      ctrl_rx = ctrl_rx,
      time_points = time_points,
      data_complete = data_complete
    )
  }
  
  return(predictions)
}