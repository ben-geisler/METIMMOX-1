# ===============================================================================
# MAIN MODEL FUNCTION - Partitioned Survival Model with PSA Support
# ===============================================================================
# This function runs the cost-effectiveness model for all strategies
# Supports both deterministic and probabilistic sensitivity analysis (PSA)
# In PSA mode, uses correlated resampled models for PFS and OS
# ===============================================================================

model_fun <- function(params, time_horizon = 520, cl = 1/52, determpsa = "det",
                      return_traces = FALSE, sim_idx = NULL) {
  # Get the strategies from our global variables
  strat_names <- strategies  # This should be defined in section 2

  # =========================================================================
  # PSA MODE: Second-order Monte Carlo simulation
  # =========================================================================
  # Each PSA iteration samples ONE individual patient from the dataset
  # That patient is evaluated across all strategies with treatment assigned
  # based on their biomarker status. This approach captures:
  #   1. Parameter uncertainty (from resampled survival models)
  #   2. Patient heterogeneity (from sampling individuals with varying covariates)
  # Over many iterations, biomarker prevalence is naturally reflected in the
  # proportion of iterations where sampled patients are biomarker-positive.
  # See GitHub Issue #70 for detailed methodology discussion.
  # =========================================================================
  if(determpsa == "psa" && !is.null(sim_idx)) {

    if (sim_idx %% 100 == 0) {
      cat("PSA iteration", sim_idx, "- second-order Monte Carlo simulation\n")
    }

    # -----------------------------------------------------------------------
    # STEP 1: Sample ONE patient from the dataset
    # -----------------------------------------------------------------------
    # Use modulo to cycle through patients if n_samples > n_patients
    patient_idx <- ((sim_idx - 1) %% nrow(data)) + 1
    sampled_patient <- data[patient_idx, ]

    # Extract treatment level references
    exp_rx <- levels(data$Rx)[2]  # Experimental treatment
    ctrl_rx <- levels(data$Rx)[1]  # Control treatment

    # -----------------------------------------------------------------------
    # STEP 2: CONTROL STRATEGY - Sampled patient receives control treatment
    # -----------------------------------------------------------------------

    tryCatch({
      # Create patient data with control treatment assigned
      patient_ctrl <- sampled_patient
      patient_ctrl$Rx <- factor(ctrl_rx, levels = levels(data$Rx))

      # Predict for this patient using resampled model #sim_idx
      os_control <- generate_sampling_predictions(
        sampling_model_list = sampling_models$control,
        outcome = "os",
        sample_idx = sim_idx,
        newdata = patient_ctrl,
        time_points = seq(0, time_horizon)
      )

      pfs_control <- generate_sampling_predictions(
        sampling_model_list = sampling_models$control,
        outcome = "pfs",
        sample_idx = sim_idx,
        newdata = patient_ctrl,
        time_points = seq(0, time_horizon)
      )

      # Check for NA values and handle
      if (any(is.na(os_control)) || any(is.na(pfs_control))) {
        warning("NA values in control survival curves for sim ", sim_idx,
                " - using base case curves")
        os_control <- params$p_os$control_OS
        pfs_control <- params$p_pfs$control_PFS
      } else {
        # Update params with this patient's predictions
        params$p_os$control_OS <- os_control
        params$p_pfs$control_PFS <- pfs_control
      }

    }, error = function(e) {
      warning("Failed to generate control curves for sim ", sim_idx, ": ",
              conditionMessage(e), " - using base case")
      # Keep base case curves in params
    })

    # -----------------------------------------------------------------------
    # STEP 3: BIOMARKER STRATEGIES - Treatment assigned by patient's biomarker
    # -----------------------------------------------------------------------
    # For each biomarker strategy, the sampled patient's biomarker status
    # determines their treatment: biomarker+ → experimental, biomarker- → control

    for(biomarker in biomarkers) {
      tryCatch({
        # Check this patient's biomarker status
        biomarker_status <- as.character(sampled_patient[[biomarker]])

        # Determine treatment based on biomarker status
        treatment_for_patient <- if(biomarker_status == "1") exp_rx else ctrl_rx

        # Create patient data with assigned treatment
        patient_biomarker <- sampled_patient
        patient_biomarker$Rx <- factor(treatment_for_patient, levels = levels(data$Rx))

        # Predict for this patient using resampled model #sim_idx
        os_pred <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "os",
          sample_idx = sim_idx,
          newdata = patient_biomarker,
          time_points = seq(0, time_horizon)
        )

        pfs_pred <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "pfs",
          sample_idx = sim_idx,
          newdata = patient_biomarker,
          time_points = seq(0, time_horizon)
        )

        # Check for NA values
        if (any(is.na(os_pred)) || any(is.na(pfs_pred))) {
          warning("NA values in ", biomarker, " survival curves for sim ", sim_idx,
                  " - using base case curves")
          # Keep base case curves
        } else {
          # For second-order Monte Carlo, we have ONE patient per iteration
          # Store their predictions in both pos/neg curves (same values)
          # and adjust prevalence to reflect this patient's status
          params$p_os[[paste0(biomarker, "_pos_OS")]] <- os_pred
          params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- pfs_pred
          params$p_os[[paste0(biomarker, "_neg_OS")]] <- os_pred
          params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- pfs_pred

          # Update prevalence for this iteration:
          # If patient is biomarker+: prevalence = 1.0
          # If patient is biomarker-: prevalence = 0.0
          # This ensures the weighted calculation gives the patient's outcome
          params[[paste0("p_", biomarker)]] <- as.numeric(biomarker_status == "1")
        }

      }, error = function(e) {
        warning("Failed to generate ", biomarker, " curves for sim ", sim_idx, ": ",
                conditionMessage(e), " - using base case")
        # Keep base case curves in params
      })
    }
  }
  
  # =========================================================================
  # INITIALIZE RESULTS STORAGE
  # =========================================================================
  
  results <- data.frame(
    Strategy = as.character(strat_names),
    Cost = numeric(length(strat_names)),
    Effect = numeric(length(strat_names)),
    stringsAsFactors = FALSE
  )
  
  traces <- list()
  if(return_traces) {
    for(strat in strat_names) {
      traces[[strat]] <- list()
    }
  }
  
  # Create discount factors over time
  v_dw_c <- 1 / (1 + params$dr_costs)^(seq(0, time_horizon) / 52)
  v_dw_e <- 1 / (1 + params$dr_effects)^(seq(0, time_horizon) / 52)
  
  # =========================================================================
  # CONTROL STRATEGY
  # =========================================================================
  
  # Get control survival curves (either base case or from resampling)
  os_control <- params$p_os[["control_OS"]]
  pfs_control <- params$p_pfs[["control_PFS"]]
  
  # Calculate state occupancy for partitioned survival model
  # Progression-free: PFS curve
  # Progressed: OS - PFS (ensures non-negative)
  # Dead: 1 - OS
  p_pf_control <- pfs_control
  p_p_control <- pmax(os_control - pfs_control, 0)
  p_d_control <- 1 - os_control
  
  # Force initial state to be 100% progression-free (all patients start alive)
  p_pf_control[1] <- 1.0
  p_p_control[1] <- 0.0
  p_d_control[1] <- 0.0
  
  # Calculate control QALYs and costs
  control_results <- calculate_outcomes(
    params = params,
    p_pf = p_pf_control,
    p_p = p_p_control,
    p_d = p_d_control,
    treatment_type = "standard",
    biomarker = NULL,
    v_dw_c = v_dw_c,
    v_dw_e = v_dw_e,
    cl = cl
  )
  
  # Store results for control
  results$Cost[results$Strategy == "control"] <- control_results$costs_total
  results$Effect[results$Strategy == "control"] <- control_results$qalys_total
  
  # Store traces for control if requested
  if(return_traces) {
    traces$control <- list(
      cycles = seq(0, time_horizon),
      p_pf = p_pf_control,
      p_p = p_p_control,
      p_d = p_d_control,
      qalys_total = control_results$qalys_total,
      costs_total = control_results$costs_total
    )
  }
  
  # =========================================================================
  # BIOMARKER STRATEGIES
  # =========================================================================
  
  # Loop through each biomarker strategy
  for(biomarker in biomarkers) {
    # Get biomarker prevalence directly from params
    biomarker_prev <- params[[paste0("p_", biomarker)]]
    
    # -----------------------------------------------------------------------
    # BIOMARKER POSITIVE SUBGROUP (receives experimental treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker positive survival curves
    os_pos <- params$p_os[[paste0(biomarker, "_pos_OS")]]
    pfs_pos <- params$p_pfs[[paste0(biomarker, "_pos_PFS")]]
    
    # Calculate positive state occupancy
    p_pf_pos <- pfs_pos
    p_p_pos <- pmax(os_pos - pfs_pos, 0)
    p_d_pos <- 1 - os_pos
    
    # Force initial state to be 100% progression-free
    p_pf_pos[1] <- 1.0
    p_p_pos[1] <- 0.0
    p_d_pos[1] <- 0.0
    
    # -----------------------------------------------------------------------
    # BIOMARKER NEGATIVE SUBGROUP (receives standard treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker negative survival curves
    os_neg <- params$p_os[[paste0(biomarker, "_neg_OS")]]
    pfs_neg <- params$p_pfs[[paste0(biomarker, "_neg_PFS")]]
    
    # Calculate negative state occupancy
    p_pf_neg <- pfs_neg
    p_p_neg <- pmax(os_neg - pfs_neg, 0)
    p_d_neg <- 1 - os_neg
    
    # Force initial state to be 100% progression-free
    p_pf_neg[1] <- 1.0
    p_p_neg[1] <- 0.0
    p_d_neg[1] <- 0.0
    
    # -----------------------------------------------------------------------
    # CALCULATE OUTCOMES FOR EACH SUBGROUP
    # -----------------------------------------------------------------------
    
    # Calculate positive subgroup outcomes
    pos_results <- calculate_outcomes(
      params = params,
      p_pf = p_pf_pos,
      p_p = p_p_pos,
      p_d = p_d_pos,
      treatment_type = "experimental",
      biomarker = biomarker,
      v_dw_c = v_dw_c,
      v_dw_e = v_dw_e,
      cl = cl
    )
    
    # Calculate negative subgroup outcomes
    neg_results <- calculate_outcomes(
      params = params,
      p_pf = p_pf_neg,
      p_p = p_p_neg,
      p_d = p_d_neg,
      treatment_type = "standard",
      biomarker = biomarker,
      v_dw_c = v_dw_c,
      v_dw_e = v_dw_e,
      cl = cl
    )
    
    # -----------------------------------------------------------------------
    # COMBINE SUBGROUPS USING BIOMARKER PREVALENCE
    # -----------------------------------------------------------------------
    
    # Weighted average of costs and effects based on biomarker prevalence
    total_cost <- biomarker_prev * pos_results$costs_total + 
      (1 - biomarker_prev) * neg_results$costs_total
    total_effect <- biomarker_prev * pos_results$qalys_total + 
      (1 - biomarker_prev) * neg_results$qalys_total
    
    # Store results for this biomarker strategy
    results$Cost[results$Strategy == biomarker] <- total_cost
    results$Effect[results$Strategy == biomarker] <- total_effect
    
    # Store traces if requested
    if(return_traces) {
      # Create weighted values for summary traces
      traces[[biomarker]] <- list(
        cycles = seq(0, time_horizon),
        biomarker_prevalence = biomarker_prev,
        p_pf = biomarker_prev * p_pf_pos + (1 - biomarker_prev) * p_pf_neg,
        p_p = biomarker_prev * p_p_pos + (1 - biomarker_prev) * p_p_neg,
        p_d = biomarker_prev * p_d_pos + (1 - biomarker_prev) * p_d_neg,
        positive = list(
          p_pf = p_pf_pos,
          p_p = p_p_pos,
          p_d = p_d_pos,
          qalys_total = pos_results$qalys_total,
          costs_total = pos_results$costs_total
        ),
        negative = list(
          p_pf = p_pf_neg,
          p_p = p_p_neg,
          p_d = p_d_neg,
          qalys_total = neg_results$qalys_total,
          costs_total = neg_results$costs_total
        ),
        qalys_total = total_effect,
        costs_total = total_cost
      )
    }
  }
  
  # =========================================================================
  # RETURN RESULTS
  # =========================================================================
  
  # Return results with or without traces
  if(return_traces) {
    return(list(
      results = results,
      traces = traces
    ))
  } else {
    return(results)
  }
}