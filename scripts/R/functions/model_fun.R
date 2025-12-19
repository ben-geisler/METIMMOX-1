# ===============================================================================
# MAIN MODEL FUNCTION - Partitioned Survival Model with PSA Support
# ===============================================================================
# This function runs the cost-effectiveness model for all strategies.
#
# Supports both deterministic and probabilistic sensitivity analysis (PSA):
#   - Deterministic: Uses base case survival curves and parameters
#   - PSA: Uses resampled survival models with subgroup population averaging
#
# PSA Methodology (Issues #73, #74, #87):
#   - Parameter uncertainty: Captured by resampled survival models (varies
#     BETWEEN PSA iterations)
#   - Patient heterogeneity: Integrated out via population averaging (averaged
#     WITHIN each iteration)
#   - Prevalence: Kept at true population values, NOT set to 0/1
#
# This separation ensures EVPI correctly measures the value of reducing
# parameter uncertainty, not patient heterogeneity.
#
# Biological Constraint (Issue #76):
#   - PFS is capped at OS to prevent invalid state occupancy (PFS > OS can
#     occur with resampled models). Without this, states can sum to >1.
#
# Fallback Tracking (Issue #79):
#   - When PSA survival curve generation fails, the function falls back to
#     base case curves. This is now tracked via a "fallback_used" attribute
#     on the returned results, allowing the PSA loop to detect and replace
#     failed iterations instead of mixing different estimation methods.
# ===============================================================================

model_fun <- function(params, time_horizon = 520, cl = 1/52, determpsa = "det",
                      return_traces = FALSE, sim_idx = NULL) {
  # Get the strategies from our global variables
  strat_names <- strategies  # This should be defined in section 2

 # Track if fallback to base case was used (Issue #79)
  fallback_used <- FALSE

  # =========================================================================
  # PSA MODE: Subgroup Population Averaging
  # =========================================================================
  # Each PSA iteration:
  #   1. Uses resampled survival model i (captures parameter uncertainty)
  #   2. Predicts for ALL patients in each subgroup, then averages
  #   3. Keeps prevalence at true population value
  #
  # This properly separates:
  #   - Parameter uncertainty (varies BETWEEN iterations via resampled models)
  #   - Patient heterogeneity (averaged WITHIN each iteration)
  #
  # See GitHub Issues #73, #74, #87 for methodology discussion.
  # =========================================================================
  if(determpsa == "psa" && !is.null(sim_idx)) {

    # Validate sim_idx (Issue #46)
    if (!is.numeric(sim_idx) || length(sim_idx) != 1 || sim_idx < 1 ||
        sim_idx != floor(sim_idx)) {
      stop("sim_idx must be a positive integer. Got: ", sim_idx)
    }
    if (sim_idx > n_samples) {
      stop("sim_idx (", sim_idx, ") exceeds n_samples (", n_samples, ")")
    }

    if (sim_idx %% 100 == 0) {
      cat("PSA iteration", sim_idx, "- subgroup population averaging\n")
    }

    # -----------------------------------------------------------------------
    # CONTROL STRATEGY: Population-averaged predictions
    # -----------------------------------------------------------------------
    # Predict for ALL patients in control arm, then average

    tryCatch({
      # Use generate_sampling_predictions with newdata = NULL
      # This predicts for all patients and returns the average
      os_control <- generate_sampling_predictions(
        sampling_model_list = sampling_models$control,
        outcome = "os",
        sample_idx = sim_idx,
        newdata = NULL,  # NULL triggers population averaging
        time_points = seq(0, time_horizon)
      )

      pfs_control <- generate_sampling_predictions(
        sampling_model_list = sampling_models$control,
        outcome = "pfs",
        sample_idx = sim_idx,
        newdata = NULL,
        time_points = seq(0, time_horizon)
      )

      # Check for NA values and handle
      if (any(is.na(os_control)) || any(is.na(pfs_control))) {
        warning("NA values in control survival curves for sim ", sim_idx,
                " - using base case curves")
        fallback_used <- TRUE  # Issue #79: Track fallback
      } else {
        # Update params with population-averaged predictions
        params$p_os$control_OS <- os_control
        params$p_pfs$control_PFS <- pfs_control
      }

    }, error = function(e) {
      warning("Failed to generate control curves for sim ", sim_idx, ": ",
              conditionMessage(e), " - using base case")
      fallback_used <<- TRUE  # Issue #79: Track fallback (<<- for enclosing scope)
    })

    # -----------------------------------------------------------------------
    # BIOMARKER STRATEGIES: Subgroup population averaging
    # -----------------------------------------------------------------------
    # For each biomarker:
    #   - Biomarker+ subgroup: ALL biomarker+ patients with experimental Rx
    #   - Biomarker- subgroup: ALL biomarker- patients with control Rx
    # Weighted combination uses TRUE prevalence (not 0/1)

    for(biomarker in biomarkers) {
      tryCatch({
        # Get population-averaged predictions for BOTH subgroups
        preds_os <- generate_psa_population_averaged_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          biomarker_name = biomarker,
          outcome = "os",
          sample_idx = sim_idx,
          data_original = data,
          time_points = seq(0, time_horizon)
        )

        preds_pfs <- generate_psa_population_averaged_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          biomarker_name = biomarker,
          outcome = "pfs",
          sample_idx = sim_idx,
          data_original = data,
          time_points = seq(0, time_horizon)
        )

        # Check for NULL or NA values
        if (is.null(preds_os) || is.null(preds_pfs) ||
            any(is.na(preds_os$positive)) || any(is.na(preds_os$negative)) ||
            any(is.na(preds_pfs$positive)) || any(is.na(preds_pfs$negative))) {
          warning("NA/NULL values in ", biomarker, " survival curves for sim ", sim_idx,
                  " - using base case curves")
          fallback_used <- TRUE  # Issue #79: Track fallback
        } else {
          # Store subgroup-specific curves
          params$p_os[[paste0(biomarker, "_pos_OS")]] <- preds_os$positive
          params$p_os[[paste0(biomarker, "_neg_OS")]] <- preds_os$negative
          params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- preds_pfs$positive
          params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- preds_pfs$negative

          # DO NOT modify prevalence - keep true value from l_params_base
          # Weighted combination happens downstream using:
          #   prevalence * biomarker+ + (1-prevalence) * biomarker-
        }

      }, error = function(e) {
        warning("Failed to generate ", biomarker, " curves for sim ", sim_idx, ": ",
                conditionMessage(e), " - using base case")
        fallback_used <<- TRUE  # Issue #79: Track fallback (<<- for enclosing scope)
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

  # Enforce biological constraint: PFS cannot exceed OS (Issue #76)
  # Without this, state occupancy can sum to >1 when curves cross
  if (any(pfs_control > os_control)) {
    n_violations <- sum(pfs_control > os_control)
    warning("PFS > OS constraint enforced at ", n_violations,
            " time points (control)",
            if (!is.null(sim_idx)) paste0(" [sim ", sim_idx, "]") else "")
  }
  pfs_control <- pmin(pfs_control, os_control)

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

    # Enforce biological constraint: PFS cannot exceed OS (Issue #76)
    if (any(pfs_pos > os_pos)) {
      n_violations <- sum(pfs_pos > os_pos)
      warning("PFS > OS constraint enforced at ", n_violations,
              " time points (", biomarker, "+)",
              if (!is.null(sim_idx)) paste0(" [sim ", sim_idx, "]") else "")
    }
    pfs_pos <- pmin(pfs_pos, os_pos)

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

    # Enforce biological constraint: PFS cannot exceed OS (Issue #76)
    if (any(pfs_neg > os_neg)) {
      n_violations <- sum(pfs_neg > os_neg)
      warning("PFS > OS constraint enforced at ", n_violations,
              " time points (", biomarker, "-)",
              if (!is.null(sim_idx)) paste0(" [sim ", sim_idx, "]") else "")
    }
    pfs_neg <- pmin(pfs_neg, os_neg)

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
  # Attach fallback_used attribute for PSA loop to detect failures (Issue #79)
  if(return_traces) {
    result <- list(
      results = results,
      traces = traces
    )
  } else {
    result <- results
  }

  # Add fallback_used attribute so PSA loop can detect and replace failed iterations
  attr(result, "fallback_used") <- fallback_used
  return(result)
}