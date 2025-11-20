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
  # PSA MODE: Generate survival curves from correlated resampled models
  # =========================================================================
  if(determpsa == "psa" && !is.null(sim_idx)) {
    # In PSA mode with sim_idx provided, we generate survival curves
    # from resampled models for this specific simulation
    # CRITICAL: PFS and OS come from the SAME resampled model to maintain correlation
    
    if (sim_idx %% 100 == 0) {
      cat("PSA iteration", sim_idx, "- generating correlated survival curves\n")
    }
    
    # Generate control strategy survival curves from resampled model #sim_idx
    tryCatch({
      os_control <- generate_sampling_predictions(
        sampling_model_list = sampling_models$control,
        outcome = "os",
        sample_idx = sim_idx,
        newdata = NULL,
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
        os_control <- params$p_os$control_OS
        pfs_control <- params$p_pfs$control_PFS
      } else {
        # Update params with these resampled curves
        params$p_os$control_OS <- os_control
        params$p_pfs$control_PFS <- pfs_control
      }
      
    }, error = function(e) {
      warning("Failed to generate control curves for sim ", sim_idx, ": ", 
              conditionMessage(e), " - using base case")
      # Keep base case curves in params
    })
    
    # Generate curves for each biomarker strategy
    for(biomarker in biomarkers) {
      tryCatch({
        # Get reference to experimental and control treatment levels
        exp_rx <- levels(data$Rx)[2]
        ctrl_rx <- levels(data$Rx)[1]
        
        # Get reference values for ALL covariates (mean/modal)
        # CRITICAL FIX: Must get reference values for ALL THREE biomarkers
        ref_age <- mean(data$Age, na.rm = TRUE)
        ref_sex <- as.numeric(names(sort(table(data$sex), decreasing = TRUE))[1])
        ref_crp <- as.numeric(names(sort(table(data$crp), decreasing = TRUE))[1])
        ref_tlr <- as.numeric(names(sort(table(data$tlr), decreasing = TRUE))[1])
        ref_tmb_braf <- as.numeric(names(sort(table(data$tmb_braf), decreasing = TRUE))[1])
        
        # Create newdata for biomarker positive + experimental treatment
        # CRITICAL FIX: Must include ALL biomarker variables (crp, tlr, tmb_braf)
        # The resampled models were fitted with all three biomarkers as covariates
        # Setting the target biomarker to 1, others to their reference values
        newdata_pos_exp <- data.frame(
          Age = ref_age,
          sex = factor(ref_sex, levels = levels(data$sex)),
          Rx = factor(exp_rx, levels = levels(data$Rx)),
          crp = factor(ifelse(biomarker == "crp", 1, ref_crp), levels = c(0, 1)),
          tlr = factor(ifelse(biomarker == "tlr", 1, ref_tlr), levels = c(0, 1)),
          tmb_braf = factor(ifelse(biomarker == "tmb_braf", 1, ref_tmb_braf), levels = c(0, 1))
        )
        
        # Create newdata for biomarker negative + control treatment
        # CRITICAL FIX: Must include ALL biomarker variables (crp, tlr, tmb_braf)
        # Setting the target biomarker to 0, others to their reference values
        newdata_neg_ctrl <- data.frame(
          Age = ref_age,
          sex = factor(ref_sex, levels = levels(data$sex)),
          Rx = factor(ctrl_rx, levels = levels(data$Rx)),
          crp = factor(ifelse(biomarker == "crp", 0, ref_crp), levels = c(0, 1)),
          tlr = factor(ifelse(biomarker == "tlr", 0, ref_tlr), levels = c(0, 1)),
          tmb_braf = factor(ifelse(biomarker == "tmb_braf", 0, ref_tmb_braf), levels = c(0, 1))
        )
        
        # Generate all four predictions for this simulation from SAME resampled model
        # This maintains correlation between PFS and OS
        os_pos_exp <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "os",
          newdata = newdata_pos_exp,
          sample_idx = sim_idx,
          time_points = seq(0, time_horizon)
        )

        pfs_pos_exp <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "pfs",
          newdata = newdata_pos_exp,
          sample_idx = sim_idx,
          time_points = seq(0, time_horizon)
        )
        
        # DIAGNOSTIC CODE - ADD TEMPORARILY
        if (sim_idx == 1 && biomarker == "crp") {
          cat("\n=== DIAGNOSTIC OUTPUT FOR CRP, SIM 1 ===\n")
          cat("os_pos_exp length:", length(os_pos_exp), "\n")
          cat("os_pos_exp range:", range(os_pos_exp, na.rm=TRUE), "\n")
          cat("pfs_pos_exp length:", length(pfs_pos_exp), "\n")
          cat("pfs_pos_exp range:", range(pfs_pos_exp, na.rm=TRUE), "\n")
          cat("os_neg_ctrl range:", range(os_neg_ctrl, na.rm=TRUE), "\n")
          cat("pfs_neg_ctrl range:", range(pfs_neg_ctrl, na.rm=TRUE), "\n")
          cat("Base case os_pos range:", range(params$p_os$crp_pos_OS, na.rm=TRUE), "\n")
          cat("Base case pfs_pos range:", range(params$p_pfs$crp_pos_PFS, na.rm=TRUE), "\n")
        }
        
        os_neg_ctrl <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "os",
          newdata = newdata_neg_ctrl,
          sample_idx = sim_idx,
          time_points = seq(0, time_horizon)
        )

        pfs_neg_ctrl <- generate_sampling_predictions(
          sampling_model_list = sampling_models[[biomarker]],
          outcome = "pfs",
          newdata = newdata_neg_ctrl,
          sample_idx = sim_idx,
          time_points = seq(0, time_horizon)
        )
        
        # Check for NA values
        if (any(is.na(os_pos_exp)) || any(is.na(pfs_pos_exp)) ||
            any(is.na(os_neg_ctrl)) || any(is.na(pfs_neg_ctrl))) {
          warning("NA values in ", biomarker, " survival curves for sim ", sim_idx,
                  " - using base case curves")
          # Keep base case curves
        } else {
          # Update params with these curves
          params$p_os[[paste0(biomarker, "_pos_OS")]] <- os_pos_exp
          params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- pfs_pos_exp
          params$p_os[[paste0(biomarker, "_neg_OS")]] <- os_neg_ctrl
          params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- pfs_neg_ctrl
          
          # Also update weighted averages
          biomarker_prev <- params[[paste0("p_", biomarker)]]
          params$p_os[[paste0(biomarker, "_weighted_OS")]] <- 
            biomarker_prev * os_pos_exp + (1 - biomarker_prev) * os_neg_ctrl
          
          params$p_pfs[[paste0(biomarker, "_weighted_PFS")]] <- 
            biomarker_prev * pfs_pos_exp + (1 - biomarker_prev) * pfs_neg_ctrl
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