# Main model function - now with PSA functionality
model_fun <- function(params, time_horizon = 520, cl = 1/52, determpsa = "det", 
                      return_traces = FALSE, sim_idx = NULL) {
  # Get the strategies from our global variables
  strat_names <- strategies  # This should be defined in section 2
  
  # Handle PSA mode
  if(determpsa == "psa" && !is.null(sim_idx)) {
    # If in PSA mode and sim_idx is provided, we need to generate survival curves
    # from the bootstrap samples for this simulation
    
    # Generate survival curves for control
    os_control <- generate_bootstrap_predictions(
      boot_models$control$os, sample_idx = sim_idx, time_points = seq(0, time_horizon)
    )
    pfs_control <- generate_bootstrap_predictions(
      boot_models$control$pfs, sample_idx = sim_idx, time_points = seq(0, time_horizon)
    )
    
    # Update params with these curves
    params$p_os$control_OS <- os_control
    params$p_pfs$control_PFS <- pfs_control
    
    # Generate curves for each biomarker strategy
    for(biomarker in biomarkers) {
      # Get reference to experimental and control treatment levels
      exp_rx <- levels(data$Rx)[2]
      ctrl_rx <- levels(data$Rx)[1]
      
      # Create newdata for each scenario
      newdata_pos_exp <- data.frame(Rx = factor(exp_rx, levels = levels(data$Rx)))
      newdata_pos_exp[[biomarker]] <- 1
      
      newdata_neg_ctrl <- data.frame(Rx = factor(ctrl_rx, levels = levels(data$Rx)))
      newdata_neg_ctrl[[biomarker]] <- 0
      
      # Generate all four predictions for this simulation
      os_pos_exp <- generate_bootstrap_predictions(
        boot_models[[biomarker]]$os, newdata = newdata_pos_exp, 
        sample_idx = sim_idx, time_points = seq(0, time_horizon)
      )
      
      pfs_pos_exp <- generate_bootstrap_predictions(
        boot_models[[biomarker]]$pfs, newdata = newdata_pos_exp, 
        sample_idx = sim_idx, time_points = seq(0, time_horizon)
      )
      
      os_neg_ctrl <- generate_bootstrap_predictions(
        boot_models[[biomarker]]$os, newdata = newdata_neg_ctrl, 
        sample_idx = sim_idx, time_points = seq(0, time_horizon)
      )
      
      pfs_neg_ctrl <- generate_bootstrap_predictions(
        boot_models[[biomarker]]$pfs, newdata = newdata_neg_ctrl, 
        sample_idx = sim_idx, time_points = seq(0, time_horizon)
      )
      
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
  }
  
  # Initialize results dataframe
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
  
  #  control strategy
  # Get control survival curves
  os_control <- params$p_os[["control_OS"]]
  pfs_control <- params$p_pfs[["control_PFS"]]
  
  # Calculate state occupancy
  p_pf_control <- pfs_control
  p_p_control <- pmax(os_control - pfs_control, 0)
  p_d_control <- 1 - os_control
  
  # Force initial state to be 100% progression-free
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
  
  # loop biomarker strategies
  for(biomarker in biomarkers) {
    # Get biomarker prevalence directly from params
    biomarker_prev <- params[[paste0("p_", biomarker)]]
    
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
    
    # Weighted average of costs and effects based on biomarker prevalence
    total_cost <- biomarker_prev * pos_results$costs_total + (1 - biomarker_prev) * neg_results$costs_total
    total_effect <- biomarker_prev * pos_results$qalys_total + (1 - biomarker_prev) * neg_results$qalys_total
    
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