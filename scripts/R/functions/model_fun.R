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
#   - Base-case ordering is guaranteed during joint distribution selection.
#   - PFS is capped at OS only for resampled PSA models. Without this safety
#     net, states can sum to >1 when a resampled pair crosses.
#
# Fallback Tracking (Issue #79):
#   - When PSA survival curve generation fails, the function falls back to
#     base case curves. This is now tracked via a "fallback_used" attribute
#     on the returned results, allowing the PSA loop to detect and replace
#     failed iterations instead of mixing different estimation methods.
#
# One Survival Model in Base Case and PSA (Issue #151):
#   - Every PSA curve (control, biomarker-positive, biomarker-negative) comes
#     from the SAME resampled joint model; the control arm is that model
#     predicted with Rx = control, exactly as in the base case. No separate
#     age/sex-only control model is used anywhere.
# ===============================================================================

# ===============================================================================
# INPUT PARAMETER VALIDATION (Issue #47)
# ===============================================================================
# Validates that all required parameters exist and have valid values.
# Provides clear error messages when parameters are misconfigured.
# ===============================================================================

validate_model_params <- function(params) {
  biomarkers <- get_biomarkers()
  biomarker_cost_params <- unique(unname(biomarker_cost_key(biomarkers)))
  prevalence_params <- unname(biomarker_prevalence_key(biomarkers))

  # Required scalar parameters - economic biomarker prevalence values
  required_scalars <- c("dr_costs", "dr_effects", "u_np", "u_p",
                        "c_drug_nivo", "c_drug_FLOX", "c_test_CT",
                        "c_test_blood", biomarker_cost_params, "c_other_visit",
                        "c_other_baseline", "c_other_follow", "c_other_last",
                        prevalence_params)

  # Required list parameters
  required_lists <- c("p_os", "p_pfs", "c_test_biomarker")

  # Required schedule vectors
  required_vectors <- c("l_nivo", "l_FLOX_exp", "l_FLOX_control",
                        "l_CT", "l_blood", "l_visit")

  # Check all required params exist
  all_required <- c(required_scalars, required_lists, required_vectors)
  missing <- setdiff(all_required, names(params))
  if (length(missing) > 0) {
    stop("Missing required parameters: ", paste(missing, collapse = ", "))
  }

  # Validate cost parameters are non-negative
  cost_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                   biomarker_cost_params, "c_other_visit", "c_other_baseline",
                   "c_other_follow", "c_other_last")
  for (p in cost_params) {
    if (params[[p]] < 0) {
      stop("Cost parameter '", p, "' cannot be negative. Got: ", params[[p]])
    }
  }

  # Validate the data-driven diagnostic-cost mapping used by each biomarker strategy.
  missing_biomarker_costs <- setdiff(biomarkers,
                                     names(params$c_test_biomarker))
  if (length(missing_biomarker_costs) > 0) {
    stop("Missing diagnostic-test costs for biomarkers: ",
         paste(missing_biomarker_costs, collapse = ", "))
  }
  invalid_biomarker_costs <- vapply(
    params$c_test_biomarker[biomarkers],
    function(x) !is.numeric(x) || length(x) != 1 || is.na(x) || x < 0,
    logical(1)
  )
  if (any(invalid_biomarker_costs)) {
    stop("Diagnostic-test costs must be non-negative numeric scalars for: ",
         paste(biomarkers[invalid_biomarker_costs], collapse = ", "))
  }

  # Validate utilities are in [0, 1]
  if (params$u_np < 0 || params$u_np > 1) {
    stop("Utility u_np must be in [0,1]. Got: ", params$u_np)
  }
  if (params$u_p < 0 || params$u_p > 1) {
    stop("Utility u_p must be in [0,1]. Got: ", params$u_p)
  }

  # Validate prevalence parameters are in [0, 1]
  for (p in prevalence_params) {
    if (params[[p]] < 0 || params[[p]] > 1) {
      stop("Prevalence '", p, "' must be in [0,1]. Got: ", params[[p]])
    }
  }

  # Validate discount rates are non-negative
  if (params$dr_costs < 0) {
    stop("Discount rate dr_costs cannot be negative. Got: ", params$dr_costs)
  }
  if (params$dr_effects < 0) {
    stop("Discount rate dr_effects cannot be negative. Got: ", params$dr_effects)
  }

  invisible(TRUE)
}

model_fun <- function(params, time_horizon = 520, cl = 1/52, determpsa = "det",
                      return_traces = FALSE, sim_idx = NULL) {
  # NOTE: In PSA mode (determpsa = "psa"), this function depends on global variables:
  #   - n_samples: Number of resampled models (from 02_setup_and_global_variables.R)
  #   - sampling_models: Resampled survival models (from 06_sampling.R)
  #   - data_complete: Full analysis cohort (from 04_parametric_survival_analysis.R)
  #   - data: Full dataset for biomarker predictions (from 03_biomarker_strategies.R)
  # These must exist in the global environment before calling model_fun in PSA mode.

  # Get economic strategies and biomarkers from central config.
  strat_names <- get_strategies()
  control_strategy <- get_control_strategy()
  biomarkers_to_run <- get_biomarkers()

  # Validate input parameters (Issue #47)
  validate_model_params(params)

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

    # Check if any resampled models for this sim_idx are failed (using originals)
    # This helps track iterations where resampled model fitting failed and the
    # original model was substituted, which reduces uncertainty estimation.
    # The control arm is predicted from the same joint bootstrap as the
    # biomarker strategies, so the biomarker components cover it (issue #151).
    biomarker_failed <- any(vapply(biomarkers_to_run, function(bm) {
      isTRUE(sampling_models[[bm]]$samples[[sim_idx]]$failed)
    }, logical(1)))
    if (biomarker_failed) {
      fallback_used <- TRUE
    }

    # -----------------------------------------------------------------------
    # CONTROL STRATEGY: Population-averaged predictions over FULL population
    # -----------------------------------------------------------------------
    # Uses the shared population-averaging helper to predict for ALL patients in
    # data_complete (both arms) with Rx forced to the control level, using the
    # resampled JOINT model. This mirrors the base case, which predicts the
    # control curve from the same joint fit with Rx = control
    # (generate_population_averaged_predictions). A separate age/sex-only
    # control model was previously bootstrapped here; that shifted the PSA
    # away from the base case and broke the control/biomarker correlation
    # (issue #151).

    tryCatch({
      joint_sampling_models <- get_joint_sampling_models(sampling_models)

      os_control <- generate_psa_population_averaged_predictions(
        sampling_model_list = joint_sampling_models,
        outcome = "os",
        sample_idx = sim_idx,
        data_original = data_complete,
        time_points = seq(0, time_horizon)
      )

      pfs_control <- generate_psa_population_averaged_predictions(
        sampling_model_list = joint_sampling_models,
        outcome = "pfs",
        sample_idx = sim_idx,
        data_original = data_complete,
        time_points = seq(0, time_horizon)
      )

      # Check for NA values and handle
      if (is.null(os_control) || is.null(pfs_control) ||
          any(is.na(os_control)) || any(is.na(pfs_control))) {
        warning("NA values in control survival curves for sim ", sim_idx,
                " - using base case curves")
        fallback_used <- TRUE  # Issue #79: Track fallback
      } else {
        # Update params with population-averaged predictions
        params$p_os[[paste0(control_strategy, "_OS")]] <- os_control
        params$p_pfs[[paste0(control_strategy, "_PFS")]] <- pfs_control
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

    for(biomarker in biomarkers_to_run) {
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
  
  expected_length <- time_horizon + 1
  # Discount on this call's cycle length, the same cl calculate_outcomes() uses
  # for QALYs, rather than params$cl (which test fixtures may omit).
  weights <- discount_weights(params, expected_length, cl)
  v_dw_c <- weights$cost
  v_dw_e <- weights$effect

  # Helper to validate survival curve length
  validate_curve_length <- function(curve, name) {
    if (length(curve) != expected_length) {
      stop("Survival curve '", name, "' has length ", length(curve),
           " but expected ", expected_length, " (time_horizon + 1)")
    }
  }

  # Ordered base-case curves are guaranteed during distribution selection.
  # Resampled PSA fits are selected upstream and may still cross, so retain the
  # clamp only there as a last-resort safety net and emit a structured warning
  # that psa_functions.R can aggregate.
  enforce_survival_ordering <- function(pfs, os, strategy, subgroup) {
    violation_idx <- which(pfs > os)
    if (length(violation_idx) == 0) {
      return(pfs)
    }

    display_name <- unname(strategy_display_name(strategy))
    curve_label <- if (subgroup == "control") {
      tolower(display_name)
    } else {
      paste0(strategy, if (subgroup == "positive") "+" else "-")
    }

    if (determpsa != "psa") {
      stop(
        "Base-case survival ordering invariant failed for ", curve_label,
        ": PFS exceeded OS at ", length(violation_idx), " time points. ",
        "Refit using ordering-constrained distribution selection."
      )
    }

    max_excess <- max(pfs[violation_idx] - os[violation_idx])
    warning(warningCondition(
      paste0(
        "PFS > OS constraint enforced at ", length(violation_idx),
        " time points (", curve_label, ")",
        if (!is.null(sim_idx)) paste0(" [sim ", sim_idx, "]") else "",
        "; max excess=", signif(max_excess, 4)
      ),
      class = "survival_ordering_warning",
      strategy = strategy,
      subgroup = subgroup,
      n_violations = length(violation_idx),
      max_excess = max_excess,
      sim_idx = sim_idx
    ))
    pmin(pfs, os)
  }

  # =========================================================================
  # CONTROL STRATEGY
  # =========================================================================

  # Get control survival curves (either base case or from resampling)
  control_os_key <- paste0(control_strategy, "_OS")
  control_pfs_key <- paste0(control_strategy, "_PFS")
  os_control <- params$p_os[[control_os_key]]
  pfs_control <- params$p_pfs[[control_pfs_key]]
  validate_curve_length(os_control, control_os_key)
  validate_curve_length(pfs_control, control_pfs_key)

  pfs_control <- enforce_survival_ordering(
    pfs_control, os_control, control_strategy, "control"
  )

  control_states <- partitioned_survival_states(os_control, pfs_control)
  p_pf_control <- control_states$p_pf
  p_p_control <- control_states$p_p
  p_d_control <- control_states$p_d
  
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
  results$Cost[results$Strategy == control_strategy] <- control_results$costs_total
  results$Effect[results$Strategy == control_strategy] <- control_results$qalys_total
  
  # Store traces for control if requested
  if(return_traces) {
    traces[[control_strategy]] <- list(
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
  for(biomarker in biomarkers_to_run) {
    # Get biomarker prevalence directly from params
    biomarker_prev <- params[[paste0("p_", biomarker)]]
    
    # -----------------------------------------------------------------------
    # BIOMARKER POSITIVE SUBGROUP (receives experimental treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker positive survival curves
    os_pos <- params$p_os[[paste0(biomarker, "_pos_OS")]]
    pfs_pos <- params$p_pfs[[paste0(biomarker, "_pos_PFS")]]
    validate_curve_length(os_pos, paste0(biomarker, "_pos_OS"))
    validate_curve_length(pfs_pos, paste0(biomarker, "_pos_PFS"))

    pfs_pos <- enforce_survival_ordering(
      pfs_pos, os_pos, biomarker, "positive"
    )

    pos_states <- partitioned_survival_states(os_pos, pfs_pos)
    p_pf_pos <- pos_states$p_pf
    p_p_pos <- pos_states$p_p
    p_d_pos <- pos_states$p_d
    
    # -----------------------------------------------------------------------
    # BIOMARKER NEGATIVE SUBGROUP (receives standard treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker negative survival curves
    os_neg <- params$p_os[[paste0(biomarker, "_neg_OS")]]
    pfs_neg <- params$p_pfs[[paste0(biomarker, "_neg_PFS")]]
    validate_curve_length(os_neg, paste0(biomarker, "_neg_OS"))
    validate_curve_length(pfs_neg, paste0(biomarker, "_neg_PFS"))

    pfs_neg <- enforce_survival_ordering(
      pfs_neg, os_neg, biomarker, "negative"
    )

    neg_states <- partitioned_survival_states(os_neg, pfs_neg)
    p_pf_neg <- neg_states$p_pf
    p_p_neg <- neg_states$p_p
    p_d_neg <- neg_states$p_d
    
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
