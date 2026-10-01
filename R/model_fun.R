# ===============================================================================
# MAIN MODEL FUNCTION - Partitioned Survival Model with PSA Support
# ===============================================================================
# This function runs the cost-effectiveness model for all strategies.
#
# Supports both deterministic and probabilistic sensitivity analysis (PSA):
#   - Deterministic: Uses base case survival curves and parameters
#   - PSA: Uses sampled survival models (multivariate-normal coefficient draws, issue #156) with subgroup population averaging
#   - Curves: Explicitly evaluates supplied curves, warning and capping PFS > OS.
#
# PSA Methodology (Issues #73, #74, #87):
#   - Parameter uncertainty: Captured by sampled survival models (varies
#     BETWEEN PSA iterations)
#   - Patient heterogeneity: Integrated out via population averaging (averaged
#     WITHIN each iteration)
#   - Population: common joint biomarker weights for every strategy (issue #166)
#
# This separation ensures EVPI correctly measures the value of reducing
# parameter uncertainty, not patient heterogeneity.
#
# Biological Constraint (Issue #76):
#   - Base-case ordering is guaranteed during joint distribution selection.
#   - PFS is capped at OS for PSA and explicit supplied-curve evaluation. Without this safety
#     net, states can sum to >1 when a sampled pair crosses.
#
# Fallback Tracking (Issue #79):
#   - When PSA survival curve generation fails, the function falls back to
#     base case curves. This is now tracked via a "fallback_used" attribute
#     on the returned results, allowing the PSA loop to detect and replace
#     failed iterations instead of mixing different estimation methods.
#
# One Survival Model in Base Case and PSA (Issue #151):
#   - Every PSA curve (control, biomarker-positive, biomarker-negative) comes
#     from the SAME sampled joint model; the control arm is that model
#     predicted with Rx = control, exactly as in the base case. No separate
#     age/sex-only control model is used anywhere.
# ===============================================================================

# ===============================================================================
# INPUT PARAMETER VALIDATION (Issue #47)
# ===============================================================================
# Validates that all required parameters exist and have valid values.
# Provides clear error messages when parameters are misconfigured.
# ===============================================================================

validate_model_params <- function(params, n_points = NULL) {
  v_biomarkers <- get_biomarkers()
  v_biomarker_cost_params <- unique(unname(biomarker_cost_key(v_biomarkers)))
  v_prevalence_params <- unname(biomarker_prevalence_key(v_biomarkers))

  # Required scalar parameters - economic biomarker prevalence values
  v_required_scalars <- c("dr_costs", "dr_effects", "u_np", "u_p",
                        "c_drug_nivo", "c_drug_FLOX", "c_test_CT",
                        "c_test_blood", v_biomarker_cost_params, "c_other_visit",
                        "c_other_baseline", "c_other_follow", "c_other_last",
                        v_prevalence_params)

  # Required list parameters
  v_required_lists <- c("p_os", "p_pfs")

  # Required schedule vectors
  v_required_vectors <- c("l_nivo", "l_FLOX_exp", "l_FLOX_control",
                        "l_CT", "l_blood", "l_visit")

  # Check all required params exist
  v_all_required <- c(v_required_scalars, v_required_lists, v_required_vectors)
  v_missing <- setdiff(v_all_required, names(params))
  if (length(v_missing) > 0) {
    stop("Missing required parameters: ", paste(v_missing, collapse = ", "))
  }

  # Validate cost parameters are non-negative
  v_cost_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                   v_biomarker_cost_params, "c_other_visit", "c_other_baseline",
                   "c_other_follow", "c_other_last")
  for (p in v_cost_params) {
    if (!is.numeric(params[[p]]) || length(params[[p]]) != 1L || !is.finite(params[[p]])) {
      stop("Cost parameter '", p, "' must be a finite numeric scalar.")
    }
    if (params[[p]] < 0) {
      stop("Cost parameter '", p, "' cannot be negative. Got: ", params[[p]])
    }
  }

  # Utilities, prevalences and discount rates must be finite numeric scalars
  # (issue #60); a range test alone admits NA, Inf or vectors.
  require_finite_scalar <- function(p) {
    if (!is.numeric(params[[p]]) || length(params[[p]]) != 1L || !is.finite(params[[p]])) {
      stop("Parameter '", p, "' must be a finite numeric scalar.")
    }
  }

  # Validate utilities are in [0, 1]
  for (p in c("u_np", "u_p")) {
    require_finite_scalar(p)
    if (params[[p]] < 0 || params[[p]] > 1) {
      stop("Utility ", p, " must be in [0,1]. Got: ", params[[p]])
    }
  }

  # Validate prevalence parameters are in [0, 1]
  for (p in v_prevalence_params) {
    require_finite_scalar(p)
    if (params[[p]] < 0 || params[[p]] > 1) {
      stop("Prevalence '", p, "' must be in [0,1]. Got: ", params[[p]])
    }
  }

  # Validate discount rates are non-negative
  for (p in c("dr_costs", "dr_effects")) {
    require_finite_scalar(p)
    if (params[[p]] < 0) {
      stop("Discount rate ", p, " cannot be negative. Got: ", params[[p]])
    }
  }

  # Schedules are 0/1 indicators on the model grid (issues #49, #172). A
  # schedule of another length means it was built for a different horizon.
  if (!is.null(n_points)) {
    for (p in v_required_vectors) {
      v <- params[[p]]
      if (!is.numeric(v) || length(v) != n_points) {
        stop("Schedule '", p, "' has length ", length(v), " but expected ", n_points,
             " (time_horizon + 1); build schedules with build_treatment_schedules().")
      }
      if (any(!is.finite(v)) || any(!v %in% c(0, 1))) {
        stop("Schedule '", p, "' must contain only 0 and 1.")
      }
    }
  }

  invisible(TRUE)
}

# ===============================================================================
# TIME SETTINGS (issue #172)
# ===============================================================================
# The parameter list owns its grid: params$time_horizon (weeks) and params$cl
# (years per cycle). An explicit model_fun() argument must agree with it; when
# params lacks a setting, the argument is required. There is no hidden default,
# so a five-year list cannot be run on the ten-year grid, nor a monthly cl
# replaced by a weekly one.
# ===============================================================================

resolve_time_setting <- function(params, name, supplied) {
  from_params <- params[[name]]
  if (is.null(from_params)) {
    if (is.null(supplied)) {
      stop("model_fun() needs '", name, "': set params$", name,
           " or pass ", name, " explicitly.")
    }
    return(supplied)
  }
  if (!is.null(supplied) &&
      !isTRUE(all.equal(supplied, from_params, tolerance = 1e-12))) {
    stop("Argument ", name, " = ", format(supplied, digits = 10),
         " conflicts with params$", name, " = ", format(from_params, digits = 10), ".")
  }
  from_params
}

model_fun <- function(params, time_horizon = NULL, cl = NULL, determpsa = "det",
                      return_traces = FALSE, sim_idx = NULL,
                      prediction_population = params$prediction_population) {
  determpsa <- match.arg(determpsa, c("det", "psa", "curves"))
  time_horizon <- resolve_time_setting(params, "time_horizon", time_horizon)
  cl <- resolve_time_setting(params, "cl", cl)
  if (!is.numeric(time_horizon) || length(time_horizon) != 1L ||
      !is.finite(time_horizon) || time_horizon < 0 ||
      time_horizon != round(time_horizon)) {
    stop("time_horizon must be a non-negative whole number of cycles. Got: ",
         paste(time_horizon, collapse = ", "))
  }
  if (!is.numeric(cl) || length(cl) != 1L || !is.finite(cl) || cl <= 0) {
    stop("cl must be a positive finite cycle length in years. Got: ",
         paste(cl, collapse = ", "))
  }
  if (determpsa == "psa" && is.null(sim_idx)) {
    stop("PSA mode requires sim_idx; use determpsa = 'curves' to evaluate supplied curves.")
  }
  if (determpsa != "psa" && !is.null(sim_idx)) {
    stop("sim_idx is only used with determpsa = 'psa'.")
  }
  # NOTE: In PSA mode (determpsa = "psa"), this function depends on global variables:
  #   - n_samples: Number of sampled coefficient draws (from 02_setup_and_global_variables.R)
  #   - sampling_models: Sampled survival models (MVN coefficient draws) (from 06_sampling.R)
  # Every prediction uses the explicit prediction_population stored in params.
  # These must exist in the global environment before calling model_fun in PSA mode.

  # Get economic strategies and biomarkers from central config.
  v_strat_names <- get_strategies()
  control_strategy <- get_control_strategy()
  v_biomarkers_to_run <- get_biomarkers()

  # Validate input parameters (Issue #47)
  validate_model_params(params, n_points = time_horizon + 1)

  # Track if fallback to base case was used (Issue #79)
  fallback_used <- FALSE

  v_population_weights <- NULL
  if (!is.null(prediction_population)) {
    v_population_weights <- model_population_weights(params, prediction_population)
    if (!isTRUE(all.equal(v_population_weights, params$population_weights, tolerance = 1e-12))) {
      if (is.null(params$population_curves)) stop("Population changes require patient-level prediction curves.")
      params <- set_population_predictions(params, standardize_population_curves(
        params$population_curves, prediction_population, v_population_weights))
    }
  }

  # =========================================================================
  # PSA MODE: Subgroup Population Averaging
  # =========================================================================
  # Each PSA iteration:
  #   1. Uses sampled survival model i (coefficient draw i) (captures parameter uncertainty)
  #   2. Predicts for ALL patients in each subgroup, then averages
  #   3. Applies the same draw-specific population weights to every strategy
  #
  # This properly separates:
  #   - Parameter uncertainty (varies BETWEEN iterations via sampled models)
  #   - Patient heterogeneity (averaged WITHIN each iteration)
  #
  # See GitHub Issues #73, #74, #87 for methodology discussion.
  # =========================================================================
  if(determpsa == "psa") {

    # Validate sim_idx (Issue #46)
    if (!is.numeric(sim_idx) || length(sim_idx) != 1 || !is.finite(sim_idx) || sim_idx < 1 ||
        sim_idx != floor(sim_idx)) {
      stop("sim_idx must be a positive integer. Got: ", sim_idx)
    }
    if (sim_idx > n_samples) {
      stop("sim_idx (", sim_idx, ") exceeds n_samples (", n_samples, ")")
    }

    if (is.null(prediction_population)) stop("PSA requires an explicit prediction_population.")

    valid_psa_curves <- function(...) {
      all(vapply(list(...), function(curve) {
        !is.null(curve) && length(curve) == time_horizon + 1 &&
          is.null(survival_curve_problem(curve))
      }, logical(1)))
    }

    if (sim_idx %% 100 == 0) {
      cat("PSA iteration", sim_idx, "- subgroup population averaging\n")
    }

    # One joint component serves the control arm and every biomarker strategy
    # (issues #151, #156). A draw that cannot be materialised, or a legacy
    # bootstrap sample flagged as failed (original fit substituted), is
    # recorded as fallback so the PSA loop replaces it.
    joint_sampling_models <- get_joint_sampling_models(sampling_models)
    joint_sample <- tryCatch(
      sampled_survival_models(joint_sampling_models, sim_idx),
      error = function(e) NULL
    )
    if (is.null(joint_sample) || isTRUE(joint_sample$failed)) {
      fallback_used <- TRUE
    }
    rm(joint_sample)

    # -----------------------------------------------------------------------
    # CONTROL STRATEGY: Population-averaged predictions over FULL population
    # -----------------------------------------------------------------------
    # Uses the shared population-averaging helper to predict for ALL patients in
    # the explicit prediction population (both arms) with Rx forced to the control level, using the
    # sampled JOINT model. This mirrors the base case, which predicts the
    # control curve from the same joint fit with Rx = control
    # (generate_population_averaged_predictions). A separate age/sex-only
    # control model was previously bootstrapped here; that shifted the PSA
    # away from the base case and broke the control/biomarker correlation
    # (issue #151).

    tryCatch({
      v_os_control <- generate_psa_population_averaged_predictions(
        sampling_model_list = joint_sampling_models,
        outcome = "os",
        sample_idx = sim_idx,
        data_original = prediction_population,
        weights = v_population_weights,
        time_points = seq(0, time_horizon)
      )

      v_pfs_control <- generate_psa_population_averaged_predictions(
        sampling_model_list = joint_sampling_models,
        outcome = "pfs",
        sample_idx = sim_idx,
        data_original = prediction_population,
        weights = v_population_weights,
        time_points = seq(0, time_horizon)
      )

      # A missing or invalid prediction (issue #60: not finite, outside
      # [0, 1], not starting at 1, or increasing) is a failed draw.
      if (!valid_psa_curves(v_os_control, v_pfs_control)) {
        warning("Missing or invalid control survival curves for sim ", sim_idx,
                " - using base case curves")
        fallback_used <- TRUE  # Issue #79: Track fallback
      } else {
        # Update params with population-averaged predictions
        params$p_os[[paste0(control_strategy, "_OS")]] <- v_os_control
        params$p_pfs[[paste0(control_strategy, "_PFS")]] <- v_pfs_control
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
    # Weighted combination uses the marginals of the common target population

    for(biomarker in v_biomarkers_to_run) {
      tryCatch({
        # Get population-averaged predictions for BOTH subgroups
        preds_os <- generate_psa_population_averaged_predictions(
          sampling_model_list = joint_sampling_models,
          biomarker_name = biomarker,
          outcome = "os",
          sample_idx = sim_idx,
          data_original = prediction_population,
          weights = v_population_weights,
          time_points = seq(0, time_horizon)
        )

        preds_pfs <- generate_psa_population_averaged_predictions(
          sampling_model_list = joint_sampling_models,
          biomarker_name = biomarker,
          outcome = "pfs",
          sample_idx = sim_idx,
          data_original = prediction_population,
          weights = v_population_weights,
          time_points = seq(0, time_horizon)
        )

        # A missing or invalid prediction is a failed draw (issue #60).
        if (is.null(preds_os) || is.null(preds_pfs) ||
            !valid_psa_curves(preds_os$positive, preds_os$negative,
                              preds_pfs$positive, preds_pfs$negative)) {
          warning("Missing or invalid ", biomarker, " survival curves for sim ", sim_idx,
                  " - using base case curves")
          fallback_used <- TRUE  # Issue #79: Track fallback
        } else {
          # Store subgroup-specific curves
          params$p_os[[paste0(biomarker, "_pos_OS")]] <- preds_os$positive
          params$p_os[[paste0(biomarker, "_neg_OS")]] <- preds_os$negative
          params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- preds_pfs$positive
          params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- preds_pfs$negative

          # Marginal prevalence comes from the common target population weights
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
  
  df_results <- data.frame(
    Strategy = as.character(v_strat_names),
    Cost = numeric(length(v_strat_names)),
    Effect = numeric(length(v_strat_names)),
    stringsAsFactors = FALSE
  )
  
  traces <- list()
  if(return_traces) {
    for(strat in v_strat_names) {
      traces[[strat]] <- list()
    }
  }
  
  expected_length <- time_horizon + 1
  # Discount on this call's cycle length, the same cl calculate_outcomes() uses
  # for QALYs, rather than params$cl (which test fixtures may omit).
  weights <- discount_weights(params, expected_length, cl)
  v_dw_c <- weights$cost
  v_dw_e <- weights$effect

  # Validate a survival curve: length first, then the survival
  # contract of calculate_outcomes.R (issue #60: finite, within [0, 1],
  # starting at 1, nonincreasing), in every mode.
  validate_curve <- function(curve, name) {
    if (length(curve) != expected_length) {
      stop("Survival curve '", name, "' has length ", length(curve),
           " but expected ", expected_length, " (time_horizon + 1)")
    }
    validate_survival_curve(curve, name)
  }

  # Ordered base-case curves are guaranteed during distribution selection.
  # Resampled PSA fits are selected upstream and may still cross, so retain the
  # clamp for PSA and explicitly supplied curves, and emit a structured warning
  # that psa_functions.R can aggregate.
  enforce_survival_ordering <- function(pfs, os, strategy, subgroup) {
    v_violation_idx <- which(pfs > os)
    if (length(v_violation_idx) == 0) {
      return(pfs)
    }

    display_name <- unname(strategy_display_name(strategy))
    curve_label <- if (subgroup == "control") {
      tolower(display_name)
    } else {
      paste0(strategy, if (subgroup == "positive") "+" else "-")
    }

    if (determpsa == "det") {
      stop(
        "Base-case survival ordering invariant failed for ", curve_label,
        ": PFS exceeded OS at ", length(v_violation_idx), " time points. ",
        "Refit using ordering-constrained distribution selection."
      )
    }

    max_excess <- max(pfs[v_violation_idx] - os[v_violation_idx])
    warning(warningCondition(
      paste0(
        "PFS > OS constraint enforced at ", length(v_violation_idx),
        " time points (", curve_label, ")",
        if (!is.null(sim_idx)) paste0(" [sim ", sim_idx, "]") else "",
        "; max excess=", signif(max_excess, 4)
      ),
      class = "survival_ordering_warning",
      strategy = strategy,
      subgroup = subgroup,
      n_violations = length(v_violation_idx),
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
  v_os_control <- params$p_os[[control_os_key]]
  v_pfs_control <- params$p_pfs[[control_pfs_key]]
  validate_curve(v_os_control, control_os_key)
  validate_curve(v_pfs_control, control_pfs_key)

  v_pfs_control <- enforce_survival_ordering(
    v_pfs_control, v_os_control, control_strategy, "control"
  )

  control_states <- partitioned_survival_states(v_os_control, v_pfs_control)
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
  df_results$Cost[df_results$Strategy == control_strategy] <- control_results$costs_total
  df_results$Effect[df_results$Strategy == control_strategy] <- control_results$qalys_total
  
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
  for(biomarker in v_biomarkers_to_run) {
    # Get biomarker prevalence directly from params
    biomarker_prev <- params[[paste0("p_", biomarker)]]
    
    # -----------------------------------------------------------------------
    # BIOMARKER POSITIVE SUBGROUP (receives experimental treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker positive survival curves
    v_os_pos <- params$p_os[[paste0(biomarker, "_pos_OS")]]
    v_pfs_pos <- params$p_pfs[[paste0(biomarker, "_pos_PFS")]]
    validate_curve(v_os_pos, paste0(biomarker, "_pos_OS"))
    validate_curve(v_pfs_pos, paste0(biomarker, "_pos_PFS"))

    v_pfs_pos <- enforce_survival_ordering(
      v_pfs_pos, v_os_pos, biomarker, "positive"
    )

    pos_states <- partitioned_survival_states(v_os_pos, v_pfs_pos)
    p_pf_pos <- pos_states$p_pf
    p_p_pos <- pos_states$p_p
    p_d_pos <- pos_states$p_d
    
    # -----------------------------------------------------------------------
    # BIOMARKER NEGATIVE SUBGROUP (receives standard treatment)
    # -----------------------------------------------------------------------
    
    # Get biomarker negative survival curves
    v_os_neg <- params$p_os[[paste0(biomarker, "_neg_OS")]]
    v_pfs_neg <- params$p_pfs[[paste0(biomarker, "_neg_PFS")]]
    validate_curve(v_os_neg, paste0(biomarker, "_neg_OS"))
    validate_curve(v_pfs_neg, paste0(biomarker, "_neg_PFS"))

    v_pfs_neg <- enforce_survival_ordering(
      v_pfs_neg, v_os_neg, biomarker, "negative"
    )

    neg_states <- partitioned_survival_states(v_os_neg, v_pfs_neg)
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
    df_results$Cost[df_results$Strategy == biomarker] <- total_cost
    df_results$Effect[df_results$Strategy == biomarker] <- total_effect
    
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
      results = df_results,
      traces = traces
    )
  } else {
    result <- df_results
  }

  # Add fallback_used attribute so PSA loop can detect and replace failed iterations
  attr(result, "fallback_used") <- fallback_used
  return(result)
}
