# Survival-curve contract (issue #60). Tolerance for numerical noise in
# population-averaged predictions; a genuine violation is far larger.
SURVIVAL_CURVE_TOL <- 1e-10

#' Describe the first violation of the survival-curve contract
#'
#' A survival curve on the model grid must be numeric, finite, within [0, 1],
#' start at 1 and be nonincreasing (all within SURVIVAL_CURVE_TOL).
#'
#' @param curve Survival probabilities at the model time points.
#' @param tol Numerical tolerance.
#' @return NULL for a valid curve, otherwise a phrase describing the problem.
survival_curve_problem <- function(curve, tol = SURVIVAL_CURVE_TOL) {
  if (!is.numeric(curve) || length(curve) == 0L) {
    return("is not a non-empty numeric vector")
  }
  if (any(!is.finite(curve))) {
    return(paste0("contains non-finite values at ", sum(!is.finite(curve)), " time points"))
  }
  if (any(curve < -tol | curve > 1 + tol)) {
    return(paste0("leaves [0, 1] (range ", signif(min(curve), 6), " to ",
                  signif(max(curve), 6), ")"))
  }
  if (abs(curve[1] - 1) > tol) {
    return(paste0("starts at ", signif(curve[1], 6), ", not 1"))
  }
  v_rises <- diff(curve) > tol
  if (any(v_rises)) {
    return(paste0("increases at ", sum(v_rises), " time points (largest rise ",
                  signif(max(diff(curve)), 6), ")"))
  }
  NULL
}

#' Stop unless a survival curve satisfies the contract
#'
#' @param curve Survival probabilities.
#' @param name Curve name used in the error message.
validate_survival_curve <- function(curve, name) {
  problem <- survival_curve_problem(curve)
  if (!is.null(problem)) stop("Survival curve '", name, "' ", problem, ".")
  invisible(TRUE)
}

#' Build partitioned-survival state occupancy vectors
#'
#' Both curves must satisfy the survival-curve contract; the first point is set
#' to exactly 1 only after it has been checked to equal 1 within tolerance
#' (issue #60), so a curve that does not start at 1 is rejected, not repaired.
#'
#' @param os Overall-survival probabilities
#' @param pfs Progression-free-survival probabilities
#' @param curve_label Label identifying the curve in an ordering-violation error.
#'   When supplied, OS >= PFS is asserted. Callers that have already run
#'   model_fun()'s enforce_survival_ordering() leave this NULL.
#' @return List containing progression-free, progressed, and dead occupancy
partitioned_survival_states <- function(os, pfs, curve_label = NULL) {
  label <- if (is.null(curve_label)) "" else paste0(curve_label, " ")
  validate_survival_curve(os, paste0(label, "OS"))
  validate_survival_curve(pfs, paste0(label, "PFS"))
  if (length(os) != length(pfs)) {
    stop("OS and PFS curves differ in length (", length(os), " vs ", length(pfs), ").")
  }
  if (!is.null(curve_label)) {
    v_violation_idx <- which(pfs > os)
    if (length(v_violation_idx) > 0) {
      stop(
        "Survival ordering invariant failed for ", curve_label,
        ": PFS exceeded OS at ", length(v_violation_idx), " time points",
        "; max excess=",
        signif(max(pfs[v_violation_idx] - os[v_violation_idx]), 4),
        ". Refit using ordering-constrained distribution selection."
      )
    }
  }
  states <- list(p_pf = pfs, p_p = pmax(os - pfs, 0), p_d = 1 - os)
  states$p_pf[1] <- 1
  states$p_p[1] <- states$p_d[1] <- 0
  states
}

#' Create cost and effect discount weights for model cycles
#'
#' @param params Parameter list supplying dr_costs and dr_effects.
#' @param n_cycles Number of model cycles (time_horizon + 1).
#' @param cl Cycle length in years; defaults to params$cl. Discounting must use
#'   the model's actual cycle length rather than assuming weekly cycles, so that
#'   changing cl rescales the horizon correctly instead of silently discounting
#'   over the wrong number of years.
#' @return List with cost and effect discount weight vectors.
discount_weights <- function(params, n_cycles, cl = params$cl) {
  if (is.null(cl) || !is.numeric(cl) || length(cl) != 1L ||
      !is.finite(cl) || cl <= 0) {
    stop("discount_weights() requires a positive numeric cycle length 'cl'; ",
         "pass it explicitly or set params$cl.")
  }
  v_years <- (seq_len(n_cycles) - 1) * cl
  list(cost = 1 / (1 + params$dr_costs)^v_years,
       effect = 1 / (1 + params$dr_effects)^v_years)
}

# Trapezoidal integration over intervals between the supplied grid points.
# A single point spans no time. Scheduled/event costs do not use these weights.
interval_weights <- function(n_points) {
  if (n_points < 2L) return(rep(0, n_points))
  v_weights <- rep(1, n_points)
  v_weights[c(1L, n_points)] <- 0.5
  v_weights
}

# Helper function to calculate costs and QALYs based on state occupancy
calculate_outcomes <- function(params, p_pf, p_p, p_d, treatment_type, biomarker, v_dw_c, v_dw_e, cl) {
  # Number of cycles
  n_cycles <- length(p_pf)
  
  # Extend treatment schedules to match time horizon
  extend_schedule <- function(schedule, length_needed) {
    if (length(schedule) < length_needed) {
      v_extended <- c(schedule, rep(0, length_needed - length(schedule)))
      return(v_extended)
    }
    return(schedule[1:length_needed])
  }
  
  # Extend all treatment schedules
  l_nivo <- extend_schedule(params$l_nivo, n_cycles)
  l_FLOX_exp <- extend_schedule(params$l_FLOX_exp, n_cycles)
  l_FLOX_control <- extend_schedule(params$l_FLOX_control, n_cycles)
  l_CT <- extend_schedule(params$l_CT, n_cycles)
  l_blood <- extend_schedule(params$l_blood, n_cycles)
  l_visit <- extend_schedule(params$l_visit, n_cycles)
  
  v_occupancy_weights <- interval_weights(n_cycles)
  
  # Calculate QALYs
  v_qalys_pf <- p_pf * params$u_np * cl * v_occupancy_weights
  v_qalys_p <- p_p * params$u_p * cl * v_occupancy_weights
  v_qalys_undiscounted <- v_qalys_pf + v_qalys_p
  v_qalys_discounted <- v_qalys_undiscounted * v_dw_e
  qalys_total <- sum(v_qalys_discounted)
  
  # Calculate costs based on treatment type
  if(treatment_type == "experimental") {
    # For biomarker positive: nivolumab + FLOX experimental regimen
    v_drug_costs <- (l_nivo * params$c_drug_nivo + l_FLOX_exp * params$c_drug_FLOX) * p_pf
  } else {
    # For control and biomarker negative: standard FLOX regimen
    v_drug_costs <- l_FLOX_control * params$c_drug_FLOX * p_pf
  }
  
  # Test costs (standard for all patients)
  v_test_costs <- l_CT * params$c_test_CT * p_pf + l_blood * params$c_test_blood * p_pf
  
  # Add biomarker test cost if applicable
  if(!is.null(biomarker)) {
    # Scalar prices are canonical, including for direct/enriched callers.
    # A legacy c_test_biomarker lookup is never used (issue #173).
    biomarker_test_cost <- params[[biomarker_cost_key(biomarker)[[1L]]]]
    if (!is.numeric(biomarker_test_cost) || length(biomarker_test_cost) != 1L ||
        !is.finite(biomarker_test_cost) || biomarker_test_cost < 0) {
      stop("No diagnostic-test cost configured for biomarker '", biomarker, "'")
    }
    v_test_costs[1] <- v_test_costs[1] + biomarker_test_cost
  }
  
  # Visit costs
  v_visit_costs <- l_visit * params$c_other_visit * p_pf
  v_visit_costs[1] <- v_visit_costs[1] + params$c_other_baseline
  
  # Ongoing progressed-state costs are rates per quarter (four per year).
  # Integrate occupancy and discounting over intervals, just as for utilities.
  v_progressed_quarters <- p_p * (4 * cl) * v_occupancy_weights
  v_follow_up_costs <- params$c_other_follow * v_progressed_quarters
  # Progressed patients are imaged once per quarter (issue #164), charged at the
  # CT unit price as a rate over progressed occupancy like the follow-up visit.
  # Progression-free patients receive their CT scans from the l_CT schedule.
  v_progressed_imaging_costs <- params$c_test_CT * v_progressed_quarters

  # Post-progression treatment costs (issue #154). Second-line systemic therapy
  # is NOT costed in the base case, where c_other_pp = 0 and the progressed
  # state accrues only the quarterly follow-up visit, the quarterly CT and the
  # end-of-life cost. The omission is differential because the
  # strategies differ in time spent progressed, so the parameter is explicit and
  # is varied in a deterministic structural scenario, as a quarterly rate.
  c_other_pp <- if (is.null(params$c_other_pp)) 0 else params$c_other_pp
  v_post_progression_costs <- c_other_pp * v_progressed_quarters

  # End-of-life costs: one-time cost applied when patients transition to death
  # At t=0: no deaths yet (everyone starts alive)
  # At t>0: incremental deaths from previous time point
  # pmax ensures non-negative values for numerical stability in PSA runs
  v_death_transitions <- c(0, diff(p_d))
  v_death_transitions <- pmax(v_death_transitions, 0)
  v_end_life_costs <- v_death_transitions * params$c_other_last
  
  # Total costs
  v_total_costs_undiscounted <- v_drug_costs + v_test_costs + v_visit_costs +
    v_follow_up_costs + v_progressed_imaging_costs + v_post_progression_costs +
    v_end_life_costs
  v_total_costs_discounted <- v_total_costs_undiscounted * v_dw_c
  costs_total <- sum(v_total_costs_discounted)
  
  return(list(
    qalys_pf = v_qalys_pf,
    qalys_p = v_qalys_p,
    qalys_undiscounted = v_qalys_undiscounted,
    qalys_discounted = v_qalys_discounted,
    qalys_total = qalys_total,
    drug_costs = v_drug_costs,
    test_costs = v_test_costs,
    visit_costs = v_visit_costs,
    follow_up_costs = v_follow_up_costs,
    progressed_imaging_costs = v_progressed_imaging_costs,
    post_progression_costs = v_post_progression_costs,
    end_life_costs = v_end_life_costs,
    total_costs_undiscounted = v_total_costs_undiscounted,
    total_costs_discounted = v_total_costs_discounted,
    costs_total = costs_total
  ))
}
