#' Build partitioned-survival state occupancy vectors
#'
#' @param os Overall-survival probabilities
#' @param pfs Progression-free-survival probabilities
#' @param curve_label Label identifying the curve in an ordering-violation error.
#'   When supplied, OS >= PFS is asserted. Callers that have already run
#'   model_fun()'s enforce_survival_ordering() leave this NULL.
#' @return List containing progression-free, progressed, and dead occupancy
partitioned_survival_states <- function(os, pfs, curve_label = NULL) {
  if (!is.null(curve_label)) {
    violation_idx <- which(pfs > os)
    if (length(violation_idx) > 0) {
      stop(
        "Survival ordering invariant failed for ", curve_label,
        ": PFS exceeded OS at ", length(violation_idx), " time points",
        "; max excess=",
        signif(max(pfs[violation_idx] - os[violation_idx]), 4),
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
  years <- (seq_len(n_cycles) - 1) * cl
  list(cost = 1 / (1 + params$dr_costs)^years,
       effect = 1 / (1 + params$dr_effects)^years)
}

# Helper function to calculate costs and QALYs based on state occupancy
calculate_outcomes <- function(params, p_pf, p_p, p_d, treatment_type, biomarker, v_dw_c, v_dw_e, cl) {
  # Number of cycles
  n_cycles <- length(p_pf)
  
  # Extend treatment schedules to match time horizon
  extend_schedule <- function(schedule, length_needed) {
    if (length(schedule) < length_needed) {
      extended <- c(schedule, rep(0, length_needed - length(schedule)))
      return(extended)
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
  
  # Define quarterly cycles for follow-up costs
  quarterly_cycles <- seq(13, n_cycles, by = 13)
  
  # Calculate QALYs
  qalys_pf <- p_pf * params$u_np * cl
  qalys_p <- p_p * params$u_p * cl
  qalys_undiscounted <- qalys_pf + qalys_p
  qalys_discounted <- qalys_undiscounted * v_dw_e
  qalys_total <- sum(qalys_discounted)
  
  # Calculate costs based on treatment type
  if(treatment_type == "experimental") {
    # For biomarker positive: nivolumab + FLOX experimental regimen
    drug_costs <- (l_nivo * params$c_drug_nivo + l_FLOX_exp * params$c_drug_FLOX) * p_pf
  } else {
    # For control and biomarker negative: standard FLOX regimen
    drug_costs <- l_FLOX_control * params$c_drug_FLOX * p_pf
  }
  
  # Test costs (standard for all patients)
  test_costs <- l_CT * params$c_test_CT * p_pf + l_blood * params$c_test_blood * p_pf
  
  # Add biomarker test cost if applicable
  if(!is.null(biomarker)) {
    biomarker_test_cost <- params$c_test_biomarker[[biomarker]]
    if (is.null(biomarker_test_cost)) {
      stop("No diagnostic-test cost configured for biomarker '", biomarker, "'")
    }
    test_costs[1] <- test_costs[1] + biomarker_test_cost
  }
  
  # Visit costs
  visit_costs <- l_visit * params$c_other_visit * p_pf
  visit_costs[1] <- visit_costs[1] + params$c_other_baseline
  
  # Follow-up costs
  follow_up_costs <- rep(0, n_cycles)
  if(length(quarterly_cycles) > 0) {
    follow_up_costs[quarterly_cycles] <- params$c_other_follow * p_p[quarterly_cycles]
  }

  # Post-progression treatment costs (issue #154). Second-line systemic therapy,
  # imaging and visits after progression are NOT costed in the base case, where
  # c_other_pp = 0 and the progressed state accrues only the quarterly follow-up
  # contact and the end-of-life cost. The omission is differential because the
  # strategies differ in time spent progressed, so the parameter is explicit and
  # is varied in a deterministic structural scenario. Charged on the same
  # quarterly cycles as the follow-up contact.
  c_other_pp <- if (is.null(params$c_other_pp)) 0 else params$c_other_pp
  post_progression_costs <- rep(0, n_cycles)
  if (c_other_pp != 0 && length(quarterly_cycles) > 0) {
    post_progression_costs[quarterly_cycles] <- c_other_pp * p_p[quarterly_cycles]
  }

  # End-of-life costs: one-time cost applied when patients transition to death
  # At t=0: no deaths yet (everyone starts alive)
  # At t>0: incremental deaths from previous time point
  # pmax ensures non-negative values for numerical stability in PSA runs
  death_transitions <- c(0, diff(p_d))
  death_transitions <- pmax(death_transitions, 0)
  end_life_costs <- death_transitions * params$c_other_last
  
  # Total costs
  total_costs_undiscounted <- drug_costs + test_costs + visit_costs +
    follow_up_costs + post_progression_costs + end_life_costs
  total_costs_discounted <- total_costs_undiscounted * v_dw_c
  costs_total <- sum(total_costs_discounted)
  
  return(list(
    qalys_pf = qalys_pf,
    qalys_p = qalys_p,
    qalys_undiscounted = qalys_undiscounted,
    qalys_discounted = qalys_discounted,
    qalys_total = qalys_total,
    drug_costs = drug_costs,
    test_costs = test_costs,
    visit_costs = visit_costs,
    follow_up_costs = follow_up_costs,
    post_progression_costs = post_progression_costs,
    end_life_costs = end_life_costs,
    total_costs_undiscounted = total_costs_undiscounted,
    total_costs_discounted = total_costs_discounted,
    costs_total = costs_total
  ))
}
