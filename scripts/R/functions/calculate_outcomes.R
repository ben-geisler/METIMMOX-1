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
  quarterly_cycles <- quarterly_cycles[quarterly_cycles <= n_cycles]
  
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
    if(biomarker == "tmb_braf") {
      test_costs[1] <- test_costs[1] + params$c_test_NGS
    } else if(biomarker == "crp") {
      test_costs[1] <- test_costs[1] + params$c_test_blood  # CRP is a blood test
    } else if(biomarker == "tlr") {
      test_costs[1] <- test_costs[1] + params$c_test_CT  # TLR requires CT measurement
    }
  }
  
  # Visit costs
  visit_costs <- l_visit * params$c_other_visit * p_pf
  visit_costs[1] <- visit_costs[1] + params$c_other_baseline
  
  # Follow-up costs
  follow_up_costs <- rep(0, n_cycles)
  if(length(quarterly_cycles) > 0) {
    follow_up_costs[quarterly_cycles] <- params$c_other_follow * p_p[quarterly_cycles]
  }
  
  # End-of-life costs
  death_transitions <- c(p_d[1], diff(p_d))
  end_life_costs <- death_transitions * params$c_other_last
  
  # Total costs
  total_costs_undiscounted <- drug_costs + test_costs + visit_costs + follow_up_costs + end_life_costs
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
    end_life_costs = end_life_costs,
    total_costs_undiscounted = total_costs_undiscounted,
    total_costs_discounted = total_costs_discounted,
    costs_total = costs_total
  ))
}