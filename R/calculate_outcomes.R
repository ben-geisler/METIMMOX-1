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

# Resource use from the trial team's answers of 15 October 2026 (issue #207,
# docs/open_questions_2026-10-15-answers.csv).
# MON3 ("FLOX 2 is two days"): a FLOX administration (bolus 5-FU/folinic acid
# on days 1 and 2) is two outpatient contacts; a nivolumab administration is
# one. l_visit carries the first contact; the second is charged here on the
# FLOX schedule of the regimen, so the 0/1 schedule contract of model_fun()
# is unchanged and dropping a FLOX administration drops both contacts.
FLOX_CONTACTS_PER_ADMINISTRATION <- 2L
# MON2 ("every 2.5 months after progression"): progressed patients have one CT
# scan and one follow-up visit every 2.5 months, 12 / 2.5 = 4.8 per
# patient-year of progressed occupancy (was one per quarter, 4 per year).
PROGRESSED_MONITORING_INTERVAL_MONTHS <- 2.5
PROGRESSED_MONITORING_PER_YEAR <- 12 / PROGRESSED_MONITORING_INTERVAL_MONTHS
# The post-progression treatment cost c_other_pp stays a rate per quarter.
QUARTERS_PER_YEAR <- 4

#' Expected resource use of the third and fourth treatment sequences (#208)
#'
#' Per patient assigned to an arm, whatever the health state: the share
#' starting sequence k (p_seq3_<arm>, and p_seq3_<arm> *
#' p_seq4_given_seq3_<arm> for the fourth) times one full sequence of the
#' arm's regimen (params$later_sequences$regimen: administration weeks from
#' the sequence start, the first sequence's pattern), placed at the observed
#' start weeks of the trial's starters (params$later_sequences$start_weeks,
#' weeks from inclusion, the OS/PFS clock). Each FLOX administration is
#' FLOX_CONTACTS_PER_ADMINISTRATION outpatient contacts and each nivolumab
#' administration one (issue #207); every administration week has one blood
#' test (MON1). Administrations beyond the grid are dropped, as in the
#' schedules. Without params$later_sequences (test fixtures) there are none.
#'
#' @param params Parameter list.
#' @param arm "control" (FLOX) or "experimental" (FLOX + nivolumab).
#' @param n_points Number of grid points.
#' @return List of per-time-point expected counts: FLOX, nivolumab, contacts,
#'   blood.
later_sequence_resource_use <- function(params, arm, n_points) {
  v_zero <- rep(0, n_points)
  out <- list(FLOX = v_zero, nivolumab = v_zero, contacts = v_zero, blood = v_zero)
  later <- params$later_sequences
  if (is.null(later)) return(out)
  v_names <- paste0(c("p_seq3_", "p_seq4_given_seq3_"), arm)
  v_values <- unlist(params[v_names])
  if (length(v_values) != 2L || any(!is.finite(v_values)) || any(v_values < 0 | v_values > 1))
    stop("Later-sequence shares ", paste(v_names, collapse = ", "),
         " must be finite values in [0, 1].")
  v_share <- c(`3` = v_values[[1]], `4` = v_values[[1]] * v_values[[2]])
  regimen <- later$regimen[[arm]]
  for (k in names(v_share)) {
    if (v_share[[k]] == 0) next
    df_start <- later$start_weeks[later$start_weeks$arm == arm &
                                    later$start_weeks$sequence == as.integer(k), , drop = FALSE]
    if (!nrow(df_start) || sum(df_start$n) <= 0)
      stop("A positive share starts sequence ", k, " in the ", arm,
           " arm but no start weeks are given.")
    v_weight <- v_share[[k]] * df_start$n / sum(df_start$n)
    for (i in seq_len(nrow(df_start))) for (drug in c("FLOX", "nivolumab")) {
      v_pos <- df_start$week[i] + regimen[[drug]] + 1L
      v_pos <- v_pos[v_pos >= 1L & v_pos <= n_points]
      out[[drug]][v_pos] <- out[[drug]][v_pos] + v_weight[i]
    }
  }
  out$contacts <- FLOX_CONTACTS_PER_ADMINISTRATION * out$FLOX + out$nivolumab
  out$blood <- out$FLOX + out$nivolumab
  out
}

#' Costs, QALYs, life-years and nivolumab administrations from state occupancy
#'
#' Outpatient contacts while progression-free are l_visit plus, for each FLOX
#' administration of the regimen, FLOX_CONTACTS_PER_ADMINISTRATION - 1 further
#' contact (issue #207). Progressed occupancy accrues the follow-up visit and
#' the CT scan PROGRESSED_MONITORING_PER_YEAR times per patient-year and the
#' post-progression treatment cost QUARTERS_PER_YEAR times. The third and
#' fourth treatment sequences (issue #208) are expected costs per patient
#' assigned to the regimen's arm, independent of the health state
#' (later_sequence_resource_use()); their nivolumab administrations are part
#' of the discounted administrations.
#'
#' Life-years (#197) are the trapezoidal area in years under the OS curve,
#' read as the alive occupancy `1 - p_d` (partitioned_survival_states() sets
#' `p_d = 1 - os`), with the interval weights of the QALYs and, when
#' discounted, the effect discount factors. They are not `p_pf + p_p`: with an
#' uncapped PFS that crosses OS, `p_pf + p_p = max(pfs, os)` would overstate
#' them. When PFS does not exceed OS the discounted life-years equal the
#' discounted QALYs at both utilities 1. The discounted nivolumab
#' administrations (#200) are sum(l_nivo * p_pf * v_dw_c) for the experimental
#' regimen and zero otherwise; they are the vector that multiplies the
#' nivolumab price in the drug cost, so they are exactly the slope of the cost
#' in that price. Second-line therapy (#209) is charged at incident
#' progression by second_line_event_costs() and returned as
#' `second_line_costs`.
#'
#' @param params Parameter list (utilities, unit costs, schedules).
#' @param p_pf,p_p,p_d Progression-free, progressed and dead occupancy.
#' @param treatment_type "experimental" (FLOX + nivolumab) or anything else
#'   (standard FLOX).
#' @param biomarker Biomarker whose diagnostic test is charged at time 0, or
#'   NULL for no test.
#' @param v_dw_c,v_dw_e Cost and effect discount factors per time point.
#' @param cl Cycle length in years.
#' @return List of per-time-point vectors (QALYs, costs, `ly_undiscounted`,
#'   `ly_discounted`, `nivo_administrations`, `visit_contacts`: scheduled
#'   outpatient contacts per point before occupancy weighting,
#'   `later_sequence_costs`) and the totals
#'   `qalys_total`,
#'   `costs_total`, `ly_disc_total`, `ly_undisc_total` and
#'   `nivo_administrations_total` (discounted).
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

  # Life-years (issue #197): the OS curve (alive occupancy 1 - p_d) integrated
  # with the QALY interval weights and, when discounted, the effect discount
  # factors. Not p_pf + p_p, which equals max(pfs, os) when an uncapped PFS
  # crosses OS. Without a crossing, the discounted life-years equal the
  # discounted QALYs at both utilities 1.
  v_ly_undiscounted <- (1 - p_d) * cl * v_occupancy_weights
  v_ly_discounted <- v_ly_undiscounted * v_dw_e

  # Calculate costs based on treatment type. The nivolumab administrations
  # (progression-free occupancy on the nivolumab schedule) are the single
  # definition of the price slope used by the value-of-information code (#200).
  if(treatment_type == "experimental") {
    # For biomarker positive: nivolumab + FLOX experimental regimen
    l_FLOX_regimen <- l_FLOX_exp
    v_nivo_administrations <- l_nivo * p_pf
    v_drug_costs <- v_nivo_administrations * params$c_drug_nivo +
      l_FLOX_exp * params$c_drug_FLOX * p_pf
  } else {
    # For control and biomarker negative: standard FLOX regimen
    l_FLOX_regimen <- l_FLOX_control
    v_nivo_administrations <- rep(0, n_cycles)
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
  
  # Third and fourth treatment sequences (issue #208): expected per patient
  # assigned to the arm, at the observed start weeks, whatever the state.
  later <- later_sequence_resource_use(
    params, if (treatment_type == "experimental") "experimental" else "control", n_cycles)
  v_later_sequence_costs <- later$FLOX * params$c_drug_FLOX +
    later$nivolumab * params$c_drug_nivo + later$contacts * params$c_other_visit +
    later$blood * params$c_test_blood
  v_nivo_administrations <- v_nivo_administrations + later$nivolumab

  # Visit costs: one contact per l_visit point plus the second day of every
  # FLOX administration of the regimen (issue #207).
  v_visit_contacts <- l_visit + (FLOX_CONTACTS_PER_ADMINISTRATION - 1L) * l_FLOX_regimen
  v_visit_costs <- v_visit_contacts * params$c_other_visit * p_pf
  v_visit_costs[1] <- v_visit_costs[1] + params$c_other_baseline

  # Ongoing progressed-state costs are rates over progressed occupancy,
  # integrated over intervals with discounting, just as for utilities.
  # Follow-up visit and CT every 2.5 months (issue #207; quarterly before).
  v_progressed_monitoring <- p_p * (PROGRESSED_MONITORING_PER_YEAR * cl) * v_occupancy_weights
  v_follow_up_costs <- params$c_other_follow * v_progressed_monitoring
  # Progressed patients are imaged with each follow-up visit (issues #164,
  # #207), charged at the CT unit price. Progression-free patients receive
  # their CT scans from the l_CT schedule.
  v_progressed_imaging_costs <- params$c_test_CT * v_progressed_monitoring
  v_progressed_quarters <- p_p * (QUARTERS_PER_YEAR * cl) * v_occupancy_weights

  # Flat post-progression treatment rate per quarter of progressed occupancy
  # (issue #154). Zero in the base case since issue #209, which costs
  # second-line therapy at incident progression instead (below); the former
  # EUR 5,000 per quarter is the Post_progression_cost structural scenario,
  # which replaces the pathway. The progressed state otherwise accrues the
  # follow-up visit and CT every 2.5 months (issue #207) and the end-of-life
  # cost.
  c_other_pp <- if (is.null(params$c_other_pp)) 0 else params$c_other_pp
  v_post_progression_costs <- c_other_pp * v_progressed_quarters

  # Second-line therapy (issue #209): the expected course cost per
  # progressing patient, charged once at incident progression (PFS decrement
  # times the progression share of PFS events) and discounted then.
  v_second_line_costs <- second_line_event_costs(p_pf, params)

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
    v_second_line_costs + v_end_life_costs + v_later_sequence_costs
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
    visit_contacts = v_visit_contacts,
    follow_up_costs = v_follow_up_costs,
    progressed_imaging_costs = v_progressed_imaging_costs,
    post_progression_costs = v_post_progression_costs,
    second_line_costs = v_second_line_costs,
    end_life_costs = v_end_life_costs,
    later_sequence_costs = v_later_sequence_costs,
    total_costs_undiscounted = v_total_costs_undiscounted,
    total_costs_discounted = v_total_costs_discounted,
    costs_total = costs_total,
    ly_undiscounted = v_ly_undiscounted,
    ly_discounted = v_ly_discounted,
    ly_disc_total = sum(v_ly_discounted),
    ly_undisc_total = sum(v_ly_undiscounted),
    nivo_administrations = v_nivo_administrations,
    nivo_administrations_total = sum(v_nivo_administrations * v_dw_c)
  ))
}

# ===============================================================================
# CELL-LEVEL OUTCOMES (issue #200, paper 3 step 1b)
# ===============================================================================
# The PSA predicts every patient of the prediction population under both arms
# and averages within the joint biomarker cells (R/prediction_population.R).
# Each cell x arm curve pair is evaluated here exactly as model_fun() evaluates
# a strategy subgroup in PSA mode: PFS capped at OS, partitioned_survival_states()
# (which applies pmax(os - pfs, 0)), then calculate_outcomes() without any
# diagnostic-test cost. Policies are combinations of these cell outcomes
# (R/voi_policies.R); because the clamps act on each cell's own curves, a cell
# combination equals a strategy's cost and QALYs exactly only when no cell
# curve crosses. Life-years depend on OS alone, which is linear in the cell
# curves and never clamped, so the cell combination equals the strategy
# life-years exactly in every draw (issue #197).

# Outcomes stored per cell and arm (the arms, CELL_ARMS, are defined in
# R/prediction_population.R).
CELL_OUTCOMES <- c("cost", "qalys", "ly_disc", "ly_undisc", "nivo_admin")

#' Outcomes of one cell under one arm
#'
#' @param params Parameter list (as for calculate_outcomes()).
#' @param os,pfs Cell-average OS and PFS curves on the model grid.
#' @param arm "control" (FLOX) or "experimental" (FLOX + nivolumab).
#' @param v_dw_c,v_dw_e Cost and effect discount factors.
#' @param cl Cycle length in years.
#' @return Named numeric vector over CELL_OUTCOMES with attribute
#'   `n_violations` (time points at which PFS exceeded OS before the cap).
cell_arm_outcomes <- function(params, os, pfs, arm, v_dw_c, v_dw_e, cl) {
  arm <- match.arg(arm, CELL_ARMS)
  n_violations <- sum(pfs > os)
  states <- partitioned_survival_states(os, pmin(pfs, os))
  outcomes <- calculate_outcomes(params, states$p_pf, states$p_p, states$p_d,
    treatment_type = if (arm == "experimental") "experimental" else "standard",
    biomarker = NULL, v_dw_c = v_dw_c, v_dw_e = v_dw_e, cl = cl)
  v_out <- c(cost = outcomes$costs_total, qalys = outcomes$qalys_total,
             ly_disc = outcomes$ly_disc_total, ly_undisc = outcomes$ly_undisc_total,
             nivo_admin = outcomes$nivo_administrations_total)
  attr(v_out, "n_violations") <- n_violations
  v_out
}

#' Outcomes of every joint biomarker cell under both arms
#'
#' @param params Parameter list (as for calculate_outcomes()).
#' @param cell_curves Cell curves from generate_psa_cell_predictions() or
#'   predict_cell_curves(): `os` and `pfs` matrices (time x cell_arm) and
#'   `counts` (patients per cell).
#' @param v_dw_c,v_dw_e Cost and effect discount factors.
#' @param cl Cycle length in years.
#' @return List with `values` (array cell x arm x outcome; NA for a cell
#'   without patients) and `violations` (matrix cell x arm of time points at
#'   which PFS exceeded OS before the cap). Stops if an observed cell curve
#'   breaks the survival-curve contract.
calculate_cell_outcomes <- function(params, cell_curves, v_dw_c, v_dw_e, cl) {
  v_cells <- names(cell_curves$counts)
  a_values <- array(NA_real_, dim = c(length(v_cells), length(CELL_ARMS), length(CELL_OUTCOMES)),
                    dimnames = list(cell = v_cells, arm = CELL_ARMS, outcome = CELL_OUTCOMES))
  m_violations <- matrix(NA_integer_, length(v_cells), length(CELL_ARMS),
                         dimnames = list(cell = v_cells, arm = CELL_ARMS))
  for (cell in v_cells[cell_curves$counts > 0]) {
    for (arm in CELL_ARMS) {
      column <- cell_curve_columns(cell, arm)
      validate_survival_curve(cell_curves$os[, column], paste0("cell ", column, " OS"))
      validate_survival_curve(cell_curves$pfs[, column], paste0("cell ", column, " PFS"))
      v_out <- cell_arm_outcomes(params, cell_curves$os[, column], cell_curves$pfs[, column],
                                 arm, v_dw_c, v_dw_e, cl)
      a_values[cell, arm, ] <- v_out[CELL_OUTCOMES]
      m_violations[cell, arm] <- attr(v_out, "n_violations")
    }
  }
  list(values = a_values, violations = m_violations)
}

# ===============================================================================
# SECOND-LINE THERAPY AT INCIDENT PROGRESSION (issue #209)
# ===============================================================================
# The costed Norwegian second-line pathway (data/external/second_line_pathway.csv,
# R/second_line_costs.R) is charged once per progression, discounted at the
# time of progression, rather than as a rate over progressed occupancy: its
# inputs are mean courses among treated patients (early stopping and death on
# treatment included), so the cost depends on how many patients progress and
# when, not on how long they survive afterwards. Incident progressions per
# interval are the PFS decrement times the share of PFS events that are
# progressions rather than deaths without progression (p_pfs_progression, a
# trial count). The cost is linear in the PFS curve, so it adds across cells
# and subgroups exactly like the other PFS-driven costs.

# Fields of the parameter list that define the second-line cost.
SECOND_LINE_MODEL_PARAMETERS <- c(
  "p_pfs_progression", "p_sl_treated", "p_sl_folfiri", "p_sl_ras_braf_wt",
  "p_sl_antiegfr_if_wt", "p_sl_bevacizumab_if_no_antiegfr",
  "n_sl_cycles_none", "n_sl_cycles_bevacizumab", "n_sl_cycles_panitumumab",
  "n_sl_contacts_FOLFIRI", "n_sl_contacts_FLIRI", "n_sl_blood_tests",
  "c_sl_drug_FOLFIRI", "c_sl_drug_FLIRI", "c_sl_drug_bevacizumab", "c_sl_drug_panitumumab")

#' Whether a parameter list costs second-line therapy
#'
#' A list with none of SECOND_LINE_MODEL_PARAMETERS (test fixtures, lists
#' built before issue #209) costs no second-line therapy; a list with only
#' some of them is an error, so a partly configured pathway cannot pass as
#' zero cost. The base case (script 05) sets all of them.
#'
#' @param params Parameter list.
#' @return TRUE when every field is present, FALSE when none is.
second_line_costed <- function(params) {
  v_present <- SECOND_LINE_MODEL_PARAMETERS %in% names(params)
  if (all(v_present)) return(TRUE)
  if (!any(v_present)) return(FALSE)
  stop("Second-line parameters incomplete; missing: ",
       paste(SECOND_LINE_MODEL_PARAMETERS[!v_present], collapse = ", "))
}

#' Expected second-line cost per progressing patient (EUR)
#'
#' Closed form of second_line_cost_at() (R/second_line_costs.R) on the scalar
#' model inputs, so the PSA and the DSA vary it through the parameter list:
#' share treated x sum over backbone b (FOLFIRI with share p_sl_folfiri, FLIRI
#' otherwise) and biologic k (panitumumab with share p_sl_ras_braf_wt x
#' p_sl_antiegfr_if_wt, bevacizumab with share (1 - that) x
#' p_sl_bevacizumab_if_no_antiegfr, none otherwise) of share_b x share_k x
#' cycles_k x (drug_b + drug_k + contacts_b x c_other_visit + blood tests x
#' c_test_blood).
#'
#' @param params Parameter list with SECOND_LINE_MODEL_PARAMETERS,
#'   `c_other_visit` and `c_test_blood`.
#' @return Numeric scalar.
second_line_cost_per_progression <- function(params) {
  v_fields <- c(setdiff(SECOND_LINE_MODEL_PARAMETERS, "p_pfs_progression"), "c_other_visit", "c_test_blood")
  for (field in v_fields) {
    x <- params[[field]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop("Second-line parameter '", field, "' must be a finite nonnegative scalar.")
    }
    if (startsWith(field, "p_sl_") && x > 1) stop("Second-line share '", field, "' must be in [0, 1].")
  }
  p_egfr <- params$p_sl_ras_braf_wt * params$p_sl_antiegfr_if_wt
  p_bev <- (1 - p_egfr) * params$p_sl_bevacizumab_if_no_antiegfr
  v_biologic_share <- c(none = 1 - p_egfr - p_bev, bevacizumab = p_bev, panitumumab = p_egfr)
  v_biologic_drug <- c(none = 0, bevacizumab = params$c_sl_drug_bevacizumab,
                       panitumumab = params$c_sl_drug_panitumumab)
  v_cycles <- c(none = params$n_sl_cycles_none, bevacizumab = params$n_sl_cycles_bevacizumab,
                panitumumab = params$n_sl_cycles_panitumumab)
  v_backbone_share <- c(FOLFIRI = params$p_sl_folfiri, FLIRI = 1 - params$p_sl_folfiri)
  v_backbone_cycle <- c(
    FOLFIRI = params$c_sl_drug_FOLFIRI + params$n_sl_contacts_FOLFIRI * params$c_other_visit,
    FLIRI = params$c_sl_drug_FLIRI + params$n_sl_contacts_FLIRI * params$c_other_visit) +
    params$n_sl_blood_tests * params$c_test_blood
  per_treated <- sum(outer(v_backbone_share, v_biologic_share) *
                       outer(v_backbone_cycle, v_biologic_drug, "+") *
                       matrix(v_cycles, 2, 3, byrow = TRUE))
  params$p_sl_treated * per_treated
}

#' Second-line costs at incident progression, per time point (undiscounted)
#'
#' The PFS decrement of the interval ending at grid point i (floored at 0)
#' times p_pfs_progression is the share of the cohort progressing in that
#' interval; it is charged second_line_cost_per_progression() at grid point
#' i, so calculate_outcomes() discounts it at the end of the interval of
#' progression. Zero at grid point 1 (no progression at time 0) and for a
#' list without second-line fields (second_line_costed()).
#'
#' @param p_pf Progression-free occupancy (the PFS curve, capped at OS where
#'   the caller caps it).
#' @param params Parameter list.
#' @return Numeric vector of the length of p_pf.
second_line_event_costs <- function(p_pf, params) {
  if (!second_line_costed(params)) return(rep(0, length(p_pf)))
  share <- params$p_pfs_progression
  if (!is.numeric(share) || length(share) != 1L || !is.finite(share) || share < 0 || share > 1) {
    stop("p_pfs_progression must be a finite scalar in [0, 1].")
  }
  v_pfs_events <- c(0, pmax(-diff(p_pf), 0))
  share * v_pfs_events * second_line_cost_per_progression(params)
}
