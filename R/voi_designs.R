# ===============================================================================
# CANDIDATE STUDY DESIGNS FOR THE EVSI OF PAPER 3 (issue #203, plan step 5)
# ===============================================================================
# Five candidate studies and their power-based versions, each turned into the
# input of the EVSI engine (evsi_design() in R/evsi_functions.R):
#
#   S1  All-comers RCT, FLOX +/- nivolumab, stratified by week-4 CRP and
#       TMB/BRAF, randomised at week 4 after CRP is known: only patients alive
#       and progression-free at week 4 are randomised (trial_design() with
#       randomisation_week = 4; the cell mix is reweighted by P(PFS > 4)).
#       Identifies Rx, Rx:crp, Rx:tmb_braf on both endpoints.
#   S2  CRP-enriched RCT: CRP-normal at week 4 (cells "10", "11"), also
#       randomised at week 4. Identifies Rx:tmb_braf and cell_effect_10.
#   S3  METIMMOX-3 as in protocol v1.1 of 1 October 2026 (EU CT
#       2026-526317-29-00; #213): non-randomised, single-arm enrichment design
#       stratified by TMB/BRAF. Every patient is enrolled at baseline and
#       starts two cycles of FLOX/FOLFOX with bevacizumab
#       (METIMMOX3_STARTUP_WEEKS) while the tumour is sequenced (OCAP NGS);
#       TMB >= 9 mut/Mb or BRAF V600E (cells "01", "11") continue with
#       FLOX/FOLFOX and bevacizumab alternating with tislelizumab (12-cycle
#       sequences), the others (cells "00", "10") continue standard-of-care
#       chemotherapy at the oncologist's discretion. Type "stratified"
#       (stratified_design_information()): experimental-arm data in the
#       positive cells from week 0, whose treatment-effect directions may
#       differ from METIMMOX's regimen (tislelizumab for nivolumab,
#       bevacizumab, optional FOLFOX) by b ~ N(0, tau_m3^2 I); control-arm data
#       in the negative cells, whose baseline log rate may differ from
#       METIMMOX's FLOX by b ~ N(0, tau_soc^2) over the whole follow-up (the
#       start-up already adds bevacizumab; afterwards FLOX or FOLFOX +/-
#       bevacizumab). Single centre, n = 77 at the planned accrual rate,
#       follow-up to the protocol's data cut-off, and on the n grid at the
#       same rate. The size is descriptive (protocol section 9.5: 31
#       biomarker-positive patients for a precise median PFS), not a power
#       calculation. Within the model's structure (no CRP x TMB/BRAF term) it
#       identifies the control log rate, the CRP main effect, the
#       experimental log rate of the TMB/BRAF-positive cells and Rx:crp
#       (stratified_targets()), but not the treatment effect in either
#       TMB/BRAF stratum, which it informs only through the prior correlation.
#   S3R Randomised counterpart of S3 (the design assumed until 12 October
#       2026), kept as a labelled comparison: FLOX +/- nivolumab in
#       TMB/BRAF-positive patients (cells "01", "11"), randomised at baseline,
#       n = 77 and the n grid, METIMMOX's multicentre accrual. Identifies
#       Rx:crp and cell_effect_01. It keeps METIMMOX's regimen (no tau_m3 or
#       tau_soc), so S3 against S3R also carries the regimen heterogeneity;
#       S3 with tau_m3 = tau_soc = 0 isolates what randomisation adds.
#   S4  METIMMOX + METIMMOX-2 pooled. METIMMOX is already in the prior, so the
#       new information is METIMMOX-2's single experimental arm (80 patients,
#       every cell). Approximation (single_arm_information()): the first
#       METIMMOX2_FULL_WEIGHT_WEEKS weeks are experimental-arm information at
#       full weight (follow-up censored at that week); later follow-up, which
#       describes METIMMOX-2's response-adapted regimen (only target lesion
#       responders continue), is left-truncated at that week and enters with a
#       trial-specific heterogeneity term: the treatment-effect coefficients
#       (Rx and both Rx interactions, per endpoint) of the later follow-up
#       differ from METIMMOX's by b ~ N(0, tau_m2^2 I), which removes
#       information in those directions by the Woodbury identity
#       (heterogeneity_adjusted_information()). Single-arm data identify the
#       experimental-arm linear predictors only; the summary statistic is the
#       experimental-arm log rate at the reference covariate mean and its two
#       biomarker shifts per endpoint (s4_targets()), a subset of what the data
#       identify (shape, Age and sex are left out), so the EVSI is a lower bound.
#   S6  Power-based versions of S1, S2 and S3R (ruling (4) of 11 October
#       2026; S3 has no randomised contrast to power): the
#       smallest n at which the engine's standard error of the design's
#       target (design_information(); power_based_n_information()) gives 80%
#       power, two-sided 5%, primary; Schoenfeld's variance criterion with the
#       design's expected events (power_based_n()) as the cross-check, both
#       reported (design_power_sizes(), column `method`). Two hypotheses
#       (column `hypothesis`): the conventional rate ratio 0.5 and the
#       decision-relevant effect, the target value at which the guided policy
#       breaks even against standard of care at the price in question
#       (decision_relevant_effect(): the interaction coefficient shifted from
#       the mean coefficients, deterministic cell-path evaluation, uniroot).
#
# Sequential value: the information of S4 plus that of Sx is summed per
# endpoint and the EVSI of the combined design is computed on the union of
# both target sets (combined_design_information()); EVSI(S4 + Sx) - EVSI(S4)
# is then the value of Sx after S4 (never a nested preposterior analysis).
#
# Accrual. METIMMOX randomised METIMMOX_ACCRUAL_PATIENTS patients at
# METIMMOX_ACCRUAL_SITES Norwegian sites between METIMMOX_ACCRUAL_START and
# METIMMOX_ACCRUAL_END; that rate of registered patients per year is the
# Norwegian reference, scaled by the eligible population of each scenario
# (annual_incidence_scenarios / annual_incidence_norway) for the Nordic and
# European rates (accrual_scenarios()). A design's randomised rate is the
# registered rate times its enrolment fraction (share of registered patients
# in the enrolled cells, times P(PFS > L) for randomisation at week L). S3 is
# a single-centre trial (protocol v1.1 lists Akershus University Hospital as
# the only site) with its own planned rate (METIMMOX3_ACCRUAL_START to
# METIMMOX3_ACCRUAL_END) in every accrual scenario; the scenarios change only
# the population that benefits (accrual_independent_designs()).
#
# ASSUMPTIONS, each a named parameter below (inputs of plan.md section 4; the
# trial team's answers of 15 October 2026 and the protocols, #213): the
# METIMMOX-2 report lag, TAU_M2_DEFAULT, TAU_SOC_DEFAULT and TAU_M3_DEFAULT
# with their ranges (elicitation pending), the two-year minimum follow-up of
# S1, S2 and S3R, METIMMOX-3's first inclusion (protocol v1.1 gives no
# recruitment dates; the 2024 synopsis's 24-month recruitment is kept and
# shifted to start after the protocol), its analysis at the protocol's data
# cut-off, the model's TMB/BRAF definition (TMB >= 9 mut/Mb on TSO500 or any
# BRAF mutation) standing in for METIMMOX-3's TMB >= 9 mut/Mb on OCAP or BRAF
# V600E, METIMMOX-3's regimen differing from METIMMOX's only through tau_m3
# (experimental cells) and tau_soc (control cells), and Schoenfeld's formula
# applied to gamma rate ratios as if they were hazard ratios. Confirmed by
# the protocol: a CRP at every visit, including cycle 3 day 1 (week 4), and
# no CRP continuation rule (continuation depends on TMB/BRAF only).
#
# Requires R/design_information.R (trial_design(), design_information(),
# design_patient_mix(), expected_gamma_information(), information_pseudo_inverse(),
# joint_target_matrix(), design_treatment_targets(), cox_interaction_heuristic(),
# gamma_information_matrix(), gamma_linear_predictor(), model_rx_levels()),
# R/evsi_functions.R (evsi_design()), R/voi_survival_parameters.R
# (survival_coefficient_map()), R/joint_survival_sampling.R
# (covariance_root()), R/pfs_os_violation_diagnostics.R
# (flexsurv_design_matrix()) and R/model_configs.R (biomarker_cost_key());
# the decision-relevant effect also needs R/prediction_functions.R
# (build_sampled_flexsurv_model(), predict_cell_curves()),
# R/prediction_population.R, R/calculate_outcomes.R (discount_weights(),
# calculate_cell_outcomes()) and R/voi_policies.R (policy_conditions_on()).
# ===============================================================================

#' Date of the research decision ("now"), the origin of every reporting year
VOI_DECISION_DATE <- as.Date("2026-10-10")

#' METIMMOX accrual: randomised patients (Ree et al. 2024, Br J Cancer)
METIMMOX_ACCRUAL_PATIENTS <- 80
#' METIMMOX accrual: Norwegian sites (NCT03388190)
METIMMOX_ACCRUAL_SITES <- 5
#' METIMMOX accrual start: registry study start (NCT03388190). ASSUMPTION
#' until I5 confirms the first inclusion date.
METIMMOX_ACCRUAL_START <- as.Date("2018-05-29")
#' METIMMOX accrual end: last enrolment 22 October 2021 (Ree et al. 2024, Br
#' J Cancer 130:1921, doi:10.1038/s41416-024-02696-6: enrolment 29 May 2018 to
#' 22 October 2021; #213, replaces the December 2021 placeholder). The
#' patients per site (M1-2) are still open.
METIMMOX_ACCRUAL_END <- as.Date("2021-10-22")

#' METIMMOX-2 (NCT05504252): patients, single experimental arm (registry
#' enrolment). Protocol v3.0 gives 110 in its synopsis and section 4.2 but 57
#' in sections 10.1 and 10.3 (v1.2 had 110 throughout); which is the actual
#' enrolment is open (M2-2).
METIMMOX2_PATIENTS <- 80
#' METIMMOX-2 accrual start (NCT05504252)
METIMMOX2_ACCRUAL_START <- as.Date("2022-10-01")
#' METIMMOX-2 accrual end: last inclusion September 2024 (trial team, 15
#' October 2026, M2-1; the day is not given, the month end is used).
METIMMOX2_ACCRUAL_END <- as.Date("2024-09-30")
#' METIMMOX-2 primary completion (NCT05504252), the analysis date
METIMMOX2_ANALYSIS_DATE <- as.Date("2028-01-31")
#' Years from METIMMOX-2's primary completion to the pooled report. ASSUMPTION.
METIMMOX2_REPORT_LAG_YEARS <- 0.5
#' Weeks of METIMMOX-2 follow-up at full weight (two FLOX and two nivolumab
#' cycles for everyone before the response-adapted continuation)
METIMMOX2_FULL_WEIGHT_WEEKS <- 8
#' SD of the difference in the treatment-effect coefficients (log rate scale)
#' between METIMMOX's regimen and METIMMOX-2's response-adapted continuation.
#' ASSUMPTION, elicitation pending (I13): 0.35 puts the rate ratio between the
#' regimens within a factor of about 2 with 95% probability.
TAU_M2_DEFAULT <- 0.35
#' Scenario range of tau_m2 (elicitation pending)
TAU_M2_RANGE <- c(low = 0.1, high = 1)

#' Minimum follow-up after the last randomisation (years). ASSUMPTION.
DESIGN_MIN_FOLLOW_UP_YEARS <- 2
#' Planned size of METIMMOX-3: at least 77 enrolled, about 40%
#' biomarker-positive (protocol v1.1, section 9.5). Also the proposal size of
#' S3R.
METIMMOX3_N <- 77
#' Biomarker-positive METIMMOX-3 patients the size is built on (protocol v1.1,
#' section 9.5): descriptive, about 28 PFS events and a 95% CI of about 14 to
#' 29 months around an assumed median PFS of 19.8 months; no power calculation
#' (trial team, 15 October 2026, M3-3)
METIMMOX3_N_POSITIVE <- 31
#' METIMMOX-3 first inclusion. ASSUMPTION: protocol v1.1 (1 October 2026)
#' gives no recruitment dates; first inclusion is taken as the first quarter
#' after the protocol, 1 January 2027 (the 2024 synopsis planned 2025 Q3).
METIMMOX3_ACCRUAL_START <- as.Date("2027-01-01")
#' METIMMOX-3 last inclusion. ASSUMPTION: the 24 months of recruitment of the
#' 2024 synopsis (one site, about 38.5 patients a year) after
#' METIMMOX3_ACCRUAL_START.
METIMMOX3_ACCRUAL_END <- as.Date("2028-12-31")
#' METIMMOX-3 data cut-off and end of study (protocol v1.1, section 4.4:
#' participation ends at the latest on 31 December 2032; section 9: analysis
#' and reporting on all data at the end of the study)
METIMMOX3_DATA_CUTOFF <- as.Date("2032-12-31")
#' METIMMOX-3 minimum follow-up after the last inclusion (years): from
#' METIMMOX3_ACCRUAL_END to the data cut-off (4 years; the 2024 synopsis had
#' 1.5). Kept fixed on the n grid.
METIMMOX3_MIN_FOLLOW_UP_YEARS <- as.numeric(difftime(METIMMOX3_DATA_CUTOFF, METIMMOX3_ACCRUAL_END,
                                                     units = "days")) / 365.25
#' Weeks of FLOX/FOLFOX with bevacizumab for every METIMMOX-3 patient while
#' TMB/BRAF is sequenced (two cycles; protocol v1.1, section 4.1 and Figure 3)
METIMMOX3_STARTUP_WEEKS <- 4
#' SD of the difference in the baseline log rate between the control regimen
#' of METIMMOX-3's TMB/BRAF-negative patients (FLOX/FOLFOX with bevacizumab in
#' the start-up, then FLOX or FOLFOX +/- bevacizumab at the oncologist's
#' discretion; protocol v1.1 and M3-5) and METIMMOX's FLOX, over the whole
#' control follow-up. ASSUMPTION, elicitation pending (I13): 0.2 puts the
#' rate ratio between the two within about 0.67 to 1.48 with 95% probability.
TAU_SOC_DEFAULT <- 0.2
#' Scenario range of tau_soc (elicitation pending)
TAU_SOC_RANGE <- c(low = 0, high = 0.5)
#' SD of the difference in the treatment-effect coefficients (Rx and both Rx
#' interactions, log rate scale, per endpoint) between METIMMOX-3's regimen
#' (FLOX/FOLFOX with bevacizumab alternating with tislelizumab, 12-cycle
#' sequences) and METIMMOX's (FLOX alternating with nivolumab). ASSUMPTION,
#' elicitation pending (I13; #213, proposal to confirm): 0.2 like tau_soc,
#' below METIMMOX-2's 0.35 because the regimen is the same alternating schedule
#' with a PD-1 antibody of the same class, not a response-adapted strategy.
TAU_M3_DEFAULT <- 0.2
#' Scenario range of tau_m3 (elicitation pending)
TAU_M3_RANGE <- c(low = 0, high = 0.5)
#' Rate ratio for the power-based sizes: the interaction for S1, the
#' conditional treatment effect for S2 and S3
POWER_RATE_RATIO <- 0.5
#' Power and two-sided significance level of the power-based sizes
POWER_LEVEL <- 0.8
POWER_ALPHA <- 0.05
#' Endpoint of the power-based sizes (METIMMOX's primary endpoint)
POWER_ENDPOINT <- "pfs"
#' Default sample-size grid of the EVSI curves
VOI_N_GRID <- c(40, 80, 120, 160, 240, 320, 480, 640, 960, 1280, 1600, 2000)
#' Maximum quadrature nodes per design_information() call (the node spacing
#' grows with the trial duration)
DESIGN_MAX_NODES <- 2000

#' Candidate studies S1, S2, S3, S3R and S4 of paper 3
#'
#' @return Named list of design specifications, each with `id`, `label`,
#'   `type` ("rct", "stratified" or "single_arm"), `enrolled_cells` (cell
#'   codes, CRP digit first), `experimental_cells` (cells on the experimental
#'   arm: the enrolled cells of a randomised or single-arm design, the
#'   positive cells of the stratifying biomarker of a stratified one),
#'   `randomisation_week`, `p_experimental` (NA for a stratified design, where
#'   it follows from the cell mix; design_allocation()),
#'   `screening_biomarkers` (tests paid for every registered patient),
#'   `n_fixed` (fixed size or NA), `n_proposal` (proposed size or NA),
#'   `power_target` ("interaction", "conditional" or NA), `power_biomarkers`,
#'   `notes` and `sample_size_basis`; a stratified design also has `accrual_per_year` (enrolled
#'   patients per year), `accrual_start`, `min_follow_up_years` and
#'   `startup_weeks`.
candidate_designs <- function() {
  spec <- function(id, label, type, cells, week, screening, n_fixed = NA_real_, n_proposal = NA_real_,
                   power_target = NA_character_, power_biomarkers = character(0), notes = "",
                   experimental_cells = cells, sample_size_basis = "", ...) {
    c(list(id = id, label = label, type = type, enrolled_cells = cells, experimental_cells = experimental_cells,
           randomisation_week = week,
           p_experimental = switch(type, single_arm = 1, rct = 0.5, stratified = NA_real_),
           screening_biomarkers = screening, n_fixed = n_fixed, n_proposal = n_proposal,
           power_target = power_target, power_biomarkers = power_biomarkers, notes = notes,
           sample_size_basis = sample_size_basis), list(...))
  }
  basis_power <- "n grid; power-based n (S6) from the engine, Schoenfeld as cross-check"
  list(
    S1 = spec("S1", "All-comers RCT stratified by week-4 CRP and TMB/BRAF, randomised at week 4", "rct",
              DESIGN_ALL_CELLS, 4, c("crp", "tmb_braf"), power_target = "interaction",
              power_biomarkers = c("crp", "tmb_braf"), sample_size_basis = basis_power,
              notes = "Only patients progression-free at week 4 are randomised; both tests for every registered patient"),
    S2 = spec("S2", "CRP-enriched RCT (CRP-normal at week 4), randomised at week 4", "rct",
              c("10", "11"), 4, "crp", power_target = "conditional", sample_size_basis = basis_power,
              notes = "CRP test for every screened patient; nothing on the CRP-raised cells"),
    S3 = spec("S3", "METIMMOX-3 (protocol v1.1): non-randomised, stratified by TMB/BRAF", "stratified",
              DESIGN_ALL_CELLS, 0, "tmb_braf", n_proposal = METIMMOX3_N, experimental_cells = c("01", "11"),
              sample_size_basis = paste0("descriptive (protocol v1.1, section 9.5): ", METIMMOX3_N_POSITIVE,
                                         " biomarker-positive patients (about 40% of ", METIMMOX3_N, " enrolled) ",
                                         "for a precise median PFS; no power calculation; also on the n grid"),
              notes = paste0("TMB >= 9 or BRAF V600E on FLOX/FOLFOX + bevacizumab alternating with tislelizumab ",
                             "(donated; heterogeneity tau_m3), negative on standard-of-care chemotherapy after ",
                             METIMMOX3_STARTUP_WEEKS, " weeks of FLOX/FOLFOX + bevacizumab (heterogeneity tau_soc); ",
                             "single centre; NGS for every patient; follow-up to the data cut-off ",
                             format(METIMMOX3_DATA_CUTOFF), "; no randomised contrast (treatment effects only ",
                             "through the prior)"),
              accrual_per_year = accrual_per_year(METIMMOX3_N, METIMMOX3_ACCRUAL_START, METIMMOX3_ACCRUAL_END),
              accrual_start = METIMMOX3_ACCRUAL_START, min_follow_up_years = METIMMOX3_MIN_FOLLOW_UP_YEARS,
              startup_weeks = METIMMOX3_STARTUP_WEEKS),
    S3R = spec("S3R", "Randomised TMB/BRAF-enriched phase 2 (hypothetical counterpart of METIMMOX-3)", "rct",
               c("01", "11"), 0, "tmb_braf", n_proposal = METIMMOX3_N, power_target = "conditional",
               sample_size_basis = basis_power,
               notes = paste0("Comparison for S3 at the same n: FLOX +/- nivolumab (METIMMOX's regimen, no ",
                              "regimen heterogeneity) in TMB/BRAF-positive patients; NGS for every screened patient")),
    S4 = spec("S4", "METIMMOX + METIMMOX-2 pooled (single experimental arm)", "single_arm", DESIGN_ALL_CELLS, 0,
              character(0), n_fixed = METIMMOX2_PATIENTS,
              sample_size_basis = paste0("fixed: registry enrolment (protocol v3.0 gives 110 in section 4.2 and 57 ",
                                         "in section 10.1: HR 0.66, one-sided 5%, about 60% power)"),
              notes = paste0("First ", METIMMOX2_FULL_WEIGHT_WEEKS, " weeks at full weight; later follow-up ",
                             "under the response-adapted regimen with heterogeneity tau_m2"))
  )
}

#' Registered patients per year of a past accrual
#'
#' @param patients Patients registered.
#' @param start,end Accrual start and end dates.
#' @return Patients per year (365.25 days).
accrual_per_year <- function(patients, start, end) {
  years <- as.numeric(difftime(as.Date(end), as.Date(start), units = "days")) / 365.25
  if (!is.finite(years) || years <= 0) stop("The accrual end must follow its start")
  patients / years
}

#' Accrual scenarios: Norway (METIMMOX reference) and the scaled regions
#'
#' @param eligible Named eligible patients per year: `norway` and the
#'   regions to scale to (e.g. c(norway = annual_incidence_norway,
#'   annual_incidence_scenarios[c("nordic", "europe")])).
#' @param reference_rate Registered patients per year in Norway (default the
#'   METIMMOX accrual).
#' @return Data frame with `accrual_scenario`, `eligible_per_year`, `scale`
#'   (eligible relative to Norway) and `registered_per_year`.
accrual_scenarios <- function(eligible, reference_rate = accrual_per_year(
                                METIMMOX_ACCRUAL_PATIENTS, METIMMOX_ACCRUAL_START, METIMMOX_ACCRUAL_END)) {
  if (is.null(names(eligible)) || !"norway" %in% names(eligible)) stop("eligible must name 'norway'")
  if (any(!is.finite(eligible)) || any(eligible <= 0)) stop("Eligible populations must be positive")
  v_scale <- eligible / eligible[["norway"]]
  data.frame(accrual_scenario = names(eligible), eligible_per_year = unname(eligible),
             scale = unname(v_scale), registered_per_year = unname(reference_rate * v_scale),
             stringsAsFactors = FALSE)
}

#' Diagnostic-test cost per registered patient of a design
#'
#' @param spec Design specification (candidate_designs()).
#' @param params Parameter list with the test prices (c_test_CRP, c_test_NGS).
#' @return EUR per registered patient.
design_screening_cost <- function(spec, params) {
  if (!length(spec$screening_biomarkers)) return(0)
  sum(vapply(spec$screening_biomarkers, function(b) {
    key <- biomarker_cost_key(b)[[1]]
    price <- params[[key]]
    if (!is.numeric(price) || length(price) != 1L || !is.finite(price) || price < 0) {
      stop("No finite diagnostic-test price ", key, " for ", b)
    }
    price
  }, numeric(1)))
}

#' Reference cell mix and the share of registered patients a design randomises
#'
#' @param spec Design specification.
#' @param reference_population Reference patients (biomarker columns).
#' @param weights Reference weights (default equal: the observed joint mix).
#' @param component Joint sampling component.
#' @param draw list(os, pfs) of the draw defining P(PFS > L).
#' @return List with `fraction` (randomised per registered patient),
#'   `cell_mix_registration` (enrolled-cell shares among enrolled patients at
#'   registration), `cell_mix_randomisation` (after the landmark selection) and
#'   `enrolled_mass` (share of all patients in the enrolled cells).
design_enrolment_fraction <- function(spec, reference_population, weights = NULL, component, draw) {
  if (is.null(weights)) weights <- rep(1, nrow(reference_population))
  v_cell <- design_cell_codes(reference_population)
  v_w <- weights / sum(weights)
  v_enrolled <- v_cell %in% spec$enrolled_cells
  v_selection <- rep(1, length(v_w))
  if (spec$randomisation_week > 0) {
    fit_pfs <- component$original_pfs
    v_eta <- gamma_linear_predictor(fit_pfs, draw$pfs, reference_population, model_rx_levels(fit_pfs)[1])
    v_selection <- stats::pgamma(spec$randomisation_week, shape = exp(draw$pfs[["shape"]]),
                                 rate = exp(v_eta), lower.tail = FALSE)
  }
  v_mix_reg <- tapply(v_w[v_enrolled], factor(v_cell[v_enrolled], spec$enrolled_cells), sum)
  v_mix_rand <- tapply((v_w * v_selection)[v_enrolled], factor(v_cell[v_enrolled], spec$enrolled_cells), sum)
  v_mix_reg[is.na(v_mix_reg)] <- 0
  v_mix_rand[is.na(v_mix_rand)] <- 0
  list(fraction = sum(v_w[v_enrolled] * v_selection[v_enrolled]),
       cell_mix_registration = v_mix_reg / sum(v_mix_reg),
       cell_mix_randomisation = v_mix_rand / sum(v_mix_rand),
       enrolled_mass = sum(v_w[v_enrolled]))
}

#' Quadrature node spacing that keeps a design within DESIGN_MAX_NODES
#'
#' @param total_weeks Follow-up plus accrual (weeks).
#' @param base Smallest spacing (weeks).
#' @return Node spacing (weeks).
design_node_step <- function(total_weeks, base = 0.5) max(base, total_weeks / DESIGN_MAX_NODES)

#' trial_design() of a randomised candidate at an accrual rate
#'
#' @param spec Design specification of type "rct".
#' @param randomised_per_year Randomised patients per year.
#' @param min_follow_up_years Minimum follow-up after the last randomisation.
#' @return trial_design() object.
as_trial_design <- function(spec, randomised_per_year, min_follow_up_years = DESIGN_MIN_FOLLOW_UP_YEARS) {
  if (spec$type != "rct") stop("as_trial_design() is for randomised designs; use single_arm_information()")
  trial_design(spec$id, enrolled_cells = spec$enrolled_cells, p_experimental = spec$p_experimental,
               min_follow_up_weeks = 52 * min_follow_up_years, accrual_rate = randomised_per_year / 52,
               randomisation_week = spec$randomisation_week)
}

#' Reporting year of a randomised design
#'
#' Continuous time from the research decision: the randomisation landmark,
#' the accrual of n randomised patients and the minimum follow-up.
#'
#' @param n Randomised patients.
#' @param randomised_per_year Randomised patients per year.
#' @param min_follow_up_years Minimum follow-up (years).
#' @param randomisation_week Randomisation landmark (weeks).
#' @return Years (vectorised in n).
design_reporting_year <- function(n, randomised_per_year, min_follow_up_years = DESIGN_MIN_FOLLOW_UP_YEARS,
                                  randomisation_week = 0) {
  if (any(!is.finite(randomised_per_year)) || any(randomised_per_year <= 0)) {
    stop("randomised_per_year must be positive")
  }
  randomisation_week / 52 + n / randomised_per_year + min_follow_up_years
}

#' Reporting year of the METIMMOX-2 pooling (S4)
#'
#' @param decision_date,analysis_date,lag_years Research decision, METIMMOX-2
#'   analysis date and reporting lag.
#' @return Years from the decision to the pooled report.
s4_reporting_year <- function(decision_date = VOI_DECISION_DATE, analysis_date = METIMMOX2_ANALYSIS_DATE,
                              lag_years = METIMMOX2_REPORT_LAG_YEARS) {
  max(0, as.numeric(difftime(analysis_date, decision_date, units = "days")) / 365.25 + lag_years)
}

#' Reporting year of a stratified design (S3)
#'
#' Continuous time from the research decision: the planned first inclusion
#' (before the decision for a trial that is already recruiting), the accrual
#' of n patients at the design's own rate and the minimum follow-up. At the
#' planned n this is the planned last inclusion plus the minimum follow-up.
#'
#' @param n Enrolled patients (vector).
#' @param spec Design specification of type "stratified".
#' @param decision_date Research decision date.
#' @return Years (vectorised in n; at least 0).
stratified_reporting_year <- function(n, spec, decision_date = VOI_DECISION_DATE) {
  if (!identical(spec$type, "stratified")) stop("stratified_reporting_year() is for stratified designs")
  pmax(0, as.numeric(difftime(spec$accrual_start, decision_date, units = "days")) / 365.25 +
         n / spec$accrual_per_year + spec$min_follow_up_years)
}

#' Designs whose information does not depend on the accrual scenario
#'
#' The single-arm pooling (fixed data) and the stratified single-centre
#' design (its own accrual rate), and their combinations ("S4+S3").
#'
#' @param specs candidate_designs().
#' @return Design ids.
accrual_independent_designs <- function(specs) {
  v_type <- vapply(specs, function(s) s$type, character(1))
  v_single <- names(specs)[v_type == "single_arm"]
  v_strat <- names(specs)[v_type == "stratified"]
  c(v_single, v_strat, as.vector(outer(v_single, v_strat, paste, sep = "+")))
}

#' Allocation of a design: experimental share and the experimental cell mix
#'
#' @param spec Design specification.
#' @param reference_population,weights Reference patients and weights.
#' @param component Joint sampling component.
#' @param draw list(os, pfs) of the draw defining P(PFS > L).
#' @return List with `p_experimental` (the allocation of a randomised design,
#'   1 for a single arm, the mass of the experimental cells among the enrolled
#'   cells of a stratified design), `experimental_mix` (enrolled-cell shares of
#'   the experimental patients) and `enrolment` (design_enrolment_fraction()).
design_allocation <- function(spec, reference_population, weights = NULL, component, draw) {
  enrol <- design_enrolment_fraction(spec, reference_population, weights, component, draw)
  if (identical(spec$type, "stratified")) {
    v_mix <- enrol$cell_mix_registration
    v_exp <- v_mix[spec$experimental_cells]
    p <- sum(v_exp)
    if (!isTRUE(p > 0 && p < 1)) {
      stop("Design ", spec$id, " needs reference patients in its experimental and its control cells")
    }
    return(list(p_experimental = p, experimental_mix = v_exp / p, enrolment = enrol))
  }
  list(p_experimental = spec$p_experimental, experimental_mix = enrol$cell_mix_randomisation, enrolment = enrol)
}

# -------------------------------------------------------------------------------
# Information of the single-arm pooling (S4) and of combined designs
# -------------------------------------------------------------------------------

#' Information after a random bias on some coefficient directions (Woodbury)
#'
#' Data with information P about theta that estimate theta + D b, with
#' b ~ N(0, tau^2 I), carry the information
#' P - P D (tau^-2 I + D' P D)^-1 D' P (valid for singular P).
#'
#' @param m_info Information matrix (p x p, named).
#' @param directions Matrix p x q of the biased directions (rows named as
#'   `m_info`).
#' @param tau Heterogeneity SD (0 returns `m_info`; Inf removes the
#'   information in the directions entirely).
#' @return Symmetric information matrix.
heterogeneity_adjusted_information <- function(m_info, directions, tau) {
  if (!is.finite(tau) && !identical(tau, Inf)) stop("tau must be non-negative")
  if (tau < 0) stop("tau must be non-negative")
  if (tau == 0) return(m_info)
  m_d <- directions[rownames(m_info), , drop = FALSE]
  m_pd <- m_info %*% m_d
  m_core <- t(m_d) %*% m_pd
  if (is.finite(tau)) {
    m_core <- m_core + diag(1 / tau^2, ncol(m_d))
    m_out <- m_info - m_pd %*% solve(m_core, t(m_pd))
  } else {
    m_out <- m_info - m_pd %*% information_pseudo_inverse(m_core)$inverse %*% t(m_pd)
  }
  m_out <- (m_out + t(m_out)) / 2
  dimnames(m_out) <- dimnames(m_info)
  m_out
}

#' Treatment-effect coefficient directions of one endpoint
#'
#' @param coef_names Coefficient names of the endpoint (rownames(fit$res)).
#' @param outcome Endpoint label.
#' @return Matrix p x q (rows `coef_names`), one unit column per Rx main
#'   effect and Rx interaction.
treatment_effect_directions <- function(coef_names, outcome) {
  df_map <- survival_coefficient_map(coef_names, outcome, design_biomarkers())
  v_cols <- df_map$coefficient[df_map$role %in% c("rx_main", "interaction")]
  m_d <- matrix(0, length(coef_names), length(v_cols), dimnames = list(coef_names, v_cols))
  m_d[cbind(match(v_cols, coef_names), seq_along(v_cols))] <- 1
  m_d
}

#' Expected information of one endpoint under a design and patient mix
#'
#' @param design trial_design() (its allocation may be 1 for a single arm).
#' @param n Patients entering the period.
#' @param draw list(os, pfs).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param outcome Endpoint.
#' @param step Quadrature spacing (weeks).
#' @return List with `information` (named p x p) and `expected_events` (data
#'   frame by arm and cell).
endpoint_design_information <- function(design, n, draw, component, reference_population, weights = NULL,
                                        outcome, step = 0.5) {
  df_mix <- design_patient_mix(design, reference_population, weights, component, draw)
  df_mix <- df_mix[df_mix$weight > 0, , drop = FALSE]
  fit <- component[[paste0("original_", outcome)]]
  assert_gamma_flexsurv(fit)
  v_par <- draw[[outcome]]
  m_x <- flexsurv_design_matrix(fit, df_mix)
  v_eta <- as.vector(v_par[["rate"]] + m_x %*% v_par[colnames(m_x)])
  m_e <- expected_gamma_information(v_par[["shape"]], v_eta, design$min_follow_up_weeks,
                                    design_accrual_weeks(design, n), design$randomisation_week, step)
  v_w <- n * df_mix$weight
  m_info <- gamma_information_matrix(cbind(rate = 1, m_x), v_w * m_e[, "aa"], v_w * m_e[, "aeta"],
                                     v_w * m_e[, "etaeta"], rownames(fit$res))
  v_stratum <- paste(df_mix$arm, df_mix$cell, sep = "|")
  df_events <- data.frame(stratum = names(tapply(v_w, v_stratum, sum)),
                          n = as.numeric(tapply(v_w, v_stratum, sum)),
                          events = as.numeric(tapply(v_w * m_e[, "events"], v_stratum, sum)),
                          stringsAsFactors = FALSE)
  df_events$arm <- sub("\\|.*$", "", df_events$stratum)
  df_events$cell <- sub("^.*\\|", "", df_events$stratum)
  list(information = m_info, expected_events = df_events)
}

#' Single-arm trial_design() (experimental or control arm only)
#'
#' trial_design() requires a randomised allocation; a single arm sets the
#' allocation to 1 (experimental) or 0 (control) after construction, so every
#' patient type of the other arm has weight zero.
#'
#' @param name Label.
#' @param cells Enrolled cells.
#' @param min_follow_up_weeks,accrual_weeks,randomisation_week As in trial_design().
#' @param arm "experimental" (default) or "control".
#' @return trial_design() object with `p_experimental` 1 or 0.
single_arm_trial_design <- function(name, cells, min_follow_up_weeks, accrual_weeks, randomisation_week = 0,
                                    arm = c("experimental", "control")) {
  arm <- match.arg(arm)
  design <- trial_design(name, enrolled_cells = cells, p_experimental = 0.5,
                         min_follow_up_weeks = min_follow_up_weeks, accrual_weeks = accrual_weeks,
                         randomisation_week = randomisation_week)
  design$p_experimental <- if (arm == "experimental") 1 else 0
  design
}

#' Expected information of the METIMMOX-2 single arm (S4), both endpoints
#'
#' Early period: every patient, censored at `full_weight_weeks`. Late period:
#' the patients progression-free at that week (experimental curves),
#' left-truncated there, followed to the analysis, with the treatment-effect
#' directions discounted by `tau_m2` (heterogeneity_adjusted_information()).
#' The early and late likelihoods partition the full likelihood, so with
#' `tau_m2 = 0` their sum is the information of the whole follow-up.
#' Approximations: the late cell mix is reweighted by P(PFS > week) under
#' control (design_patient_mix()), and the late OS likelihood conditions on
#' OS > week, not on PFS > week.
#'
#' @param n Patients (METIMMOX2_PATIENTS).
#' @param draw list(os, pfs).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param tau_m2 Heterogeneity SD of the late follow-up.
#' @param accrual_weeks,analysis_weeks Accrual duration and time from the
#'   last entry to the analysis (weeks; defaults from the METIMMOX-2 dates).
#' @param full_weight_weeks Weeks at full weight.
#' @return List with `endpoints` (per endpoint `information`, `early`,
#'   `late`, `expected_events`), `n`, `n_late` and the settings.
single_arm_information <- function(n = METIMMOX2_PATIENTS, draw, component, reference_population, weights = NULL,
                                   tau_m2 = TAU_M2_DEFAULT,
                                   accrual_weeks = as.numeric(difftime(METIMMOX2_ACCRUAL_END, METIMMOX2_ACCRUAL_START,
                                                                       units = "weeks")),
                                   analysis_weeks = as.numeric(difftime(METIMMOX2_ANALYSIS_DATE, METIMMOX2_ACCRUAL_END,
                                                                        units = "weeks")),
                                   full_weight_weeks = METIMMOX2_FULL_WEIGHT_WEEKS) {
  draw <- split_joint_draw(draw, component)
  if (analysis_weeks <= full_weight_weeks) stop("The analysis must follow the full-weight period")
  early <- single_arm_trial_design("S4_early", DESIGN_ALL_CELLS, full_weight_weeks, 0)
  late <- single_arm_trial_design("S4_late", DESIGN_ALL_CELLS, analysis_weeks - full_weight_weeks,
                                  accrual_weeks, full_weight_weeks)
  # Share of the arm still progression-free at the landmark (experimental curves).
  fit_pfs <- component$original_pfs
  df_ref <- reference_population
  v_w <- if (is.null(weights)) rep(1, nrow(df_ref)) else weights
  v_eta <- gamma_linear_predictor(fit_pfs, draw$pfs, df_ref, model_rx_levels(fit_pfs)[2])
  p_late <- sum(v_w * stats::pgamma(full_weight_weeks, shape = exp(draw$pfs[["shape"]]), rate = exp(v_eta),
                                    lower.tail = FALSE)) / sum(v_w)
  n_late <- n * p_late
  step <- design_node_step(analysis_weeks + accrual_weeks)
  endpoints <- lapply(c(os = "os", pfs = "pfs"), function(outcome) {
    e1 <- endpoint_design_information(early, n, draw, component, reference_population, weights, outcome, 0.25)
    e2 <- endpoint_design_information(late, n_late, draw, component, reference_population, weights, outcome, step)
    v_names <- rownames(component[[paste0("original_", outcome)]]$res)
    m_late <- heterogeneity_adjusted_information(e2$information, treatment_effect_directions(v_names, outcome),
                                                 tau_m2)
    list(information = e1$information + m_late, early = e1$information, late = m_late,
         late_unadjusted = e2$information,
         expected_events = rbind(cbind(period = "early", e1$expected_events),
                                 cbind(period = "late", e2$expected_events)))
  })
  list(endpoints = endpoints, n = n, n_late = n_late, tau_m2 = tau_m2, accrual_weeks = accrual_weeks,
       analysis_weeks = analysis_weeks, full_weight_weeks = full_weight_weeks)
}

#' Summary statistic of the single-arm pooling: experimental-arm log rates
#'
#' Per endpoint: `exp_base` (rate + Rx + the baseline covariates at the
#' reference mean, the experimental log rate of a biomarker-negative patient
#' of average age and sex), `exp_<biomarker>` (the biomarker main effect plus
#' its Rx interaction, the biomarker shift of the experimental log rate).
#'
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @return Joint target matrix (rows "<outcome>::exp_*", columns
#'   "<outcome>::<coef>").
s4_targets <- function(component, reference_population, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, nrow(reference_population))
  v_w <- weights / sum(weights)
  l_targets <- lapply(c("os", "pfs"), function(outcome) {
    fit <- component[[paste0("original_", outcome)]]
    v_names <- rownames(fit$res)
    df_map <- survival_coefficient_map(v_names, outcome, design_biomarkers())
    df_ref <- reference_population
    df_ref$Rx <- factor(model_rx_levels(fit)[2], levels = model_rx_levels(fit))
    m_x <- flexsurv_design_matrix(fit, df_ref)
    v_mean <- colSums(v_w * m_x)
    v_base <- stats::setNames(numeric(length(v_names)), v_names)
    v_base[["rate"]] <- 1
    v_cov_baseline <- intersect(df_map$coefficient[df_map$role %in% c("baseline", "rx_main")], colnames(m_x))
    v_base[v_cov_baseline] <- v_mean[v_cov_baseline]
    v_base[df_map$coefficient[df_map$role == "rx_main"]] <- 1
    rows <- list(exp_base = v_base)
    for (b in design_biomarkers()) {
      v <- stats::setNames(numeric(length(v_names)), v_names)
      v[df_map$coefficient[df_map$role %in% c("prognostic", "interaction") & df_map$biomarker %in% b]] <- 1
      rows[[paste0("exp_", b)]] <- v
    }
    m <- do.call(rbind, rows)
    dimnames(m) <- list(paste0(outcome, "::", rownames(m)), paste0(outcome, "::", v_names))
    m
  })
  joint_target_matrix(l_targets, c(paste0("os::", rownames(component$original_os$res)),
                                   paste0("pfs::", rownames(component$original_pfs$res))))
}

#' Linearly independent union of target matrices, per endpoint
#'
#' @param l_targets List of joint target matrices with the same columns.
#' @param tol Rank tolerance.
#' @return Joint target matrix: the rows of each matrix in turn, dropping
#'   rows in the span of those already kept.
combined_target_matrix <- function(l_targets, tol = 1e-8) {
  m_all <- do.call(rbind, l_targets)
  v_keep <- logical(nrow(m_all))
  for (i in seq_len(nrow(m_all))) {
    m_try <- m_all[c(which(v_keep), i), , drop = FALSE]
    v_keep[i] <- qr(t(m_try), tol = tol)$rank == nrow(m_try)
  }
  m_out <- m_all[v_keep, , drop = FALSE]
  if (anyDuplicated(rownames(m_out))) stop("Combined targets repeat a row name")
  m_out
}

#' Covariance of the targets of a design from its information matrices
#'
#' The design_information() logic for given information: per endpoint, the
#' pseudo-inverse of the information, a stop on targets in the null space,
#' and (cross_endpoint = "recolour") the whitened cross-correlation of the
#' targets under the joint coefficient covariance recoloured with the design's
#' marginal covariances.
#'
#' @param l_info list(os = , pfs = ) information matrices named as
#'   rownames(fit$res).
#' @param targets Joint target matrix.
#' @param component Joint sampling component (joint covariance for the
#'   recolouring).
#' @param cross_endpoint "recolour" or "independent".
#' @param label Design label for messages.
#' @return Target covariance matrix (q x q, named by the targets).
information_target_covariance <- function(l_info, targets, component, cross_endpoint = c("recolour", "independent"),
                                          label = "design") {
  cross_endpoint <- match.arg(cross_endpoint)
  v_all_cols <- c(paste0("os::", rownames(component$original_os$res)),
                  paste0("pfs::", rownames(component$original_pfs$res)))
  m_targets <- joint_target_matrix(list(targets), v_all_cols)
  v_outcome <- vapply(rownames(m_targets), function(r) {
    v_out <- unique(sub("::.*$", "", v_all_cols[m_targets[r, ] != 0]))
    if (length(v_out) != 1) stop("Target ", r, " must load on exactly one endpoint")
    v_out
  }, character(1))
  m_cov <- matrix(0, nrow(m_targets), nrow(m_targets), dimnames = list(rownames(m_targets), rownames(m_targets)))
  v_unidentified <- character(0)
  rows <- list()
  for (outcome in c("os", "pfs")) {
    v_rows <- which(v_outcome == outcome)
    if (!length(v_rows)) next
    v_names <- rownames(component[[paste0("original_", outcome)]]$res)
    m_a <- m_targets[v_rows, paste0(outcome, "::", v_names), drop = FALSE]
    pinv <- information_pseudo_inverse(l_info[[outcome]][v_names, v_names])
    if (ncol(pinv$null_space)) {
      v_proj <- apply(abs(m_a %*% pinv$null_space), 1, max) / pmax(1, sqrt(rowSums(m_a^2)))
      v_unidentified <- c(v_unidentified, rownames(m_a)[v_proj >= 1e-6])
    }
    m_cov[v_rows, v_rows] <- m_a %*% pinv$inverse %*% t(m_a)
    rows[[outcome]] <- v_rows
  }
  if (length(v_unidentified)) stop("Design ", label, " does not identify: ", paste(v_unidentified, collapse = ", "))
  if (cross_endpoint == "recolour" && length(rows) == 2 && !is.null(component$joint_covariance$covariance)) {
    m_joint <- component$joint_covariance$covariance
    m_prior_t <- m_targets[, colnames(m_joint), drop = FALSE] %*% m_joint %*%
      t(m_targets[, colnames(m_joint), drop = FALSE])
    i_os <- rows$os; i_pfs <- rows$pfs
    m_r <- covariance_root(m_prior_t[i_os, i_os, drop = FALSE], inverse = TRUE) %*%
      m_prior_t[i_os, i_pfs, drop = FALSE] %*% covariance_root(m_prior_t[i_pfs, i_pfs, drop = FALSE], inverse = TRUE)
    m_cross <- covariance_root(m_cov[i_os, i_os, drop = FALSE]) %*% m_r %*%
      covariance_root(m_cov[i_pfs, i_pfs, drop = FALSE])
    m_cov[i_os, i_pfs] <- m_cross
    m_cov[i_pfs, i_os] <- t(m_cross)
  }
  (m_cov + t(m_cov)) / 2
}

#' Information and targets of a combined design (S4 plus a randomised design)
#'
#' @param info_s4 single_arm_information() result.
#' @param info_x design_information() result of the randomised design.
#' @param targets_s4 s4_targets().
#' @param component Joint sampling component.
#' @param label Label of the combined design.
#' @return List with `information` (per endpoint, the sums), `targets` (union)
#'   and `covariance`.
combined_design_information <- function(info_s4, info_x, targets_s4, component, label) {
  l_info <- lapply(c(os = "os", pfs = "pfs"), function(outcome) {
    info_s4$endpoints[[outcome]]$information + info_x$endpoints[[outcome]]$information
  })
  m_targets <- combined_target_matrix(list(info_x$targets, targets_s4))
  v_order <- order(sub("::.*$", "", rownames(m_targets)) != "os")
  m_targets <- m_targets[v_order, , drop = FALSE]
  list(information = l_info, targets = m_targets,
       covariance = information_target_covariance(l_info, m_targets, component, label = label))
}

# -------------------------------------------------------------------------------
# Information of the biomarker-stratified, non-randomised METIMMOX-3 (S3)
# -------------------------------------------------------------------------------

#' Biomarker whose positive cells are a stratified design's experimental cells
#'
#' @param spec Design specification of type "stratified".
#' @param biomarkers Biomarker identifiers (cell-code order).
#' @return Biomarker identifier; stops unless exactly one biomarker is 1 in
#'   every experimental cell and 0 in every control cell.
stratifying_biomarker <- function(spec, biomarkers = design_biomarkers()) {
  v_ctl <- setdiff(spec$enrolled_cells, spec$experimental_cells)
  v_hit <- vapply(seq_along(biomarkers), function(j) {
    all(substr(spec$experimental_cells, j, j) == "1") && all(substr(v_ctl, j, j) == "0")
  }, logical(1))
  if (sum(v_hit) != 1 || !length(v_ctl)) {
    stop("Design ", spec$id, ": the experimental cells must be the positive cells of exactly one biomarker")
  }
  biomarkers[v_hit]
}

#' Coefficient direction of a change of control regimen: the baseline log rate
#'
#' A standard of care other than FLOX shifts the log rate of every control
#' patient by the same amount (heterogeneity_adjusted_information()).
#'
#' @param coef_names Coefficient names of one endpoint (rownames(fit$res)).
#' @return Matrix p x 1 (rows `coef_names`), the unit vector of "rate".
control_regimen_directions <- function(coef_names) {
  if (!"rate" %in% coef_names) stop("The model has no 'rate' parameter")
  m_d <- matrix(0, length(coef_names), 1, dimnames = list(coef_names, "rate"))
  m_d["rate", 1] <- 1
  m_d
}

#' Summary statistic of a stratified design: what its two arms identify
#'
#' Per endpoint, with s the stratifying biomarker (TMB/BRAF) and o each other
#' biomarker (CRP): `ctl_base` (rate plus the baseline covariates at the
#' reference mean: the control log rate of a patient negative for every
#' biomarker), `ctl_<o>` (the prognostic coefficient of o, from the control
#' cells), `exp_<s>` (`ctl_base` plus Rx, the prognostic coefficient of s and
#' Rx:s: the experimental log rate of the s-positive, o-negative cell) and
#' `rx_x_<o>` (Rx:o: the o shift of the experimental cells minus that of the
#' control cells, identified because the model has no o x s term). The
#' treatment effect in either s stratum is not identified.
#'
#' @param spec Design specification of type "stratified".
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param biomarkers Biomarker identifiers.
#' @return Joint target matrix (rows "<outcome>::<target>", columns
#'   "<outcome>::<coef>").
stratified_targets <- function(spec, component, reference_population, weights = NULL,
                               biomarkers = design_biomarkers()) {
  s <- stratifying_biomarker(spec, biomarkers)
  v_other <- setdiff(biomarkers, s)
  if (is.null(weights)) weights <- rep(1, nrow(reference_population))
  v_w <- weights / sum(weights)
  l_targets <- lapply(c("os", "pfs"), function(outcome) {
    fit <- component[[paste0("original_", outcome)]]
    v_names <- rownames(fit$res)
    df_map <- survival_coefficient_map(v_names, outcome, biomarkers)
    coef_of <- function(role, b = NULL) {
      v <- df_map$coefficient[df_map$role == role & (is.null(b) | df_map$biomarker %in% b)]
      if (length(v) != 1) stop("No unique ", outcome, " coefficient with role ", role, " ", b)
      v
    }
    unit <- function(coef) { v <- stats::setNames(numeric(length(v_names)), v_names); v[coef] <- 1; v }
    df_ref <- reference_population
    df_ref$Rx <- factor(model_rx_levels(fit)[1], levels = model_rx_levels(fit))
    m_x <- flexsurv_design_matrix(fit, df_ref)
    v_mean <- colSums(v_w * m_x)
    v_baseline <- intersect(df_map$coefficient[df_map$role == "baseline"], colnames(m_x))
    v_ctl <- unit("rate")
    v_ctl[v_baseline] <- v_mean[v_baseline]
    rows <- list(ctl_base = v_ctl)
    for (b in v_other) rows[[paste0("ctl_", b)]] <- unit(coef_of("prognostic", b))
    v_exp <- v_ctl
    v_exp[c(coef_of("rx_main"), coef_of("prognostic", s), coef_of("interaction", s))] <- 1
    rows[[paste0("exp_", s)]] <- v_exp
    for (b in v_other) rows[[paste0("rx_x_", b)]] <- unit(coef_of("interaction", b))
    m <- do.call(rbind, rows)
    dimnames(m) <- list(paste0(outcome, "::", rownames(m)), paste0(outcome, "::", v_names))
    m
  })
  joint_target_matrix(l_targets, c(paste0("os::", rownames(component$original_os$res)),
                                   paste0("pfs::", rownames(component$original_pfs$res))))
}

#' Expected information of a stratified design (S3), both endpoints
#'
#' n patients enrolled uniformly over n / `accrual_per_year` years and
#' followed to `min_follow_up_years` after the last inclusion. The
#' experimental cells contribute experimental-arm data from week 0 (the
#' regimen starts with the same two chemotherapy cycles); their
#' treatment-effect directions (Rx and both Rx interactions, per endpoint;
#' treatment_effect_directions()) may differ from METIMMOX's regimen by
#' b ~ N(0, tau_m3^2 I) (METIMMOX-3: tislelizumab for nivolumab, bevacizumab,
#' optional FOLFOX; heterogeneity_adjusted_information()). The control cells
#' contribute control-arm data in two periods that partition their
#' likelihood: the start-up (every patient, censored at `startup_weeks`)
#' and, for the patients still at risk on the endpoint at that week (control
#' curves; the cell mix reweighted by the endpoint's own survival), the
#' standard of care after it, left-truncated there. The control regimen
#' (FLOX/FOLFOX with bevacizumab in the start-up, then FLOX or FOLFOX +/-
#' bevacizumab) may differ from METIMMOX's FLOX in its baseline log rate by
#' one b ~ N(0, tau_soc^2) shared by both periods
#' (heterogeneity_adjusted_information() of their sum on
#' control_regimen_directions()). With `tau_soc = 0` the two periods sum to
#' the information of the whole control follow-up. The two biases are
#' independent of each other and between endpoints.
#'
#' @param spec Design specification of type "stratified".
#' @param n Enrolled patients (scalar).
#' @param draw Draw defining the expected information (list(os, pfs)).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param tau_soc Heterogeneity SD of the control regimen (whole control
#'   follow-up).
#' @param tau_m3 Heterogeneity SD of the experimental regimen's treatment
#'   effects.
#' @param cross_endpoint "recolour" or "independent".
#' @return List with `design`, `n`, `n_experimental`, `n_control`,
#'   `accrual_weeks`, `min_follow_up_weeks`, `startup_weeks`, `tau_soc`,
#'   `tau_m3`, `endpoints` (per endpoint `information`, the parts
#'   `experimental` (adjusted by tau_m3), `experimental_unadjusted`, `control`
#'   (start-up plus standard of care, adjusted by tau_soc), `control_startup`
#'   and `control_soc` (both unadjusted), `n_at_risk_soc` and
#'   `expected_events` by period, arm and cell), `targets`
#'   (stratified_targets()) and `covariance` (V_d(n) of the targets).
stratified_design_information <- function(spec, n, draw, component, reference_population, weights = NULL,
                                          tau_soc = TAU_SOC_DEFAULT, tau_m3 = TAU_M3_DEFAULT,
                                          cross_endpoint = c("recolour", "independent")) {
  if (!identical(spec$type, "stratified")) stop("stratified_design_information() is for stratified designs")
  cross_endpoint <- match.arg(cross_endpoint)
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0) stop("n must be one positive number")
  draw <- split_joint_draw(draw, component)
  alloc <- design_allocation(spec, reference_population, weights, component, draw)
  v_ctl_cells <- setdiff(spec$enrolled_cells, spec$experimental_cells)
  accrual <- 52 * n / spec$accrual_per_year
  follow <- 52 * spec$min_follow_up_years
  startup <- spec$startup_weeks
  if (follow <= startup) stop("The minimum follow-up must exceed the start-up period")
  n_exp <- n * alloc$p_experimental
  n_ctl <- n - n_exp
  step <- design_node_step(follow + accrual)
  d_exp <- single_arm_trial_design(paste0(spec$id, "_experimental"), spec$experimental_cells, follow, accrual)
  d_startup <- single_arm_trial_design(paste0(spec$id, "_control_startup"), v_ctl_cells, startup, 0, arm = "control")
  d_soc <- single_arm_trial_design(paste0(spec$id, "_control_soc"), v_ctl_cells, follow - startup, accrual, startup,
                                   arm = "control")
  v_w <- if (is.null(weights)) rep(1, nrow(reference_population)) else weights
  # Survival at the end of the start-up under control, per reference patient.
  survival_at_startup <- function(outcome) {
    fit <- component[[paste0("original_", outcome)]]
    v_eta <- gamma_linear_predictor(fit, draw[[outcome]], reference_population, model_rx_levels(fit)[1])
    stats::pgamma(startup, shape = exp(draw[[outcome]][["shape"]]), rate = exp(v_eta), lower.tail = FALSE)
  }
  v_s_pfs <- survival_at_startup("pfs")
  v_ctl <- design_cell_codes(reference_population) %in% v_ctl_cells
  endpoints <- lapply(c(os = "os", pfs = "pfs"), function(outcome) {
    v_names <- rownames(component[[paste0("original_", outcome)]]$res)
    v_s <- survival_at_startup(outcome)
    p_at_risk <- sum((v_w * v_s)[v_ctl]) / sum(v_w[v_ctl])
    # design_patient_mix() reweights a landmark design by P(PFS > L); dividing
    # by it here leaves the endpoint's own survival as the selection.
    v_w_soc <- v_w * v_s / v_s_pfs
    e_exp <- endpoint_design_information(d_exp, n_exp, draw, component, reference_population, weights, outcome, step)
    e_startup <- endpoint_design_information(d_startup, n_ctl, draw, component, reference_population, weights,
                                             outcome, 0.25)
    e_soc <- endpoint_design_information(d_soc, n_ctl * p_at_risk, draw, component, reference_population, v_w_soc,
                                         outcome, step)
    m_ctl <- heterogeneity_adjusted_information(e_startup$information + e_soc$information,
                                                control_regimen_directions(v_names), tau_soc)
    m_exp <- heterogeneity_adjusted_information(e_exp$information, treatment_effect_directions(v_names, outcome),
                                                tau_m3)
    list(information = m_exp + m_ctl, experimental = m_exp, experimental_unadjusted = e_exp$information,
         control = m_ctl, control_startup = e_startup$information, control_soc = e_soc$information,
         n_at_risk_soc = n_ctl * p_at_risk,
         expected_events = rbind(cbind(period = "whole", e_exp$expected_events),
                                 cbind(period = "startup", e_startup$expected_events),
                                 cbind(period = "soc", e_soc$expected_events)))
  })
  m_targets <- stratified_targets(spec, component, reference_population, weights)
  list(design = spec$id, n = n, n_experimental = n_exp, n_control = n_ctl, accrual_weeks = accrual,
       min_follow_up_weeks = follow, startup_weeks = startup, tau_soc = tau_soc, tau_m3 = tau_m3,
       endpoints = endpoints,
       targets = m_targets,
       covariance = information_target_covariance(lapply(endpoints, function(e) e$information), m_targets,
                                                  component, cross_endpoint, label = spec$id))
}

# -------------------------------------------------------------------------------
# Power-based sample sizes (S6)
# -------------------------------------------------------------------------------

#' Variance of a log effect that gives the requested power (Schoenfeld)
#'
#' A two-sided test at level alpha of a log effect beta has power
#' `power` when Var(beta_hat) = beta^2 / (z_(1 - alpha/2) + z_power)^2.
#'
#' @param rate_ratio Effect to detect (rate or hazard ratio).
#' @param power,alpha Power and two-sided level.
#' @return Target variance.
power_variance_target <- function(rate_ratio = POWER_RATE_RATIO, power = POWER_LEVEL, alpha = POWER_ALPHA) {
  if (!is.finite(rate_ratio) || rate_ratio <= 0 || rate_ratio == 1) stop("rate_ratio must be positive and not 1")
  log(rate_ratio)^2 / (stats::qnorm(1 - alpha / 2) + stats::qnorm(power))^2
}

#' Events required with balanced strata (Schoenfeld's closed forms)
#'
#' @param rate_ratio,power,alpha As in power_variance_target().
#' @param type "main" (two arms, 1:1) or "interaction" (two arms x two
#'   biomarker levels with equal events in the four strata).
#' @return Total events: 4 z^2 / log(RR)^2 for a main effect and 16 z^2 /
#'   log(RR)^2 for an interaction (four times as many).
schoenfeld_events <- function(rate_ratio = POWER_RATE_RATIO, power = POWER_LEVEL, alpha = POWER_ALPHA,
                              type = c("main", "interaction")) {
  type <- match.arg(type)
  (if (type == "main") 4 else 16) / power_variance_target(rate_ratio, power, alpha)
}

#' Heuristic variance of a design's power target from its expected events
#'
#' @param df_events Expected events by `arm` and `cell` (one endpoint).
#' @param spec Design specification (`power_target`, `power_biomarkers`).
#' @return For an interaction: the largest sum(1/d) over the arm x biomarker
#'   strata of the `power_biomarkers` (cox_interaction_heuristic()); for a
#'   conditional effect: sum over arms of 1 / (events in that arm).
power_target_variance <- function(df_events, spec) {
  if (identical(spec$power_target, "interaction")) {
    info <- list(endpoints = list(x = list(expected_events = df_events)), design = list(biomarkers = design_biomarkers()))
    max(vapply(spec$power_biomarkers, function(b) cox_interaction_heuristic(info, "x", b), numeric(1)))
  } else if (identical(spec$power_target, "conditional")) {
    v_d <- tapply(df_events$events, df_events$arm, sum)
    if (length(v_d) != 2) stop("A conditional effect needs events in both arms")
    sum(1 / v_d)
  } else {
    stop("Design ", spec$id, " has no power target")
  }
}

#' Expected events of a randomised design with n patients, one endpoint
#'
#' @param design trial_design().
#' @param n Randomised patients.
#' @param draw,component,reference_population,weights As in design_information().
#' @param outcome Endpoint.
#' @return Data frame of expected events by `arm` and `cell`.
design_expected_events <- function(design, n, draw, component, reference_population, weights = NULL,
                                   outcome = POWER_ENDPOINT) {
  draw <- split_joint_draw(draw, component)
  step <- design_node_step(design$min_follow_up_weeks + design_accrual_weeks(design, n))
  endpoint_design_information(design, n, draw, component, reference_population, weights, outcome, step)$expected_events
}

#' Power-based sample size of a randomised candidate from Schoenfeld's events
#'
#' The cross-check of power_based_n_information() (ruling (4) of 11 October
#' 2026). The smallest n whose expected events, through the draw's own curves with
#' the design's accrual and follow-up, give the power target the variance of
#' power_variance_target(). Schoenfeld's formula is applied to the gamma rate
#' ratio as if it were a hazard ratio (stated approximation).
#'
#' @param spec Design specification of type "rct".
#' @param randomised_per_year Randomised patients per year.
#' @param draw,component,reference_population,weights Inputs of the expected events.
#' @param rate_ratio,power,alpha Power settings.
#' @param outcome Endpoint.
#' @param n_max Upper end of the search.
#' @return List with `n` (integer, NA when above n_max), `events`,
#'   `variance`, `target_variance`.
power_based_n <- function(spec, randomised_per_year, draw, component, reference_population, weights = NULL,
                          rate_ratio = POWER_RATE_RATIO, power = POWER_LEVEL, alpha = POWER_ALPHA,
                          outcome = POWER_ENDPOINT, n_max = 1e5) {
  design <- as_trial_design(spec, randomised_per_year)
  target <- power_variance_target(rate_ratio, power, alpha)
  variance_at <- function(n) power_target_variance(design_expected_events(design, n, draw, component,
                                                                          reference_population, weights, outcome), spec)
  if (variance_at(n_max) > target) return(list(n = NA_integer_, events = NA_real_, variance = NA_real_,
                                               target_variance = target))
  lo <- 2; hi <- n_max
  if (variance_at(lo) <= target) hi <- lo
  # Bisection on integers: variance(n) is decreasing in n.
  while (hi - lo > 1) {
    mid <- floor(sqrt(lo * hi))
    if (mid <= lo) mid <- lo + 1
    if (mid >= hi) mid <- hi - 1
    if (variance_at(mid) <= target) hi <- mid else lo <- mid
  }
  df_events <- design_expected_events(design, hi, draw, component, reference_population, weights, outcome)
  list(n = as.integer(hi), events = sum(df_events$events), variance = power_target_variance(df_events, spec),
       target_variance = target)
}


# -------------------------------------------------------------------------------
# Power from the engine's standard error and the decision-relevant effect
# (ruling (4) of 11 October 2026)
# -------------------------------------------------------------------------------

#' Hypotheses of the power-based sizes
POWER_HYPOTHESES <- c("conventional", "decision_relevant")
#' Methods of the power-based sizes: the engine's standard error (primary)
#' and Schoenfeld's event formula (cross-check)
POWER_METHODS <- c("engine", "schoenfeld")
#' Bracket of the shift of the target coefficient searched by
#' decision_relevant_effect() (log rate scale: rate ratios 0.05 to 20)
DECISION_EFFECT_BRACKET <- c(-3, 3)

#' Policy pattern that treats the positive cells of one biomarker
#'
#' @param biomarker Biomarker identifier.
#' @param cells Cell codes (biomarker digits in design_biomarkers() order).
#' @param biomarkers Biomarker identifiers.
#' @return Pattern string, one digit per cell (1 = nivolumab), e.g. "0011"
#'   for CRP and "0101" for TMB/BRAF.
guided_policy_pattern <- function(biomarker, cells = DESIGN_ALL_CELLS, biomarkers = design_biomarkers()) {
  j <- match(biomarker, biomarkers)
  if (is.na(j)) stop("Unknown biomarker ", biomarker)
  paste0(substr(cells, j, j), collapse = "")
}

#' Power targets of a randomised candidate, one per biomarker it is powered for
#'
#' S1 is powered for each treatment-by-biomarker interaction (target: the
#' interaction coefficient); an enriched design (S2, S3) for the conditional
#' treatment effect in its enrolled cells, the randomisation-mix-weighted mean
#' of the enrolled cells' log rate ratios (Rx + c Rx:crp + m Rx:tmb_braf). In
#' both cases the target loads with weight one on the interaction coefficient
#' of the biomarker (every enrolled cell of an enriched design is positive for
#' it), so shifting that coefficient by delta shifts the target by delta and
#' leaves the biomarker-negative cells unchanged; the corresponding guided
#' policy treats the biomarker-positive cells (guided_policy_pattern()).
#'
#' @param spec Design specification (`power_target`, `power_biomarkers`,
#'   `enrolled_cells`).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param draw Draw defining the randomisation mix (list(os, pfs)).
#' @param outcome Endpoint of the power calculation.
#' @param biomarkers Biomarker identifiers.
#' @return List over biomarkers of list(`biomarker`, `target` (1 x p matrix
#'   over "<outcome>::<coefficient>" with row name "<outcome>::power_<b>"),
#'   `target_mean` (value at `draw`), `shift_coefs` (named by endpoint: the
#'   interaction coefficient of the biomarker), `policy` (pattern),
#'   `power_target`).
design_power_targets <- function(spec, component, reference_population, weights = NULL, draw,
                                 outcome = POWER_ENDPOINT, biomarkers = design_biomarkers()) {
  draw <- split_joint_draw(draw, component)
  interaction_of <- function(o, b) {
    df_map <- survival_coefficient_map(rownames(component[[paste0("original_", o)]]$res), o, biomarkers)
    v <- df_map$coefficient[df_map$role == "interaction" & df_map$biomarker == b]
    if (length(v) != 1) stop("No unique ", o, " interaction for ", b)
    v
  }
  coef_names <- rownames(component[[paste0("original_", outcome)]]$res)
  df_map <- survival_coefficient_map(coef_names, outcome, biomarkers)
  rx <- df_map$coefficient[df_map$role == "rx_main"]
  if (identical(spec$power_target, "interaction")) {
    v_biomarkers <- spec$power_biomarkers
  } else if (identical(spec$power_target, "conditional")) {
    v_pos <- vapply(seq_along(biomarkers), function(j) all(substr(spec$enrolled_cells, j, j) == "1"), logical(1))
    if (sum(v_pos) != 1) stop("Design ", spec$id, ": the enrolled cells are not positive for exactly one biomarker")
    v_biomarkers <- biomarkers[v_pos]
    enrol <- design_enrolment_fraction(spec, reference_population, weights, component, draw)
    v_mix <- enrol$cell_mix_randomisation
  } else {
    stop("Design ", spec$id, " has no power target")
  }
  lapply(stats::setNames(v_biomarkers, v_biomarkers), function(b) {
    v_row <- stats::setNames(numeric(length(coef_names)), coef_names)
    if (identical(spec$power_target, "interaction")) {
      v_row[interaction_of(outcome, b)] <- 1
    } else {
      for (cell in names(v_mix)) {
        v_digits <- as.numeric(strsplit(cell, "")[[1]])
        v_row[rx] <- v_row[rx] + v_mix[[cell]]
        for (j in seq_along(biomarkers)) {
          v_int <- interaction_of(outcome, biomarkers[j])
          v_row[v_int] <- v_row[v_int] + v_mix[[cell]] * v_digits[j]
        }
      }
    }
    m_target <- matrix(v_row, nrow = 1, dimnames = list(paste0(outcome, "::power_", b),
                                                        paste0(outcome, "::", coef_names)))
    list(biomarker = b, target = m_target, target_mean = sum(v_row * draw[[outcome]][coef_names]),
         shift_coefs = c(os = interaction_of("os", b), pfs = interaction_of("pfs", b)),
         policy = guided_policy_pattern(b, biomarkers = biomarkers), power_target = spec$power_target)
  })
}

#' Power-based sample size from the engine's standard error
#'
#' The smallest n at which the standard error of the design's target from
#' design_information() (the expected information of the gamma models, with
#' shape, age, sex and the other coefficients estimated) gives power `power`
#' for a true log effect `effect` in a two-sided test at level `alpha`:
#' Var(target_hat)(n) <= effect^2 / (z_(1 - alpha/2) + z_power)^2. With
#' several target rows every row must reach the power (the largest variance
#' counts). Doubling from n = 32 brackets the size, then integer bisection,
#' since the variance falls with n (accrual and
#' follow-up change with n through the design's accrual rate).
#'
#' @param design trial_design().
#' @param target Target matrix (rows: targets; columns joint
#'   "<outcome>::<coefficient>") or a named numeric vector (one target).
#' @param effect True log effect (log rate ratio) to detect; not 0.
#' @param power,alpha Power and two-sided level.
#' @param draw,component,reference_population,weights Inputs of
#'   design_information().
#' @param n_max Upper end of the search.
#' @return List with `n` (integer; NA when above n_max), `variance` and `se`
#'   at n, `target_variance`, `events` (expected events at n on the target's
#'   endpoint) and `method` ("engine").
power_based_n_information <- function(design, target, effect, power = POWER_LEVEL, alpha = POWER_ALPHA, draw,
                                      component, reference_population, weights = NULL, n_max = 1e5) {
  if (!is.matrix(target)) target <- matrix(target, nrow = 1, dimnames = list("target", names(target)))
  if (!is.finite(effect) || effect == 0) stop("effect must be a finite, non-zero log effect")
  v_outcome <- unique(sub("::.*$", "", colnames(target)[colSums(abs(target)) > 0]))
  if (length(v_outcome) != 1) stop("The power target must load on one endpoint")
  target_var <- power_variance_target(exp(effect), power, alpha)
  info_at <- function(n) {
    design_information(design, n, draw, component, reference_population, weights, targets = target,
                       step = design_node_step(design$min_follow_up_weeks + design_accrual_weeks(design, n)))
  }
  variance_at <- function(n) max(diag(info_at(n)$covariance))
  empty <- list(n = NA_integer_, variance = NA_real_, se = NA_real_, target_variance = target_var,
                events = NA_real_, method = "engine")
  # Doubling from a small n brackets the size without evaluating the
  # information at absurd accrual durations (n_max patients).
  lo <- 1; hi <- 32
  while (variance_at(hi) > target_var) {
    if (hi >= n_max) return(empty)
    lo <- hi
    hi <- min(2 * hi, n_max)
  }
  while (hi - lo > 1) {
    mid <- floor(sqrt(lo * hi))
    if (mid <= lo) mid <- lo + 1
    if (mid >= hi) mid <- hi - 1
    if (variance_at(mid) <= target_var) hi <- mid else lo <- mid
  }
  info <- info_at(hi)
  v <- max(diag(info$covariance))
  list(n = as.integer(hi), variance = v, se = sqrt(v), target_variance = target_var,
       events = sum(info$endpoints[[v_outcome]]$expected_events$events), method = "engine")
}

#' Expected incremental NMB of a guided policy versus standard of care
#'
#' Deterministic evaluation at the coefficient vector `par` (list(os, pfs) of
#' full parameter vectors): the flexsurvreg objects are rebuilt with
#' build_sampled_flexsurv_model(), the cell curves predicted on the model
#' grid (predict_cell_curves(), the cell path of the PSA) and the cell
#' outcomes computed (calculate_cell_outcomes(); PFS capped at OS per cell).
#' Costs are linear in the nivolumab price with slope the discounted
#' administrations, so one evaluation serves every price.
#'
#' @param par list(os, pfs) of full coefficient vectors.
#' @param component Joint sampling component (the fits to rebuild).
#' @param params Parameter list with the time grid, discounting, utilities,
#'   unit costs, schedules and test prices (l_params_base).
#' @param population Prediction population (Age, sex, Rx, biomarkers).
#' @param weights Patient weights (NULL: equal).
#' @return List with `values` (cell x arm x outcome), `mass` (cell weights),
#'   `cells` (cell definition) and `price_ref` (params$c_drug_nivo).
guided_policy_cell_outcomes <- function(par, component, params, population, weights = NULL) {
  model_os <- build_sampled_flexsurv_model(component$original_os, par$os)
  model_pfs <- build_sampled_flexsurv_model(component$original_pfs, par$pfs)
  v_time <- seq(0, params$time_horizon)
  cell_curves <- predict_cell_curves(model_os, model_pfs, population, v_time, weights)
  dw <- discount_weights(params, length(v_time), params$cl)
  out <- calculate_cell_outcomes(params, cell_curves, dw$cost, dw$effect, params$cl)
  list(values = out$values, mass = cell_curves$mass, cells = cell_curves$cells, price_ref = params$c_drug_nivo)
}

#' Incremental NMB of a policy pattern from cell outcomes
#'
#' @param cell_outcomes guided_policy_cell_outcomes() result.
#' @param policy Pattern ("0011") or policy id ("nivo_0011").
#' @param price Nivolumab price (vector).
#' @param wtp Willingness-to-pay threshold.
#' @param params Parameter list with the test prices.
#' @return Incremental NMB versus treating nobody, per price, net of the
#'   tests the policy needs (policy_conditions_on()).
policy_incremental_nmb <- function(cell_outcomes, policy, price, wtp, params) {
  pattern <- sub("^nivo_", "", policy)
  v_cells <- rownames(cell_outcomes$cells)
  v_assign <- as.integer(strsplit(pattern, "")[[1]])
  if (length(v_assign) != length(v_cells)) stop("Policy ", policy, " does not have one digit per cell")
  a <- cell_outcomes$values
  v_w <- cell_outcomes$mass[v_cells]
  v_treat <- v_cells[v_assign == 1 & v_w > 0]
  d_q <- sum(v_w[v_treat] * (a[v_treat, "experimental", "qalys"] - a[v_treat, "control", "qalys"]))
  d_c <- sum(v_w[v_treat] * (a[v_treat, "experimental", "cost"] - a[v_treat, "control", "cost"]))
  d_admin <- sum(v_w[v_treat] * (a[v_treat, "experimental", "nivo_admin"] - a[v_treat, "control", "nivo_admin"]))
  df_cells <- cell_outcomes$cells
  v_tested <- names(df_cells)[vapply(names(df_cells), function(b) policy_conditions_on(v_assign, df_cells, b),
                                     logical(1))]
  test_cost <- sum(vapply(v_tested, function(b) params[[biomarker_cost_key(b)[[1]]]], numeric(1)))
  wtp * d_q - (d_c + (price - cell_outcomes$price_ref) * d_admin + test_cost)
}

#' Decision-relevant effect: the target effect at which a guided policy breaks even
#'
#' The value of the design's target (the interaction for S1, the conditional
#' effect for S2 and S3; design_power_targets()) at which the expected
#' incremental NMB of the guided policy versus standard of care is zero at
#' the nivolumab price and threshold. The target coefficient (the biomarker's
#' interaction, `shift_coefs`) is shifted by delta on the endpoints in
#' `endpoints` from the mean coefficients and the model evaluated
#' deterministically (guided_policy_cell_outcomes()); uniroot() finds delta
#' in `bracket`. Evaluations are shared across prices. If the incremental
#' NMB has the same sign at both ends of the bracket there is no root: the
#' effect is NA with the reason (the policy is not cost-effective for any
#' effect in the bracket, or is cost-effective throughout).
#'
#' @param policy Pattern ("0011") or policy id ("nivo_0011").
#' @param price Nivolumab price(s).
#' @param wtp Willingness-to-pay threshold.
#' @param component Joint sampling component (mean coefficients: its fits).
#' @param params Parameter list (l_params_base: grid, costs, schedules, tests).
#' @param population Prediction population.
#' @param shift_coefs Coefficient shifted, named by endpoint (e.g.
#'   c(os = "RxExperimental arm:crp1", pfs = "RxExperimental arm:crp1")).
#' @param weights Patient weights (NULL: equal).
#' @param target_mean Value of the design's target at the mean coefficients on
#'   `outcome` (default: the shifted coefficient's value); the reported
#'   effect is target_mean + delta.
#' @param endpoints Endpoints whose coefficient is shifted (default both).
#' @param outcome Endpoint of the reported effect.
#' @param bracket Search interval of delta (log rate scale).
#' @param tol uniroot() tolerance.
#' @return Data frame, one row per price: `policy`, `price`, `wtp`,
#'   `inmb_at_mean`, `shift` (delta at the root), `effect` (target value at
#'   the root on `outcome`), `effect_at_mean`, `rate_ratio` (exp(effect)),
#'   `method` ("uniroot" or "no_root") and `reason`.
decision_relevant_effect <- function(policy, price, wtp, component, params, population, shift_coefs,
                                     weights = NULL, target_mean = NULL, endpoints = c("os", "pfs"),
                                     outcome = POWER_ENDPOINT, bracket = DECISION_EFFECT_BRACKET, tol = 1e-4) {
  base <- component_mean_draw(component)
  if (is.null(target_mean)) target_mean <- base[[outcome]][[shift_coefs[[outcome]]]]
  memo <- new.env()
  outcomes_at <- function(delta) {
    key <- format(delta, digits = 15)
    if (!is.null(memo[[key]])) return(memo[[key]])
    par <- base
    for (o in endpoints) par[[o]][shift_coefs[[o]]] <- par[[o]][shift_coefs[[o]]] + delta
    memo[[key]] <- guided_policy_cell_outcomes(par, component, params, population, weights)
    memo[[key]]
  }
  rows <- lapply(price, function(p) {
    inmb_at <- function(delta) policy_incremental_nmb(outcomes_at(delta), policy, p, wtp, params)
    f0 <- inmb_at(0)
    f_lo <- inmb_at(bracket[1]); f_hi <- inmb_at(bracket[2])
    out <- data.frame(policy = paste0("nivo_", sub("^nivo_", "", policy)), price = p, wtp = wtp,
                      inmb_at_mean = f0, shift = NA_real_, effect = NA_real_, effect_at_mean = target_mean,
                      rate_ratio = NA_real_, method = "no_root", reason = "", stringsAsFactors = FALSE)
    if (sign(f_lo) == sign(f_hi)) {
      out$reason <- if (f_lo < 0) {
        sprintf("not cost-effective for any shift in [%g, %g] (incremental NMB %.0f at the most favourable end)",
                bracket[1], bracket[2], max(f_lo, f_hi))
      } else {
        sprintf("cost-effective for every shift in [%g, %g]", bracket[1], bracket[2])
      }
      return(out)
    }
    root <- stats::uniroot(inmb_at, bracket, f.lower = f_lo, f.upper = f_hi, tol = tol)$root
    out$shift <- root
    out$effect <- target_mean + root
    out$rate_ratio <- exp(out$effect)
    out$method <- "uniroot"
    out$reason <- "incremental NMB versus standard of care is zero"
    out
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

#' Decision-relevant effects of every randomised candidate's power targets
#'
#' One root search per biomarker (S1's interaction and the enriched design of
#' the same biomarker shift the same coefficient for the same policy, so the
#' shift is shared and only the reported target value differs).
#'
#' @param specs candidate_designs().
#' @param prices Nivolumab prices.
#' @param wtp Willingness-to-pay threshold.
#' @param component Joint sampling component (mean coefficients).
#' @param params Parameter list (l_params_base).
#' @param population,weights Prediction population and weights.
#' @param outcome Endpoint of the reported effect.
#' @return Data frame per design, biomarker and price: decision_relevant_effect()
#'   columns plus `design`, `biomarker`, `power_target`.
design_decision_effects <- function(specs, prices, wtp, component, params, population, weights = NULL,
                                    outcome = POWER_ENDPOINT) {
  draw <- component_mean_draw(component)
  l_targets <- list()
  for (spec in specs) {
    if (spec$type != "rct" || is.na(spec$power_target)) next
    for (tg in design_power_targets(spec, component, population, weights, draw, outcome)) {
      l_targets[[length(l_targets) + 1]] <- c(tg, design = spec$id)
    }
  }
  v_biomarkers <- unique(vapply(l_targets, function(t) t$biomarker, character(1)))
  l_roots <- lapply(stats::setNames(v_biomarkers, v_biomarkers), function(b) {
    tg <- l_targets[[which(vapply(l_targets, function(t) t$biomarker, character(1)) == b)[1]]]
    decision_relevant_effect(tg$policy, prices, wtp, component, params, population, tg$shift_coefs, weights,
                             target_mean = 0, outcome = outcome)
  })
  rows <- lapply(l_targets, function(tg) {
    df <- l_roots[[tg$biomarker]]
    df$effect_at_mean <- tg$target_mean
    df$effect <- tg$target_mean + df$shift
    df$rate_ratio <- exp(df$effect)
    cbind(design = tg$design, biomarker = tg$biomarker, power_target = tg$power_target, df,
          stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

#' Power-based sizes of one randomised candidate by method and hypothesis
#'
#' For each biomarker the design is powered for: the conventional hypothesis
#' (rate ratio POWER_RATE_RATIO) and, when `decision_effects` are given, the
#' decision-relevant effect at each price (design_decision_effects()); each
#' by the engine's standard error (power_based_n_information(), primary) and
#' by Schoenfeld's events (power_based_n(), cross-check). A decision-relevant
#' effect that is NA or not a benefit (>= 0) gives n = NA with the reason.
#'
#' @param spec Design specification of type "rct".
#' @param randomised_per_year Randomised patients per year.
#' @param draw,component,reference_population,weights Inputs of the information.
#' @param decision_effects design_decision_effects() rows or NULL.
#' @param methods Subset of POWER_METHODS.
#' @param power,alpha Power and two-sided level.
#' @param outcome Endpoint.
#' @param n_max Upper end of the search.
#' @return Data frame: `design`, `biomarker`, `hypothesis`, `price` (NA for the
#'   conventional hypothesis), `effect`, `rate_ratio`, `method`, `n`, `events`,
#'   `se` (engine), `target_variance`, `reason`.
design_power_sizes <- function(spec, randomised_per_year, draw, component, reference_population, weights = NULL,
                               decision_effects = NULL, methods = POWER_METHODS, power = POWER_LEVEL,
                               alpha = POWER_ALPHA, outcome = POWER_ENDPOINT, n_max = 1e5) {
  draw <- split_joint_draw(draw, component)
  design <- as_trial_design(spec, randomised_per_year)
  targets <- design_power_targets(spec, component, reference_population, weights, draw, outcome)
  rows <- list()
  for (tg in targets) {
    df_h <- data.frame(hypothesis = "conventional", price = NA_real_, effect = log(POWER_RATE_RATIO),
                       reason = "", stringsAsFactors = FALSE)
    if (!is.null(decision_effects)) {
      d <- decision_effects[decision_effects$design == spec$id & decision_effects$biomarker == tg$biomarker, ,
                            drop = FALSE]
      if (nrow(d)) df_h <- rbind(df_h, data.frame(hypothesis = "decision_relevant", price = d$price,
                                                  effect = d$effect, reason = ifelse(is.na(d$effect), d$reason, ""),
                                                  stringsAsFactors = FALSE))
    }
    spec_b <- spec
    if (identical(spec$power_target, "interaction")) spec_b$power_biomarkers <- tg$biomarker
    for (h in seq_len(nrow(df_h))) for (method in methods) {
      eff <- df_h$effect[h]
      res <- list(n = NA_integer_, events = NA_real_, se = NA_real_, target_variance = NA_real_)
      reason <- df_h$reason[h]
      if (!is.finite(eff)) {
        if (!nzchar(reason)) reason <- "no decision-relevant effect"
      } else if (eff >= 0) {
        reason <- "the decision-relevant effect is not a benefit (rate ratio >= 1)"
      } else if (method == "engine") {
        res <- power_based_n_information(design, tg$target, eff, power, alpha, draw, component,
                                         reference_population, weights, n_max)
      } else {
        sch <- power_based_n(spec_b, randomised_per_year, draw, component, reference_population, weights,
                             rate_ratio = exp(eff), power = power, alpha = alpha, outcome = outcome, n_max = n_max)
        res <- list(n = sch$n, events = sch$events, se = NA_real_, target_variance = sch$target_variance)
      }
      if (is.finite(eff) && eff < 0 && is.na(res$n)) reason <- paste0("above n_max = ", n_max)
      rows[[length(rows) + 1]] <- data.frame(
        design = spec$id, biomarker = tg$biomarker, power_target = spec$power_target,
        hypothesis = df_h$hypothesis[h], price = df_h$price[h], effect = eff, rate_ratio = exp(eff),
        method = method, n = as.numeric(res$n), events = res$events, se = res$se,
        target_variance = res$target_variance, reason = reason, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

# -------------------------------------------------------------------------------
# EVSI design objects for the engine
# -------------------------------------------------------------------------------

#' EVSI design objects of the candidates for one accrual scenario
#'
#' @param specs candidate_designs().
#' @param registered_per_year Registered patients per year of the scenario.
#' @param n_grid Sample sizes of the randomised designs.
#' @param draw Draw defining the expected information
#'   (component_mean_draw(component)).
#' @param component Joint sampling component.
#' @param reference_population,weights Reference patients and weights.
#' @param max_reporting_year Designs reporting later carry no population value;
#'   their n are not evaluated (EVSI NA, method "reports_after_horizon").
#' @param tau_m2 Heterogeneity SD of S4's late follow-up.
#' @param include_power Add the power-based n to each randomised grid.
#' @param decision_effects design_decision_effects() rows for the
#'   decision-relevant hypothesis, or NULL (conventional hypothesis only).
#' @param power_method Method of the primary power-based n added to the grid
#'   ("engine", the default, or "schoenfeld"); both are reported.
#' @param tau_soc Heterogeneity SD of the control regimen of a stratified
#'   design's control cells.
#' @param tau_m3 Heterogeneity SD of the experimental regimen of a stratified
#'   design's experimental cells.
#' @return List with `designs` (named evsi_design() objects: every candidate
#'   and "S4+<id>" combined), `meta` (data frame per design and n: role,
#'   n_role, reporting year, enrolment fraction, randomised and registered
#'   rates and counts; a stratified design's rate is its own, its randomised
#'   rate NA), `power` (primary_power_size() per randomised design) and
#'   `power_sizes` (design_power_sizes() rows of every randomised design).
build_evsi_designs <- function(specs, registered_per_year, n_grid, draw, component, reference_population,
                               weights = NULL, max_reporting_year = Inf, tau_m2 = TAU_M2_DEFAULT,
                               include_power = TRUE, decision_effects = NULL, power_method = "engine",
                               tau_soc = TAU_SOC_DEFAULT, tau_m3 = TAU_M3_DEFAULT) {
  draw <- split_joint_draw(draw, component)
  power_method <- match.arg(power_method, POWER_METHODS)
  designs <- list()
  meta <- list()
  power <- list()
  power_sizes <- list()
  information <- list()
  s4 <- specs[vapply(specs, function(s) s$type == "single_arm", logical(1))]
  info_s4 <- NULL
  if (length(s4)) {
    spec4 <- s4[[1]]
    info_s4 <- single_arm_information(spec4$n_fixed, draw, component, reference_population, weights, tau_m2)
    targets_s4 <- s4_targets(component, reference_population, weights)
    cov_s4 <- information_target_covariance(lapply(info_s4$endpoints, function(e) e$information), targets_s4,
                                            component, label = spec4$id)
    designs[[spec4$id]] <- evsi_design(spec4$id, targets_s4, stats::setNames(list(cov_s4), spec4$n_fixed))
    meta[[spec4$id]] <- data.frame(design = spec4$id, base_design = spec4$id, role = "standalone",
                                   n = spec4$n_fixed, n_role = "fixed", reporting_year = s4_reporting_year(),
                                   enrolment_fraction = 1, randomised_per_year = NA_real_,
                                   registered_per_year = NA_real_, n_registered = spec4$n_fixed,
                                   evaluated = s4_reporting_year() <= max_reporting_year, stringsAsFactors = FALSE)
  }
  for (spec in specs[vapply(specs, function(s) s$type %in% c("rct", "stratified"), logical(1))]) {
    v_n <- sort(unique(c(n_grid, spec$n_proposal[is.finite(spec$n_proposal)])))
    if (spec$type == "stratified") {
      # Single centre at its own planned rate; enrolled = registered.
      v_report <- stratified_reporting_year(v_n, spec)
      v_eval <- v_report <= max_reporting_year
      meta[[spec$id]] <- data.frame(design = spec$id, base_design = spec$id, role = "standalone", n = v_n,
                                    n_role = ifelse(v_n %in% spec$n_proposal, "proposal", "grid"),
                                    reporting_year = v_report, enrolment_fraction = 1, randomised_per_year = NA_real_,
                                    registered_per_year = spec$accrual_per_year, n_registered = v_n,
                                    evaluated = v_eval, stringsAsFactors = FALSE)
      if (!any(v_eval)) next
      infos <- lapply(v_n[v_eval], function(n) {
        stratified_design_information(spec, n, draw, component, reference_population, weights, tau_soc, tau_m3)
      })
    } else {
      enrol <- design_enrolment_fraction(spec, reference_population, weights, component, draw)
      rate <- registered_per_year * enrol$fraction
      pw <- NULL
      if (include_power && !is.na(spec$power_target)) {
        df_ps <- design_power_sizes(spec, rate, draw, component, reference_population, weights,
                                    decision_effects = decision_effects)
        df_ps$randomised_per_year <- rate
        power_sizes[[spec$id]] <- df_ps
        pw <- primary_power_size(df_ps, power_method)
      }
      power[[spec$id]] <- pw
      if (!is.null(pw) && is.finite(pw$n)) v_n <- sort(unique(c(v_n, pw$n)))
      v_report <- design_reporting_year(v_n, rate, randomisation_week = spec$randomisation_week)
      v_eval <- v_report <= max_reporting_year
      v_role <- ifelse(!is.null(pw) & v_n %in% pw$n, "power",
                       ifelse(v_n %in% spec$n_proposal, "proposal", "grid"))
      meta[[spec$id]] <- data.frame(design = spec$id, base_design = spec$id, role = "standalone", n = v_n,
                                    n_role = v_role, reporting_year = v_report,
                                    enrolment_fraction = enrol$fraction, randomised_per_year = rate,
                                    registered_per_year = registered_per_year, n_registered = v_n / enrol$fraction,
                                    evaluated = v_eval, stringsAsFactors = FALSE)
      if (!any(v_eval)) next
      trial <- as_trial_design(spec, rate)
      infos <- lapply(v_n[v_eval], function(n) {
        design_information(trial, n, draw, component, reference_population, weights,
                           step = design_node_step(trial$min_follow_up_weeks + design_accrual_weeks(trial, n)))
      })
    }
    names(infos) <- as.character(v_n[v_eval])
    information[[spec$id]] <- infos
    designs[[spec$id]] <- evsi_design(spec$id, infos[[1]]$targets, lapply(infos, function(x) x$covariance))
    if (!is.null(info_s4)) {
      id_c <- paste0(s4[[1]]$id, "+", spec$id)
      l_comb <- lapply(infos, function(info) combined_design_information(info_s4, info, targets_s4, component, id_c))
      designs[[id_c]] <- evsi_design(id_c, l_comb[[1]]$targets, lapply(l_comb, function(x) x$covariance))
      df_c <- meta[[spec$id]]
      df_c$design <- id_c
      df_c$role <- "combined"
      df_c$reporting_year <- pmax(df_c$reporting_year, s4_reporting_year())
      df_c$evaluated <- df_c$evaluated & df_c$reporting_year <= max_reporting_year
      meta[[id_c]] <- df_c
    }
  }
  df_meta <- do.call(rbind, meta)
  rownames(df_meta) <- NULL
  df_power_sizes <- if (length(power_sizes)) do.call(rbind, power_sizes) else NULL
  if (!is.null(df_power_sizes)) rownames(df_power_sizes) <- NULL
  list(designs = designs, meta = df_meta, power = power, power_sizes = df_power_sizes, information = information,
       info_s4 = info_s4)
}

#' Primary power-based n of a design from its power-size table
#'
#' The conventional hypothesis; with several biomarkers (S1) the largest n,
#' since the trial must be powered for every interaction.
#'
#' @param df_sizes design_power_sizes() rows of one design.
#' @param method "engine" or "schoenfeld".
#' @return List with `n` (NA when any biomarker's n is NA), `events` (at the
#'   largest n), `method`, `n_engine` and `n_schoenfeld`.
primary_power_size <- function(df_sizes, method = "engine") {
  pick_method <- function(m) {
    d <- df_sizes[df_sizes$hypothesis == "conventional" & df_sizes$method == m, , drop = FALSE]
    if (!nrow(d) || anyNA(d$n)) return(list(n = NA_real_, events = NA_real_))
    k <- which.max(d$n)
    list(n = d$n[k], events = d$events[k])
  }
  primary <- pick_method(method)
  list(n = primary$n, events = primary$events, method = method, n_engine = pick_method("engine")$n,
       n_schoenfeld = pick_method("schoenfeld")$n)
}

#' Design table (voi_designs.csv)
#'
#' @param specs candidate_designs().
#' @param df_accrual accrual_scenarios().
#' @param params Parameter list with the test prices.
#' @param reference_population,weights Reference patients and weights.
#' @param component Joint sampling component.
#' @param draw Draw defining P(PFS > L) and the expected events.
#' @param power List over accrual scenarios of build_evsi_designs()$power, or
#'   NULL.
#' @param tau_m2 Heterogeneity SD of S4.
#' @param tau_soc Heterogeneity SD of a stratified design's control regimen.
#' @param tau_m3 Heterogeneity SD of a stratified design's experimental regimen.
#' @return Data frame with one row per design (the candidates and the S6
#'   power-based versions) and accrual scenario.
design_table <- function(specs, df_accrual, params, reference_population, weights = NULL, component, draw,
                         power = NULL, tau_m2 = TAU_M2_DEFAULT, tau_soc = TAU_SOC_DEFAULT,
                         tau_m3 = TAU_M3_DEFAULT) {
  draw <- split_joint_draw(draw, component)
  v_target_names <- function(spec) {
    if (spec$type == "single_arm") return(paste(rownames(s4_targets(component, reference_population, weights)),
                                                collapse = "; "))
    if (spec$type == "stratified") {
      return(paste(rownames(stratified_targets(spec, component, reference_population, weights)), collapse = "; "))
    }
    m_t <- joint_target_matrix(list(
      design_treatment_targets(rownames(component$original_os$res), "os", spec$enrolled_cells),
      design_treatment_targets(rownames(component$original_pfs$res), "pfs", spec$enrolled_cells)),
      c(paste0("os::", rownames(component$original_os$res)), paste0("pfs::", rownames(component$original_pfs$res))))
    paste(rownames(m_t), collapse = "; ")
  }
  rows <- list()
  for (k in seq_len(nrow(df_accrual))) {
    sc <- df_accrual[k, ]
    for (spec in specs) {
      alloc <- design_allocation(spec, reference_population, weights, component, draw)
      enrol <- alloc$enrolment
      single <- spec$type == "single_arm"
      strat <- spec$type == "stratified"
      rate <- if (single || strat) NA_real_ else sc$registered_per_year * enrol$fraction
      pw <- if (!is.null(power)) power[[sc$accrual_scenario]][[spec$id]] else NULL
      base <- data.frame(
        design = spec$id, label = spec$label, accrual_scenario = sc$accrual_scenario, type = spec$type,
        enrolled_cells = paste(spec$enrolled_cells, collapse = ";"),
        experimental_cells = paste(spec$experimental_cells, collapse = ";"),
        randomisation_week = spec$randomisation_week, p_experimental = alloc$p_experimental,
        cell_mix_registration = paste(sprintf("%s=%.4f", names(enrol$cell_mix_registration),
                                              enrol$cell_mix_registration), collapse = ";"),
        cell_mix_randomisation = paste(sprintf("%s=%.4f", names(enrol$cell_mix_randomisation),
                                               enrol$cell_mix_randomisation), collapse = ";"),
        enrolment_fraction = if (single || strat) 1 else enrol$fraction,
        registered_per_year = if (single) NA_real_ else if (strat) spec$accrual_per_year else sc$registered_per_year,
        randomised_per_year = rate,
        min_follow_up_years = if (single) NA_real_ else if (strat) spec$min_follow_up_years else
          DESIGN_MIN_FOLLOW_UP_YEARS,
        identified_targets = v_target_names(spec),
        screening_tests = paste(spec$screening_biomarkers, collapse = ";"),
        screening_cost_per_registered = design_screening_cost(spec, params),
        n_fixed = spec$n_fixed, n_proposal = spec$n_proposal,
        reporting_year_proposal = if (strat) stratified_reporting_year(spec$n_proposal, spec) else
          if (!single && is.finite(spec$n_proposal))
            design_reporting_year(spec$n_proposal, rate, randomisation_week = spec$randomisation_week) else
              if (single) s4_reporting_year() else NA_real_,
        power_target = spec$power_target, n_power = NA_real_, events_power = NA_real_,
        power_method = NA_character_, n_power_schoenfeld = NA_real_,
        reporting_year_power = NA_real_, tau_m2 = if (single) tau_m2 else NA_real_,
        tau_soc = if (strat) tau_soc else NA_real_, tau_m3 = if (strat) tau_m3 else NA_real_,
        sample_size_basis = if (is.null(spec$sample_size_basis)) "" else spec$sample_size_basis,
        notes = spec$notes, stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- base
      if (!is.null(pw)) {
        s6 <- base
        s6$design <- paste0("S6_", spec$id)
        s6$label <- paste0("Power-based ", spec$id, " (rate ratio ", POWER_RATE_RATIO, ", power ", POWER_LEVEL,
                           ", two-sided ", POWER_ALPHA, ", ", toupper(POWER_ENDPOINT), ")")
        s6$n_power <- pw$n
        s6$events_power <- pw$events
        s6$power_method <- if (is.null(pw$method)) "schoenfeld" else pw$method
        s6$n_power_schoenfeld <- if (is.null(pw$n_schoenfeld)) NA_real_ else pw$n_schoenfeld
        s6$reporting_year_power <- if (is.finite(pw$n))
          design_reporting_year(pw$n, rate, randomisation_week = spec$randomisation_week) else NA_real_
        rows[[length(rows) + 1]] <- s6
      }
    }
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}
