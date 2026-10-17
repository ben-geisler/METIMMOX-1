# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# Create treatment schedules.
#
# SECOND TREATMENT SEQUENCE (issue #154): positions 25-39 (modeled weeks 24-38)
# repeat the first eight-cycle sequence. Both arms receive it, and the model
# gives it to EVERY patient still progression-free at that point. In METIMMOX
# the second sequence was started on progression during the treatment break, so
# this is a simplifying assumption that charges second-sequence drug and visit
# costs to progression-free patients who would not have been re-treated. A
# partitioned survival model has no on-treatment substate, so re-treatment
# cannot be triggered on progression without restructuring the model; the
# assumption is instead bounded by the Second_sequence structural scenario in
# 09_DSA.R, which removes the second sequence and holds survival at the trial
# estimate (a cost-side bound).
#
# The positions and monitoring rules live in build_treatment_schedules()
# (parameter_distributions.R, issues #49, #172), shared with the time-horizon
# scenario: nivolumab at modeled weeks 4, 6, 12, 14, 28, 30, 36, 38; FLOX
# (experimental) at weeks 0, 2, 8, 10, 24, 26, 32, 34; FLOX (control) at the
# union; CT at baseline and every 8 weeks; blood tests at every administration
# week and with every CT; visits at baseline, every administration and every CT
# scan (issue #164: a surveillance visit with each post-treatment CT), and a
# second contact on day 2 of every FLOX administration, charged in
# calculate_outcomes() (issue #207: trial team's answers MON1-MON3 of 15
# October 2026). Positions beyond the horizon are dropped, so every vector has
# exactly time_points_length points.
schedules <- build_treatment_schedules(time_points_length)
l_nivo <- schedules$l_nivo
l_FLOX_exp <- schedules$l_FLOX_exp
l_FLOX_control <- schedules$l_FLOX_control
l_CT <- schedules$l_CT
l_blood <- schedules$l_blood
l_visit <- schedules$l_visit

# Map configured strategy predictions into the curve keys consumed by model_fun().
control_strategy <- get_control_strategy()
p_os <- setNames(list(predictions[[control_strategy]]$os),
                 paste0(control_strategy, "_OS"))
p_pfs <- setNames(list(predictions[[control_strategy]]$pfs),
                  paste0(control_strategy, "_PFS"))
for (biomarker in get_biomarkers()) {
  prediction <- predictions[[biomarker]]
  p_os[[paste0(biomarker, "_pos_OS")]] <- prediction$biomarker_positive$os
  p_os[[paste0(biomarker, "_neg_OS")]] <- prediction$biomarker_negative$os
  p_os[[paste0(biomarker, "_weighted_OS")]] <- prediction$os
  p_pfs[[paste0(biomarker, "_pos_PFS")]] <- prediction$biomarker_positive$pfs
  p_pfs[[paste0(biomarker, "_neg_PFS")]] <- prediction$biomarker_negative$pfs
  p_pfs[[paste0(biomarker, "_weighted_PFS")]] <- prediction$pfs
}

# Biomarker diagnostic costs are distinct from routine monitoring costs.
c_test_CRP <- 16
c_test_NGS <- 2518

# IPD utilities (UTILITY_SOURCE = 0, issue #188): derived in the pipeline from
# the METIMMOX EQ-5D-5L responses with the EQ5D_VALUE_SET value set, with
# patient-clustered bootstrap standard errors for the PSA. The literature base
# case (UTILITY_SOURCE = 1, CORRECT trial) reads no utility data. There is no
# fallback to constants: a missing EQ-5D file stops.
if (UTILITY_SOURCE == 0) {
  source(here::here("R", "eq5d5l_utility.R"))
  source(here::here("R", "trial_utilities.R"))
  ipd_utilities <- trial_eq5d_utility_uncertainty(
    read_trial_eq5d(), readRDS(here::here("data", "tidy", "METIMMOX.rds")),
    value_set = EQ5D_VALUE_SET, seed = analysis_seed)
}

# Compile all parameters into a list for the model function
l_params_base <- list(
  # Time parameters
  cl = cl,
  time_horizon = time_horizon,
  
  # Discount rates
  dr_costs = dr,
  dr_effects = dr,
  
  # Utilities - determined by UTILITY_SOURCE switch (set in 02_setup_and_global_variables.R)
  # IPD-derived (issue #188): trial_eq5d_utility_uncertainty() above
  # CORRECT trial: u_np = 0.73, u_p = 0.59 (Gourzoulidis et al. 2018)
  u_np = if (UTILITY_SOURCE == 0) ipd_utilities$u_np else 0.73,
  u_p  = if (UTILITY_SOURCE == 0) ipd_utilities$u_p else 0.59,
  
  # Drug costs
  # c_drug_nivo: ASSUMPTION, not a published tariff (issue #154). EUR 1,392 per
  # administration is the assumed Norwegian hospital acquisition cost of one
  # nivolumab administration in this regimen, already about 40% below list;
  # eight administrations give about EUR 11,100 per treated patient. (Before
  # #210 the value was EUR 13,923, ten times the intended price.) Norwegian
  # hospital prices are set by confidential LIS tender and cannot be cited. The
  # value is a major cost driver, so it is varied in the one-way DSA (+/-20%)
  # and in the biosimilar scenario (NIVO_PRICE_BIOSIMILAR in
  # scenario_analysis.R, one third of the base case) rather than sampled in the
  # PSA.
  c_drug_nivo = 1392,   # assumed cost of nivolumab per administration
  c_drug_FLOX = 427,    # cost of FLOX per administration
  
  # Test costs
  c_test_CT = 386,      # cost of CT scan
  c_test_blood = 16,    # routine CBC and chemistry monitoring; 8.77 NOKs per parameter except basic chemistry panel which is 4.40 NOKs per parameter;
  # assuming that this covers 40% of the actual lab costs; 193 Norwegian Krone equals 16,41 Euro
  c_test_CRP = c_test_CRP, # one-time CRP biomarker test (independent of routine blood monitoring)
  c_test_NGS = c_test_NGS, # cost of next-generation sequencing (for TMB/BRAF), now updated to reflect Pia's paper
  
  # Other costs
  # Outpatient visits (issue #36): Norwegian DRG (ISF) weight 0.047 for an
  # outpatient contact times the 2023 somatic-care unit price of NOK 49,484
  # (EUR 4,541.60), i.e. NOK 2,325.75 = EUR 213.46 per visit. The weight is the
  # same for every visit type, so the baseline, treatment/surveillance and
  # progressed-state follow-up visits share one value; the three parameters are kept
  # separate so that each can be varied on its own. (Before #36: EUR 33, 530
  # and 33 without a source.)
  c_other_visit = 213.46,    # outpatient visit (DRG weight 0.047)
  c_other_baseline = 213.46, # baseline visit (DRG weight 0.047)
  c_other_follow = 213.46,   # follow-up visit in the progressed state (one every 2.5 months, #207)
  # Flat post-progression treatment cost per quarter in the progressed state
  # (issue #154). Zero in the base case: since issue #209 second-line therapy
  # is costed once at incident progression from the Norwegian pathway (the
  # second-line fields appended below), and the progressed state otherwise
  # accrues the follow-up visit and the CT scan every 2.5 months (issues
  # #164, #207) and the end-of-life cost. The former EUR 5,000 per quarter is
  # the Post_progression_cost structural scenario in 09_DSA.R, which replaces
  # the pathway.
  c_other_pp = 0,
  c_other_last = 13803,     # cost of end-of-life care
  
  # Treatment schedules
  l_nivo = l_nivo,
  l_FLOX_exp = l_FLOX_exp,
  l_FLOX_control = l_FLOX_control,
  l_CT = l_CT,
  l_blood = l_blood,
  l_visit = l_visit,
  
  # Survival curves
  p_os = p_os,
  p_pfs = p_pfs
  
  # Biomarker prevalence parameters are appended below from configured IDs.
)

# PSA utility decrement (issue #154). u_p is not sampled independently of u_np;
# the PSA draws this non-negative decrement and sets u_p = u_np - u_decrement,
# so no draw can put the progressed utility above the progression-free utility.
# Derived here so the two utilities and the decrement can never disagree.
l_params_base$u_decrement <- l_params_base$u_np - l_params_base$u_p
if (l_params_base$u_decrement <= 0) {
  stop("The base-case progressed utility must be below the progression-free ",
       "utility for the PSA decrement parameterisation.")
}

# IPD mode only: bootstrap SEs replace the assumed CV 0.15 of u_np and
# u_decrement in create_parameter_distributions(), and the EQ-5D input (value
# set, file digest, counts) enters l_params_base and hence the PSA fingerprint.
if (UTILITY_SOURCE == 0) {
  l_params_base$utility_se <- ipd_utilities$se[c("u_np", "u_decrement")]
  l_params_base$utility_input <- list(
    source = "METIMMOX EQ-5D-5L IPD", value_set = ipd_utilities$value_set,
    eq5d_md5 = unname(tools::md5sum(here::here(TRIAL_EQ5D_PATH))),
    n_responses = ipd_utilities$n_responses, n_patients = ipd_utilities$n_patients,
    n_boot = ipd_utilities$n_boot, seed = ipd_utilities$seed)
}

for (biomarker in get_biomarkers()) {
  prevalence_key <- biomarker_prevalence_key(biomarker)
  l_params_base[[prevalence_key]] <-
    strategies_df$prevalence[strategies_df$id == biomarker]
}

# Shared prediction population and joint biomarker distribution (issue #166).
l_params_base <- set_population_predictions(l_params_base, predictions)
joint_population <- joint_biomarker_population(data_complete)
l_params_base$joint_counts <- joint_population$counts
for (cell in names(joint_population$probabilities))
  l_params_base[[paste0("p_joint_", cell)]] <- unname(joint_population$probabilities[cell])

# Third and fourth treatment sequences (issue #208). Trial team, 15 October
# 2026 (SEQ1-SEQ5): later sequences were also started on progression during a
# break, so they follow the trial's own data, not a schedule. Script 01 derives
# per patient from the export's administration block whether a third and a
# fourth sequence were started and when (weeks from inclusion). Here, in the
# economic prediction population, the shares starting each sequence by arm
# (with their counts, for the beta distributions of the PSA) and the observed
# start weeks of the starters enter l_params_base; calculate_outcomes() charges
# share x one full sequence of the arm's regimen at those weeks, whatever the
# health state (later_sequence_resource_use()). Sequence 2 is unchanged (weeks
# 24-38 for every progression-free patient); no fifth or sixth sequence is
# modelled.
source(here::here("R", "later_sequences.R"))
later_sequence_data <- later_sequence_inputs(
  data_complete, readRDS(here::here("data", "tidy", "METIMMOX.rds")))
l_params_base$later_sequences <- later_sequence_data$later_sequences
l_params_base[names(later_sequence_data$shares)] <- later_sequence_data$shares
rm(later_sequence_data)

# Second-line therapy at incident progression (issue #209; user decision of
# 17 October 2026: the costed Norwegian pathway is the base case). The
# pathway table data/external/second_line_pathway.csv (R/second_line_costs.R)
# gives the shares, cycles, contacts and blood tests, and the drug cost per
# cycle of each backbone and biologic at the table's doses and public list
# prices. The progression share of PFS events is counted in the economic
# cohort under the #181 PFS rule, pooled over both arms: the trial team
# expects no difference by arm (SL3), and by arm one count is zero, which
# would give a degenerate beta. calculate_outcomes() charges the expected
# course cost per progressing patient once per progression
# (second_line_event_costs()).
source(here::here("R", "second_line_costs.R"))
second_line_inputs <- second_line_model_inputs(read_second_line_pathway())
progression_share <- pfs_progression_share(data_complete)
for (field in names(second_line_inputs$values)) {
  l_params_base[[field]] <- unname(second_line_inputs$values[[field]])
}
l_params_base$p_pfs_progression <- progression_share$value
# PSA standard errors (create_parameter_distributions()) and sourced one-way
# bounds (build_dsa_ranges()) of the second-line parameters.
l_params_base$psa_se <- c(p_pfs_progression = progression_share$se, second_line_inputs$se)
v_ranged <- !is.na(second_line_inputs$lower)
l_params_base$dsa_bounds <- data.frame(
  pars = c("p_pfs_progression", names(second_line_inputs$lower)[v_ranged]),
  min = c(progression_share$lower, unname(second_line_inputs$lower[v_ranged])),
  max = c(progression_share$upper, unname(second_line_inputs$upper[v_ranged])),
  stringsAsFactors = FALSE)
rm(field, v_ranged)
