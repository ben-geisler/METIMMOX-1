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
l_nivo <- rep(0, time_points_length)
l_nivo[c(5, 7, 13, 15, 29, 31, 37, 39)] <- 1  

l_FLOX_exp <- rep(0, time_points_length)
l_FLOX_exp[c(1, 3, 9, 11, 25, 27, 33, 35)] <- 1

l_FLOX_control <- rep(0, time_points_length)
l_FLOX_control[c(1, 3, 5, 7, 9, 11, 13, 15, 25, 27, 29, 31, 33, 35, 37, 39)] <- 1

l_CT <- rep(0, time_points_length)
l_CT[1] <- 1  # baseline
l_CT[seq(13, time_points_length, by=12)] <- 1  # every 12 weeks

l_blood <- rep(0, time_points_length)
l_blood[1] <- 1  # baseline
l_blood[seq(5, time_points_length, by=4)] <- 1  # every 4 weeks

l_visit <- rep(0, time_points_length)
l_visit[1] <- 1  # baseline
l_visit[which(l_nivo == 1 | l_FLOX_exp == 1 | l_FLOX_control == 1)] <- 1

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

# Compile all parameters into a list for the model function
l_params_base <- list(
  # Time parameters
  cl = cl,
  time_horizon = time_horizon,
  
  # Discount rates
  dr_costs = dr,
  dr_effects = dr,
  
  # Utilities - determined by UTILITY_SOURCE switch (set in 02_setup_and_global_variables.R)
  # IPD-derived: u_np = 0.9077, u_p = 0.9005 (from METIMMOX trial)
  # CORRECT trial: u_np = 0.73, u_p = 0.59 (Gourzoulidis et al. 2018)
  u_np = if (UTILITY_SOURCE == 0) 0.9077 else 0.73,
  u_p  = if (UTILITY_SOURCE == 0) 0.9005 else 0.59,
  
  # Drug costs
  # c_drug_nivo: ASSUMPTION, not a published tariff (issue #154). EUR 13,923 per
  # administration is the assumed Norwegian hospital acquisition cost of one
  # nivolumab administration in this regimen; eight administrations give about
  # EUR 111,000 per treated patient. Norwegian hospital prices are set by
  # confidential LIS tender and cannot be cited; published list prices are
  # higher (a US list price is about USD 8,100 per 240 mg vial). The value is
  # the largest single cost driver in the model, so it is varied in the one-way
  # DSA (+/-20%) and in the biosimilar scenario (EUR 4,641 per administration,
  # one third of the base case) rather than sampled in the PSA.
  c_drug_nivo = 13923,   # assumed cost of nivolumab per administration
  c_drug_FLOX = 427,    # cost of FLOX per administration
  
  # Test costs
  c_test_CT = 386,      # cost of CT scan
  c_test_blood = 16,    # routine CBC and chemistry monitoring; 8.77 NOKs per parameter except basic chemistry panel which is 4.40 NOKs per parameter;
  # assuming that this covers 40% of the actual lab costs; 193 Norwegian Krone equals 16,41 Euro
  c_test_CRP = c_test_CRP, # one-time CRP biomarker test (independent of routine blood monitoring)
  c_test_NGS = c_test_NGS, # cost of next-generation sequencing (for TMB/BRAF), now updated to reflect Pia's paper
  
  # Other costs
  c_other_visit = 33,     # cost of standard outpatient visit
  c_other_baseline = 530,  # cost of comprehensive baseline visit
  c_other_follow = 33,    # cost of follow-up visits (quarterly)
  # Post-progression treatment cost per quarter in the progressed state
  # (issue #154). Zero in the base case: second-line systemic therapy, imaging
  # and visits after progression are outside the modeled scope, so the
  # progressed state accrues only the quarterly follow-up contact and the
  # end-of-life cost. The strategies differ in time spent progressed, so the
  # omission is differential; the parameter makes it explicit and is varied in
  # the Post_progression_cost structural scenario in 09_DSA.R.
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
