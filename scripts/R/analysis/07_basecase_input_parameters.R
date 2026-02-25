# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# Create treatment schedules
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

# Map survival curves directly into a list structure
# All three biomarker strategies (CRP, TLR, TMB/BRAF) are included in all models
p_os <- list(
  control_OS = predictions$control$os,
  crp_pos_OS = predictions$crp$biomarker_positive$os,
  crp_neg_OS = predictions$crp$biomarker_negative$os,
  tlr_pos_OS = predictions$tlr$biomarker_positive$os,
  tlr_neg_OS = predictions$tlr$biomarker_negative$os,
  tmb_braf_pos_OS = predictions$tmb_braf$biomarker_positive$os,
  tmb_braf_neg_OS = predictions$tmb_braf$biomarker_negative$os,
  crp_weighted_OS = predictions$crp$os,      # Population-marginalized
  tlr_weighted_OS = predictions$tlr$os,      # Population-marginalized
  tmb_braf_weighted_OS = predictions$tmb_braf$os  # Population-marginalized
)

p_pfs <- list(
  control_PFS = predictions$control$pfs,
  crp_pos_PFS = predictions$crp$biomarker_positive$pfs,
  crp_neg_PFS = predictions$crp$biomarker_negative$pfs,
  tlr_pos_PFS = predictions$tlr$biomarker_positive$pfs,
  tlr_neg_PFS = predictions$tlr$biomarker_negative$pfs,
  tmb_braf_pos_PFS = predictions$tmb_braf$biomarker_positive$pfs,
  tmb_braf_neg_PFS = predictions$tmb_braf$biomarker_negative$pfs,
  crp_weighted_PFS = predictions$crp$pfs,    # Population-marginalized
  tlr_weighted_PFS = predictions$tlr$pfs,    # Population-marginalized
  tmb_braf_weighted_PFS = predictions$tmb_braf$pfs  # Population-marginalized
)


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
  c_drug_nivo = 13923,   # cost of nivolumab per administration
  c_drug_FLOX = 427,    # cost of FLOX per administration
  
  # Test costs
  c_test_CT = 386,      # cost of CT scan
  c_test_blood = 16,    # assuming CRP, CBC, and chem-7; 8.77 NOKs per parameter except basic chemistry panel which is 4.40 NOKs per parameter; 
  # assuming that this covers 40% of the actual lab costs; 193 Norwegian Krone equals 16,41 Euro
  c_test_NGS = 2518,    # cost of next-generation sequencing (for TMB/BRAF), now updated to reflect Pia's paper
  
  # Other costs
  c_other_visit = 33,     # cost of standard outpatient visit
  c_other_baseline = 530,  # cost of comprehensive baseline visit
  c_other_follow = 33,    # cost of follow-up visits (quarterly)
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
  p_pfs = p_pfs,
  
  # Biomarker prevalence - using values from strategies_df
  # All three biomarkers are included in all models
  p_crp = strategies_df$prevalence[strategies_df$id == "crp"],
  p_tlr = strategies_df$prevalence[strategies_df$id == "tlr"],
  p_tmb_braf = strategies_df$prevalence[strategies_df$id == "tmb_braf"]
)