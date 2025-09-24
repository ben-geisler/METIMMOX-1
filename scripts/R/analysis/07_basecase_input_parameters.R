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
p_os <- list(
  control_OS = predictions$control$os,  
  crp_pos_OS = predictions$crp$os_pos_exp,
  crp_neg_OS = predictions$crp$os_neg_ctrl,
  tlr_pos_OS = predictions$tlr$os_pos_exp,
  tlr_neg_OS = predictions$tlr$os_neg_ctrl,
  tmb_braf_pos_OS = predictions$tmb_braf$os_pos_exp,
  tmb_braf_neg_OS = predictions$tmb_braf$os_neg_ctrl,
  crp_weighted_OS = predictions$crp$os_weighted,
  tlr_weighted_OS = predictions$tlr$os_weighted,
  tmb_braf_weighted_OS = predictions$tmb_braf$os_weighted
)

p_pfs <- list(
  control_PFS = predictions$control$pfs, 
  crp_pos_PFS = predictions$crp$pfs_pos_exp,
  crp_neg_PFS = predictions$crp$pfs_neg_ctrl,
  tlr_pos_PFS = predictions$tlr$pfs_pos_exp,
  tlr_neg_PFS = predictions$tlr$pfs_neg_ctrl,
  tmb_braf_pos_PFS = predictions$tmb_braf$pfs_pos_exp,
  tmb_braf_neg_PFS = predictions$tmb_braf$pfs_neg_ctrl,
  crp_weighted_PFS = predictions$crp$pfs_weighted,
  tlr_weighted_PFS = predictions$tlr$pfs_weighted,
  tmb_braf_weighted_PFS = predictions$tmb_braf$pfs_weighted
)


# Compile all parameters into a list for the model function
l_params_base <- list(
  # Time parameters
  cl = cl,
  time_horizon = time_horizon,
  
  # Discount rates
  dr_costs = dr,
  dr_effects = dr,
  
  # Utilities
  u_np =  0.9077, # utility in non-progressed state
  u_p =  0.9005,  # utility in progressed state
  
  # Drug costs
  c_drug_nivo = 13923,   # cost of nivolumab per administration
  c_drug_FLOX = 427,    # cost of FLOX per administration
  
  # Test costs
  c_test_CT = 386,      # cost of CT scan
  c_test_blood = 0,    # cost of blood tests
  c_test_NGS = 1439,    # cost of next-generation sequencing (for TMB/BRAF)
  
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
  
  # Biomarker prevalence - using values from data
  p_crp = mean(data$crp, na.rm = TRUE),
  p_tlr = mean(data$tlr, na.rm = TRUE),
  p_tmb_braf = mean(data$tmb_braf, na.rm = TRUE)
)