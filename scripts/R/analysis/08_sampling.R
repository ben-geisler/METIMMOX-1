# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

source("scripts/R/functions/bootstrap_survival_model.R")

# Initialize list to store bootstrapped models
boot_models <- list()

# Bootstrap control models
boot_models$control <- list(
  os = bootstrap_survival_model(Surv(OSwk, Death) ~ 1, data = data_control),
  pfs = bootstrap_survival_model(Surv(PFSwk, Progression) ~ 1, data = data_control)
)

# Bootstrap biomarker models
for (biomarker in biomarkers) {
  # OS model with biomarker, treatment, and interaction
  os_formula <- as.formula(paste0("Surv(OSwk, Death) ~ ", biomarker, " + Rx + ", biomarker, ":Rx"))
  boot_models[[biomarker]]$os <- bootstrap_survival_model(os_formula, data = data)
  
  # PFS model with biomarker, treatment, and interaction
  pfs_formula <- as.formula(paste0("Surv(PFSwk, Progression) ~ ", biomarker, " + Rx + ", biomarker, ":Rx"))
  boot_models[[biomarker]]$pfs <- bootstrap_survival_model(pfs_formula, data = data)
}

# Function to generate survival predictions from normboot samples
generate_bootstrap_predictions <- function(boot_model, newdata = NULL, sample_idx = 1, time_points) {
  # Make sure time_points passed in has correct length
  if (length(time_points) != time_points_length) {
    warning("time_points length mismatch in bootstrap predictions")
    time_points <- time_points[1:min(length(time_points), time_points_length)]
  }
  # Get the sampling coefficients for this sample
  # Format is different from boot() - we need to extract from transposed matrix
  boot_coeffs <- boot_model$bootstrap[, sample_idx]
  
  # Create a temporary model with sample coefficients
  temp_model <- boot_model$original
  temp_model$coefficients <- boot_coeffs
  
  # Generate predictions
  if (is.null(newdata)) {
    pred <- predict(temp_model, type = "survival", times = time_points)
  } else {
    pred <- predict(temp_model, newdata = newdata, type = "survival", times = time_points)
  }
  
  # Unnest and return survival probabilities
  pred_unnested <- unnest(pred, .pred)
  return(pred_unnested$.pred_survival)
}

# Create parameter distributions for other model parameters
# Cost parameters (assume 20% coefficient of variation for costs)
cv_costs <- 0.2

# Function to calculate gamma parameters from mean and CV
get_gamma_params <- function(mean_val, cv) {
  shape <- 1/cv^2
  rate <- shape/mean_val
  return(list(shape = shape, rate = rate))
}

# Drug costs - gamma distributions
dist_c_drug_nivo <- c(list(dist = "gamma"), 
                      get_gamma_params(l_params_base$c_drug_nivo, cv_costs))

dist_c_drug_FLOX <- c(list(dist = "gamma"), 
                      get_gamma_params(l_params_base$c_drug_FLOX, cv_costs))

# Test costs - gamma distributions
dist_c_test_CT <- c(list(dist = "gamma"), 
                    get_gamma_params(l_params_base$c_test_CT, cv_costs))

dist_c_test_blood <- c(list(dist = "gamma"), 
                       get_gamma_params(l_params_base$c_test_blood, cv_costs))

dist_c_test_NGS <- c(list(dist = "gamma"), 
                     get_gamma_params(l_params_base$c_test_NGS, cv_costs))

# Other costs - gamma distributions
dist_c_other_visit <- c(list(dist = "gamma"), 
                        get_gamma_params(l_params_base$c_other_visit, cv_costs))

dist_c_other_baseline <- c(list(dist = "gamma"), 
                           get_gamma_params(l_params_base$c_other_baseline, cv_costs))

dist_c_other_follow <- c(list(dist = "gamma"), 
                         get_gamma_params(l_params_base$c_other_follow, cv_costs))

dist_c_other_last <- c(list(dist = "gamma"), 
                       get_gamma_params(l_params_base$c_other_last, cv_costs))


# Utility parameters - beta distributions (bounded between 0 and 1)
# Using method of moments to get alpha and beta parameters
get_beta_params <- function(mean_val, cv) {
  var_val <- (mean_val * cv)^2
  alpha <- mean_val * (mean_val * (1 - mean_val) / var_val - 1)
  beta <- (1 - mean_val) * (mean_val * (1 - mean_val) / var_val - 1)
  return(list(shape1 = alpha, shape2 = beta))
}

cv_utilities <- 0.15  # 15% CV for utilities
dist_u_np <- c(list(dist = "beta"), get_beta_params(l_params_base$u_np, cv_utilities))
dist_u_p <- c(list(dist = "beta"), get_beta_params(l_params_base$u_p, cv_utilities))

# Create comprehensive parameter distributions list
param_distributions <- list(
  # Cost parameters
  c_drug_nivo = dist_c_drug_nivo,
  c_drug_FLOX = dist_c_drug_FLOX,
  c_test_CT = dist_c_test_CT,
  c_test_blood = dist_c_test_blood,
  c_test_NGS = dist_c_test_NGS,
  c_other_visit = dist_c_other_visit,
  c_other_baseline = dist_c_other_baseline,
  c_other_follow = dist_c_other_follow,
  c_other_last = dist_c_other_last,
  
  # Utility parameters
  u_np = dist_u_np,
  u_p = dist_u_p
)