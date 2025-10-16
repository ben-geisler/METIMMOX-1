# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/psa_functions.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}

# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Generate PSA samples
cat("Generating PSA samples for", n_sim, "simulations\n")
psa_params <- generate_psa_samples(param_distributions, n_sim)

# Run the PSA
cat("Starting PSA with", n_sim, "simulations\n")
psa_results <- run_psa_analysis(
  psa_params = psa_params,
  l_params_base = l_params_base,
  param_distributions = param_distributions,
  strategies = strategies,
  time_horizon = time_horizon,
  cl = cl,
  n_sim = n_sim
)

# Create the PSA object using dampack's make_psa_obj function
psa_obj <- dampack::make_psa_obj(
  cost = as.data.frame(psa_results$cost),
  effect = as.data.frame(psa_results$effect),
  strategies = strategies,
  currency = "€"
)

# Verify the PSA object structure
cat("PSA object created with:\n")
cat(" - Simulations:", psa_obj$n_sim, "\n")
cat(" - Strategies:", psa_obj$n_strategies, "\n")
cat(" - Cost matrix dimensions:", dim(psa_obj$cost), "\n")
cat(" - Effect matrix dimensions:", dim(psa_obj$effect), "\n")

psa_summary <- summary(psa_obj)
print(psa_summary)

psa_summary_c95ci <- summary(psa_obj, calc_sds = TRUE)
print(psa_summary_c95ci)

# Calculate ICERs using the correct arguments
cea_psa <- calculate_icers_psa(psa_obj, uncertainty = TRUE)
print(cea_psa)

# ICE scatter plot
plot(psa_obj, frontier = TRUE, points = TRUE)

# Create cost-effectiveness acceptability curves
c_wtp <- seq(from = 0, to = 2e5, by = 5e3)
ceac_obj <- ceac(c_wtp, psa_obj)
ceac_sum <- summary(ceac_obj)
print(ceac_sum)
plot(ceac_obj, frontier = TRUE, points = TRUE, currency = "€")