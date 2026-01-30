# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/evppi_functions.R"))

# Ensure consistent time indexing
if (!exists("time_points_length")) {
  time_points_length <- length(time_points)
}

# Validate time_points consistency
if (length(time_points) != time_points_length) {
  stop("time_points length inconsistency detected")
}

# Check if PSA object exists
if (!exists("psa_obj")) {
  stop("PSA object not found. Please run 12_PSA.R first to create the PSA object.")
}

# ===============================================================================
# DIAGNOSTIC CHECKS
# ===============================================================================

cat("=== PSA OBJECT DIAGNOSTICS ===\n")
cat("PSA object structure:\n")
cat("  - Number of simulations:", psa_obj$n_sim, "\n")
cat("  - Number of strategies:", psa_obj$n_strategies, "\n")
cat("  - Strategy names:", paste(psa_obj$strategies, collapse = ", "), "\n")

# Check cost and effect matrices
cost_matrix <- as.matrix(psa_obj$cost)
effect_matrix <- as.matrix(psa_obj$effect)

cat("\nCost matrix summary:\n")
print(summary(cost_matrix))
cat("\nEffect matrix summary:\n")
print(summary(effect_matrix))

# Calculate NMB and check for variation
cat("\n=== NMB ANALYSIS ===\n")
nmb_matrix <- effect_matrix * WTP - cost_matrix
cat("NMB matrix summary:\n")
print(summary(nmb_matrix))

# Calculate EVPI manually for verification
max_nmb_per_sim <- apply(nmb_matrix, 1, max, na.rm = TRUE)
expected_max_nmb <- mean(max_nmb_per_sim, na.rm = TRUE)
max_expected_nmb <- max(colMeans(nmb_matrix, na.rm = TRUE))
evpi_manual <- expected_max_nmb - max_expected_nmb

cat("EVPI calculation:\n")
cat("  - Expected value of max NMB per simulation:", round(expected_max_nmb, 2), "\n")
cat("  - Max expected NMB across strategies:", round(max_expected_nmb, 2), "\n")
cat("  - EVPI:", round(evpi_manual, 2), "\n")

# ===============================================================================
# RUN EVPPI ANALYSIS
# ===============================================================================

# Define parameters for EVPPI analysis
# Model B (focused) excludes TLR, so p_tlr is only included for Models A and C
evppi_params_base <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS",
                       "c_test_CT", "u_np", "u_p", "c_other_last",
                       "c_test_blood", "c_other_visit", "c_other_baseline", "c_other_follow",
                       "p_crp", "p_tmb_braf")
if (MODEL_STRUCTURE != 1) {
  evppi_params <- c(evppi_params_base, "p_tlr")
} else {
  evppi_params <- evppi_params_base
}

# Run EVPPI analysis (param_groups defined in 08_sampling.R)
evppi_results <- run_evppi_analysis(
  psa_obj = psa_obj,
  psa_params = psa_params,
  wtp = WTP,
  evppi_params = evppi_params,
  param_groups = param_groups
)

# ===============================================================================
# SUMMARY AND OUTPUT
# ===============================================================================

cat("\n=== EVPPI ANALYSIS SUMMARY ===\n")
cat("Willingness-to-pay threshold:", WTP, "\n")
cat("Number of PSA simulations:", psa_obj$n_sim, "\n")
cat("Total EVPI:", round(evpi_manual, 4), "\n\n")

if (nrow(evppi_results) > 0) {
  cat("Individual Parameter EVPPIs (ranked by importance):\n")
  print(evppi_results[, c("parameter", "evppi", "evppi_percent_of_evpi", "error")])
  
  # Identify high-priority parameters
  high_priority_threshold <- 0.01  # 1% of EVPI
  high_priority_params <- evppi_results[evppi_results$evppi_percent_of_evpi > high_priority_threshold * 100, ]
  
  if (nrow(high_priority_params) > 0) {
    cat("\nHigh-priority parameters for future research (>", high_priority_threshold * 100, "% of EVPI):\n")
    for (i in 1:nrow(high_priority_params)) {
      cat("  ", high_priority_params$parameter[i], ": €", 
          round(high_priority_params$evppi[i], 4), "\n")
    }
  } else {
    cat("\nNo parameters exceed the high-priority threshold\n")
  }
  
  # Create visualization if we have meaningful results
  meaningful_results <- evppi_results[evppi_results$evppi > 0, ]
  if (nrow(meaningful_results) > 0) {
    par(mar = c(5, 8, 4, 2))
    barplot(meaningful_results$evppi, 
            names.arg = meaningful_results$parameter,
            horiz = TRUE,
            las = 1,
            main = "Expected Value of Partially Perfect Information (EVPPI)",
            xlab = paste("EVPPI (€) at WTP =", WTP),
            col = "steelblue")
  }
  
} else {
  cat("No EVPPI results available\n")
}

# Save results
save(evppi_results, evpi_manual, file = here::here("data/tidy/evppi_results.RData"))
cat("\nEVPPI analysis complete. Results saved to data/tidy/evppi_results.RData\n")