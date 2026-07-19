# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))
source(here::here("scripts/R/functions/evppi_functions.R"))

# ===============================================================================
# EXTRACT INTERACTION COEFFICIENTS FROM SAMPLING MODELS
# ===============================================================================
# Treatment-biomarker interaction coefficients are extracted from the cached
# sampling_models and appended to psa_params. This enables EVPPI analysis for
# the interaction parameters that drive biomarker-guided treatment decisions.
# No cache regeneration required - extraction uses existing sampling_models.
# ===============================================================================

cat("\n=== Extracting interaction coefficients for EVPPI ===\n")

if (!exists("sampling_models") || is.null(sampling_models)) {
  warning("sampling_models not available - interaction EVPPI will be skipped")
  interaction_params_available <- FALSE
} else {
  interaction_coefs_all <- extract_interaction_coefficients(
    sampling_models = sampling_models,
    n_sim = n_sim
  )

  retained_sim_ids <- if ("sim" %in% names(psa_params)) {
    psa_params$sim
  } else {
    seq_len(nrow(psa_params))
  }
  valid_sim_ids <-
    length(retained_sim_ids) == nrow(psa_params) &&
    all(is.finite(retained_sim_ids)) &&
    all(retained_sim_ids == as.integer(retained_sim_ids)) &&
    all(retained_sim_ids >= 1L & retained_sim_ids <= n_sim)

  if (!is.null(interaction_coefs_all) &&
      nrow(interaction_coefs_all) == n_sim && valid_sim_ids) {
    interaction_coefs <- interaction_coefs_all[retained_sim_ids, , drop = FALSE]
    psa_params <- cbind(psa_params, interaction_coefs)
    attr(psa_params, "seed") <- analysis_seed
    cat("Interaction coefficients appended to psa_params\n")
    cat("psa_params now has", ncol(psa_params), "columns\n")
    interaction_params_available <- TRUE
  } else {
    warning("Failed to extract interaction coefficients - skipping")
    interaction_params_available <- FALSE
  }
}

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

# Total EVPI via dampack. This equals the population-NMB definition used above
# (mean(max NMB per draw) - max(mean NMB per strategy)); verified equal to the
# hand-rolled value to < 1e-6. Named evpi_manual for cache/report compatibility.
evpi_manual <- dampack::calc_evpi(psa = psa_obj, wtp = WTP)$EVPI[1]

cat("EVPI calculation:\n")
cat("  - EVPI (dampack::calc_evpi):", round(evpi_manual, 2), "\n")

# ===============================================================================
# RUN EVPPI ANALYSIS
# ===============================================================================

# Define parameters for EVPPI analysis
evppi_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CRP", "c_test_NGS",
                  "c_test_CT", "u_np", "u_p", "c_other_last",
                  "c_test_blood", "c_other_visit", "c_other_baseline",
                  "c_other_follow", "p_crp", "p_tmb_braf")

# Add interaction parameters if available
if (interaction_params_available) {
  interaction_evppi_params <- get_interaction_evppi_params()
  interaction_evppi_params <- intersect(
    interaction_evppi_params, colnames(psa_params)
  )
  evppi_params <- c(evppi_params, interaction_evppi_params)
  cat("EVPPI parameters include", length(interaction_evppi_params),
      "interaction coefficients\n")
}

# Verify param_groups exists (defined in 08_sampling.R)
if (!exists("param_groups")) {
  stop("'param_groups' not found. Run 08_sampling.R first to define parameter groups.")
}

# Add interaction parameter groups (if available)
if (interaction_params_available) {
  param_groups <- add_interaction_param_groups(param_groups)
  cat("Added interaction parameter groups to param_groups\n")
  cat("Total param_groups:", length(param_groups), "\n")
}

# Run EVPPI analysis
evppi_results <- run_evppi_analysis(
  psa_obj = psa_obj,
  psa_params = psa_params,
  wtp = WTP,
  evppi_params = evppi_params,
  param_groups = param_groups,
  seed = analysis_seed
)

# ===============================================================================
# ADD POPULATION-LEVEL EVPPI
# ===============================================================================

# Add population-level EVPPI columns
evppi_results <- evppi_results %>%
  mutate(
    evppi_population_millions = calculate_population_evppi(
      evppi,
      annual_incidence = annual_incidence_norway,
      research_horizon = research_horizon_years,
      discount_rate = discount_rate_research
    ),
    evpi_population_millions = calculate_population_evppi(
      evpi,
      annual_incidence = annual_incidence_norway,
      research_horizon = research_horizon_years,
      discount_rate = discount_rate_research
    )
  )

# Print summary of population-level values
cat("\n=== Population-Level EVPPI Summary ===\n")
cat(sprintf("Annual incidence: %d patients/year\n", annual_incidence_norway))
cat(sprintf("Research horizon: %d years\n", research_horizon_years))
cat(sprintf("Discount rate: %.1f%%\n", discount_rate_research * 100))
cat(sprintf("Total eligible population (discounted): %.0f patient-years\n",
            annual_incidence_norway * sum(1 / (1 + discount_rate_research)^(1:research_horizon_years))))
cat("\nTop 5 parameters by population-level EVPPI:\n")
print(evppi_results %>%
        arrange(desc(evppi_population_millions)) %>%
        select(parameter, evppi_population_millions, evppi_percent_of_evpi) %>%
        head(5) %>%
        mutate(evppi_population_millions = sprintf("€%.2fM", evppi_population_millions)))

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
evppi_cache_file <- evppi_path()
save(evppi_results, evpi_manual, file = evppi_cache_file)
cat("\nEVPPI analysis complete. Results saved to", evppi_cache_file, "\n")
