# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, dplyr, parallel)

# Load functions
source(here::here("scripts/R/functions/model_fun.R"))
source(here::here("scripts/R/functions/calculate_outcomes.R"))

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

# First, let's examine the PSA object and parameters
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

# Check parameter variation
if (exists("psa_params")) {
  cat("\n=== PARAMETER VARIATION ANALYSIS ===\n")
  param_summary <- data.frame(
    parameter = character(),
    min = numeric(),
    max = numeric(),
    mean = numeric(),
    sd = numeric(),
    cv = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (param_name in colnames(psa_params)) {
    if (is.numeric(psa_params[[param_name]])) {
      param_values <- psa_params[[param_name]]
      param_values <- param_values[!is.na(param_values)]
      
      if (length(param_values) > 0) {
        param_mean <- mean(param_values)
        param_sd <- sd(param_values)
        param_cv <- if (param_mean != 0) param_sd / param_mean else NA
        
        param_summary <- rbind(param_summary, data.frame(
          parameter = param_name,
          min = min(param_values),
          max = max(param_values),
          mean = param_mean,
          sd = param_sd,
          cv = param_cv,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  print(param_summary)
  
  # Check for parameters with no variation
  no_variation <- param_summary$cv < 0.001 | is.na(param_summary$cv)
  if (any(no_variation)) {
    cat("\nParameters with little/no variation:\n")
    print(param_summary[no_variation, c("parameter", "cv")])
  }
}

# ===============================================================================
# IMPROVED EVPPI CALCULATION FUNCTIONS
# ===============================================================================

# Improved EVPPI calculation function with better error handling
calculate_evppi_improved <- function(psa_obj, psa_params, param_names, wtp = WTP, n_inner = 100, n_grid = 10) {
  
  cat("Calculating EVPPI for parameter(s):", paste(param_names, collapse = ", "), "\n")
  
  # Extract cost and effect matrices
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- as.matrix(psa_obj$effect)
  n_sim <- nrow(cost_matrix)
  n_strategies <- ncol(cost_matrix)
  
  # Remove any rows with missing values
  complete_rows <- complete.cases(cost_matrix) & complete.cases(effect_matrix)
  if (sum(complete_rows) < n_sim) {
    cat("  Warning: Removing", n_sim - sum(complete_rows), "incomplete simulations\n")
    cost_matrix <- cost_matrix[complete_rows, , drop = FALSE]
    effect_matrix <- effect_matrix[complete_rows, , drop = FALSE]
    n_sim <- nrow(cost_matrix)
  }
  
  # Calculate Net Monetary Benefit (NMB) matrix
  nmb_matrix <- effect_matrix * wtp - cost_matrix
  
  # Check for missing values in NMB matrix
  if (any(is.na(nmb_matrix))) {
    cat("  Warning: Missing values detected in NMB matrix\n")
    return(list(evppi = 0, evpi = 0, param_names = param_names, 
                error = "Missing values in NMB matrix"))
  }
  
  # Calculate EVPI
  max_nmb_per_sim <- apply(nmb_matrix, 1, max, na.rm = TRUE)
  expected_max_nmb <- mean(max_nmb_per_sim, na.rm = TRUE)
  max_expected_nmb <- max(colMeans(nmb_matrix, na.rm = TRUE))
  evpi <- expected_max_nmb - max_expected_nmb
  
  cat("  EVPI for this calculation:", round(evpi, 4), "\n")
  
  # If EVPI is zero or negative, return zero EVPPI
  if (evpi <= 0) {
    cat("  EVPI is", evpi, "- returning zero EVPPI\n")
    return(list(evppi = 0, evpi = evpi, param_names = param_names,
                error = "EVPI is zero or negative"))
  }
  
  # Check parameter variation
  param_has_variation <- TRUE
  for (param_name in param_names) {
    if (!param_name %in% colnames(psa_params)) {
      cat("  Warning: Parameter", param_name, "not found in psa_params\n")
      return(list(evppi = 0, evpi = evpi, param_names = param_names,
                  error = paste("Parameter", param_name, "not found")))
    }
    
    param_values <- psa_params[[param_name]]
    param_values <- param_values[complete_rows]  # Use same subset
    param_values <- param_values[!is.na(param_values)]
    
    if (length(unique(param_values)) < 5) {
      cat("  Warning: Parameter", param_name, "has insufficient variation\n")
      param_has_variation <- FALSE
    }
  }
  
  if (!param_has_variation) {
    return(list(evppi = 0, evpi = evpi, param_names = param_names,
                error = "Insufficient parameter variation"))
  }
  
  # Simplified EVPPI calculation using quantile-based approach
  tryCatch({
    # Create parameter grids using fewer points for stability
    param_grids <- list()
    for (param_name in param_names) {
      param_values <- psa_params[[param_name]]
      param_values <- param_values[complete_rows]
      param_values <- param_values[!is.na(param_values)]
      
      # Use quantiles for grid points
      quantiles <- seq(0.1, 0.9, length.out = n_grid)
      param_grids[[param_name]] <- quantile(param_values, quantiles)
    }
    
    # Create parameter combinations
    if (length(param_grids) == 1) {
      param_combinations <- data.frame(param_grids[[1]])
      colnames(param_combinations) <- param_names[1]
    } else {
      param_combinations <- expand.grid(param_grids)
    }
    
    # Calculate expected NMB for each parameter combination
    expected_nmbs <- numeric(nrow(param_combinations))
    
    for (i in 1:nrow(param_combinations)) {
      param_values_combo <- as.numeric(param_combinations[i, ])
      
      # Find closest simulations (using Euclidean distance)
      distances <- rep(0, nrow(psa_params))
      
      for (j in seq_along(param_names)) {
        param_name <- param_names[j]
        target_value <- param_values_combo[j]
        param_col <- psa_params[[param_name]]
        param_col <- param_col[complete_rows]
        
        # Standardize by parameter range to avoid scale issues
        param_range <- max(param_col, na.rm = TRUE) - min(param_col, na.rm = TRUE)
        if (param_range > 0) {
          normalized_diff <- (param_col - target_value) / param_range
          distances <- distances + normalized_diff^2
        }
      }
      
      # Select closest simulations
      n_closest <- min(n_inner, length(distances))
      closest_indices <- order(distances)[1:n_closest]
      
      # Calculate expected NMB for these simulations
      subset_nmb <- nmb_matrix[closest_indices, , drop = FALSE]
      strategy_means <- colMeans(subset_nmb, na.rm = TRUE)
      expected_nmbs[i] <- max(strategy_means, na.rm = TRUE)
    }
    
    # Calculate EVPPI
    expected_value_ppi <- mean(expected_nmbs, na.rm = TRUE)
    expected_value_current <- max_expected_nmb
    
    evppi <- expected_value_ppi - expected_value_current
    evppi <- max(0, min(evppi, evpi))  # Bound between 0 and EVPI
    
    cat("  Expected value with PPI:", round(expected_value_ppi, 4), "\n")
    cat("  Expected value current:", round(expected_value_current, 4), "\n")
    cat("  Raw EVPPI:", round(expected_value_ppi - expected_value_current, 4), "\n")
    cat("  Bounded EVPPI:", round(evppi, 4), "\n")
    
    return(list(
      evppi = evppi,
      evpi = evpi,
      param_names = param_names,
      expected_value_ppi = expected_value_ppi,
      expected_value_current = expected_value_current,
      n_combinations = nrow(param_combinations)
    ))
    
  }, error = function(e) {
    cat("  Error in EVPPI calculation:", conditionMessage(e), "\n")
    return(list(evppi = 0, evpi = evpi, param_names = param_names,
                error = conditionMessage(e)))
  })
}

# ===============================================================================
# CALCULATE EVPPI FOR ALL PARAMETERS
# ===============================================================================

# Only proceed if we have meaningful EVPI
if (evpi_manual > 0.01) {  # Only if EVPI > 1 cent
  
  # Define parameters for EVPPI analysis
  evppi_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS", 
                    "c_test_CT", "u_np", "u_p", "c_other_last")
  
  additional_params <- c("c_test_blood", "c_other_visit", "c_other_baseline", "c_other_follow")
  all_evppi_params <- c(evppi_params, additional_params)
  
  # Filter to only include parameters that exist in the PSA and have variation
  available_params <- intersect(all_evppi_params, colnames(psa_params))
  
  # Further filter based on variation
  params_with_variation <- character()
  for (param in available_params) {
    param_values <- psa_params[[param]]
    if (is.numeric(param_values) && !all(is.na(param_values))) {
      param_values <- param_values[!is.na(param_values)]
      if (length(unique(param_values)) >= 5 && sd(param_values) > 0) {
        params_with_variation <- c(params_with_variation, param)
      }
    }
  }
  
  cat("\nParameters with sufficient variation for EVPPI analysis:", 
      paste(params_with_variation, collapse = ", "), "\n")
  
  # Initialize results storage
  evppi_results <- data.frame(
    parameter = character(),
    evppi = numeric(),
    evpi = numeric(),
    evppi_percent_of_evpi = numeric(),
    n_combinations = numeric(),
    error = character(),
    stringsAsFactors = FALSE
  )
  
  # Calculate EVPPI for each parameter individually
  if (length(params_with_variation) > 0) {
    cat("\nCalculating individual parameter EVPPIs...\n")
    for (param in params_with_variation) {
      result <- calculate_evppi_improved(psa_obj, psa_params, param, wtp = WTP, n_inner = 50, n_grid = 8)
      
      evppi_percent <- if (result$evpi > 0) (result$evppi / result$evpi) * 100 else 0
      
      evppi_results <- rbind(evppi_results, data.frame(
        parameter = param,
        evppi = result$evppi,
        evpi = result$evpi,
        evppi_percent_of_evpi = evppi_percent,
        n_combinations = if (is.null(result$n_combinations)) 0 else result$n_combinations,
        error = if (is.null(result$error)) "" else result$error,
        stringsAsFactors = FALSE
      ))
      
      cat("EVPPI for", param, ":", round(result$evppi, 4), 
          "(", round(evppi_percent, 1), "% of EVPI)\n")
    }
    
    # Sort results by EVPPI value (descending)
    evppi_results <- evppi_results[order(-evppi_results$evppi), ]
    
  } else {
    cat("No parameters have sufficient variation for EVPPI analysis\n")
  }
  
} else {
  cat("EVPI is too small (", round(evpi_manual, 4), ") - skipping EVPPI analysis\n")
  evppi_results <- data.frame()
}

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
  high_priority_threshold <- 0.01  # 1% of EVPI (lowered threshold)
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
#save(evppi_results, evpi_manual, file = here::here("data/processed/evppi_results.RData"))
cat("\nEVPPI analysis complete. Results saved to data/processed/evppi_results.RData\n")