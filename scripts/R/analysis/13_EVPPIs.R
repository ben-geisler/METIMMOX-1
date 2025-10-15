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
# EVPPI CALCULATION FUNCTIONS
# ===============================================================================

# Function to calculate EVPPI for a single parameter or group of parameters
calculate_evppi <- function(psa_obj, psa_params, param_names, wtp = WTP, n_inner = 100, use_parallel = TRUE) {
  
  cat("Calculating EVPPI for parameter(s):", paste(param_names, collapse = ", "), "\n")
  
  # Extract cost and effect matrices
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- as.matrix(psa_obj$effect)
  n_sim <- nrow(cost_matrix)
  n_strategies <- ncol(cost_matrix)
  
  # Calculate Net Monetary Benefit (NMB) matrix
  nmb_matrix <- effect_matrix * wtp - cost_matrix
  
  # Calculate expected value of perfect information (EVPI) as baseline
  max_nmb_per_sim <- apply(nmb_matrix, 1, max)
  evpi <- mean(max_nmb_per_sim) - max(colMeans(nmb_matrix))
  
  # Function to calculate expected NMB given partial perfect information
  calculate_expected_nmb_ppi <- function(param_values, cost_matrix, effect_matrix, psa_params, param_names, wtp, n_inner) {
    
    # Create a subset of simulations where the parameters of interest match the given values
    # For continuous parameters, we'll use a tolerance-based approach
    tolerance <- 0.01  # 1% tolerance for matching
    
    matching_sims <- rep(TRUE, nrow(psa_params))
    
    for (i in seq_along(param_names)) {
      param_name <- param_names[i]
      target_value <- param_values[i]
      
      if (param_name %in% colnames(psa_params)) {
        param_values_col <- psa_params[[param_name]]
        # Use relative tolerance for matching
        relative_diff <- abs(param_values_col - target_value) / target_value
        matching_sims <- matching_sims & (relative_diff <= tolerance)
      }
    }
    
    # If no exact matches found, use closest matches
    if (sum(matching_sims) < 10) {
      # Calculate distance for each simulation
      distances <- rep(0, nrow(psa_params))
      for (i in seq_along(param_names)) {
        param_name <- param_names[i]
        if (param_name %in% colnames(psa_params)) {
          target_value <- param_values[i]
          param_values_col <- psa_params[[param_name]]
          # Normalize by parameter value to get relative distance
          distances <- distances + ((param_values_col - target_value) / target_value)^2
        }
      }
      
      # Select the n_inner closest simulations
      closest_indices <- order(distances)[1:min(n_inner, length(distances))]
      matching_sims <- rep(FALSE, nrow(psa_params))
      matching_sims[closest_indices] <- TRUE
    }
    
    # Extract matching simulations
    if (sum(matching_sims) > 0) {
      subset_cost <- cost_matrix[matching_sims, , drop = FALSE]
      subset_effect <- effect_matrix[matching_sims, , drop = FALSE]
      subset_nmb <- subset_effect * wtp - subset_cost
      
      # Calculate expected NMB for each strategy given this parameter value
      strategy_nmb <- colMeans(subset_nmb)
      return(max(strategy_nmb))
    } else {
      # If no matching simulations, return overall expected max NMB
      overall_nmb <- effect_matrix * wtp - cost_matrix
      return(max(colMeans(overall_nmb)))
    }
  }
  
  # Generate a grid of parameter values to integrate over
  # Use quantiles of the parameter distributions
  n_grid_points <- 20  # Number of points for numerical integration
  
  param_grids <- list()
  for (param_name in param_names) {
    if (param_name %in% colnames(psa_params)) {
      param_values <- psa_params[[param_name]]
      # Create grid based on quantiles
      quantiles <- seq(0.05, 0.95, length.out = n_grid_points)
      param_grids[[param_name]] <- quantile(param_values, quantiles)
    }
  }
  
  # Create all combinations of parameter values
  if (length(param_grids) == 1) {
    param_combinations <- data.frame(param_grids[[1]])
    colnames(param_combinations) <- param_names[1]
  } else {
    param_combinations <- expand.grid(param_grids)
  }
  
  # Calculate expected NMB for each parameter combination
  # Check if we can use parallel processing (only on non-Windows systems)
  can_use_parallel <- use_parallel && .Platform$OS.type != "windows" && require(parallel, quietly = TRUE)
  
  if (can_use_parallel) {
    # Use parallel processing on Unix-like systems
    expected_nmbs <- mclapply(1:nrow(param_combinations), function(i) {
      param_values <- as.numeric(param_combinations[i, ])
      calculate_expected_nmb_ppi(param_values, cost_matrix, effect_matrix, 
                                 psa_params, param_names, wtp, n_inner)
    }, mc.cores = detectCores() - 1)
    expected_nmbs <- unlist(expected_nmbs)
  } else {
    # Sequential processing (for Windows or when parallel is not available)
    expected_nmbs <- sapply(1:nrow(param_combinations), function(i) {
      param_values <- as.numeric(param_combinations[i, ])
      calculate_expected_nmb_ppi(param_values, cost_matrix, effect_matrix, 
                                 psa_params, param_names, wtp, n_inner)
    })
  }
  
  # Calculate EVPPI as the difference between expected value with partial perfect information
  # and the expected value with current information
  expected_value_ppi <- mean(expected_nmbs, na.rm = TRUE)
  expected_value_current <- max(colMeans(nmb_matrix))
  
  evppi <- expected_value_ppi - expected_value_current
  
  # Ensure EVPPI is non-negative and not greater than EVPI
  evppi <- max(0, min(evppi, evpi))
  
  return(list(
    evppi = evppi,
    evpi = evpi,
    param_names = param_names,
    expected_value_ppi = expected_value_ppi,
    expected_value_current = expected_value_current
  ))
}

# ===============================================================================
# CALCULATE EVPPI FOR ALL PARAMETERS
# ===============================================================================

# Define parameters for EVPPI analysis (matching those used in DSA)
evppi_params <- c("c_drug_nivo", "c_drug_FLOX", "c_test_NGS", 
                  "c_test_CT", "u_np", "u_p", "c_other_last")

# Also include additional parameters from the PSA
additional_params <- c("c_test_blood", "c_other_visit", "c_other_baseline", "c_other_follow")
all_evppi_params <- c(evppi_params, additional_params)

# Filter to only include parameters that exist in the PSA
available_params <- intersect(all_evppi_params, colnames(psa_params))

cat("Available parameters for EVPPI analysis:", paste(available_params, collapse = ", "), "\n")

# Check operating system and inform user about parallel processing
if (.Platform$OS.type == "windows") {
  cat("Running on Windows - using sequential processing (parallel processing not supported)\n")
  use_parallel_processing <- FALSE
} else {
  cat("Running on Unix-like system - parallel processing available\n")
  use_parallel_processing <- TRUE
}

# Initialize results storage
evppi_results <- data.frame(
  parameter = character(),
  evppi = numeric(),
  evpi = numeric(),
  evppi_percent_of_evpi = numeric(),
  stringsAsFactors = FALSE
)

# Calculate EVPPI for each parameter individually
cat("Calculating individual parameter EVPPIs...\n")
for (param in available_params) {
  tryCatch({
    result <- calculate_evppi(psa_obj, psa_params, param, wtp = WTP, n_inner = 50, use_parallel = use_parallel_processing)
    
    evppi_results <- rbind(evppi_results, data.frame(
      parameter = param,
      evppi = result$evppi,
      evpi = result$evpi,
      evppi_percent_of_evpi = (result$evppi / result$evpi) * 100,
      stringsAsFactors = FALSE
    ))
    
    cat("EVPPI for", param, ":", round(result$evppi, 2), 
        "(", round((result$evppi / result$evpi) * 100, 1), "% of EVPI)\n")
    
  }, error = function(e) {
    cat("Error calculating EVPPI for", param, ":", conditionMessage(e), "\n")
  })
}

# Sort results by EVPPI value (descending)
evppi_results <- evppi_results[order(-evppi_results$evppi), ]

# ===============================================================================
# PARAMETER GROUP EVPPI ANALYSIS
# ===============================================================================

# Calculate EVPPI for parameter groups
cat("\nCalculating parameter group EVPPIs...\n")

# Group 1: Cost parameters
cost_params <- available_params[grepl("^c_", available_params)]
if (length(cost_params) > 1) {
  tryCatch({
    cost_evppi <- calculate_evppi(psa_obj, psa_params, cost_params, wtp = WTP, n_inner = 30, use_parallel = use_parallel_processing)
    cat("EVPPI for all cost parameters:", round(cost_evppi$evppi, 2), 
        "(", round((cost_evppi$evppi / cost_evppi$evpi) * 100, 1), "% of EVPI)\n")
  }, error = function(e) {
    cat("Error calculating EVPPI for cost parameters:", conditionMessage(e), "\n")
  })
}

# Group 2: Utility parameters
utility_params <- available_params[grepl("^u_", available_params)]
if (length(utility_params) > 1) {
  tryCatch({
    utility_evppi <- calculate_evppi(psa_obj, psa_params, utility_params, wtp = WTP, n_inner = 30, use_parallel = use_parallel_processing)
    cat("EVPPI for all utility parameters:", round(utility_evppi$evppi, 2), 
        "(", round((utility_evppi$evppi / utility_evppi$evpi) * 100, 1), "% of EVPI)\n")
  }, error = function(e) {
    cat("Error calculating EVPPI for utility parameters:", conditionMessage(e), "\n")
  })
}

# Group 3: Drug cost parameters
drug_params <- available_params[grepl("c_drug", available_params)]
if (length(drug_params) > 1) {
  tryCatch({
    drug_evppi <- calculate_evppi(psa_obj, psa_params, drug_params, wtp = WTP, n_inner = 30, use_parallel = use_parallel_processing)
    cat("EVPPI for drug cost parameters:", round(drug_evppi$evppi, 2), 
        "(", round((drug_evppi$evppi / drug_evppi$evpi) * 100, 1), "% of EVPI)\n")
  }, error = function(e) {
    cat("Error calculating EVPPI for drug cost parameters:", conditionMessage(e), "\n")
  })
}

# Group 4: Test cost parameters
test_params <- available_params[grepl("c_test", available_params)]
if (length(test_params) > 1) {
  tryCatch({
    test_evppi <- calculate_evppi(psa_obj, psa_params, test_params, wtp = WTP, n_inner = 30, use_parallel = use_parallel_processing)
    cat("EVPPI for test cost parameters:", round(test_evppi$evppi, 2), 
        "(", round((test_evppi$evppi / test_evppi$evpi) * 100, 1), "% of EVPI)\n")
  }, error = function(e) {
    cat("Error calculating EVPPI for test cost parameters:", conditionMessage(e), "\n")
  })
}

# ===============================================================================
# SUMMARY AND OUTPUT
# ===============================================================================

# Print summary table
cat("\n=== EVPPI ANALYSIS SUMMARY ===\n")
cat("Willingness-to-pay threshold:", WTP, "\n")
cat("Number of PSA simulations:", psa_obj$n_sim, "\n")
if (nrow(evppi_results) > 0) {
  cat("Total EVPI:", round(evppi_results$evpi[1], 2), "\n\n")
} else {
  cat("No EVPPI results calculated\n\n")
}

cat("Individual Parameter EVPPIs (ranked by importance):\n")
print(evppi_results)

# Identify high-priority parameters for future research
high_priority_threshold <- 0.05  # 5% of EVPI
if (nrow(evppi_results) > 0) {
  high_priority_params <- evppi_results[evppi_results$evppi_percent_of_evpi > high_priority_threshold * 100, ]
  
  if (nrow(high_priority_params) > 0) {
    cat("\nHigh-priority parameters for future research (>", high_priority_threshold * 100, "% of EVPI):\n")
    for (i in 1:nrow(high_priority_params)) {
      cat("  ", high_priority_params$parameter[i], ": €", 
          round(high_priority_params$evppi[i], 0), "\n")
    }
  } else {
    cat("\nNo parameters exceed the high-priority threshold of", high_priority_threshold * 100, "% of EVPI\n")
  }
}

# Create a simple visualization of results
if (nrow(evppi_results) > 0) {
  # Create a horizontal bar plot
  par(mar = c(5, 8, 4, 2))
  barplot(evppi_results$evppi, 
          names.arg = evppi_results$parameter,
          horiz = TRUE,
          las = 1,
          main = "Expected Value of Partially Perfect Information (EVPPI)",
          xlab = paste("EVPPI (€) at WTP =", WTP),
          col = "steelblue")
  
  # Add percentage labels
  text(evppi_results$evppi + max(evppi_results$evppi) * 0.02, 
       1:nrow(evppi_results) * 1.2 - 0.5,
       paste0(round(evppi_results$evppi_percent_of_evpi, 1), "%"),
       pos = 4, cex = 0.8)
}

# Save results for use in Quarto report
save(evppi_results, file = here::here("data/processed/evppi_results.RData"))

cat("\nEVPPI analysis complete. Results saved to data/processed/evppi_results.RData\n")