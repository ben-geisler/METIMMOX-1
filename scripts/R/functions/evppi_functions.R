# EVPPI calculation functions for value of information analysis

#' Calculate EVPPI for specified parameters with improved methodology
#' 
#' @param psa_obj PSA object created by dampack::make_psa_obj
#' @param psa_params Data frame with PSA parameter samples
#' @param param_names Vector of parameter names to calculate EVPPI for
#' @param wtp Willingness-to-pay threshold
#' @param n_inner Number of inner loop simulations
#' @param n_grid Number of grid points for parameter values
#' @return List with EVPPI results and diagnostics
calculate_evppi_improved <- function(psa_obj, psa_params, param_names, wtp, 
                                     n_inner = 1000, n_grid = 500) {
  
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
    # Create parameter grids using quantiles
    param_grids <- list()
    for (param_name in param_names) {
      param_values <- psa_params[[param_name]]
      param_values <- param_values[complete_rows]
      param_values <- param_values[!is.na(param_values)]
      
      # Use quantiles for grid points (500 points from 0.001 to 0.999)
      quantiles <- seq(0.001, 0.999, length.out = n_grid)
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
      
      # Select closest simulations (1000 nearest neighbors)
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
      n_combinations = nrow(param_combinations),
      n_inner = n_inner,
      n_grid = n_grid
    ))
    
  }, error = function(e) {
    cat("  Error in EVPPI calculation:", conditionMessage(e), "\n")
    return(list(evppi = 0, evpi = evpi, param_names = param_names,
                error = conditionMessage(e)))
  })
}

#' Run complete EVPPI analysis for all parameters
#'
#' @param psa_obj PSA object from dampack
#' @param psa_params Data frame with PSA parameter samples
#' @param wtp Willingness-to-pay threshold
#' @param evppi_params Vector of parameter names to analyze
#' @param param_groups Optional named list of parameter groups for joint EVPPI
#' @return Data frame with EVPPI results for all parameters and groups
run_evppi_analysis <- function(psa_obj, psa_params, wtp, evppi_params,
                               param_groups = NULL) {
  
  # Calculate NMB and EVPI for diagnostics
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- as.matrix(psa_obj$effect)
  nmb_matrix <- effect_matrix * wtp - cost_matrix
  
  # Calculate EVPI manually for verification
  max_nmb_per_sim <- apply(nmb_matrix, 1, max, na.rm = TRUE)
  expected_max_nmb <- mean(max_nmb_per_sim, na.rm = TRUE)
  max_expected_nmb <- max(colMeans(nmb_matrix, na.rm = TRUE))
  evpi_manual <- expected_max_nmb - max_expected_nmb
  
  cat("\n=== EVPPI Analysis at WTP =", wtp, "===\n")
  cat("Total EVPI:", round(evpi_manual, 4), "\n")
  cat("Using n_inner = 1000 and n_grid = 500 for EVPPI calculations\n\n")
  
  # Only proceed if we have meaningful EVPI
  if (evpi_manual <= 0.01) {
    cat("EVPI is too small (", round(evpi_manual, 4), ") - skipping EVPPI analysis\n")
    return(data.frame())
  }
  
  # Filter to only include parameters that exist and have variation
  available_params <- intersect(evppi_params, colnames(psa_params))
  
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
  
  cat("Parameters with sufficient variation:", 
      paste(params_with_variation, collapse = ", "), "\n\n")
  
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
    for (param in params_with_variation) {
      result <- calculate_evppi_improved(psa_obj, psa_params, param, 
                                         wtp = wtp, n_inner = 1000, n_grid = 500)
      
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
  }

  # Calculate EVPPI for parameter groups (if provided)
  if (!is.null(param_groups) && length(param_groups) > 0) {
    cat("\n=== Grouped Parameter EVPPI ===\n")

    for (group_name in names(param_groups)) {
      group_params <- param_groups[[group_name]]
      # Filter to params with variation
      group_params <- intersect(group_params, params_with_variation)

      if (length(group_params) >= 2) {
        cat("\nCalculating EVPPI for group:", group_name,
            "(", paste(group_params, collapse = ", "), ")\n")

        result <- calculate_evppi_improved(psa_obj, psa_params, group_params,
                                           wtp = wtp, n_inner = 1000, n_grid = 100)

        evppi_percent <- if (result$evpi > 0) {
          (result$evppi / result$evpi) * 100
        } else {
          0
        }

        evppi_results <- rbind(evppi_results, data.frame(
          parameter = paste0("[GROUP] ", group_name),
          evppi = result$evppi,
          evpi = result$evpi,
          evppi_percent_of_evpi = evppi_percent,
          n_combinations = if (is.null(result$n_combinations)) 0
                           else result$n_combinations,
          error = if (is.null(result$error)) "" else result$error,
          stringsAsFactors = FALSE
        ))

        cat("EVPPI for group", group_name, ":", round(result$evppi, 4),
            "(", round(evppi_percent, 1), "% of EVPI)\n")
      }
    }
  }

  return(evppi_results)
}

#' Calculate Population-Level EVPPI
#'
#' Scales per-patient EVPPI values to the eligible population level
#' using annual incidence, research horizon, and discounting.
#'
#' @param evppi_per_patient Numeric. Per-patient EVPPI value in EUR.
#' @param annual_incidence Numeric. Annual eligible patients (default: 1500).
#' @param research_horizon Numeric. Years of research value (default: 10).
#' @param discount_rate Numeric. Annual discount rate for research benefits (default: 0.035).
#' @return Numeric. Population-level EVPPI in millions of EUR.
#'
#' @details
#' The calculation accounts for the time value of research by discounting
#' future benefits. The formula is:
#'   Population EVPPI = Per-patient EVPPI × Annual incidence ×
#'                      Sum of discount factors over research horizon
#'
#' Example: For EUR 1,000 per-patient EVPPI with 1,500 patients/year over 10 years:
#'   Sum of discount factors ≈ 8.317 (using 3.5% discount rate)
#'   Population EVPPI = EUR 1,000 × 1,500 × 8.317 = EUR 12.48 million
#'
#' @export
calculate_population_evppi <- function(evppi_per_patient,
                                       annual_incidence = 1500,
                                       research_horizon = 10,
                                       discount_rate = 0.035) {
  # Calculate sum of discount factors over research horizon
  # (Present value of an annuity of 1 unit per year for n years)
  years <- 1:research_horizon
  discount_factors <- 1 / (1 + discount_rate)^years
  sum_discount_factors <- sum(discount_factors)

  # Scale per-patient value to population
  population_evppi_eur <- evppi_per_patient * annual_incidence * sum_discount_factors

  # Convert to millions
  population_evppi_millions <- population_evppi_eur / 1e6

  return(population_evppi_millions)
}