# EVPPI calculation functions for value of information analysis

#' Calculate EVPPI for specified parameters with improved methodology
#' 
#' @param psa_obj PSA object created by dampack::make_psa_obj
#' @param psa_params Data frame with PSA parameter samples
#' @param param_names Vector of parameter names to calculate EVPPI for
#' @param wtp Willingness-to-pay threshold
#' @param n_inner Number of inner loop simulations
#' @param n_grid Number of grid points for parameter values
#' @param seed RNG seed used before subsampling multi-parameter combinations
#' @return List with EVPPI results and diagnostics
calculate_evppi_improved <- function(psa_obj, psa_params, param_names, wtp, 
                                     n_inner = 1000, n_grid = 500,
                                     seed = 123L) {
  
  cat("Calculating EVPPI for parameter(s):", paste(param_names, collapse = ", "), "\n")
  
  # Extract cost and effect matrices
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- as.matrix(psa_obj$effect)
  n_sim <- nrow(cost_matrix)
  n_strategies <- ncol(cost_matrix)
  
  # Remove any rows with missing values
  # IMPORTANT: Filter psa_params alongside cost/effect matrices to maintain
  # consistent row alignment for kNN distance calculations
  complete_rows <- complete.cases(cost_matrix) & complete.cases(effect_matrix)
  if (sum(complete_rows) < n_sim) {
    cat("  Warning: Removing", n_sim - sum(complete_rows), "incomplete simulations\n")
    cost_matrix <- cost_matrix[complete_rows, , drop = FALSE]
    effect_matrix <- effect_matrix[complete_rows, , drop = FALSE]
    psa_params <- psa_params[complete_rows, , drop = FALSE]
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
      # For multi-parameter groups, expand.grid creates n_grid^d points which is
      # computationally infeasible (e.g., 500^4 = 62.5 billion for 4 params).
      # Instead, use a random subsample of the actual PSA parameter combinations
      # which naturally cover the joint distribution.
      max_combos <- n_grid  # Cap at n_grid evaluation points
      all_param_values <- psa_params[, param_names, drop = FALSE]
      if (nrow(all_param_values) > max_combos) {
        # Make the subset independent of prior EVPPI groups and scenarios.
        set.seed(seed)
        sample_idx <- sample(nrow(all_param_values), max_combos)
        param_combinations <- all_param_values[sample_idx, , drop = FALSE]
      } else {
        param_combinations <- all_param_values
      }
      rownames(param_combinations) <- NULL
    }
    
    # Calculate expected NMB for each parameter combination
    expected_nmbs <- numeric(nrow(param_combinations))
    
    for (i in seq_len(nrow(param_combinations))) {
      param_values_combo <- as.numeric(param_combinations[i, ])
      
      # Find closest simulations (using Euclidean distance)
      # psa_params is already filtered by complete_rows, so nrow matches nmb_matrix
      distances <- rep(0, nrow(psa_params))

      for (j in seq_along(param_names)) {
        param_name <- param_names[j]
        target_value <- param_values_combo[j]
        param_col <- psa_params[[param_name]]
        
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
#' @param seed RNG seed used for grouped-parameter subsampling
#' @return Data frame with EVPPI results for all parameters and groups
run_evppi_analysis <- function(psa_obj, psa_params, wtp, evppi_params,
                               param_groups = NULL, seed = 123L) {
  
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
                                         wtp = wtp, n_inner = 1000, n_grid = 500,
                                         seed = seed)
      
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
                                           wtp = wtp, n_inner = 1000, n_grid = 100,
                                           seed = seed)

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


#' Extract treatment-biomarker interaction coefficients from sampling models
#'
#' Extracts the biomarker:Rx interaction coefficients from resampled survival
#' models and returns them as additional columns aligned with PSA iterations.
#' These can be appended to psa_params for EVPPI analysis of treatment effect
#' modification parameters.
#'
#' @param sampling_models List of resampled models (global variable from 06_sampling.R)
#' @param n_sim Number of PSA iterations (must match length of sampling_models$*$samples)
#' @return Data frame with n_sim rows and columns for each extracted interaction
#'   coefficient (b_<biomarker>_rx_<outcome>). Returns NULL if extraction fails.
extract_interaction_coefficients <- function(sampling_models, n_sim) {
  if (is.null(sampling_models)) {
    warning("sampling_models is NULL - cannot extract interaction coefficients")
    return(NULL)
  }

  biomarkers <- if (exists("get_biomarkers")) {
    get_biomarkers()
  } else {
    c("crp", "tmb_braf")
  }
  biomarkers <- intersect(biomarkers, c("crp", "tmb_braf"))
  outcomes <- c("os", "pfs")

  # Coefficient name patterns for each biomarker
  # flexsurvreg encodes factor levels: crp -> crp1, tmb_braf -> tmb_braf1
  coef_patterns <- list(
    crp = "crp.*:Rx",
    tmb_braf = "tmb_braf.*:Rx"
  )

  result <- data.frame(row.names = seq_len(n_sim))

  cat("\n=== Extracting interaction coefficients from sampling models ===\n")

  for (biomarker in biomarkers) {
    if (!biomarker %in% names(sampling_models)) {
      cat("  Warning:", biomarker, "not found in sampling_models - skipping\n")
      next
    }

    strategy_models <- sampling_models[[biomarker]]

    for (outcome in outcomes) {
      col_name <- paste0("b_", biomarker, "_rx_", outcome)
      values <- numeric(n_sim)
      n_extracted <- 0
      n_failed <- 0
      n_missing <- 0

      for (i in seq_len(n_sim)) {
        sample_i <- strategy_models$samples[[i]]

        if (is.null(sample_i)) {
          values[i] <- NA
          n_failed <- n_failed + 1
          next
        }

        coefs <- sample_i[[outcome]]$coefficients

        if (is.null(coefs)) {
          values[i] <- NA
          n_failed <- n_failed + 1
          next
        }

        # Find the interaction coefficient by pattern matching
        pattern <- coef_patterns[[biomarker]]
        matching_names <- grep(pattern, names(coefs), value = TRUE,
                               ignore.case = TRUE)

        if (length(matching_names) > 0) {
          values[i] <- coefs[matching_names[1]]
          n_extracted <- n_extracted + 1
        } else {
          values[i] <- NA
          n_missing <- n_missing + 1
        }
      }

      result[[col_name]] <- values

      cat(sprintf("  %s: extracted=%d, failed=%d, missing=%d\n",
                  col_name, n_extracted, n_failed, n_missing))

      if (n_extracted / n_sim < 0.9) {
        warning("Only ", round(n_extracted / n_sim * 100, 1),
                "% of iterations have valid ", col_name, " values")
      }
    }
  }

  # Remove columns that are entirely NA (coefficient not in model formula)
  all_na_cols <- vapply(result, function(x) all(is.na(x)), logical(1))
  if (any(all_na_cols)) {
    cat("  Removing all-NA columns:",
        paste(names(result)[all_na_cols], collapse = ", "), "\n")
    result <- result[, !all_na_cols, drop = FALSE]
  }

  # Report summary statistics
  cat("\nExtracted coefficient summary:\n")
  for (col in names(result)) {
    vals <- result[[col]]
    valid_vals <- vals[!is.na(vals)]
    if (length(valid_vals) > 0) {
      cat(sprintf("  %s: mean=%.4f, sd=%.4f, range=[%.4f, %.4f], n_valid=%d\n",
                  col, mean(valid_vals), sd(valid_vals),
                  min(valid_vals), max(valid_vals), length(valid_vals)))
    }
  }

  if (ncol(result) == 0) {
    warning("No interaction coefficients could be extracted")
    return(NULL)
  }

  return(result)
}


#' Get interaction EVPPI parameter names for the economic model
#'
#' @return Character vector of parameter names (e.g., "b_crp_rx_os")
get_interaction_evppi_params <- function() {
  biomarkers <- if (exists("get_biomarkers")) {
    get_biomarkers()
  } else {
    c("crp", "tmb_braf")
  }
  biomarkers <- intersect(biomarkers, c("crp", "tmb_braf"))
  outcomes <- c("os", "pfs")

  params <- character(0)
  for (biomarker in biomarkers) {
    for (outcome in outcomes) {
      params <- c(params, paste0("b_", biomarker, "_rx_", outcome))
    }
  }

  return(params)
}


#' Add interaction parameter groups for EVPPI analysis
#'
#' Extends existing param_groups with per-biomarker interaction groups.
#' Only adds groups of 2 parameters (OS + PFS per biomarker) because the
#' kNN EVPPI method uses expand.grid() making larger groups infeasible.
#'
#' @param existing_groups Named list of existing parameter groups
#' @return Updated named list with additional interaction groups
add_interaction_param_groups <- function(existing_groups) {
  interaction_params <- get_interaction_evppi_params()

  if (length(interaction_params) == 0) {
    return(existing_groups)
  }

  # Per-biomarker groups (OS + PFS together, 2 params each - feasible)
  biomarkers <- if (exists("get_biomarkers")) {
    get_biomarkers()
  } else {
    c("crp", "tmb_braf")
  }
  biomarkers <- intersect(biomarkers, c("crp", "tmb_braf"))

  for (biomarker in biomarkers) {
    bio_params <- grep(paste0("^b_", biomarker, "_rx_"),
                       interaction_params, value = TRUE)
    if (length(bio_params) >= 2) {
      existing_groups[[paste0("interaction_", biomarker)]] <- bio_params
    }
  }

  return(existing_groups)
}
