# PSA-related functions for probabilistic sensitivity analysis

#' Generate PSA parameter samples from distributions
#' 
#' @param param_distributions List of parameter distributions with dist, shape, and rate/scale info
#' @param n_sim Number of simulations
#' @return Data frame with sampled parameter values
generate_psa_samples <- function(param_distributions, n_sim) {
  samples <- data.frame(sim = 1:n_sim)
  
  # Sample from each parameter distribution
  for (param_name in names(param_distributions)) {
    dist_info <- param_distributions[[param_name]]
    
    if (dist_info$dist == "lnorm") {
      samples[[param_name]] <- rlnorm(n_sim, 
                                      meanlog = dist_info$meanlog, 
                                      sdlog = dist_info$sdlog)
    } else if (dist_info$dist == "beta") {
      samples[[param_name]] <- rbeta(n_sim, 
                                     shape1 = dist_info$shape1, 
                                     shape2 = dist_info$shape2)
    } else if (dist_info$dist == "gamma") {
      samples[[param_name]] <- rgamma(n_sim, 
                                      shape = dist_info$shape, 
                                      rate = dist_info$rate)
    } else if (dist_info$dist == "unif") {
      samples[[param_name]] <- runif(n_sim, 
                                     min = dist_info$min, 
                                     max = dist_info$max)
    }
  }
  
  return(samples)
}

#' Run PSA analysis across all simulations
#'
#' @param psa_params Data frame with PSA parameter samples
#' @param l_params_base Base case parameter list
#' @param param_distributions Parameter distributions list
#' @param strategies Vector of strategy names
#' @param time_horizon Time horizon for analysis
#' @param cl Cycle length (1/52 weeks)
#' @param n_sim Number of simulations
#' @param fallback_threshold Maximum permitted proportion of initial PSA
#'   iterations using fallback before the analysis stops
#' @return List with cost and effect matrices and fallback diagnostics
run_psa_analysis <- function(psa_params, l_params_base, param_distributions,
                             strategies, time_horizon, cl, n_sim,
                             fallback_threshold = 0.02) {
  if (length(fallback_threshold) != 1L || !is.finite(fallback_threshold) ||
      fallback_threshold < 0 || fallback_threshold > 1) {
    stop("fallback_threshold must be a single finite value between 0 and 1")
  }

  # Create empty matrices to store results
  cost_matrix <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                        dimnames = list(NULL, strategies))
  effect_matrix <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                          dimnames = list(NULL, strategies))

  # Track PFS > OS constraint violations (Issue #76)
  # Instead of printing thousands of warnings, we aggregate and summarize
  pfs_os_violations <- list(
    control = integer(0),
    crp_pos = integer(0),
    crp_neg = integer(0),
    tmb_braf_pos = integer(0),
    tmb_braf_neg = integer(0)
  )

  # Track iterations that used fallback to base case (Issue #79)
  # These will be replaced with additional successful simulations
  fallback_iterations <- integer(0)

  # Progress reporting
  cat("Running PSA with", n_sim, "simulations...\n")
  start_time <- Sys.time()
  
  # Run the model for each simulation
  for (i in 1:n_sim) {
    # Create parameter set for this simulation
    sim_params <- l_params_base
    
    # Update parameters with PSA sample values (costs and utilities)
    for (param_name in names(param_distributions)) {
      sim_params[[param_name]] <- psa_params[[param_name]][i]
    }
    # Keep the data-driven diagnostic-cost lookup aligned with sampled scalars.
    sim_params$c_test_biomarker$crp <- sim_params$c_test_CRP
    sim_params$c_test_biomarker$tmb_braf <- sim_params$c_test_NGS
    
    # Use withCallingHandlers to capture warnings, tryCatch for errors
    # This allows us to aggregate PFS > OS warnings instead of printing thousands
    tryCatch({
      withCallingHandlers({
        sim_results <- model_fun(
          params = sim_params,
          time_horizon = time_horizon,
          cl = cl,
          determpsa = "psa",
          return_traces = FALSE,
          sim_idx = i  # THIS IS THE CRITICAL FIX - passes resampled model index
        )
      }, warning = function(w) {
        # Parse PFS > OS constraint warnings and track them (Issue #76)
        msg <- conditionMessage(w)
        if (grepl("PFS > OS constraint enforced", msg)) {
          # Extract which subgroup triggered the warning
          if (grepl("\\(control\\)", msg)) {
            pfs_os_violations$control <<- c(pfs_os_violations$control, i)
          } else if (grepl("\\(crp\\+\\)", msg)) {
            pfs_os_violations$crp_pos <<- c(pfs_os_violations$crp_pos, i)
          } else if (grepl("\\(crp-\\)", msg)) {
            pfs_os_violations$crp_neg <<- c(pfs_os_violations$crp_neg, i)
          } else if (grepl("\\(tmb_braf\\+\\)", msg)) {
            pfs_os_violations$tmb_braf_pos <<- c(pfs_os_violations$tmb_braf_pos, i)
          } else if (grepl("\\(tmb_braf-\\)", msg)) {
            pfs_os_violations$tmb_braf_neg <<- c(pfs_os_violations$tmb_braf_neg, i)
          }
          # Suppress the warning (don't print it)
          invokeRestart("muffleWarning")
        }
        # Let other warnings through
      })

      # Check if fallback was used (Issue #79)
      if (isTRUE(attr(sim_results, "fallback_used"))) {
        fallback_iterations <- c(fallback_iterations, i)
      }

      # Store results
      for (strat in strategies) {
        strat_row <- which(sim_results$Strategy == strat)
        if (length(strat_row) > 0) {
          cost_matrix[i, strat] <- sim_results$Cost[strat_row]
          effect_matrix[i, strat] <- sim_results$Effect[strat_row]
        } else {
          # If strategy not found, flag for replacement rather than silently
          # using running mean which would introduce bias
          warning("Strategy '", strat, "' not found in sim ", i, " results - flagging for replacement")
          cost_matrix[i, strat] <- NA
          effect_matrix[i, strat] <- NA
          if (!(i %in% fallback_iterations)) {
            fallback_iterations <- c(fallback_iterations, i)
          }
        }
      }
    }, error = function(e) {
      cat("Error in simulation", i, ":", conditionMessage(e), "\n")
      # Track this as a fallback iteration (Issue #79) - will be replaced
      fallback_iterations <<- c(fallback_iterations, i)
      # Temporarily store NAs - will be replaced in replacement phase
      for (strat in strategies) {
        cost_matrix[i, strat] <<- NA
        effect_matrix[i, strat] <<- NA
      }
    })
    
    # Progress reporting
    if (i %% 100 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- i / elapsed
      remaining <- (n_sim - i) / rate
      cat(sprintf("Completed %d/%d PSA iterations (%.1f%%) - ETA: %.1f min\n", 
                  i, n_sim, 100*i/n_sim, remaining))
    }
  }
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  cat(sprintf("PSA initial pass complete: %d simulations in %.1f minutes\n", n_sim, total_time))

  # Report and bound initial-pass fallback use before attempting replacements.
  # A high rate indicates systematic survival-prediction failures and would
  # otherwise silently reduce the uncertainty represented by the PSA.
  fallback_iterations <- unique(fallback_iterations)
  fallback_count <- length(fallback_iterations)
  fallback_rate <- fallback_count / n_sim
  cat(sprintf(
    "PSA initial-pass fallback rate: %d/%d (%.2f%%; maximum permitted %.2f%%)\n",
    fallback_count, n_sim, 100 * fallback_rate, 100 * fallback_threshold
  ))

  if (fallback_rate > fallback_threshold) {
    stop(sprintf(
      paste0(
        "PSA initial-pass fallback rate %.2f%% (%d/%d) exceeds the ",
        "permitted %.2f%% threshold; investigate survival-prediction failures."
      ),
      100 * fallback_rate, fallback_count, n_sim, 100 * fallback_threshold
    ))
  }

  # =========================================================================
  # REPLACE FALLBACK ITERATIONS (Issue #79)
  # =========================================================================
  # Instead of mixing base case fallbacks with resampled model results,
  # run additional simulations to replace failed iterations

  if (length(fallback_iterations) > 0) {
    cat(sprintf("\n--- Replacing %d fallback iterations (Issue #79) ---\n",
                length(fallback_iterations)))

    replacement_idx <- n_sim + 1  # Start with models beyond initial set
    replaced_count <- 0
    max_attempts <- length(fallback_iterations) * 10  # Safety limit
    attempt <- 0

    for (failed_i in fallback_iterations) {
      success <- FALSE

      while (!success && attempt < max_attempts) {
        attempt <- attempt + 1

        # Use modulo to wrap around if we exceed available models
        # This ensures we always have a valid model index
        actual_sim_idx <- ((replacement_idx - 1) %% n_sim) + 1

        # Create parameter set (use same PSA samples as original iteration)
        sim_params <- l_params_base
        for (param_name in names(param_distributions)) {
          sim_params[[param_name]] <- psa_params[[param_name]][failed_i]
        }
        sim_params$c_test_biomarker$crp <- sim_params$c_test_CRP
        sim_params$c_test_biomarker$tmb_braf <- sim_params$c_test_NGS

        # Try to run with a different resampled model
        tryCatch({
          suppressWarnings({
            sim_results <- model_fun(
              params = sim_params,
              time_horizon = time_horizon,
              cl = cl,
              determpsa = "psa",
              return_traces = FALSE,
              sim_idx = actual_sim_idx
            )
          })

          # Check if this one succeeded (no fallback)
          if (!isTRUE(attr(sim_results, "fallback_used"))) {
            # Replace the failed iteration's results
            for (strat in strategies) {
              strat_row <- which(sim_results$Strategy == strat)
              if (length(strat_row) > 0) {
                cost_matrix[failed_i, strat] <- sim_results$Cost[strat_row]
                effect_matrix[failed_i, strat] <- sim_results$Effect[strat_row]
              }
            }
            success <- TRUE
            replaced_count <- replaced_count + 1
          }
        }, error = function(e) {
          # This model also failed, try next
        })

        replacement_idx <- replacement_idx + 1
      }

      if (!success) {
        warning("Could not find valid replacement for iteration ", failed_i,
                " after ", max_attempts, " attempts")
      }
    }

    cat(sprintf("Successfully replaced %d of %d fallback iterations\n",
                replaced_count, length(fallback_iterations)))
    cat(sprintf("Models tried: %d (wrapped around %d times)\n",
                replacement_idx - n_sim - 1,
                (replacement_idx - n_sim - 1) %/% n_sim))
  }

  # Replace any remaining NAs with column means (safety fallback)
  na_count <- sum(is.na(cost_matrix))
  if (na_count > 0) {
    warning("Still have ", na_count, " NA values after replacement - using column means")
    for (strat in strategies) {
      cost_matrix[is.na(cost_matrix[, strat]), strat] <- mean(cost_matrix[, strat], na.rm = TRUE)
      effect_matrix[is.na(effect_matrix[, strat]), strat] <- mean(effect_matrix[, strat], na.rm = TRUE)
    }
  }

  # Print summary of PFS > OS constraint violations (Issue #76)
  total_violations <- sum(sapply(pfs_os_violations, length))
  if (total_violations > 0) {
    # Count unique iterations with any violation
    all_violation_iters <- unique(unlist(pfs_os_violations))
    n_iters_with_violations <- length(all_violation_iters)

    cat("\n--- PFS > OS Constraint Summary (Issue #76) ---\n")
    cat(sprintf("Violations occurred in %d of %d iterations (%.1f%%)\n",
                n_iters_with_violations, n_sim, 100 * n_iters_with_violations / n_sim))
    cat("Breakdown by subgroup:\n")

    # Only print subgroups that had violations
    subgroup_names <- c(
      control = "  Control",
      crp_pos = "  CRP+",
      crp_neg = "  CRP-",
      tmb_braf_pos = "  TMB/BRAF+",
      tmb_braf_neg = "  TMB/BRAF-"
    )

    for (subgroup in names(pfs_os_violations)) {
      n_violations <- length(pfs_os_violations[[subgroup]])
      if (n_violations > 0) {
        cat(sprintf("%s: %d iterations (%.1f%%)\n",
                    subgroup_names[subgroup], n_violations, 100 * n_violations / n_sim))
      }
    }
    cat("Constraint enforced: PFS capped at OS to ensure valid state occupancy.\n")
    cat("------------------------------------------------\n\n")
  }

  # Return the cost and effect matrices separately
  return(list(
    cost = cost_matrix,
    effect = effect_matrix,
    fallback_count = fallback_count,
    fallback_rate = fallback_rate,
    fallback_threshold = fallback_threshold
  ))
}
