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
#' @return List with cost and effect matrices
run_psa_analysis <- function(psa_params, l_params_base, param_distributions,
                             strategies, time_horizon, cl, n_sim) {
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
    tlr_pos = integer(0),
    tlr_neg = integer(0),
    tmb_braf_pos = integer(0),
    tmb_braf_neg = integer(0)
  )

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
          } else if (grepl("\\(tlr\\+\\)", msg)) {
            pfs_os_violations$tlr_pos <<- c(pfs_os_violations$tlr_pos, i)
          } else if (grepl("\\(tlr-\\)", msg)) {
            pfs_os_violations$tlr_neg <<- c(pfs_os_violations$tlr_neg, i)
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

      # Store results
      for (strat in strategies) {
        strat_row <- which(sim_results$Strategy == strat)
        if (length(strat_row) > 0) {
          cost_matrix[i, strat] <- sim_results$Cost[strat_row]
          effect_matrix[i, strat] <- sim_results$Effect[strat_row]
        } else {
          # If strategy not found, use mean or base case
          if (i > 1) {
            cost_matrix[i, strat] <- mean(cost_matrix[1:(i-1), strat], na.rm = TRUE)
            effect_matrix[i, strat] <- mean(effect_matrix[1:(i-1), strat], na.rm = TRUE)
          } else {
            # For first simulation, use NA and fix later
            cost_matrix[i, strat] <- NA
            effect_matrix[i, strat] <- NA
          }
        }
      }
    }, error = function(e) {
      cat("Error in simulation", i, ":", conditionMessage(e), "\n")
      # Use mean values for failed simulations
      for (strat in strategies) {
        if (i > 1) {
          cost_matrix[i, strat] <<- mean(cost_matrix[1:(i-1), strat], na.rm = TRUE)
          effect_matrix[i, strat] <<- mean(effect_matrix[1:(i-1), strat], na.rm = TRUE)
        } else {
          # For first simulation, use base case
          base_results <- model_fun(
            params = l_params_base,
            time_horizon = time_horizon,
            cl = cl,
            determpsa = "det"
          )
          strat_row <- which(base_results$Strategy == strat)
          if (length(strat_row) > 0) {
            cost_matrix[i, strat] <<- base_results$Cost[strat_row]
            effect_matrix[i, strat] <<- base_results$Effect[strat_row]
          } else {
            cost_matrix[i, strat] <<- NA
            effect_matrix[i, strat] <<- NA
          }
        }
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
  cat(sprintf("PSA complete: %d simulations in %.1f minutes\n", n_sim, total_time))
  
  # Replace any remaining NAs with column means
  for (strat in strategies) {
    cost_matrix[is.na(cost_matrix[, strat]), strat] <- mean(cost_matrix[, strat], na.rm = TRUE)
    effect_matrix[is.na(effect_matrix[, strat]), strat] <- mean(effect_matrix[, strat], na.rm = TRUE)
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
      tlr_pos = "  TLR+",
      tlr_neg = "  TLR-",
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
    effect = effect_matrix
  ))
}