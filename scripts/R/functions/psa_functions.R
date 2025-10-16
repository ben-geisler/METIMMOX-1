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
#' @param cl Cluster object for parallel processing
#' @param n_sim Number of simulations
#' @return List with cost and effect matrices
run_psa_analysis <- function(psa_params, l_params_base, param_distributions, 
                             strategies, time_horizon, cl, n_sim) {
  # Create empty matrices to store results
  cost_matrix <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                        dimnames = list(NULL, strategies))
  effect_matrix <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                          dimnames = list(NULL, strategies))
  
  # Run the model for each simulation
  for (i in 1:n_sim) {
    # Create parameter set for this simulation
    sim_params <- l_params_base
    
    # Update parameters with PSA sample values
    for (param_name in names(param_distributions)) {
      sim_params[[param_name]] <- psa_params[[param_name]][i]
    }
    
    # Use tryCatch to handle any errors in individual simulations
    tryCatch({
      sim_results <- model_fun(
        params = sim_params,
        time_horizon = time_horizon,
        cl = cl,
        determpsa = "psa",
        return_traces = FALSE
      )
      
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
  }
  
  # Replace any remaining NAs with column means
  for (strat in strategies) {
    cost_matrix[is.na(cost_matrix[, strat]), strat] <- mean(cost_matrix[, strat], na.rm = TRUE)
    effect_matrix[is.na(effect_matrix[, strat]), strat] <- mean(effect_matrix[, strat], na.rm = TRUE)
  }
  
  # Return the cost and effect matrices separately
  return(list(
    cost = cost_matrix,
    effect = effect_matrix
  ))
}