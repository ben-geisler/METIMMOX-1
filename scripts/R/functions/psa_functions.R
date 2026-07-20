# PSA-related functions for probabilistic sensitivity analysis

#' Generate PSA parameter samples from distributions
#' 
#' @param param_distributions List of parameter distributions with dist, shape, and rate/scale info
#' @param n_sim Number of simulations
#' @param seed RNG seed used immediately before drawing PSA parameters
#' @return Data frame with sampled parameter values
generate_psa_samples <- function(param_distributions, n_sim, seed = 123L) {
  samples <- data.frame(sim = 1:n_sim)

  # Make draws independent of RNG use in earlier analysis or cache branches.
  set.seed(seed)

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

  attr(samples, "seed") <- seed
  return(samples)
}

#' Check whether cached PSA samples were generated with the requested seed
#'
#' @param psa_params Data frame returned by generate_psa_samples()
#' @param seed Expected RNG seed
#' @return Logical scalar
psa_samples_seed_matches <- function(psa_params, seed = 123L) {
  identical(attr(psa_params, "seed"), seed)
}

#' Return the failed-draw handling policy used by PSA caches
#'
#' @return Character policy identifier
psa_failed_draw_policy <- function() {
  "cyclic_other_models_drop_unreplaced_v2"
}

#' List replacement survival-model indices for a failed PSA draw
#'
#' The cache contains exactly `n_sim` resampled survival models, all of which
#' are assigned during the initial PSA pass. There is therefore no unused model
#' bank for replacements. A failed non-survival parameter draw is instead
#' re-paired with up to 10 *other* cached survival models, each used at most
#' once for that draw, in deterministic cyclic order starting after its
#' original model. If every candidate fails, the complete PSA draw is dropped.
#'
#' @param failed_i Index of the failed PSA draw
#' @param n_sim Number of cached resampled survival models
#' @param max_attempts Maximum distinct alternative models to return
#' @return Integer vector of candidate model indices
psa_replacement_candidates <- function(failed_i, n_sim, max_attempts = 10L) {
  if (n_sim <= 1L) {
    return(integer(0))
  }

  candidate_count <- min(as.integer(max_attempts), n_sim - 1L)
  as.integer(((failed_i + seq_len(candidate_count) - 1L) %% n_sim) + 1L)
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
#' @return List with cost and effect matrices, retained iteration indices, and
#'   fallback/drop diagnostics
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
  for (i in seq_len(n_sim)) {
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

  replaced_iterations <- integer(0)

  if (length(fallback_iterations) > 0) {
    cat(sprintf("\n--- Replacing %d fallback iterations (Issue #79) ---\n",
                length(fallback_iterations)))

    replaced_count <- 0
    total_attempts <- 0L

    for (failed_i in fallback_iterations) {
      success <- FALSE

      # There are no models beyond the initial n_sim cached draws. Re-pair the
      # failed economic-parameter draw with up to 10 other cached survival
      # models, each at most once. Start after its original index to avoid
      # favoring low indices and bound the work for an unrecoverable draw.
      replacement_candidates <- psa_replacement_candidates(failed_i, n_sim)
      max_attempts <- length(replacement_candidates)

      # Each failed draw receives its own complete candidate scan. A shared
      # counter would let early failures exhaust the budget for later failures.
      attempt <- 0L

      while (!success && attempt < max_attempts) {
        attempt <- attempt + 1
        total_attempts <- total_attempts + 1L
        actual_sim_idx <- replacement_candidates[attempt]

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

          # A replacement is valid only when it did not use fallback and
          # returned exactly one complete result for every strategy.
          strategy_rows <- match(strategies, sim_results$Strategy)
          replacement_complete <-
            all(!is.na(strategy_rows)) &&
            all(is.finite(sim_results$Cost[strategy_rows])) &&
            all(is.finite(sim_results$Effect[strategy_rows]))

          if (!isTRUE(attr(sim_results, "fallback_used")) &&
              replacement_complete) {
            # Replace the failed iteration's results
            cost_matrix[failed_i, ] <- sim_results$Cost[strategy_rows]
            effect_matrix[failed_i, ] <- sim_results$Effect[strategy_rows]
            success <- TRUE
            replaced_count <- replaced_count + 1
            replaced_iterations <- c(replaced_iterations, failed_i)
          }
        }, error = function(e) {
          # This model also failed, try next
        })
      }

      if (!success) {
        warning("Could not find valid replacement for iteration ", failed_i,
                " after ", max_attempts, " attempts")
      }
    }

    cat(sprintf("Successfully replaced %d of %d fallback iterations\n",
                replaced_count, length(fallback_iterations)))
    cat(sprintf(
      paste0(
        "Replacement model attempts: %d ",
        "(policy: up to 10 other cached models, each at most once ",
        "per failed draw)\n"
      ),
      total_attempts
    ))
  }

  # Drop entire draws that could not be replaced. Keeping base-case fallback
  # values or mean-imputing individual cells would artificially reduce PSA
  # variance and bias CEAC/EVPI/EVPPI estimates.
  unreplaced_iterations <- setdiff(fallback_iterations, replaced_iterations)
  incomplete_iterations <- which(
    !apply(cost_matrix, 1L, function(x) all(is.finite(x))) |
      !apply(effect_matrix, 1L, function(x) all(is.finite(x)))
  )
  dropped_iterations <- sort(unique(c(
    unreplaced_iterations,
    incomplete_iterations
  )))
  dropped_count <- length(dropped_iterations)
  retained_iterations <- setdiff(seq_len(n_sim), dropped_iterations)

  if (dropped_count > 0L) {
    warning(
      "Dropping ", dropped_count,
      " PSA iteration(s) that remained invalid after replacement; ",
      length(retained_iterations), " simulation(s) remain.",
      call. = FALSE
    )
    cost_matrix <- cost_matrix[retained_iterations, , drop = FALSE]
    effect_matrix <- effect_matrix[retained_iterations, , drop = FALSE]
  }

  if (length(retained_iterations) == 0L) {
    stop("No valid PSA iterations remain after dropping failed draws.")
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
    dropped_count = dropped_count,
    dropped_iterations = dropped_iterations,
    retained_iterations = retained_iterations,
    n_sim = length(retained_iterations),
    failed_draw_policy = psa_failed_draw_policy(),
    replacement_model_policy = "cyclic_next_10_other_cached_models_v1",
    fallback_rate = fallback_rate,
    fallback_threshold = fallback_threshold
  ))
}
