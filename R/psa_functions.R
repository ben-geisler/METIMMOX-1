# PSA-related functions for probabilistic sensitivity analysis

#' Generate PSA parameter samples from distributions
#' 
#' @param param_distributions List of parameter distributions with dist, shape, and rate/scale info
#' @param n_sim Number of simulations
#' @param seed RNG seed used immediately before drawing PSA parameters
#' @return Data frame with sampled parameter values
generate_psa_samples <- function(param_distributions, n_sim, seed = 123L) {
  df_samples <- data.frame(sim = 1:n_sim)

  # Make draws independent of RNG use in earlier analysis or cache branches.
  set.seed(seed)
  m_joint_draws <- NULL

  # Sample from each parameter distribution
  for (param_name in names(param_distributions)) {
    l_dist_info <- param_distributions[[param_name]]

    if (l_dist_info$dist == "derived") {
      # Filled in below from the parameters it is derived from (issue #154).
      next
    } else if (l_dist_info$dist == "dirichlet") {
      if (is.null(m_joint_draws)) {
        # Bayesian-bootstrap masses of the observed joint biomarker cells:
        # Dirichlet(n_00, n_01, n_10, n_11). Marginals are derived together.
        m_joint_draws <- vapply(l_dist_info$alpha, function(shape)
          rgamma(n_sim, shape = shape, rate = 1), numeric(n_sim))
        m_joint_draws <- m_joint_draws / rowSums(m_joint_draws)
      }
      df_samples[[param_name]] <- m_joint_draws[, l_dist_info$component]
    } else if (l_dist_info$dist == "lnorm") {
      df_samples[[param_name]] <- rlnorm(n_sim, 
                                      meanlog = l_dist_info$meanlog, 
                                      sdlog = l_dist_info$sdlog)
    } else if (l_dist_info$dist == "beta") {
      df_samples[[param_name]] <- rbeta(n_sim, 
                                     shape1 = l_dist_info$shape1, 
                                     shape2 = l_dist_info$shape2)
    } else if (l_dist_info$dist == "gamma") {
      df_samples[[param_name]] <- rgamma(n_sim, 
                                      shape = l_dist_info$shape, 
                                      rate = l_dist_info$rate)
    } else if (l_dist_info$dist == "unif") {
      df_samples[[param_name]] <- runif(n_sim, 
                                     min = l_dist_info$min, 
                                     max = l_dist_info$max)
    }
  }

  # Deterministic functions of the sampled parameters (currently u_p).
  df_samples <- apply_derived_psa_parameters(df_samples)

  v_missing_columns <- setdiff(names(param_distributions), names(df_samples))
  if (length(v_missing_columns) > 0) {
    stop("PSA sampling produced no column for: ",
         paste(v_missing_columns, collapse = ", "))
  }

  attr(df_samples, "seed") <- seed
  return(df_samples)
}

#' Fraction of PSA draws with an implausible utility ordering
#'
#' Reported alongside the PSA results so the ordering assumption is visible
#' (issue #154). The decrement parameterisation makes this zero by
#' construction; independent beta draws for u_np and u_p did not.
#'
#' @param psa_params Data frame of PSA parameter draws.
#' @param margin Excess of u_p over u_np counted as a material reversal.
#' @return List with n, n_reversed, fraction, and fraction_material, or NULL
#'   when the utility columns are absent.
utility_reversal_fraction <- function(psa_params, margin = 0.05) {
  if (!all(c("u_np", "u_p") %in% names(psa_params))) return(NULL)
  v_excess <- psa_params$u_p - psa_params$u_np
  v_excess <- v_excess[is.finite(v_excess)]
  if (length(v_excess) == 0) return(NULL)
  list(
    n = length(v_excess),
    n_reversed = sum(v_excess > 0),
    fraction = mean(v_excess > 0),
    fraction_material = mean(v_excess > margin),
    margin = margin
  )
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
  "finite_outcomes_cyclic_other_models_drop_unreplaced_v3"
}

#' Get keys used to aggregate survival-ordering warnings
#' @return Named list of empty integer vectors.
pfs_os_violation_groups <- function() {
  v_biomarkers <- get_biomarkers()
  v_subgroup_keys <- unlist(lapply(
    v_biomarkers,
    function(x) paste(x, c("positive", "negative"), sep = "_")
  ))
  setNames(replicate(1L + length(v_subgroup_keys), integer(0), simplify = FALSE),
           c(get_control_strategy(), v_subgroup_keys))
}

#' Convert a structured survival-ordering warning to an aggregation key
#' @param warning_condition A `survival_ordering_warning` condition.
#' @return Character scalar.
pfs_os_violation_key <- function(warning_condition) {
  if (identical(warning_condition$subgroup, "control")) {
    return(warning_condition$strategy)
  }
  paste(warning_condition$strategy, warning_condition$subgroup, sep = "_")
}

#' Create display labels for survival-ordering aggregation groups
#' @return Named character vector.
pfs_os_violation_labels <- function() {
  control_strategy <- get_control_strategy()
  v_biomarkers <- get_biomarkers()
  v_labels <- setNames(
    paste0("  ", unname(strategy_display_name(control_strategy))),
    control_strategy
  )
  for (biomarker in v_biomarkers) {
    display_name <- unname(strategy_display_name(biomarker))
    v_labels[[paste0(biomarker, "_positive")]] <- paste0("  ", display_name, "+")
    v_labels[[paste0(biomarker, "_negative")]] <- paste0("  ", display_name, "-")
  }
  v_labels
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

# Outcome columns every PSA draw must return, one finite value per strategy:
# cost, QALYs and discounted/undiscounted life-years (issue #197).
PSA_OUTCOME_COLUMNS <- c("Cost", "Effect", "LY_disc", "LY_undisc")

#' Check that a model_fun() result is a complete PSA outcome row set
#'
#' Both initial draws and replacements must return exactly one finite cost,
#' effect and discounted and undiscounted life-year value per strategy (the
#' result columns PSA_OUTCOME_COLUMNS). Duplicates, missing columns and extra
#' rows are failures.
#'
#' @param results Data frame returned by model_fun().
#' @param strategies Expected strategy identifiers.
#' @return TRUE or FALSE.
psa_outcomes_complete <- function(results, strategies) {
  is.data.frame(results) &&
    all(c("Strategy", PSA_OUTCOME_COLUMNS) %in% names(results)) &&
    nrow(results) == length(strategies) &&
    !anyNA(results$Strategy) && !anyDuplicated(results$Strategy) &&
    setequal(as.character(results$Strategy), strategies) &&
    all(vapply(PSA_OUTCOME_COLUMNS, function(column) {
      is.numeric(results[[column]]) && all(is.finite(results[[column]]))
    }, logical(1)))
}

#' Check the outcomes of one draw of the extended PSA
#'
#' A draw of the extended PSA (issues #197, #200) must carry the complete
#' outcome columns of psa_outcomes_complete() (cost, QALYs and both life-year
#' columns for every strategy) and the cell outcomes of model_fun(return_cells
#' = TRUE) (attribute `cells`): finite values for every observed cell, both
#' arms and every outcome in CELL_OUTCOMES, and joint cell weights that are
#' finite, non-negative and sum to one.
#'
#' @param results Data frame returned by model_fun().
#' @param strategies Strategy identifiers.
#' @param cells Expected cell keys (e.g. "00", "01", "10", "11").
#' @return Logical scalar.
psa_extended_outcomes_complete <- function(results, strategies, cells) {
  if (!psa_outcomes_complete(results, strategies)) return(FALSE)
  cell_results <- attr(results, "cells")
  if (is.null(cell_results) || !is.array(cell_results$values)) return(FALSE)
  v_dim <- dimnames(cell_results$values)
  if (!identical(v_dim[[1]], cells) || !identical(v_dim[[2]], CELL_ARMS) ||
      !identical(v_dim[[3]], CELL_OUTCOMES)) return(FALSE)
  v_mass <- cell_results$mass
  if (!identical(names(v_mass), cells) || any(!is.finite(v_mass)) || any(v_mass < 0) ||
      abs(sum(v_mass) - 1) > 1e-10) return(FALSE)
  v_observed <- cell_results$counts[cells] > 0
  all(is.finite(cell_results$values[v_observed, , ]))
}

#' Run PSA analysis across all simulations
#'
#' Every draw must return one finite cost, QALY and discounted and
#' undiscounted life-year value per strategy (the result columns of
#' model_fun(); issue #197). With extended outcomes (issue #200) every draw
#' also returns the outcomes of the four joint biomarker cells under both arms
#' (model_fun(return_cells = TRUE)). A draw without them is a failed draw and
#' is replaced like any other.
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
#' @param extended_outcomes Collect the cell-level outcomes. NULL (default)
#'   uses them whenever model_fun() has a `return_cells` argument, as the
#'   production model does; test stubs without it run without cells (they
#'   must still return the four outcome columns).
#' @return List with cost, effect and life-year (`ly_disc`, `ly_undisc`)
#'   matrices (retained draws x strategies), retained iteration indices, and
#'   fallback/drop diagnostics; with extended outcomes also `cells` (list of
#'   `values`, draw x cell x arm x outcome; `mass`, draw x cell; `violations`,
#'   draw x cell x arm; `cells`, the cell keys).
run_psa_analysis <- function(psa_params, l_params_base, param_distributions,
                             strategies, time_horizon, cl, n_sim,
                             fallback_threshold = 0.02,
                             extended_outcomes = NULL) {
  if (length(fallback_threshold) != 1L || !is.finite(fallback_threshold) ||
      fallback_threshold < 0 || fallback_threshold > 1) {
    stop("fallback_threshold must be a single finite value between 0 and 1")
  }
  if (is.null(extended_outcomes)) {
    extended_outcomes <- "return_cells" %in% names(formals(model_fun))
  }

  # Create empty matrices to store results
  m_cost <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                        dimnames = list(NULL, strategies))
  m_effect <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                          dimnames = list(NULL, strategies))
  # Life-years of the same draws (issue #197), filled together with cost and
  # effect so a replaced draw carries the outcomes of the model actually used.
  m_ly_disc <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                      dimnames = list(NULL, strategies))
  m_ly_undisc <- matrix(NA, nrow = n_sim, ncol = length(strategies),
                        dimnames = list(NULL, strategies))

  # Cell-level outcomes (#200), filled per draw.
  v_cells <- NULL
  if (isTRUE(extended_outcomes)) {
    if (is.null(l_params_base$prediction_population)) {
      stop("Extended PSA outcomes need l_params_base$prediction_population.")
    }
    v_cells <- rownames(joint_biomarker_population(l_params_base$prediction_population)$cells)
    a_cell_values <- array(NA_real_, dim = c(n_sim, length(v_cells), length(CELL_ARMS),
                                             length(CELL_OUTCOMES)),
      dimnames = list(draw = NULL, cell = v_cells, arm = CELL_ARMS, outcome = CELL_OUTCOMES))
    m_cell_mass <- matrix(NA_real_, n_sim, length(v_cells), dimnames = list(NULL, v_cells))
    a_cell_violations <- array(NA_integer_, dim = c(n_sim, length(v_cells), length(CELL_ARMS)),
      dimnames = list(draw = NULL, cell = v_cells, arm = CELL_ARMS))
  }

  # A draw is usable when it did not fall back and returned every outcome
  # column (and, with extended outcomes, the cells).
  draw_complete <- function(df_sim_results) {
    !isTRUE(attr(df_sim_results, "fallback_used")) &&
      psa_outcomes_complete(df_sim_results, strategies) &&
      (!isTRUE(extended_outcomes) ||
         psa_extended_outcomes_complete(df_sim_results, strategies, v_cells))
  }
  # Store a complete draw's outcomes in row i; the life-years come from the
  # result columns LY_disc and LY_undisc (issue #197).
  store_draw <- function(i, df_sim_results) {
    v_strategy_rows <- match(strategies, df_sim_results$Strategy)
    m_cost[i, ] <<- df_sim_results$Cost[v_strategy_rows]
    m_effect[i, ] <<- df_sim_results$Effect[v_strategy_rows]
    m_ly_disc[i, ] <<- df_sim_results$LY_disc[v_strategy_rows]
    m_ly_undisc[i, ] <<- df_sim_results$LY_undisc[v_strategy_rows]
    if (isTRUE(extended_outcomes)) {
      cell_results <- attr(df_sim_results, "cells")
      a_cell_values[i, , , ] <<- cell_results$values
      m_cell_mass[i, ] <<- cell_results$mass[v_cells]
      a_cell_violations[i, , ] <<- cell_results$violations
    }
  }
  # model_fun() arguments common to the initial pass and the replacements.
  run_draw <- function(sim_params, model_index) {
    v_args <- list(params = sim_params, time_horizon = time_horizon, cl = cl,
                   determpsa = "psa", return_traces = FALSE, sim_idx = model_index)
    if (isTRUE(extended_outcomes)) v_args$return_cells <- TRUE
    do.call(model_fun, v_args)
  }

  # Track PFS > OS constraint violations (Issue #76)
  # Instead of printing thousands of warnings, we aggregate and summarize
  pfs_os_violations <- pfs_os_violation_groups()

  # Track iterations that used fallback to base case (Issue #79)
  # These will be replaced with additional successful simulations
  v_fallback_iterations <- integer(0)

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
    
    # Use withCallingHandlers to capture warnings, tryCatch for errors
    # This allows us to aggregate PFS > OS warnings instead of printing thousands
    tryCatch({
      withCallingHandlers({
        # sim_idx = i passes the resampled model index (the critical fix of #79).
        df_sim_results <- run_draw(sim_params, i)
      }, warning = function(w) {
        # Route structured PFS > OS warnings without parsing their message text.
        if (inherits(w, "survival_ordering_warning")) {
          violation_key <- pfs_os_violation_key(w)
          if (!violation_key %in% names(pfs_os_violations)) {
            stop("Unknown survival-ordering warning subgroup: ", violation_key)
          }
          pfs_os_violations[[violation_key]] <<-
            c(pfs_os_violations[[violation_key]], i)
          invokeRestart("muffleWarning")
        }
        # Let other warnings through
      })

      # Validate immediately, before the threshold and replacement phase (#171).
      if (!draw_complete(df_sim_results)) {
        v_fallback_iterations <- c(v_fallback_iterations, i)
      } else {
        store_draw(i, df_sim_results)
      }
    }, error = function(e) {
      cat("Error in simulation", i, ":", conditionMessage(e), "\n")
      # Track this as a fallback iteration (Issue #79) - will be replaced
      v_fallback_iterations <<- c(v_fallback_iterations, i)
      # Temporarily store NAs - will be replaced in replacement phase
      for (strat in strategies) {
        m_cost[i, strat] <<- NA
        m_effect[i, strat] <<- NA
        m_ly_disc[i, strat] <<- NA
        m_ly_undisc[i, strat] <<- NA
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
  v_fallback_iterations <- unique(v_fallback_iterations)
  fallback_count <- length(v_fallback_iterations)
  fallback_rate <- fallback_count / n_sim
  cat(sprintf(
    "PSA initial-pass fallback rate: %d/%d (%.2f%%; maximum permitted %.2f%%)\n",
    fallback_count, n_sim, 100 * fallback_rate, 100 * fallback_threshold
  ))

  if (fallback_rate > fallback_threshold) {
    stop(sprintf(
      paste0(
        "PSA initial-pass fallback rate %.2f%% (%d/%d) exceeds the ",
        "permitted %.2f%% threshold; investigate prediction failures and invalid outcomes."
      ),
      100 * fallback_rate, fallback_count, n_sim, 100 * fallback_threshold
    ))
  }

  # =========================================================================
  # REPLACE FALLBACK ITERATIONS (Issue #79)
  # =========================================================================
  # Instead of mixing base case fallbacks with resampled model results,
  # run additional simulations to replace failed iterations

  v_replaced_iterations <- integer(0)
  # Cached survival model actually used by each draw (issue #156). A replaced
  # draw keeps its economic-parameter row but its outcomes come from another
  # cached model, so downstream code that pairs outcomes with model-derived
  # quantities (the interaction coefficients) must index by model_idx.
  v_model_idx <- seq_len(n_sim)

  if (length(v_fallback_iterations) > 0) {
    cat(sprintf("\n--- Replacing %d fallback iterations (Issue #79) ---\n",
                length(v_fallback_iterations)))

    replaced_count <- 0
    total_attempts <- 0L

    for (failed_i in v_fallback_iterations) {
      success <- FALSE

      # There are no models beyond the initial n_sim cached draws. Re-pair the
      # failed economic-parameter draw with up to 10 other cached survival
      # models, each at most once. Start after its original index to avoid
      # favoring low indices and bound the work for an unrecoverable draw.
      v_replacement_candidates <- psa_replacement_candidates(failed_i, n_sim)
      max_attempts <- length(v_replacement_candidates)

      # Each failed draw receives its own complete candidate scan. A shared
      # counter would let early failures exhaust the budget for later failures.
      attempt <- 0L

      while (!success && attempt < max_attempts) {
        attempt <- attempt + 1
        total_attempts <- total_attempts + 1L
        actual_sim_idx <- v_replacement_candidates[attempt]

        # Create parameter set (use same PSA samples as original iteration)
        sim_params <- l_params_base
        for (param_name in names(param_distributions)) {
          sim_params[[param_name]] <- psa_params[[param_name]][failed_i]
        }

        # Try to run with a different resampled model
        tryCatch({
          suppressWarnings({
            df_sim_results <- run_draw(sim_params, actual_sim_idx)
          })

          # A replacement is valid only when it did not use fallback and
          # returned exactly one complete result for every strategy (and,
          # with extended outcomes, complete life-years and cell outcomes).
          if (draw_complete(df_sim_results)) {
            # Replace the failed iteration's results
            store_draw(failed_i, df_sim_results)
            success <- TRUE
            replaced_count <- replaced_count + 1
            v_replaced_iterations <- c(v_replaced_iterations, failed_i)
            v_model_idx[failed_i] <- actual_sim_idx
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
                replaced_count, length(v_fallback_iterations)))
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
  v_unreplaced_iterations <- setdiff(v_fallback_iterations, v_replaced_iterations)
  v_incomplete_iterations <- which(
    !apply(m_cost, 1L, function(x) all(is.finite(x))) |
      !apply(m_effect, 1L, function(x) all(is.finite(x))) |
      !apply(m_ly_disc, 1L, function(x) all(is.finite(x))) |
      !apply(m_ly_undisc, 1L, function(x) all(is.finite(x)))
  )
  v_dropped_iterations <- sort(unique(c(
    v_unreplaced_iterations,
    v_incomplete_iterations
  )))
  dropped_count <- length(v_dropped_iterations)
  v_retained_iterations <- setdiff(seq_len(n_sim), v_dropped_iterations)

  if (dropped_count > 0L) {
    warning(
      "Dropping ", dropped_count,
      " PSA iteration(s) that remained invalid after replacement; ",
      length(v_retained_iterations), " simulation(s) remain.",
      call. = FALSE
    )
    m_cost <- m_cost[v_retained_iterations, , drop = FALSE]
    m_effect <- m_effect[v_retained_iterations, , drop = FALSE]
    m_ly_disc <- m_ly_disc[v_retained_iterations, , drop = FALSE]
    m_ly_undisc <- m_ly_undisc[v_retained_iterations, , drop = FALSE]
  }
  # The cell arrays were not subset above; keep the retained draws only.
  cells <- NULL
  if (isTRUE(extended_outcomes)) {
    cells <- list(values = a_cell_values[v_retained_iterations, , , , drop = FALSE],
                  mass = m_cell_mass[v_retained_iterations, , drop = FALSE],
                  violations = a_cell_violations[v_retained_iterations, , , drop = FALSE],
                  cells = v_cells)
  }

  if (length(v_retained_iterations) == 0L) {
    stop("No valid PSA iterations remain after dropping failed draws.")
  }

  # Print summary of PFS > OS constraint violations (Issue #76)
  total_violations <- sum(sapply(pfs_os_violations, length))
  if (total_violations > 0) {
    # Count unique iterations with any violation
    v_all_violation_iters <- unique(unlist(pfs_os_violations))
    n_iters_with_violations <- length(v_all_violation_iters)

    cat("\n--- PFS > OS Constraint Summary (Issue #76) ---\n")
    cat(sprintf("Violations occurred in %d of %d iterations (%.1f%%)\n",
                n_iters_with_violations, n_sim, 100 * n_iters_with_violations / n_sim))
    cat("Breakdown by subgroup:\n")

    # Only print subgroups that had violations
    v_subgroup_names <- pfs_os_violation_labels()

    for (subgroup in names(pfs_os_violations)) {
      n_violations <- length(pfs_os_violations[[subgroup]])
      if (n_violations > 0) {
        cat(sprintf("%s: %d iterations (%.1f%%)\n",
                    v_subgroup_names[subgroup], n_violations, 100 * n_violations / n_sim))
      }
    }
    cat("Constraint enforced: PFS capped at OS to ensure valid state occupancy.\n")
    cat("------------------------------------------------\n\n")
  }

  # Return the cost, effect and life-year matrices separately
  return(list(
    cost = m_cost,
    effect = m_effect,
    ly_disc = m_ly_disc,
    ly_undisc = m_ly_undisc,
    fallback_count = fallback_count,
    dropped_count = dropped_count,
    dropped_iterations = v_dropped_iterations,
    retained_iterations = v_retained_iterations,
    model_idx = v_model_idx[v_retained_iterations],
    n_sim = length(v_retained_iterations),
    failed_draw_policy = psa_failed_draw_policy(),
    replacement_model_policy = "cyclic_next_10_other_cached_models_v1",
    fallback_rate = fallback_rate,
    fallback_threshold = fallback_threshold,
    pfs_os_violations = pfs_os_violations,
    extended_outcomes = isTRUE(extended_outcomes),
    cells = cells
  ))
}

# Schema of the cell-level component of a PSA object (issue #200).
PSA_CELLS_SCHEMA <- "psa_cells_v1"

#' Assemble the cell-level component stored in a PSA object
#'
#' @param psa_results Output of run_psa_analysis() with extended outcomes.
#' @param l_params_base Base-case parameter list of the PSA (prediction
#'   population, nivolumab price, diagnostic-test prices).
#' @param df_psa_params Retained PSA parameter rows (row-aligned with the
#'   outcomes); a sampled `c_drug_nivo` column is used per draw.
#' @return List with `schema`, `cells` (keys), `cell_definition` (biomarker
#'   values per cell), `biomarkers`, `arms`, `outcomes`, `values` (draw x cell
#'   x arm x outcome), `weights` (draw x cell joint cell masses),
#'   `violations` (draw x cell x arm), `nivolumab_price` (per draw: the price
#'   at which the cell costs were evaluated) and `test_costs` (per biomarker).
build_psa_cell_component <- function(psa_results, l_params_base, df_psa_params) {
  v_biomarkers <- get_biomarkers()
  joint <- joint_biomarker_population(l_params_base$prediction_population, v_biomarkers)
  n_draws <- nrow(psa_results$cells$mass)
  v_price <- if ("c_drug_nivo" %in% names(df_psa_params)) df_psa_params$c_drug_nivo else
    rep(l_params_base$c_drug_nivo, n_draws)
  list(schema = PSA_CELLS_SCHEMA, cells = psa_results$cells$cells,
       cell_definition = joint$cells, biomarkers = v_biomarkers, arms = CELL_ARMS,
       outcomes = CELL_OUTCOMES, values = psa_results$cells$values,
       weights = psa_results$cells$mass, violations = psa_results$cells$violations,
       nivolumab_price = as.numeric(v_price),
       test_costs = vapply(v_biomarkers, function(b) l_params_base[[biomarker_cost_key(b)[[1]]]],
                           numeric(1)))
}


#' Generate PSA draws, run the model, and create a row-aligned PSA object
#'
#' @param l_params_base Base-case parameter list
#' @param param_distributions Named parameter distributions
#' @param strategies Strategy identifiers
#' @param time_horizon Model time horizon
#' @param cl Cycle length
#' @param n_sim Requested simulations
#' @param seed RNG seed
#' @param currency Currency label passed to `dampack::make_psa_obj()`
#' @param additional_params Optional draw-aligned columns, such as interaction coefficients
#' @param fallback_threshold Maximum initial fallback rate
#' @param extended_outcomes See run_psa_analysis().
#' @return List containing `psa_obj` (with the life-year draws as
#'   `psa_obj$ly_disc` and `psa_obj$ly_undisc`, n_sim x strategies data frames
#'   aligned with dampack's `cost` and `effectiveness`; issue #197), retained
#'   `psa_params`, and `psa_results`. With extended outcomes the PSA object
#'   also holds `cells` (build_psa_cell_component(); issue #200).
build_psa_obj <- function(l_params_base, param_distributions, strategies,
                          time_horizon, cl, n_sim, seed = 123L,
                          currency = "€", additional_params = NULL,
                          fallback_threshold = 0.02, extended_outcomes = NULL) {
  df_psa_params <- generate_psa_samples(param_distributions, n_sim, seed)
  if (!is.null(additional_params) && nrow(additional_params) != n_sim) {
    stop("additional_params must contain exactly n_sim rows")
  }

  psa_results <- run_psa_analysis(
    df_psa_params, l_params_base, param_distributions, strategies,
    time_horizon, cl, n_sim, fallback_threshold, extended_outcomes
  )
  df_psa_params <- df_psa_params[psa_results$retained_iterations, , drop = FALSE]
  # The cached survival model behind each retained row. Replaced draws keep
  # their economic-parameter row but were run with another cached model, so
  # model-derived columns (interaction coefficients) are aligned by model_idx,
  # not by the draw number (issue #156).
  df_psa_params$model_idx <- psa_results$model_idx
  if (!is.null(additional_params)) {
    df_psa_params <- cbind(
      df_psa_params,
      additional_params[psa_results$model_idx, , drop = FALSE]
    )
  }
  rownames(df_psa_params) <- NULL
  attr(df_psa_params, "seed") <- seed

  psa_obj <- dampack::make_psa_obj(
    cost = as.data.frame(psa_results$cost),
    effect = as.data.frame(psa_results$effect),
    strategies = strategies,
    currency = currency
  )
  # Life-years of the same retained draws (issue #197), stored like the dampack
  # cost and effectiveness tables: one column per strategy, one row per draw.
  psa_obj$ly_disc <- as.data.frame(psa_results$ly_disc)
  psa_obj$ly_undisc <- as.data.frame(psa_results$ly_undisc)
  v_metadata <- c(
    requested_n_sim = "requested_n_sim", fallback_count = "fallback_count",
    model_idx = "model_idx",
    dropped_count = "dropped_count", dropped_iterations = "dropped_iterations",
    failed_draw_policy = "failed_draw_policy",
    replacement_model_policy = "replacement_model_policy"
  )
  for (field in names(v_metadata)) {
    psa_obj[[field]] <- if (field == "requested_n_sim") n_sim else psa_results[[v_metadata[[field]]]]
  }
  if (isTRUE(psa_results$extended_outcomes)) {
    psa_obj$cells <- build_psa_cell_component(psa_results, l_params_base, df_psa_params)
  }

  list(psa_obj = psa_obj, psa_params = df_psa_params, psa_results = psa_results)
}
