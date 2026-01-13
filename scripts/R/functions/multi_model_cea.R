# ===============================================================================
# MULTI-MODEL CEA HELPER FUNCTIONS
# ===============================================================================
# Functions for running and comparing cost-effectiveness analysis across multiple
# model structures (Model A/B/C - joint/focused/separate).
#
# Used by CEA.qmd to generate multi-model comparison tables and figures.
#
# Prerequisites: Before sourcing this file, ensure:
#   1. Global variables are set (source 02_setup_and_global_variables.R)
#   2. Biomarker strategies are defined (source 03_biomarker_strategies.R)
#   3. model_fun.R and calculate_outcomes.R are loaded
# ===============================================================================

# Load required packages (don't reload here - CEA.qmd handles this)
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dplyr, ggplot2, scales, dampack)

# ===============================================================================
# MODEL CONFIGURATION
# ===============================================================================

#' Get model configuration definitions
#'
#' @return Named list of model configurations
get_model_configs <- function() {
  list(
    Model_A = list(
      structure = 0,
      label = "Model A: Joint",
      short_label = "A",
      description = "All biomarkers + all treatment interactions in one model"
    ),
    Model_B = list(
      structure = 1,
      label = "Model B: Focused",
      short_label = "B",
      description = "All biomarkers as main effects + one interaction per model"
    ),
    Model_C = list(
      structure = 2,
      label = "Model C: Separate",
      short_label = "C",
      description = "One biomarker + its interaction only per model"
    )
  )
}

# ===============================================================================
# BASE CASE ANALYSIS FUNCTIONS
# ===============================================================================

#' Run base case analysis for a specific model structure
#'
#' @param model_structure Integer (0, 1, or 2) indicating model structure
#' @param verbose Logical, print progress messages
#' @return List containing base_results and icer_obj
run_basecase_for_structure <- function(model_structure, verbose = TRUE) {

  if (verbose) {
    model_label <- c("joint", "focused", "separate")[model_structure + 1]
    cat("Running base case for MODEL_STRUCTURE =", model_structure,
        "(", model_label, ")...\n")
  }

  # Set MODEL_STRUCTURE globally
  MODEL_STRUCTURE <<- model_structure

  # Re-source survival analysis to update predictions
  suppressMessages(suppressWarnings(
    source(here::here("scripts/R/analysis/06_parametric_survival_analysis.R"))
  ))

  # Re-source input parameters to update l_params_base
  suppressMessages(suppressWarnings(
    source(here::here("scripts/R/analysis/07_basecase_input_parameters.R"))
  ))

  # Run base case model
  base_results <- model_fun(l_params_base)

  # Calculate ICERs
  icer_obj <- dampack::calculate_icers(
    cost = base_results$Cost,
    effect = base_results$Effect,
    strategies = base_results$Strategy
  )

  if (verbose) {
    cat("  Completed. Strategies:", paste(base_results$Strategy, collapse = ", "), "\n")
  }

  return(list(
    base_results = base_results,
    icer_obj = icer_obj,
    model_structure = model_structure
  ))
}

#' Run base case analysis for all model structures
#'
#' @param verbose Logical, print progress messages
#' @return Named list of results for each model structure
run_all_basecase_analyses <- function(verbose = TRUE) {

  model_configs <- get_model_configs()
  all_results <- list()

  for (model_name in names(model_configs)) {
    config <- model_configs[[model_name]]

    if (verbose) {
      cat("\n--- Processing", config$label, "---\n")
    }

    all_results[[model_name]] <- run_basecase_for_structure(
      model_structure = config$structure,
      verbose = verbose
    )
    all_results[[model_name]]$config <- config
  }

  return(all_results)
}

# ===============================================================================
# PSA CACHE LOADING FUNCTIONS
# ===============================================================================

#' Load PSA cache for a specific model structure
#'
#' @param model_structure Integer (0, 1, or 2) indicating model structure
#' @param n_samples Number of samples (for cache filename)
#' @param use_both_models USE_BOTH_MODELS setting (for cache filename)
#' @return PSA object or NULL if not found
load_psa_cache_for_structure <- function(model_structure,
                                          n_samples = 5000,
                                          use_both_models = 0) {

  model_label <- c("joint", "focused", "separate")[model_structure + 1]

  cache_file <- here::here("data", "tidy",
                           paste0("psa_obj_", model_label, ".rds"))

  if (file.exists(cache_file)) {
    psa_obj <- readRDS(cache_file)
    cat("Loaded PSA cache for", model_label, "model\n")
    return(psa_obj)
  } else {
    warning("PSA cache not found: ", cache_file)
    return(NULL)
  }
}

#' Load PSA caches for all model structures
#'
#' @param n_samples Number of samples (for cache filename)
#' @param use_both_models USE_BOTH_MODELS setting (for cache filename)
#' @return Named list of PSA objects
load_all_psa_caches <- function(n_samples = 5000, use_both_models = 0) {

  model_configs <- get_model_configs()
  all_psa <- list()

  for (model_name in names(model_configs)) {
    config <- model_configs[[model_name]]

    psa_obj <- load_psa_cache_for_structure(
      model_structure = config$structure,
      n_samples = n_samples,
      use_both_models = use_both_models
    )

    if (!is.null(psa_obj)) {
      all_psa[[model_name]] <- psa_obj
    }
  }

  return(all_psa)
}

# ===============================================================================
# COMPARISON TABLE FUNCTIONS
# ===============================================================================

#' Format multi-model comparison table for base case results
#'
#' @param all_results List of base case results from run_all_basecase_analyses()
#' @return Data frame with side-by-side comparison
format_multimodel_comparison_table <- function(all_results) {

  model_configs <- get_model_configs()

  # Get strategy names from first result
  strategies <- all_results[[1]]$base_results$Strategy

  # Initialize comparison data frame
  comparison_df <- data.frame(Strategy = strategies)

  for (model_name in names(all_results)) {
    result <- all_results[[model_name]]
    short_label <- model_configs[[model_name]]$short_label

    # Add cost and effect columns
    comparison_df[[paste0("Cost_", short_label)]] <- result$base_results$Cost
    comparison_df[[paste0("QALY_", short_label)]] <- result$base_results$Effect

    # Add ICER from icer_obj
    icer_df <- as.data.frame(result$icer_obj)
    comparison_df[[paste0("ICER_", short_label)]] <- icer_df$ICER
  }

  return(comparison_df)
}

#' Format comparison table for display with kable
#'
#' @param comparison_df Data frame from format_multimodel_comparison_table()
#' @return Formatted data frame ready for kable
format_comparison_for_display <- function(comparison_df) {

  display_df <- comparison_df

  # Format cost columns
  cost_cols <- grep("^Cost_", names(display_df), value = TRUE)
  for (col in cost_cols) {
    display_df[[col]] <- scales::dollar(display_df[[col]], accuracy = 1)
  }

  # Format QALY columns
  qaly_cols <- grep("^QALY_", names(display_df), value = TRUE)
  for (col in qaly_cols) {
    display_df[[col]] <- sprintf("%.3f", display_df[[col]])
  }

  # Format ICER columns
  icer_cols <- grep("^ICER_", names(display_df), value = TRUE)
  for (col in icer_cols) {
    display_df[[col]] <- ifelse(
      is.na(display_df[[col]]),
      "Reference",
      scales::dollar(display_df[[col]], accuracy = 1)
    )
  }

  return(display_df)
}

#' Create optimal strategy comparison table
#'
#' @param all_results List of base case results
#' @param wtp Willingness-to-pay threshold
#' @return Data frame showing optimal strategy per model
create_optimal_strategy_table <- function(all_results, wtp = 51000) {

  model_configs <- get_model_configs()

  optimal_df <- data.frame(
    Model = character(),
    Optimal_Strategy = character(),
    NMB = numeric(),
    stringsAsFactors = FALSE
  )

  for (model_name in names(all_results)) {
    result <- all_results[[model_name]]
    config <- model_configs[[model_name]]

    # Calculate NMB for each strategy
    nmb <- result$base_results$Effect * wtp - result$base_results$Cost

    # Find optimal strategy
    optimal_idx <- which.max(nmb)
    optimal_strategy <- result$base_results$Strategy[optimal_idx]
    optimal_nmb <- nmb[optimal_idx]

    optimal_df <- rbind(optimal_df, data.frame(
      Model = config$label,
      Optimal_Strategy = optimal_strategy,
      NMB = optimal_nmb,
      stringsAsFactors = FALSE
    ))
  }

  return(optimal_df)
}

# ===============================================================================
# PSA COMPARISON FUNCTIONS
# ===============================================================================

#' Create multi-model CEAC plot
#'
#' @param all_psa_results Named list of PSA objects
#' @param wtp_range Vector of WTP thresholds
#' @return ggplot object
create_multimodel_ceac_plot <- function(all_psa_results,
                                         wtp_range = seq(0, 100000, by = 2000)) {

  model_configs <- get_model_configs()

  # Calculate CEAC for each model
  ceac_list <- list()

  for (model_name in names(all_psa_results)) {
    psa_obj <- all_psa_results[[model_name]]
    config <- model_configs[[model_name]]

    if (is.null(psa_obj)) next

    # Calculate CEAC using dampack
    ceac_obj <- tryCatch({
      dampack::ceac(wtp = wtp_range, psa = psa_obj)
    }, error = function(e) {
      warning("CEAC calculation failed for ", model_name, ": ", e$message)
      NULL
    })

    if (is.null(ceac_obj)) next

    # Add to combined data
    for (i in seq_along(wtp_range)) {
      for (strategy in psa_obj$strategies) {
        prob <- tryCatch({
          ceac_obj$ceac[[strategy]][i]
        }, error = function(e) {
          NA_real_
        })

        if (is.null(prob) || length(prob) == 0) {
          prob <- NA_real_
        }

        ceac_list[[length(ceac_list) + 1]] <- data.frame(
          WTP = wtp_range[i],
          Strategy = strategy,
          Model = config$label,
          Probability = prob,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Combine all rows
  if (length(ceac_list) > 0) {
    ceac_data <- do.call(rbind, ceac_list)
  } else {
    # Return empty plot if no data
    ceac_data <- data.frame(
      WTP = numeric(),
      Strategy = character(),
      Model = character(),
      Probability = numeric(),
      stringsAsFactors = FALSE
    )
  }

  # Remove NA values for plotting
  ceac_data <- ceac_data[!is.na(ceac_data$Probability), ]

  # Create plot
  if (nrow(ceac_data) > 0) {
    p <- ggplot(ceac_data, aes(x = WTP, y = Probability,
                                color = Strategy, linetype = Model)) +
      geom_line(linewidth = 0.8) +
      scale_x_continuous(labels = scales::dollar_format(scale = 0.001,
                                                         suffix = "K")) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(
        x = "Willingness-to-Pay Threshold",
        y = "Probability of Being Cost-Effective",
        color = "Strategy",
        linetype = "Model Structure"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.minor = element_blank()
      )
  } else {
    # Return a placeholder plot if no data
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = "No CEAC data available", size = 5) +
      theme_void()
  }

  return(p)
}

#' Create PSA summary table across models
#'
#' @param all_psa_results Named list of PSA objects
#' @param wtp WTP threshold for probability calculation
#' @return Data frame with PSA summary
create_psa_summary_table <- function(all_psa_results, wtp = 51000) {

  model_configs <- get_model_configs()

  summary_list <- list()

  for (model_name in names(all_psa_results)) {
    psa_obj <- all_psa_results[[model_name]]
    config <- model_configs[[model_name]]

    if (is.null(psa_obj)) next

    # Calculate CEAC at specified WTP
    # Use a small range around WTP to ensure proper CEAC calculation
    ceac_obj <- dampack::ceac(wtp = c(wtp - 1, wtp, wtp + 1), psa = psa_obj)

    # Get mean costs and effects
    mean_costs <- colMeans(psa_obj$cost)
    mean_effects <- colMeans(psa_obj$effect)

    for (i in seq_along(psa_obj$strategies)) {
      strategy <- psa_obj$strategies[i]

      # Get probability at the middle WTP value (index 2)
      prob_ce <- tryCatch({
        ceac_obj$ceac[[strategy]][2]
      }, error = function(e) {
        NA_real_
      })

      # If that failed, try index 1
      if (is.null(prob_ce) || length(prob_ce) == 0) {
        prob_ce <- tryCatch({
          ceac_obj$ceac[[strategy]][1]
        }, error = function(e) {
          NA_real_
        })
      }

      # Final fallback
      if (is.null(prob_ce) || length(prob_ce) == 0) {
        prob_ce <- NA_real_
      }

      summary_list[[length(summary_list) + 1]] <- data.frame(
        Model = config$short_label,
        Strategy = strategy,
        Mean_Cost = mean_costs[i],
        Mean_QALY = mean_effects[i],
        Prob_CE = prob_ce,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all rows
  if (length(summary_list) > 0) {
    summary_df <- do.call(rbind, summary_list)
  } else {
    summary_df <- data.frame(
      Model = character(),
      Strategy = character(),
      Mean_Cost = numeric(),
      Mean_QALY = numeric(),
      Prob_CE = numeric(),
      stringsAsFactors = FALSE
    )
  }

  return(summary_df)
}

message("Multi-model CEA helper functions loaded successfully.")
