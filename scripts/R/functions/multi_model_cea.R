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
# Model configurations are defined in scripts/R/functions/model_configs.R
# which provides the single source of truth for all model structures.
# The get_model_configs() function is imported from that file.

# Source model configurations if not already loaded
if (!exists("get_model_configs") || !exists("get_strategies")) {
  model_configs_path <- here::here("scripts/R/functions/model_configs.R")
  if (file.exists(model_configs_path)) {
    source(model_configs_path)
  } else {
    stop("model_configs.R not found. Please ensure scripts/R/functions/model_configs.R exists.")
  }
}

# Note: get_model_configs() is now provided by model_configs.R
# All models include all three biomarker strategies (CRP, TLR, TMB/BRAF)

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
#' @return Data frame with side-by-side comparison, sorted by cost (Model A)
format_multimodel_comparison_table <- function(all_results) {

  model_configs <- get_model_configs()

  # Get strategy names from first result
  first_icer_df <- as.data.frame(all_results[[1]]$icer_obj)
  strategies <- first_icer_df$Strategy

  # Initialize comparison data frame
  comparison_df <- data.frame(Strategy = strategies, stringsAsFactors = FALSE)

  for (model_name in names(all_results)) {
    result <- all_results[[model_name]]
    short_label <- model_configs[[model_name]]$short_label

    # Get icer_df which contains all the data we need
    icer_df <- as.data.frame(result$icer_obj)

    # Match strategies by name to handle potential ordering differences
    match_idx <- match(comparison_df$Strategy, icer_df$Strategy)

    # Add cost and effect columns (from icer_df to ensure consistency)
    comparison_df[[paste0("Cost_", short_label)]] <- icer_df$Cost[match_idx]
    comparison_df[[paste0("QALY_", short_label)]] <- icer_df$Effect[match_idx]

    # Add ICER with proper handling of dominated strategies
    # dampack marks dominated strategies with Status "D" or "ED"
    icer_values <- icer_df$ICER[match_idx]
    status_values <- icer_df$Status[match_idx]

    # For dominated strategies, show "Dominated" instead of ICER
    icer_display <- ifelse(
      status_values %in% c("D", "ED"),
      NA,  # Will be converted to "Dominated" in format_comparison_for_display
      icer_values
    )

    comparison_df[[paste0("ICER_", short_label)]] <- icer_display
    comparison_df[[paste0("Status_", short_label)]] <- status_values
  }

  # Sort by Cost from first model (Model A) in ascending order
  comparison_df <- comparison_df[order(comparison_df$Cost_A), ]
  rownames(comparison_df) <- NULL

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

  # Format ICER columns with proper handling of reference and dominated strategies
  icer_cols <- grep("^ICER_", names(display_df), value = TRUE)
  status_cols <- grep("^Status_", names(display_df), value = TRUE)

  for (i in seq_along(icer_cols)) {
    icer_col <- icer_cols[i]
    # Extract model label (e.g., "A" from "ICER_A")
    model_label <- gsub("^ICER_", "", icer_col)
    status_col <- paste0("Status_", model_label)

    if (status_col %in% names(display_df)) {
      display_df[[icer_col]] <- ifelse(
        display_df[[status_col]] %in% c("D", "ED"),
        "Dominated",
        ifelse(
          is.na(display_df[[icer_col]]),
          "Reference",
          scales::dollar(display_df[[icer_col]], accuracy = 1)
        )
      )
    } else {
      # Fallback if no status column
      display_df[[icer_col]] <- ifelse(
        is.na(display_df[[icer_col]]),
        "Reference",
        scales::dollar(display_df[[icer_col]], accuracy = 1)
      )
    }
  }

  # Remove Status columns from display
  display_df <- display_df[, !grepl("^Status_", names(display_df)), drop = FALSE]

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

    # dampack::ceac() returns a data frame with columns:
    # WTP, Strategy, Proportion, On_Frontier
    # Add Model column and append to list
    ceac_obj$Model <- config$label
    ceac_list[[length(ceac_list) + 1]] <- ceac_obj
  }

  # Combine all rows
  if (length(ceac_list) > 0) {
    ceac_data <- do.call(rbind, ceac_list)
    # Rename Proportion to Probability for consistency with plot labels
    names(ceac_data)[names(ceac_data) == "Proportion"] <- "Probability"
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

  # Create plot with faceted panels (one per model structure)
  if (nrow(ceac_data) > 0) {
    p <- ggplot(ceac_data, aes(x = WTP, y = Probability, color = Strategy)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~Model, ncol = 1) +
      geom_vline(xintercept = 51000, linetype = "dotted", color = "gray40", linewidth = 0.5) +
      scale_x_continuous(labels = scales::dollar_format(prefix = "EUR", scale = 0.001, suffix = "K")) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
      labs(
        x = "Willingness-to-Pay Threshold (EUR/QALY)",
        y = "Probability of Being Cost-Effective",
        color = "Strategy"
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 11, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(nrow = 1))
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

      # dampack::ceac() returns data frame with columns: WTP, Strategy, Proportion, On_Frontier
      # Filter for the target WTP and strategy to get the probability
      prob_ce <- tryCatch({
        ceac_row <- ceac_obj[ceac_obj$WTP == wtp & ceac_obj$Strategy == strategy, ]
        if (nrow(ceac_row) > 0) {
          ceac_row$Proportion[1]
        } else {
          NA_real_
        }
      }, error = function(e) {
        NA_real_
      })

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

# ===============================================================================
# EVPPI HELPER FUNCTIONS
# ===============================================================================

#' Load scenario EVPPI results
#'
#' @param results_file Path to the .rds file containing scenario EVPPI results
#' @return List containing all_scenario_results, all_model_scenario_results,
#'         evppi_all_scenarios, and scenarios
load_scenario_evppi_results <- function(results_file = NULL) {

  if (is.null(results_file)) {
    results_file <- here::here("data/tidy/scenario_evppi_results.rds")
  }

  if (!file.exists(results_file)) {
    warning("EVPPI results file not found: ", results_file)
    return(NULL)
  }

  # Load results from .rds file (stored as list)
  result <- readRDS(results_file)

  return(result)
}

#' Check if EVPPI results have multi-model data
#'
#' @param evppi_data Data frame of EVPPI results (evppi_all_scenarios)
#' @return Logical indicating whether model_structure column exists
has_multimodel_evppi <- function(evppi_data) {
  if (is.null(evppi_data) || nrow(evppi_data) == 0) {
    return(FALSE)
  }
  return("model_structure" %in% names(evppi_data))
}

#' Get unique model structures from EVPPI results
#'
#' @param evppi_data Data frame of EVPPI results
#' @return Vector of unique model labels, or "joint" if single model
get_evppi_model_labels <- function(evppi_data) {
  if (!has_multimodel_evppi(evppi_data)) {
    return("joint")
  }
  return(unique(evppi_data$model_label))
}

#' Create multi-model EVPPI comparison table
#'
#' @param evppi_data Data frame of EVPPI results with model_structure column
#' @param scenario_id Scenario ID to filter (e.g., "base")
#' @param top_n Number of top parameters to include
#' @return Data frame with side-by-side EVPPI comparison across models
create_multimodel_evppi_table <- function(evppi_data, scenario_id = "base", top_n = 10) {

  model_configs <- get_model_configs()

  if (!has_multimodel_evppi(evppi_data)) {
    # Single model - return simple table
    scenario_data <- evppi_data %>%
      filter(scenario_id == !!scenario_id) %>%
      arrange(desc(evppi)) %>%
      slice(1:top_n)

    return(scenario_data %>%
             select(Parameter = parameter, EVPPI = evppi, `% of EVPI` = evppi_percent_of_evpi))
  }

  # Multi-model comparison
  comparison_list <- list()

  for (model_name in names(model_configs)) {
    config <- model_configs[[model_name]]
    short_label <- config$short_label
    model_label <- c("joint", "focused", "separate")[config$structure + 1]

    model_data <- evppi_data %>%
      filter(scenario_id == !!scenario_id,
             model_label == !!model_label) %>%
      arrange(desc(evppi)) %>%
      slice(1:top_n) %>%
      select(parameter, evppi, evppi_percent_of_evpi)

    if (nrow(model_data) > 0) {
      comparison_list[[model_name]] <- model_data %>%
        rename(
          !!paste0("EVPPI_", short_label) := evppi,
          !!paste0("Pct_", short_label) := evppi_percent_of_evpi
        )
    }
  }

  # Merge all models by parameter
  if (length(comparison_list) == 0) {
    return(data.frame())
  }

  result <- comparison_list[[1]]
  for (i in 2:length(comparison_list)) {
    result <- full_join(result, comparison_list[[i]], by = "parameter")
  }

  # Rename parameter column and sort
  result <- result %>%
    rename(Parameter = parameter) %>%
    arrange(desc(rowMeans(select(., starts_with("EVPPI")), na.rm = TRUE)))

  return(result)
}

#' Create multi-model EVPI comparison table
#'
#' @param evppi_data Data frame of EVPPI results with model_structure column
#' @param scenario_id Scenario ID to filter (e.g., "base"), or NULL for all scenarios
#' @return Data frame with EVPI values across models
create_multimodel_evpi_table <- function(evppi_data, scenario_id = NULL) {

  if (!has_multimodel_evppi(evppi_data)) {
    # Single model
    evpi_data <- evppi_data %>%
      select(scenario_name, wtp, evpi) %>%
      distinct()

    if (!is.null(scenario_id)) {
      evpi_data <- evpi_data %>% filter(scenario_id == !!scenario_id)
    }

    return(evpi_data)
  }

  # Multi-model EVPI comparison
  evpi_data <- evppi_data %>%
    select(scenario_id, scenario_name, wtp, model_label, model_name, evpi) %>%
    distinct()

  if (!is.null(scenario_id)) {
    evpi_data <- evpi_data %>% filter(scenario_id == !!scenario_id)
  }

  # Pivot to wide format
  evpi_wide <- evpi_data %>%
    select(scenario_name, wtp, model_label, evpi) %>%
    pivot_wider(
      names_from = model_label,
      values_from = evpi,
      names_prefix = "EVPI_"
    )

  return(evpi_wide)
}

#' Format EVPPI table for display with kable
#'
#' @param evppi_table Data frame from create_multimodel_evppi_table()
#' @return Formatted data frame ready for kable
format_evppi_for_display <- function(evppi_table) {

  display_df <- evppi_table

  # Format EVPPI columns (currency)
  evppi_cols <- grep("^EVPPI", names(display_df), value = TRUE)
  for (col in evppi_cols) {
    display_df[[col]] <- ifelse(
      is.na(display_df[[col]]),
      "-",
      paste0("EUR", format(round(display_df[[col]], 2), big.mark = ","))
    )
  }

  # Format percentage columns
  pct_cols <- grep("^Pct", names(display_df), value = TRUE)
  for (col in pct_cols) {
    display_df[[col]] <- ifelse(
      is.na(display_df[[col]]),
      "-",
      paste0(round(display_df[[col]], 1), "%")
    )
  }

  return(display_df)
}

#' Create faceted EVPPI bar plot for multi-model comparison
#'
#' @param evppi_data Data frame of EVPPI results
#' @param scenario_id Scenario ID to filter
#' @param top_n Number of top parameters to show per model
#' @return ggplot object
create_multimodel_evppi_plot <- function(evppi_data, scenario_id = "base", top_n = 8) {

  # Filter to scenario
  plot_data <- evppi_data %>%
    filter(scenario_id == !!scenario_id, evppi > 0)

  if (nrow(plot_data) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No EVPPI data available for this scenario") +
             theme_void())
  }

  # Get top parameters per model
  if (has_multimodel_evppi(evppi_data)) {
    top_params <- plot_data %>%
      group_by(model_label) %>%
      arrange(desc(evppi)) %>%
      slice(1:top_n) %>%
      ungroup() %>%
      pull(parameter) %>%
      unique()

    plot_data <- plot_data %>%
      filter(parameter %in% top_params) %>%
      mutate(parameter = factor(parameter, levels = rev(top_params)))

    p <- ggplot(plot_data, aes(x = evppi, y = parameter, fill = model_name)) +
      geom_col(position = "dodge", alpha = 0.8) +
      facet_wrap(~model_name, ncol = 3) +
      scale_x_continuous(labels = scales::dollar_format(prefix = "EUR")) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        x = "EVPPI (EUR)",
        y = NULL,
        fill = "Model Structure"
      ) +
      theme_minimal() +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))

  } else {
    # Single model
    top_params <- plot_data %>%
      arrange(desc(evppi)) %>%
      slice(1:top_n) %>%
      pull(parameter)

    plot_data <- plot_data %>%
      filter(parameter %in% top_params) %>%
      mutate(parameter = factor(parameter, levels = rev(top_params)))

    p <- ggplot(plot_data, aes(x = evppi, y = parameter)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      scale_x_continuous(labels = scales::dollar_format(prefix = "EUR")) +
      labs(
        x = "EVPPI (EUR)",
        y = NULL
      ) +
      theme_minimal()
  }

  return(p)
}

message("Multi-model CEA helper functions loaded successfully.")
