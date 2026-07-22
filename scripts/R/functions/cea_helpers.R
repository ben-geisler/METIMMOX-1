# ===============================================================================
# SINGLE-MODEL CEA HELPER FUNCTIONS
# ===============================================================================
# Functions for running and summarizing the single joint economic model.
#
# Prerequisites before calling run_basecase():
#   1. Global variables are set (source 02_setup_and_global_variables.R)
#   2. Biomarker strategies are defined (source 03_biomarker_strategies.R)
#   3. Survival predictions and l_params_base are available (source 06/07)
#   4. model_fun.R and calculate_outcomes.R are loaded
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dplyr, ggplot2, scales, dampack)

# Source model configuration if needed. The configuration file defines the single
# joint economic model and the economic strategies included in CEA.
if (!exists("get_current_model_config") || !exists("get_strategies")) {
  model_configs_path <- here::here("scripts/R/functions/model_configs.R")
  if (file.exists(model_configs_path)) {
    source(model_configs_path)
  } else {
    stop("model_configs.R not found. Please ensure scripts/R/functions/model_configs.R exists.")
  }
}

# Source cache-path helpers if needed (utility label + cache file locations).
if (!exists("psa_obj_path") || !exists("resolve_util_label")) {
  cache_paths_path <- here::here("scripts/R/functions/cache_paths.R")
  if (file.exists(cache_paths_path)) {
    source(cache_paths_path)
  } else {
    stop("cache_paths.R not found. Please ensure scripts/R/functions/cache_paths.R exists.")
  }
}

#' Run base case analysis for the single joint economic model
#'
#' @param params Parameter list. Defaults to global l_params_base.
#' @param verbose Logical, print progress messages.
#' @return List containing base_results, icer_obj, and model_config.
run_basecase <- function(params = NULL, verbose = TRUE) {

  if (is.null(params)) {
    if (!exists("l_params_base", inherits = TRUE)) {
      stop("l_params_base not found. Source 06/07 before calling run_basecase().")
    }
    params <- get("l_params_base", inherits = TRUE)
  }

  if (!exists("model_fun", mode = "function", inherits = TRUE)) {
    stop("model_fun() not found. Source scripts/R/functions/model_fun.R first.")
  }

  if (verbose) {
    config <- get_current_model_config()
    cat("Running base case for:", config$label, "\n")
    cat("Strategies:", paste(get_strategies(), collapse = ", "), "\n")
  }

  base_results <- model_fun(params)

  icer_obj <- dampack::calculate_icers(
    cost = base_results$Cost,
    effect = base_results$Effect,
    strategies = base_results$Strategy
  )

  if (verbose) {
    cat("Base case complete. Strategies:",
        paste(base_results$Strategy, collapse = ", "), "\n")
  }

  list(
    base_results = base_results,
    icer_obj = icer_obj,
    model_config = get_current_model_config()
  )
}

#' Calculate pairwise ICERs against the control strategy
#'
#' @param results Data frame with Strategy, Cost, and Effect columns
#' @return Results with incremental outcomes, ICER, and pairwise status
calculate_pairwise_icers <- function(results) {
  control <- results[results$Strategy == get_control_strategy(), ]
  if (nrow(control) != 1) stop("Control strategy must appear exactly once.")
  out <- transform(
    results,
    Inc_Cost = Cost - control$Cost,
    Inc_Effect = Effect - control$Effect
  )
  out$ICER <- out$Inc_Cost / out$Inc_Effect
  out$ICER[out$Strategy == get_control_strategy()] <- NA_real_
  out$Status <- "Non-dominated"
  out$Status[out$Inc_Cost > 0 & out$Inc_Effect <= 0] <- "Dominated"
  out$Status[out$Inc_Cost < 0 & out$Inc_Effect > 0] <- "Cost-saving"
  out$Status[out$Strategy == get_control_strategy()] <- "Reference"
  out
}

#' Load the single-model PSA cache
#'
#' @param util_label Utility-source label. Defaults to global utility_source_label.
#' @param directory Directory containing PSA cache files.
#' @param verbose Logical, print progress messages.
#' @return dampack PSA object, or NULL if the cache is not available.
load_psa_cache <- function(util_label = NULL,
                           directory = cache_dir(),
                           verbose = TRUE) {

  cache_file <- psa_obj_path(util_label, directory)

  if (!file.exists(cache_file)) {
    warning("PSA cache not found: ", cache_file)
    return(NULL)
  }

  psa_obj <- readRDS(cache_file)

  if (verbose) {
    cat("Loaded PSA cache:", cache_file, "\n")
    if (!is.null(psa_obj$strategies)) {
      cat("Strategies:", paste(psa_obj$strategies, collapse = ", "), "\n")
    }
  }

  psa_obj
}

#' Load the single-model PSA parameter cache
#'
#' @param util_label Utility-source label. Defaults to global utility_source_label.
#' @param directory Directory containing PSA cache files.
#' @param verbose Logical, print progress messages.
#' @return Data frame of PSA parameters, or NULL if the cache is not available.
load_psa_params_cache <- function(util_label = NULL,
                                  directory = cache_dir(),
                                  verbose = TRUE) {

  cache_file <- psa_params_path(util_label, directory)

  if (!file.exists(cache_file)) {
    warning("PSA parameter cache not found: ", cache_file)
    return(NULL)
  }

  psa_params <- readRDS(cache_file)

  if (verbose) {
    cat("Loaded PSA parameter cache:", cache_file, "\n")
  }

  psa_params
}

#' Create a single-model CEAC plot
#'
#' @param psa_obj dampack PSA object.
#' @param wtp_range Vector of willingness-to-pay thresholds.
#' @param wtp_line Reference willingness-to-pay threshold.
#' @return ggplot object.
create_ceac_plot <- function(psa_obj,
                             wtp_range = seq(0, 100000, by = 2000),
                             wtp_line = NULL) {

  if (is.null(wtp_line)) {
    wtp_line <- if (exists("WTP", inherits = TRUE)) get("WTP", inherits = TRUE) else 51000
  }

  if (is.null(psa_obj)) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No CEAC data available", size = 5) +
             theme_void())
  }

  ceac_data <- dampack::ceac(wtp = wtp_range, psa = psa_obj)
  names(ceac_data)[names(ceac_data) == "Proportion"] <- "Probability"

  ggplot(ceac_data, aes(x = WTP, y = Probability, color = Strategy)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = wtp_line, linetype = "dotted",
               color = "gray40", linewidth = 0.5) +
    scale_x_continuous(labels = scales::dollar_format(prefix = "EUR", scale = 0.001, suffix = "K")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    labs(
      x = "Willingness-to-Pay Threshold (EUR/QALY)",
      y = "Probability of Being Cost-Effective",
      color = "Strategy"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

#' Create a PSA summary table for the single model
#'
#' @param psa_obj dampack PSA object.
#' @param wtp WTP threshold for probability of cost-effectiveness.
#' @return Data frame with mean costs/effects and probability cost-effective.
create_psa_summary_table <- function(psa_obj, wtp = NULL) {

  if (is.null(wtp)) {
    wtp <- if (exists("WTP", inherits = TRUE)) get("WTP", inherits = TRUE) else 51000
  }

  if (is.null(psa_obj)) {
    return(data.frame(
      Strategy = character(),
      Mean_Cost = numeric(),
      Mean_QALY = numeric(),
      Prob_CE = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  ceac_obj <- dampack::ceac(wtp = c(wtp - 1, wtp, wtp + 1), psa = psa_obj)
  cost_matrix <- as.matrix(psa_obj$cost)
  effect_matrix <- if (!is.null(psa_obj$effect)) {
    as.matrix(psa_obj$effect)
  } else {
    as.matrix(psa_obj$effectiveness)
  }

  summary_list <- lapply(seq_along(psa_obj$strategies), function(i) {
    strategy <- psa_obj$strategies[i]
    ceac_row <- ceac_obj[ceac_obj$WTP == wtp & ceac_obj$Strategy == strategy, ]
    prob_ce <- if (nrow(ceac_row) > 0) ceac_row$Proportion[1] else NA_real_

    data.frame(
      Strategy = strategy,
      Mean_Cost = mean(cost_matrix[, i], na.rm = TRUE),
      Mean_QALY = mean(effect_matrix[, i], na.rm = TRUE),
      Prob_CE = prob_ce,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, summary_list)
}

#' Load scenario EVPPI results
#'
#' @param results_file Path to the .rds file containing scenario EVPPI results.
#' @return Scenario EVPPI result object, or NULL if unavailable.
load_scenario_evppi_results <- function(results_file = NULL) {

  if (is.null(results_file)) {
    results_file <- scenario_evppi_path()
  }

  if (!file.exists(results_file)) {
    warning("Scenario EVPPI results file not found: ", results_file)
    return(NULL)
  }

  readRDS(results_file)
}

message("Single-model CEA helper functions loaded successfully.")
