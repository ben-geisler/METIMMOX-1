# ===============================================================================
# PARAMETRIC SURVIVAL ANALYSIS WITH MODEL SELECTION
# ===============================================================================
# This script performs parametric survival analysis for the single economic
# survival model:
#   OS:  Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
#   PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
#
# TLR is not included in the economic survival model. It remains in the prepared
# dataset for clinical effectiveness, DAG, and biomarker distribution analyses.
#
# Testing and validation outputs are handled in tests/para_models.R
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr)

# Load parametric model fitting functions
para_model_fit_path <- here::here("R", "para_model_fit.R")
if (file.exists(para_model_fit_path)) {
  source(para_model_fit_path)
} else {
  message("File not found: ", para_model_fit_path)
}
rm(para_model_fit_path)

# Load prediction functions
prediction_functions_path <- here::here("R", "prediction_functions.R")
if (file.exists(prediction_functions_path)) {
  source(prediction_functions_path)
} else {
  message("File not found: ", prediction_functions_path)
}
rm(prediction_functions_path)

# ===============================================================================
# DATA PREPARATION
# ===============================================================================

# Create time points sequence for prediction (from setup file)
time_points <- seq(0, time_horizon, by = 1)

# Ensure categorical variables are properly coded as factors
# Note: sex is converted to factor in script 03 before subsets are created
# This prevents type mismatch in resampled models (issue #72)
if ("tlr" %in% names(data) && !is.null(data$tlr)) {
  data$tlr <- as.factor(data$tlr)
}
for (biomarker in get_biomarkers()) {
  if (biomarker %in% names(data) && !is.null(data[[biomarker]])) {
    data[[biomarker]] <- as.factor(data[[biomarker]])
  }
}

# Remove rows with missing values required for economic survival modeling.
# TLR is intentionally not part of this complete-case filter.
v_required_model_vars <- c(
  "Age", "sex", "Rx", get_biomarkers(),
  "OSwk", "Death", "PFSwk", "Progression"
)
data_complete <- economic_prediction_population(data)
rm(v_required_model_vars)

# ===============================================================================
# MODEL FORMULA DEFINITIONS
# ===============================================================================
# Model formulas are defined in R/model_configs.R, the single
# source of truth for economic survival strategies and formulas.
# ===============================================================================

# Source model configurations if not already loaded
if (!exists("get_model_configs")) {
  source(here::here("R/model_configs.R"))
}

# Get model configuration from central source
model_config <- get_current_model_config()
model_formulas <- get_model_formulas()

# Print configuration summary
cat("Configuration:", model_config$label, "\n")
cat("Description:", model_config$description, "\n")
cat("Economic strategies:", paste(model_config$strategies, collapse = ", "), "\n")
cat("Economic biomarkers:", paste(model_config$biomarkers, collapse = ", "), "\n")

# Print which formula sets will be fitted
cat("Formula sets to be fitted:", paste(names(model_formulas), collapse = ", "), "\n")

# Define parametric distributions to test
v_distributions_to_test <- c(
  "exponential", "weibull", "weibullph", "llogis",
  "lognormal", "gamma", "gompertz", "gengamma", "genf"
)

# ===============================================================================
# MODEL FITTING
# ===============================================================================

# Initialize storage structure:
# models$[formula_set]$[outcome]$[distribution] = fitted flexsurvreg object
# models$[formula_set]$[outcome]_ic = information criteria table
# models$best_fit$[outcome] = best model overall
models <- list()

# Initialize model formula lists based on configuration
for (formula_set in names(model_formulas)) {
  models[[formula_set]] <- list()
}

# Fit models for each selected formula set
for (formula_set in names(model_formulas)) {

  cat("Fitting models for formula set:", formula_set, "\n")

  # Fit OS models using fit_all_direct function
  os_results <- fit_all_direct(
    fit_data = data_complete,
    fit_formula = model_formulas[[formula_set]]$os,
    fit_dists = v_distributions_to_test
  )

  # Fit PFS models using fit_all_direct function
  pfs_results <- fit_all_direct(
    fit_data = data_complete,
    fit_formula = model_formulas[[formula_set]]$pfs,
    fit_dists = v_distributions_to_test
  )

  # Store fitted models in organized structure
  models[[formula_set]]$os <- os_results$fitted_models
  models[[formula_set]]$pfs <- pfs_results$fitted_models

  # Store Kaplan-Meier fits for reference
  # Warnings from candidate fits (per distribution; NULL when none)
  models[[formula_set]]$os_fit_warnings <- os_results$fit_warnings
  models[[formula_set]]$pfs_fit_warnings <- pfs_results$fit_warnings

  models[[formula_set]]$km_os <- os_results$km_all
  models[[formula_set]]$km_pfs <- pfs_results$km_all

  # Extract and store information criteria (AIC/BIC)
  models[[formula_set]]$os_ic <- extract_ic_single(os_results)
  models[[formula_set]]$pfs_ic <- extract_ic_single(pfs_results)

  # Identify best models within this formula set
  best_os <- find_best_model(models[[formula_set]]$os_ic, criterion = "AIC")
  best_pfs <- find_best_model(models[[formula_set]]$pfs_ic, criterion = "AIC")

  if (!is.null(best_os)) {
    models[[formula_set]]$best_os_dist <- best_os$distribution
    models[[formula_set]]$best_os_aic <- best_os$criterion_value
  }

  if (!is.null(best_pfs)) {
    models[[formula_set]]$best_pfs_dist <- best_pfs$distribution
    models[[formula_set]]$best_pfs_aic <- best_pfs$criterion_value
  }
}

# ===============================================================================
# ORDERING-CONSTRAINED JOINT MODEL SELECTION
# ===============================================================================

#' Collect the information-criteria tables of one outcome across formula sets
#'
#' Reads `models[[formula_set]][[paste0(outcome, "_ic")]]` from the script-level
#' `models` object for every formula set in `model_formulas` and stacks them.
#'
#' @param outcome Outcome prefix, "os" or "pfs".
#' @return Data frame with a `formula_set` column followed by the columns of
#'   `extract_ic_single()` (one row per fitted distribution), or NULL when no
#'   formula set has an information-criteria table for `outcome`.
collect_ic <- function(outcome) {
  pieces <- lapply(names(model_formulas), function(formula_set) {
    ic_name <- paste0(outcome, "_ic")
    df_ic <- models[[formula_set]][[ic_name]]

    if (is.null(df_ic)) {
      return(NULL)
    }

    data.frame(formula_set = formula_set, df_ic, stringsAsFactors = FALSE)
  })

  pieces <- Filter(Negate(is.null), pieces)

  if (length(pieces) == 0) {
    return(NULL)
  }

  do.call(rbind, pieces)
}

df_all_os_aic <- collect_ic("os")
df_all_pfs_aic <- collect_ic("pfs")

#' Build the candidate table for ordering-constrained joint model selection
#'
#' Drops candidates without an AIC, looks up each fitted model in the
#' script-level `models` object and predicts its population-averaged endpoint
#' curves for `data_complete` on `time_points` (control and every economic
#' biomarker subgroup). A prediction error is stored as the condition object
#' rather than raised, so `ordering_check()` can report it per pair.
#'
#' @param ic_data Output of `collect_ic()` for one outcome (columns
#'   `formula_set`, `Distribution`, `AIC`, ...).
#' @param outcome Outcome name in `models[[formula_set]]`, "os" or "pfs".
#' @return Data frame with columns `formula_set`, `distribution`, `AIC` and the
#'   list column `ordering_data` (curves from
#'   `generate_population_averaged_endpoint_curves()` or an error condition),
#'   as expected by `find_best_ordered_model_pair()`.
make_candidates <- function(ic_data, outcome) {
  ic_data <- ic_data[!is.na(ic_data$AIC), , drop = FALSE]
  candidate_models <- Map(
    function(formula_set, distribution) {
      models[[formula_set]][[outcome]][[distribution]]
    },
    ic_data$formula_set,
    ic_data$Distribution
  )
  ordering_data <- lapply(candidate_models, function(candidate_model) {
    tryCatch(
      generate_population_averaged_endpoint_curves(
        model = candidate_model,
        data_complete = data_complete,
        time_points = time_points
      ),
      error = identity
    )
  })

  data.frame(
    formula_set = ic_data$formula_set,
    distribution = ic_data$Distribution,
    AIC = ic_data$AIC,
    ordering_data = I(ordering_data),
    stringsAsFactors = FALSE
  )
}

if (!is.null(df_all_os_aic) && !is.null(df_all_pfs_aic)) {
  os_candidates <- make_candidates(df_all_os_aic, "os")
  pfs_candidates <- make_candidates(df_all_pfs_aic, "pfs")

  if (nrow(os_candidates) == 0 || nrow(pfs_candidates) == 0) {
    stop("No valid fitted OS or PFS candidates are available for joint selection.")
  }

  #' Check OS >= PFS for one candidate OS/PFS pair
  #'
  #' Re-raises a prediction error stored by `make_candidates()`, otherwise
  #' compares the two sets of population-averaged curves.
  #'
  #' @param os_curves `ordering_data` entry of an OS candidate (curves or an
  #'   error condition).
  #' @param pfs_curves `ordering_data` entry of a PFS candidate (curves or an
  #'   error condition).
  #' @return Result of `check_endpoint_curve_ordering()` for the pair.
  ordering_check <- function(os_curves, pfs_curves) {
    if (inherits(os_curves, "condition")) stop(conditionMessage(os_curves))
    if (inherits(pfs_curves, "condition")) stop(conditionMessage(pfs_curves))
    check_endpoint_curve_ordering(os_curves, pfs_curves)
  }

  ordered_selection <- find_best_ordered_model_pair(
    os_candidates = os_candidates,
    pfs_candidates = pfs_candidates,
    ordering_check = ordering_check,
    criterion = "AIC"
  )

  # Preserve the full audit trail, including rejected pairs and AIC penalty.
  models$ordered_selection <- ordered_selection

  if (is.null(ordered_selection$selected)) {
    stop(
      "No fitted OS/PFS distribution pair satisfies OS >= PFS over the full ",
      "time horizon for control and all biomarker subgroups. A PFS + ",
      "post-progression-survival decomposition is required."
    )
  }

  df_selected <- ordered_selection$selected[1, ]
  models$best_fit$os <- models[[df_selected$os_formula_set]]$os[[df_selected$os_distribution]]
  models$best_fit$pfs <- models[[df_selected$pfs_formula_set]]$pfs[[df_selected$pfs_distribution]]
  models$best_fit$os_structure <- df_selected$os_formula_set
  models$best_fit$pfs_structure <- df_selected$pfs_formula_set
  models$best_fit$os_distribution <- df_selected$os_distribution
  models$best_fit$pfs_distribution <- df_selected$pfs_distribution
  models$best_fit$os_aic <- df_selected$os_ic
  models$best_fit$pfs_aic <- df_selected$pfs_ic
  models$best_fit$combined_aic <- df_selected$combined_ic
  models$best_fit$ordering_aic_penalty <- ordered_selection$ic_penalty

  # Candidate-fit warnings are recorded, not raised; a warning from a selected
  # family still reaches the caller.
  v_selected_fit_warnings <- c(
    models[[df_selected$os_formula_set]]$os_fit_warnings[[df_selected$os_distribution]],
    models[[df_selected$pfs_formula_set]]$pfs_fit_warnings[[df_selected$pfs_distribution]])
  for (msg in unique(v_selected_fit_warnings)) {
    warning("Selected survival model fit: ", msg, call. = FALSE)
  }

  cat("Ordering-constrained model pair:", df_selected$os_distribution, "(OS),",
      df_selected$pfs_distribution, "(PFS); combined AIC:",
      round(df_selected$combined_ic, 2), "; constraint penalty:",
      round(ordered_selection$ic_penalty, 2), "\n")
}

# ===============================================================================
# POPULATION-LEVEL PREDICTIONS
# ===============================================================================
# Uses population averaging to account for heterogeneity in patient
# characteristics (age, sex, biomarkers). Predictions are generated for all
# patients using their actual covariate values, then averaged within relevant
# subgroups. This is the methodologically rigorous approach for CEA.
# See GitHub Issue #69 for methodology discussion.
# ===============================================================================

if (!is.null(models$best_fit$os) && !is.null(models$best_fit$pfs)) {
  best_models <- list(
    os = models$best_fit$os,
    pfs = models$best_fit$pfs
  )

  predictions <- generate_population_averaged_predictions(
    models = best_models,
    strategies_df = strategies_df,
    data_complete = data_complete,
    time_points = time_points,
    prevalences = setNames(strategies_df$prevalence, strategies_df$id)
  )

  basecase_ordering_check <- check_population_survival_ordering(predictions)
  if (!basecase_ordering_check$ordered) {
    stop("Internal error: selected base-case curves violate OS >= PFS.")
  }
} else {
  predictions <- list()
  warning("Model fitting failed - predictions object is empty")
}

# ===============================================================================
# FINAL STORAGE STRUCTURE DOCUMENTATION
# ===============================================================================
#
# CONFIGURATION:
# - Single economic survival model:
#     OS:  Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
#     PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
# - TLR is excluded from economic survival modeling.
#
# MODELS OBJECT STRUCTURE:
# - models$full
#     - $os: OS models by distribution
#     - $pfs: PFS models by distribution
#     - $os_ic: AIC/BIC table for OS models
#     - $pfs_ic: AIC/BIC table for PFS models
#     - $km_os: Kaplan-Meier fit for OS
#     - $km_pfs: Kaplan-Meier fit for PFS
# - models$best_fit
#     - $os: Best OS model object
#     - $os_structure: Formula set name ("full")
#     - $os_distribution: Best OS distribution
#     - $os_aic: AIC value
#     - equivalent fields for PFS
#     - $combined_aic: Combined AIC for the selected ordered pair
#     - $ordering_aic_penalty: Combined-AIC cost versus independent selection
# - models$ordered_selection
#     - $pairs: Audit table for every candidate OS/PFS distribution pair
#     - $selected: Minimum-combined-AIC pair satisfying OS >= PFS
#
# PREDICTIONS OBJECT STRUCTURE:
# - predictions$control: Standard of care strategy
#     - $os: Population OS curve
#     - $pfs: Population PFS curve
# - predictions$crp: CRP-guided strategy
#     - $os/$pfs: Population-weighted curves
#     - $biomarker_positive: CRP+ subgroup predictions
#     - $biomarker_negative: CRP- subgroup predictions
#     - $prevalence: CRP prevalence used for weighting
# - predictions$tmb_braf: TMB/BRAF-guided strategy
#     - same structure as CRP
# ===============================================================================
