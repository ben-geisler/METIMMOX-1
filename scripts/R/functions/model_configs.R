# ===============================================================================
# MODEL CONFIGURATION - SINGLE SOURCE OF TRUTH
# ===============================================================================
# This file defines all model structures, strategies, biomarkers, and formulas.
# All other files should source this and use its helper functions.
#
# Key design principles:
# 1. All three biomarker strategies (CRP, TLR, TMB/BRAF) are available in ALL models
# 2. Model B uses per-strategy formulas: TLR gets its own formula
# 3. Helper functions provide MODEL_STRUCTURE-aware access to configurations
#
# Usage:
#   source(here::here("scripts/R/functions/model_configs.R"))
#   strategies <- get_strategies()
#   biomarkers <- get_biomarkers()
#   formula <- get_strategy_formula("crp", "os")
# ===============================================================================

# Load required packages
if (!require("survival")) {
  stop("survival package is required for model_configs.R")
}

# ===============================================================================
# CORE CONSTANTS (NOT MODEL_STRUCTURE dependent)
# ===============================================================================

#' All possible strategies (always includes all biomarkers + control)
ALL_STRATEGIES <- c("control", "crp", "tlr", "tmb_braf")

#' All possible biomarkers
ALL_BIOMARKERS <- c("crp", "tlr", "tmb_braf")

# ===============================================================================
# MODEL STRUCTURE DEFINITIONS
# ===============================================================================

#' Get complete model configurations
#'
#' Returns a named list defining all three model structures with their
#' formulas, strategies, and biomarkers. Model B uses per-strategy formulas
#' to allow different modeling approaches for TLR vs CRP/TMB_BRAF.
#'
#' @return Named list with Model_A, Model_B, Model_C configurations
#' @export
get_model_configs <- function() {
  list(
    # =========================================================================
    # MODEL A: JOINT
    # All biomarkers + all treatment interactions in ONE model
    # =========================================================================
    Model_A = list(
      structure = 0,
      label = "Model A: Joint",
      short_label = "A",
      model_type_label = "joint",
      description = "All biomarkers + all treatment interactions in one model",
      strategies = ALL_STRATEGIES,
      biomarkers = ALL_BIOMARKERS,
      # Single formula used for all biomarker strategies
      formulas = list(
        shared = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tlr*Rx + tmb_braf*Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tlr*Rx + tmb_braf*Rx
        )
      ),
      control_formulas = list(
        os = Surv(OSwk, Death) ~ Age + sex,
        pfs = Surv(PFSwk, Progression) ~ Age + sex
      )
    ),

    # =========================================================================
    # MODEL B: FOCUSED (Per-Strategy Formulas)
    # CRP and TMB/BRAF share a formula with both interactions
    # TLR uses its own formula with TLR-only interaction (adjusting for other biomarkers)
    # =========================================================================
    Model_B = list(
      structure = 1,
      label = "Model B: Focused",
      short_label = "B",
      model_type_label = "focused",
      description = "CRP/TMB_BRAF with both interactions; TLR with TLR-only interaction (adjusts for other biomarkers)",
      strategies = ALL_STRATEGIES,
      biomarkers = ALL_BIOMARKERS,
      # Per-strategy formulas
      formulas = list(
        # CRP strategy: includes CRP and TMB/BRAF as main effects + BOTH interactions
        crp = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx
        ),
        # TLR strategy: includes CRP and TMB/BRAF as main effects but only TLR:Rx interaction
        tlr = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tmb_braf + tlr:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tmb_braf + tlr:Rx
        ),
        # TMB/BRAF strategy: same as CRP (shared focused formula)
        tmb_braf = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx
        )
      ),
      control_formulas = list(
        os = Surv(OSwk, Death) ~ Age + sex,
        pfs = Surv(PFSwk, Progression) ~ Age + sex
      )
    ),

    # =========================================================================
    # MODEL C: SEPARATE
    # Only ONE biomarker + its interaction per model (most parsimonious)
    # =========================================================================
    Model_C = list(
      structure = 2,
      label = "Model C: Separate",
      short_label = "C",
      model_type_label = "separate",
      description = "One biomarker + its interaction only per model",
      strategies = ALL_STRATEGIES,
      biomarkers = ALL_BIOMARKERS,
      # Per-strategy formulas (one biomarker + its interaction only)
      formulas = list(
        crp = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp:Rx
        ),
        tlr = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + tlr:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tlr:Rx
        ),
        tmb_braf = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + tmb_braf:Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + tmb_braf:Rx
        )
      ),
      control_formulas = list(
        os = Surv(OSwk, Death) ~ Age + sex,
        pfs = Surv(PFSwk, Progression) ~ Age + sex
      )
    )
  )
}

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

#' Get the current model config based on MODEL_STRUCTURE
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global MODEL_STRUCTURE
#' @return Model configuration list
#' @export
get_current_model_config <- function(model_structure = NULL) {
  if (is.null(model_structure)) {
    if (!exists("MODEL_STRUCTURE", envir = .GlobalEnv)) {
      stop("MODEL_STRUCTURE not defined. Source 02_setup_and_global_variables.R first.")
    }
    model_structure <- get("MODEL_STRUCTURE", envir = .GlobalEnv)
  }

  validate_model_structure(model_structure)

  configs <- get_model_configs()
  model_name <- c("Model_A", "Model_B", "Model_C")[model_structure + 1]
  return(configs[[model_name]])
}

#' Get strategies for current MODEL_STRUCTURE
#'
#' Always returns all 4 strategies (control, crp, tlr, tmb_braf) since
#' all strategies are available in all model structures.
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Character vector of strategy IDs
#' @export
get_strategies <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)
  return(config$strategies)
}

#' Get biomarkers for current MODEL_STRUCTURE
#'
#' Always returns all 3 biomarkers (crp, tlr, tmb_braf) since all
#' biomarkers are available in all model structures.
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Character vector of biomarker names
#' @export
get_biomarkers <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)
  return(config$biomarkers)
}

#' Get formula for a specific strategy and outcome
#'
#' Returns the appropriate survival formula for the given biomarker strategy
#' and outcome (os or pfs). For Model A (joint), all strategies use the same
#' formula. For Models B and C, each strategy may have its own formula.
#'
#' @param strategy Strategy ID ("crp", "tlr", or "tmb_braf")
#' @param outcome "os" or "pfs"
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Formula object
#' @export
get_strategy_formula <- function(strategy, outcome, model_structure = NULL) {
  config <- get_current_model_config(model_structure)

  # Validate outcome
  if (!outcome %in% c("os", "pfs")) {
    stop("Outcome must be 'os' or 'pfs', got: ", outcome)
  }

  # Check if model uses shared formula (Model A) or per-strategy (Models B, C)
  if ("shared" %in% names(config$formulas)) {
    return(config$formulas$shared[[outcome]])
  } else {
    if (!strategy %in% names(config$formulas)) {
      stop("Strategy '", strategy, "' not found in model ", config$short_label,
           ". Available strategies: ", paste(names(config$formulas), collapse = ", "))
    }
    return(config$formulas[[strategy]][[outcome]])
  }
}

#' Get control formulas
#'
#' @param outcome "os" or "pfs"
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Formula object
#' @export
get_control_formula <- function(outcome, model_structure = NULL) {
  config <- get_current_model_config(model_structure)

  if (!outcome %in% c("os", "pfs")) {
    stop("Outcome must be 'os' or 'pfs', got: ", outcome)
  }

  return(config$control_formulas[[outcome]])
}

#' Check if model uses per-strategy formulas
#'
#' Returns TRUE for Models B and C (which have different formulas for
#' different biomarker strategies), FALSE for Model A (which uses one
#' shared formula for all strategies).
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Logical
#' @export
uses_per_strategy_formulas <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)
  return(!"shared" %in% names(config$formulas))
}

#' Get model type label for cache filenames
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Character: "joint", "focused", or "separate"
#' @export
get_model_type_label <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)
  return(config$model_type_label)
}

#' Get model short label (A, B, or C)
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Character: "A", "B", or "C"
#' @export
get_model_short_label <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)
  return(config$short_label)
}

#' Validate MODEL_STRUCTURE value
#'
#' @param model_structure Value to validate
#' @return TRUE if valid, throws error otherwise
#' @export
validate_model_structure <- function(model_structure) {
  if (!model_structure %in% c(0, 1, 2)) {
    stop("Invalid MODEL_STRUCTURE value: ", model_structure,
         ". Must be 0 (joint/A), 1 (focused/B), or 2 (separate/C)")
  }
  invisible(TRUE)
}

#' Get model formulas list for survival analysis
#'
#' Returns a list of formulas suitable for use in 06_parametric_survival_analysis.R.
#' For Model A, returns list(full = list(os = ..., pfs = ...)).
#' For Models B and C, returns list(crp = list(os = ..., pfs = ...), tlr = ..., tmb_braf = ...).
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @return Named list of formula lists
#' @export
get_model_formulas <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)

  if ("shared" %in% names(config$formulas)) {
    # Model A: return as "full" structure for backward compatibility
    return(list(
      full = config$formulas$shared
    ))
  } else {
    # Models B and C: return per-strategy formulas
    return(config$formulas)
  }
}

#' Print current model configuration summary
#'
#' @param model_structure Integer (0, 1, or 2) or NULL to use global
#' @export
print_model_config <- function(model_structure = NULL) {
  config <- get_current_model_config(model_structure)

  cat("================================================================================\n")
  cat("MODEL CONFIGURATION:", config$label, "\n")
  cat("================================================================================\n")
  cat("Structure ID:", config$structure, "\n")
  cat("Type label:", config$model_type_label, "\n")
  cat("Description:", config$description, "\n")
  cat("Strategies:", paste(config$strategies, collapse = ", "), "\n")
  cat("Biomarkers:", paste(config$biomarkers, collapse = ", "), "\n")
  cat("Uses per-strategy formulas:", uses_per_strategy_formulas(config$structure), "\n")
  cat("================================================================================\n")

  invisible(config)
}

# ===============================================================================
# INITIALIZATION MESSAGE
# ===============================================================================

message("Model configuration functions loaded from model_configs.R")
message("  - All models include all 3 biomarker strategies (crp, tlr, tmb_braf)")
message("  - Use get_strategies(), get_biomarkers(), get_strategy_formula() for access")
