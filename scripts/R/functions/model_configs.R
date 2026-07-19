# ===============================================================================
# MODEL CONFIGURATION - SINGLE SOURCE OF TRUTH
# ===============================================================================
# This file defines the single economic survival model, strategies, biomarkers,
# and formulas used by the cost-effectiveness analysis.
#
# Economic analyses include standard of care plus two pre-treatment biomarker
# strategies: CRP and TMB/BRAF. TLR is intentionally excluded here because it is a
# post-randomization mediator, not a baseline treatment-selection biomarker.
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
# CORE CONSTANTS
# ===============================================================================

#' Economic strategies: standard of care plus pre-treatment biomarker strategies
ALL_STRATEGIES <- c("control", "crp", "tmb_braf")

#' Pre-treatment biomarkers included in the economic model
ALL_BIOMARKERS <- c("crp", "tmb_braf")

# ===============================================================================
# MODEL DEFINITION
# ===============================================================================

#' Get complete model configuration
#'
#' Returns the single joint economic model configuration. The same shared formula
#' is used for CRP-guided and TMB/BRAF-guided strategies.
#'
#' @return Named list with one joint economic model configuration
#' @export
get_model_configs <- function() {
  list(
    joint = list(
      label = "Joint economic survival model",
      description = "CRP and TMB/BRAF treatment interactions in one model",
      strategies = ALL_STRATEGIES,
      biomarkers = ALL_BIOMARKERS,
      formulas = list(
        shared = list(
          os = Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx,
          pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
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

.validate_outcome <- function(outcome) {
  if (!outcome %in% c("os", "pfs")) {
    stop("Outcome must be 'os' or 'pfs', got: ", outcome)
  }
}

#' Get the current economic model config
#'
#' @return Model configuration list
#' @export
get_current_model_config <- function() {
  get_model_configs()$joint
}

#' Get economic strategies
#'
#' @return Character vector of strategy IDs
#' @export
get_strategies <- function() {
  get_current_model_config()$strategies
}

#' Get economic biomarkers
#'
#' @return Character vector of biomarker names
#' @export
get_biomarkers <- function() {
  get_current_model_config()$biomarkers
}

#' Get shared formula for a biomarker strategy and outcome
#'
#' @param strategy Strategy ID ("crp" or "tmb_braf")
#' @param outcome "os" or "pfs"
#' @return Formula object
#' @export
get_strategy_formula <- function(strategy, outcome) {
  config <- get_current_model_config()
  .validate_outcome(outcome)

  if (!strategy %in% config$biomarkers) {
    stop("Strategy '", strategy, "' not found. Available biomarker strategies: ",
         paste(config$biomarkers, collapse = ", "))
  }

  config$formulas$shared[[outcome]]
}

#' Get control formula
#'
#' @param outcome "os" or "pfs"
#' @return Formula object
#' @export
get_control_formula <- function(outcome) {
  config <- get_current_model_config()
  .validate_outcome(outcome)

  config$control_formulas[[outcome]]
}

#' Get model formulas list for survival analysis
#'
#' Returns a list of formulas suitable for use in
#' 06_parametric_survival_analysis.R.
#'
#' @return Named list of formula lists
#' @export
get_model_formulas <- function() {
  config <- get_current_model_config()

  list(
    full = config$formulas$shared
  )
}

#' Print current model configuration summary
#'
#' @export
print_model_config <- function() {
  config <- get_current_model_config()

  cat("================================================================================\n")
  cat("MODEL CONFIGURATION:", config$label, "\n")
  cat("================================================================================\n")
  cat("Description:", config$description, "\n")
  cat("Strategies:", paste(config$strategies, collapse = ", "), "\n")
  cat("Biomarkers:", paste(config$biomarkers, collapse = ", "), "\n")
  cat("OS formula:", deparse(config$formulas$shared$os), "\n")
  cat("PFS formula:", deparse(config$formulas$shared$pfs), "\n")
  cat("================================================================================\n")

  invisible(config)
}

# ===============================================================================
# INITIALIZATION MESSAGE
# ===============================================================================

message("Model configuration functions loaded from model_configs.R")
message("  - Economic strategies: control, crp, tmb_braf")
message("  - Shared model: Age + sex + Rx + crp*Rx + tmb_braf*Rx")
message("  - TLR is excluded from economic strategy modeling")
