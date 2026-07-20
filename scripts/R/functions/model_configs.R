# ===============================================================================
# MODEL CONFIGURATION - SINGLE SOURCE OF TRUTH
# ===============================================================================
# Defines the single economic survival model, its strategies, biomarkers, and
# formulas. The same shared formula is used for the CRP-guided and TMB/BRAF-guided
# strategies; the control arm is age/sex-adjusted only.
#
# Economic analyses include standard of care plus two pre-treatment biomarker
# strategies (CRP and TMB/BRAF). TLR is intentionally excluded here because it is a
# post-randomization mediator, not a baseline treatment-selection biomarker.
#
# Usage:
#   source(here::here("scripts/R/functions/model_configs.R"))
#   strategies <- get_strategies()
#   formula    <- get_strategy_formula("crp", "os")
# ===============================================================================

if (!require("survival")) {
  stop("survival package is required for model_configs.R")
}

#' Economic strategies: standard of care plus pre-treatment biomarker strategies
ALL_STRATEGIES <- c("control", "crp", "tmb_braf")

#' Pre-treatment biomarkers included in the economic model
ALL_BIOMARKERS <- c("crp", "tmb_braf")

#' Get the economic model configuration
#'
#' Single source of truth for the joint economic survival model. Returns a list
#' with label/description, the strategy and biomarker sets, the shared biomarker
#' formulas, and the control formulas.
#'
#' @return Model configuration list.
#' @export
get_current_model_config <- function() {
  list(
    label = "Joint economic survival model",
    description = "CRP and TMB/BRAF treatment interactions in one model",
    strategies = ALL_STRATEGIES,
    biomarkers = ALL_BIOMARKERS,
    formulas = list(
      shared = list(
        os  = Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx,
        pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
      )
    ),
    control_formulas = list(
      os  = Surv(OSwk, Death) ~ Age + sex,
      pfs = Surv(PFSwk, Progression) ~ Age + sex
    )
  )
}

#' Backward-compatible accessor
#'
#' Historically returned a named list of alternative model configurations (the
#' retired Model A/B/C design). Only one economic model remains; this shim keeps
#' the single-key shape for any legacy caller.
#'
#' @return Named list with a single \code{joint} configuration.
#' @export
get_model_configs <- function() {
  list(joint = get_current_model_config())
}

.validate_outcome <- function(outcome) {
  if (!outcome %in% c("os", "pfs")) {
    stop("Outcome must be 'os' or 'pfs', got: ", outcome)
  }
}

#' Get economic strategies
#' @return Character vector of strategy IDs.
#' @export
get_strategies <- function() {
  get_current_model_config()$strategies
}

#' Get economic biomarkers
#' @return Character vector of biomarker names.
#' @export
get_biomarkers <- function() {
  get_current_model_config()$biomarkers
}

#' Get the shared survival formula for a biomarker strategy and outcome
#'
#' @param strategy Biomarker strategy ID ("crp" or "tmb_braf").
#' @param outcome "os" or "pfs".
#' @return Formula object.
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

#' Get the control-arm survival formula for an outcome
#'
#' @param outcome "os" or "pfs".
#' @return Formula object.
#' @export
get_control_formula <- function(outcome) {
  .validate_outcome(outcome)
  get_current_model_config()$control_formulas[[outcome]]
}

#' Get model formulas for survival analysis
#'
#' Returns the shared formulas under a single \code{full} formula set, as consumed
#' by 04_parametric_survival_analysis.R and para_model_fit_table.R.
#'
#' @return Named list with one \code{full} element.
#' @export
get_model_formulas <- function() {
  list(full = get_current_model_config()$formulas$shared)
}

#' Print a summary of the current model configuration
#' @export
print_model_config <- function() {
  config <- get_current_model_config()
  cat(strrep("=", 80), "\n", sep = "")
  cat("MODEL CONFIGURATION:", config$label, "\n")
  cat(strrep("=", 80), "\n", sep = "")
  cat("Description:", config$description, "\n")
  cat("Strategies:", paste(config$strategies, collapse = ", "), "\n")
  cat("Biomarkers:", paste(config$biomarkers, collapse = ", "), "\n")
  cat("OS formula:", deparse(config$formulas$shared$os), "\n")
  cat("PFS formula:", deparse(config$formulas$shared$pfs), "\n")
  cat(strrep("=", 80), "\n", sep = "")
  invisible(config)
}

message("Model configuration functions loaded from model_configs.R")
message("  - Economic strategies: control, crp, tmb_braf")
message("  - Shared model: Age + sex + Rx + crp*Rx + tmb_braf*Rx")
message("  - TLR is excluded from economic strategy modeling")
