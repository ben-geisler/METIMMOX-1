# ===============================================================================
# MODEL CONFIGURATION - SINGLE SOURCE OF TRUTH
# ===============================================================================
# Defines the single economic survival model, its strategies, biomarkers,
# display metadata, parameter keys, and formulas. The same shared formula is
# used for the CRP-guided and TMB/BRAF-guided strategies. The control arm is
# NOT a separate model: it is the same joint fit predicted with Rx forced to
# the control level, in both the base case and the PSA (issue #151).
#
# Economic analyses include standard of care plus two pre-immunotherapy biomarker
# strategies (CRP and TMB/BRAF), i.e. biomarkers that are available before the
# decision to add nivolumab is made. TMB/BRAF comes from baseline NGS. CRP is
# the week-4 value (cycle 3 day 1, trial visit 3; variable CRP1cat), measured
# after the two FLOX cycles that both arms receive and before the first
# nivolumab dose (issue #150). TLR is intentionally excluded here because it is
# a post-randomization mediator measured on treatment, not a treatment-selection
# biomarker.
#
# Usage:
#   source(here::here("R/model_configs.R"))
#   strategies <- get_strategies()
#   formula    <- get_strategy_formula("crp", "os")
# ===============================================================================

if (!require("survival")) {
  stop("survival package is required for model_configs.R")
}

#' Standard-of-care strategy ID
CONTROL_STRATEGY <- "control"

#' Economic strategies: standard of care plus pre-immunotherapy biomarker strategies
ALL_STRATEGIES <- c(CONTROL_STRATEGY, "crp", "tmb_braf")

#' Pre-immunotherapy biomarkers included in the economic model
ALL_BIOMARKERS <- c("crp", "tmb_braf")

#' Keyed display metadata for economic strategies
#'
#' Row order is deliberately independent of ALL_STRATEGIES. Accessors match by
#' strategy ID so changing strategy order cannot silently relabel output.
STRATEGY_METADATA <- data.frame(
  id = c("control", "crp", "tmb_braf"),
  short_name = c("Control", "CRP", "TMB/BRAF"),
  report_label = c("Standard of Care", "CRP-guided", "TMB/BRAF-guided"),
  trace_label = c("Standard of Care", "CRP Strategy", "TMB/BRAF Strategy"),
  name = c(
    "Standard of care: FLOX chemotherapy only",
    "Biomarker-guided: C-reactive protein",
    "Biomarker-guided: tumor mutation burden or BRAF mutation"
  ),
  description = c(
    "Standard of care - All patients receive only FLOX chemotherapy",
    paste0(
      "C-reactive protein with cut-off of <5 mg/L for biomarker-positive status, ",
      "measured at week 4 (cycle 3 day 1) after two FLOX cycles common to both ",
      "arms and before the first nivolumab dose. ",
      "If CRP-positive: alternating two cycles each of FLOX (chemotherapy) ",
      "and nivolumab (anti-PD1 immunotherapy); if CRP-negative: chemotherapy only"
    ),
    paste0(
      "Combined biomarker: either Tumor Mutation Burden >= 9 or BRAF V600 ",
      "mutation positive (both from next-generation sequencing). If ",
      "TMB/BRAF-positive: alternating two cycles each of FLOX (chemotherapy) ",
      "and nivolumab (anti-PD1 immunotherapy); if TMB/BRAF-negative: ",
      "chemotherapy only"
    )
  ),
  stringsAsFactors = FALSE
)

#' Keyed parameter metadata for economic biomarkers
BIOMARKER_METADATA <- data.frame(
  id = c("crp", "tmb_braf"),
  cost_key = c("c_test_CRP", "c_test_NGS"),
  stringsAsFactors = FALSE
)

#' Get the economic model configuration
#'
#' Single source of truth for the joint economic survival model. Returns a list
#' with label/description, the strategy and biomarker sets, and the shared
#' joint formulas (the control arm is these formulas predicted with
#' Rx = control; it has no formula of its own).
#'
#' @return Model configuration list.
#' @export
get_current_model_config <- function() {
  list(
    label = "Joint economic survival model",
    description = "CRP and TMB/BRAF treatment interactions in one model",
    control_strategy = CONTROL_STRATEGY,
    strategies = ALL_STRATEGIES,
    biomarkers = ALL_BIOMARKERS,
    strategy_metadata = STRATEGY_METADATA,
    biomarker_metadata = BIOMARKER_METADATA,
    formulas = list(
      shared = list(
        os  = Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx,
        pfs = Surv(PFSwk, Progression) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
      )
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

#' Get the standard-of-care strategy ID
#' @return Character scalar.
#' @export
get_control_strategy <- function() {
  get_current_model_config()$control_strategy
}

#' Get economic biomarkers
#' @return Character vector of biomarker names.
#' @export
get_biomarkers <- function() {
  get_current_model_config()$biomarkers
}

#' Get keyed economic-strategy metadata
#'
#' @param strategies Strategy IDs. Defaults to the configured strategy order.
#' @return Data frame ordered to match `strategies`.
#' @export
get_strategy_metadata <- function(strategies = get_strategies()) {
  metadata <- get_current_model_config()$strategy_metadata
  unknown <- setdiff(strategies, metadata$id)
  if (length(unknown) > 0L) {
    stop("Missing strategy metadata for: ", paste(unknown, collapse = ", "))
  }
  metadata[match(strategies, metadata$id), , drop = FALSE]
}

#' Get short display names for economic strategies
#'
#' @param strategies Strategy IDs. Defaults to the configured strategy order.
#' @return Named character vector.
#' @export
strategy_display_name <- function(strategies = get_strategies()) {
  metadata <- get_strategy_metadata(strategies)
  setNames(metadata$short_name, metadata$id)
}

#' Get scalar diagnostic-test parameter keys for biomarkers
#'
#' @param biomarkers Biomarker IDs. Defaults to all configured biomarkers.
#' @return Named character vector mapping biomarker IDs to parameter keys.
#' @export
biomarker_cost_key <- function(biomarkers = get_biomarkers()) {
  metadata <- get_current_model_config()$biomarker_metadata
  unknown <- setdiff(biomarkers, metadata$id)
  if (length(unknown) > 0L) {
    stop("Missing diagnostic-test cost key for biomarkers: ",
         paste(unknown, collapse = ", "))
  }
  setNames(metadata$cost_key[match(biomarkers, metadata$id)], biomarkers)
}

#' Get prevalence parameter keys for biomarkers
#'
#' @param biomarkers Biomarker IDs. Defaults to all configured biomarkers.
#' @return Named character vector mapping biomarker IDs to parameter keys.
#' @export
biomarker_prevalence_key <- function(biomarkers = get_biomarkers()) {
  unknown <- setdiff(biomarkers, get_biomarkers())
  if (length(unknown) > 0L) {
    stop("Unknown economic biomarkers: ", paste(unknown, collapse = ", "))
  }
  setNames(paste0("p_", biomarkers), biomarkers)
}

#' Synchronize the diagnostic-cost lookup with scalar model parameters
#'
#' @param params Model parameter list.
#' @return Updated parameter list.
#' @export
sync_biomarker_test_costs <- function(params) {
  biomarkers <- get_biomarkers()
  cost_keys <- biomarker_cost_key(biomarkers)
  missing <- setdiff(unname(cost_keys), names(params))
  if (length(missing) > 0L) {
    stop("Missing scalar diagnostic-test cost parameters: ",
         paste(missing, collapse = ", "))
  }
  if (is.null(params$c_test_biomarker)) {
    params$c_test_biomarker <- list()
  }
  for (biomarker in biomarkers) {
    params$c_test_biomarker[[biomarker]] <- params[[cost_keys[[biomarker]]]]
  }
  params
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
