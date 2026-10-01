# ===============================================================================
# SHARED REPORT LABELS AND FORMATTERS
# ===============================================================================
# Canonical presentation helpers used by the economic reports and publication
# vignettes. Keep display wording and missing-value handling consistent here.
# ===============================================================================

#' Canonical economic-strategy display labels
strategy_labels <- with(get_strategy_metadata(), setNames(report_label, id))

#' Convert economic-strategy identifiers to display labels
#'
#' Unknown identifiers are returned unchanged.
strategy_label <- function(x) {
  v_ids <- as.character(x)
  v_out <- v_ids
  v_known <- !is.na(v_ids) & v_ids %in% names(strategy_labels)
  v_out[v_known] <- unname(strategy_labels[v_ids[v_known]])
  v_out
}

#' Short labels for narrow tables; guided strategies retain their canonical names.
strategy_short_label <- function(x) {
  v_out <- strategy_label(x)
  v_out[!is.na(x) & as.character(x) == get_control_strategy()] <- "SoC"
  v_out
}

#' Convert biomarker identifiers to display labels
#'
#' Unknown identifiers are returned unchanged.
biomarker_label <- function(x) {
  v_labels <- c(strategy_display_name(get_biomarkers()), tlr = "TLR")
  v_ids <- as.character(x)
  v_out <- v_ids
  v_known <- !is.na(v_ids) & v_ids %in% names(v_labels)
  v_out[v_known] <- unname(v_labels[v_ids[v_known]])
  v_out
}

#' Convert EVPPI parameter-group identifiers to display labels
#'
#' Accepts identifiers with or without the "[GROUP] " prefix used in the EVPPI
#' result tables (e.g. "[GROUP] drug_costs" or "drug_costs"). Unknown
#' identifiers are returned without the prefix.
evppi_group_label <- function(x) {
  v_biomarkers <- get_biomarkers()
  v_labels <- c(
    # Drug and test unit prices are fixed in the PSA since issue #154, so these
    # groups no longer appear; the labels are kept for older caches.
    drug_costs = "Drug costs",
    test_costs = "Test costs",
    other_costs = "Resource-use costs",
    all_costs = "All costs",
    utilities = "Utilities",
    prevalence = "Biomarker prevalence",
    stats::setNames(paste0(biomarker_label(v_biomarkers), "-treatment interaction"),
                    paste0("interaction_", v_biomarkers)),
    interaction_all = "Biomarker-treatment interaction"
  )
  v_ids <- sub("^\\[GROUP\\] ", "", as.character(x))
  v_out <- v_ids
  v_known <- !is.na(v_ids) & v_ids %in% names(v_labels)
  v_out[v_known] <- unname(v_labels[v_ids[v_known]])
  v_out
}

#' Format values as euros, guarding missing and non-finite values
format_eur <- function(x, accuracy = 1) {
  v_out <- rep("--", length(x))
  v_valid <- !is.na(x) & is.finite(x)
  v_out[v_valid] <- scales::dollar(x[v_valid], prefix = "EUR ", accuracy = accuracy)
  v_out
}

# Backward-compatible short name used by input and poster reports.
fmt_eur <- format_eur

#' Compact display of probabilistic means and 95% percentile intervals (#180).
#' The numeric summary remains the publication CSV; this is presentation only.
format_psa_summary <- function(x) {
  interval <- function(mean, lower, upper, money = FALSE) {
    fmt <- if (money) format_eur else function(v) sprintf("%.4f", v)
    paste0(fmt(mean), " (", fmt(lower), "; ", fmt(upper), ")")
  }
  data.frame(
    Strategy = strategy_label(x$Strategy),
    Cost = interval(x$Mean_Cost, x$Cost_Lower, x$Cost_Upper, TRUE),
    QALYs = interval(x$Mean_QALY, x$QALY_Lower, x$QALY_Upper),
    Inc_Cost = interval(x$Mean_Inc_Cost, x$Inc_Cost_Lower, x$Inc_Cost_Upper, TRUE),
    Inc_QALYs = interval(x$Mean_Inc_QALY, x$Inc_QALY_Lower, x$Inc_QALY_Upper),
    ICER = format_icer(x$ICER, x$Status),
    Status = x$Status,
    Prob_CE = scales::percent(x$Prob_CE, accuracy = 0.1)
  )
}

#' Format an ICER using its frontier or pairwise status
#'
#' Accepts the dampack frontier codes ("ND", "D", "ED") and the pairwise
#' statuses of `calculate_pairwise_icers()` ("Reference", "Pairwise ICER vs
#' SoC", "Dominated by SoC", "Cost-saving vs SoC"; issue #153). Reference
#' strategies are shown as "--", dominated or extendedly dominated strategies
#' as "Dominated", and cost-saving (dominant) strategies as "Cost-saving": their
#' ratio is negative because the incremental cost is negative, and printing it
#' would read as a cost per QALY (issue #157). A negative ratio with no
#' informative status is also "--", because without the signs of the increments
#' a cost-saving strategy cannot be told from a dominated one. Missing and
#' infinite estimates are "--".
format_icer <- function(x, status = NULL) {
  if (length(x) == 0) return(character())

  status <- if (is.null(status)) {
    rep(NA_character_, length(x))
  } else {
    rep_len(as.character(status), length(x))
  }

  v_out <- format_eur(x)
  v_reference <- status %in% c("Reference", "ref")
  v_dominated <- status %in% c("D", "ED", "Dominated") |
    grepl("^Dominated", status)
  v_cost_saving <- status %in% c("Dominant", "Cost-saving") |
    grepl("^Cost-saving", status) | grepl("^Dominant", status)
  v_out[!is.na(v_dominated) & v_dominated] <- "Dominated"
  v_out[!is.na(v_cost_saving) & v_cost_saving] <- "Cost-saving"
  v_out[v_reference] <- "--"
  v_unlabelled_negative <- is.finite(x) & x < 0 & !v_reference &
    !(!is.na(v_dominated) & v_dominated) & !(!is.na(v_cost_saving) & v_cost_saving)
  v_out[v_unlabelled_negative] <- "--"
  v_out
}

#' Restricted mean survival time of a curve on the model's weekly grid
#'
#' Trapezoidal integration of a survival curve evaluated at equally spaced time
#' points (the weekly model cycle), returned in years. The mean is restricted to
#' the horizon spanned by the curve: 520 weeks (10 years) for the base-case
#' curves in `l_params_base`. The earlier left Riemann sum `sum(y) * cl`
#' included the t = 0 point as a full cycle and so extended the integral by one
#' cycle beyond the horizon (issue #157).
#'
#' @param y Survival probabilities at consecutive grid points, starting at t = 0.
#' @param cycle_length Grid spacing in years (defaults to the global `cl`).
#' @return Restricted mean survival time in years.
restricted_mean_survival <- function(y, cycle_length = cl) {
  y <- as.numeric(y)
  n <- length(y)
  if (n < 2 || anyNA(y)) return(NA_real_)
  (sum(y) - (y[1] + y[n]) / 2) * cycle_length
}
