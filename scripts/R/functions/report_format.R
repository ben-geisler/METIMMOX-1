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
  x_chr <- as.character(x)
  out <- x_chr
  known <- !is.na(x_chr) & x_chr %in% names(strategy_labels)
  out[known] <- unname(strategy_labels[x_chr[known]])
  out
}

#' Convert biomarker identifiers to display labels
#'
#' Unknown identifiers are returned unchanged.
biomarker_label <- function(x) {
  labels <- c(strategy_display_name(get_biomarkers()), tlr = "TLR")
  x_chr <- as.character(x)
  out <- x_chr
  known <- !is.na(x_chr) & x_chr %in% names(labels)
  out[known] <- unname(labels[x_chr[known]])
  out
}

#' Convert EVPPI parameter-group identifiers to display labels
#'
#' Accepts identifiers with or without the "[GROUP] " prefix used in the EVPPI
#' result tables (e.g. "[GROUP] drug_costs" or "drug_costs"). Unknown
#' identifiers are returned without the prefix.
evppi_group_label <- function(x) {
  biomarkers <- get_biomarkers()
  labels <- c(
    drug_costs = "Drug costs",
    test_costs = "Test costs",
    other_costs = "Other costs",
    all_costs = "All costs",
    utilities = "Utilities",
    prevalence = "Biomarker prevalence",
    stats::setNames(paste0(biomarker_label(biomarkers), "-treatment interaction"),
                    paste0("interaction_", biomarkers)),
    interaction_all = "Biomarker-treatment interaction"
  )
  x_chr <- sub("^\\[GROUP\\] ", "", as.character(x))
  out <- x_chr
  known <- !is.na(x_chr) & x_chr %in% names(labels)
  out[known] <- unname(labels[x_chr[known]])
  out
}

#' Format values as euros, guarding missing and non-finite values
format_eur <- function(x, accuracy = 1) {
  out <- rep("--", length(x))
  valid <- !is.na(x) & is.finite(x)
  out[valid] <- scales::dollar(x[valid], prefix = "EUR ", accuracy = accuracy)
  out
}

# Backward-compatible short name used by input and poster reports.
fmt_eur <- format_eur

#' Format an ICER using its frontier status
#'
#' Reference strategies are shown as "--" and dominated or extended-dominated
#' strategies as "Dominated". Missing and infinite estimates are also "--".
format_icer <- function(x, status = NULL) {
  if (length(x) == 0) return(character())

  status <- if (is.null(status)) {
    rep(NA_character_, length(x))
  } else {
    rep_len(as.character(status), length(x))
  }

  out <- format_eur(x)
  out[status %in% c("D", "ED", "Dominated")] <- "Dominated"
  out[status %in% c("Reference", "ref")] <- "--"
  out
}
