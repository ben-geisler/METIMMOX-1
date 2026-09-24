# PFS endpoint, single source of truth (issues #149, #181).
# Protocol: CT every 8 weeks; progression on active therapy or death, with
# censoring at the last imaging assessment. Apply a two-interval death window.
# TTPwk is the primary assessment proxy; LastEvalwk is the alternative.
# SAM OS / SAM PFS remain unread pending reconciliation with the data provider.

#' Derive PFS from preserved progression-exit variables
#' @param data Data frame containing Progression, PFSwk, Death and OSwk.
#' @param death_window_weeks Nonnegative window; equality counts as an event.
#'   Inf recovers the pre-181 rule for observed assessment times.
#' @param last_assessment Assessment-time column, TTPwk or LastEvalwk.
#' @param verbose Print aggregate rule counts only.
#' @return Data with raw ProgressionExit/TTPwk preserved, derived Progression,
#'   PFSwk, PFS_rule and a pfs_endpoint provenance attribute. Missing assessment
#'   times for non-progressors yield missing endpoints, even with window Inf.
derive_pfs_endpoint <- function(data, death_window_weeks = 16,
                                last_assessment = c("TTPwk", "LastEvalwk"),
                                verbose = TRUE) {
  fail <- function(...) stop("derive_pfs_endpoint(): ", ..., call. = FALSE)
  last_assessment <- tryCatch(match.arg(last_assessment),
                              error = function(e) fail(conditionMessage(e)))
  if (!is.numeric(death_window_weeks) || length(death_window_weeks) != 1L ||
      is.na(death_window_weeks) || death_window_weeks < 0) {
    fail("death_window_weeks must be a nonnegative numeric scalar (Inf allowed).")
  }
  required <- c("Progression", "PFSwk", "Death", "OSwk",
                if (last_assessment == "LastEvalwk") "LastEvalwk")
  missing <- setdiff(required, names(data))
  if (length(missing)) fail("missing column(s): ", paste(missing, collapse = ", "),
                            "; re-run analysis/01_data_prep.R.")
  if (!"ProgressionExit" %in% names(data)) data$ProgressionExit <- data$Progression
  if (!"TTPwk" %in% names(data)) data$TTPwk <- data$PFSwk
  for (nm in c("ProgressionExit", "Death")) {
    if (anyNA(data[[nm]]) || !all(data[[nm]] %in% c(0, 1)))
      fail(nm, " flags must be 0/1 without NA.")
  }
  for (nm in unique(c("TTPwk", "OSwk", last_assessment))) {
    if (!is.numeric(data[[nm]]) || any(!is.na(data[[nm]]) &
        (!is.finite(data[[nm]]) | data[[nm]] < 0)))
      fail(nm, " must contain nonnegative finite times or NA.")
    if (any(data[[nm]] > data$OSwk, na.rm = TRUE)) fail(nm, " exceeds OS time.")
  }
  progressed <- data$ProgressionExit == 1
  died <- data$Death == 1
  anchor <- data[[last_assessment]]
  anchor_missing <- !progressed & is.na(anchor)
  gap <- data$OSwk - anchor
  death_event <- died & !progressed & !is.na(anchor) &
    (is.infinite(death_window_weeks) | (!is.na(gap) & gap <= death_window_weeks))
  data$Progression <- as.numeric(progressed | death_event)
  data$PFSwk <- ifelse(progressed, data$TTPwk, ifelse(death_event, data$OSwk, anchor))
  data$Progression[anchor_missing] <- NA_real_
  rule <- ifelse(progressed, "progression", ifelse(death_event, "death_within_window",
    ifelse(died, "death_censored_at_assessment", "censored_alive")))
  rule[anchor_missing] <- "anchor_missing"
  data$PFS_rule <- factor(rule, levels = c("progression", "death_within_window",
    "death_censored_at_assessment", "censored_alive", "anchor_missing"))
  counts <- table(data$PFS_rule)
  attr(data, "pfs_endpoint") <- list(rule = "two_interval_v1",
    death_window_weeks = death_window_weeks, last_assessment = last_assessment,
    counts = counts)
  if (verbose) message("PFS endpoint: ", paste(names(counts), counts, collapse = "; "))
  data
}
