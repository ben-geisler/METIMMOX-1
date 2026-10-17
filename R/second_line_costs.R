# ===============================================================================
# SECOND-LINE THERAPY COST AFTER PROGRESSION (issue #209, part 1: inputs only)
# ===============================================================================
# data/external/second_line_pathway.csv holds the components of the Norwegian
# second-line pathway after first-line FLOX (trial team, 15 October 2026:
# "FOLFIRI + FLIRI, some bevacizumab or panitumumab", about 80% treated, no
# difference by arm): pathway shares, patient size, doses, public Norwegian list
# prices, vial sharing, outpatient contacts and blood tests per cycle and the
# number of cycles from the pivotal second-line trials, each with value, bounds,
# source and retrieval date. This file reads the table, checks its contract and
# computes the expected second-line cost per progressing patient. Contacts and
# blood tests are costed at the model's own unit costs (`c_other_visit`,
# `c_test_blood` of `l_params_base`), which the caller passes in; they are not
# repeated in the table.
#
# Part 2 of #209 wires the pathway into the base case: script 05 sources this
# file, turns the table into the scalar model inputs of l_params_base
# (second_line_model_inputs(): pathway shares, cycles, contacts, blood tests
# and the drug cost per cycle of each backbone and biologic) and adds the
# share of PFS events that are progressions (pfs_progression_share(), a count
# from the trial under the #181 PFS rule). The model itself evaluates the cost
# per progressing patient in closed form from those scalars
# (second_line_cost_per_progression() in calculate_outcomes.R), so the PSA
# and the DSA can vary them. The derivation is documented in the private
# repository (docs/manuscripts/2_cost_effectiveness/inputs/second_line_derivation.md).
# Only published aggregate figures and public list prices; no patient data.
# ===============================================================================

# Repository-relative path of the second-line pathway table.
SECOND_LINE_PATHWAY_PATH <- "data/external/second_line_pathway.csv"

# Columns the table must carry.
SECOND_LINE_COLUMNS <- c("component", "parameter", "regimen", "drug", "value", "lower", "upper",
                         "unit", "vial_mg", "administrations_per_cycle", "cost_parameter",
                         "source", "doi_or_url", "retrieved", "note")

# Chemotherapy backbones and the biologic add-ons of the pathway.
SECOND_LINE_BACKBONES <- c("FOLFIRI", "FLIRI")
SECOND_LINE_BIOLOGICS <- c("none", "bevacizumab", "panitumumab")

# Scalar parameters the table must define (component -> parameters).
SECOND_LINE_REQUIRED <- list(
  pathway = c("share_treated", "share_folfiri", "share_ras_braf_wt",
              "share_antiegfr_if_wt", "share_bevacizumab_if_no_antiegfr"),
  schedule = "cycle_length_weeks",
  patient = c("bsa", "weight"),
  dose_intensity = "relative_dose_intensity",
  currency = "nok_per_eur",
  wastage = "vial_sharing",
  resource = c("contacts_per_cycle_FOLFIRI", "contacts_per_cycle_FLIRI", "blood_tests_per_cycle"),
  duration = c("cycles_none", "cycles_bevacizumab", "cycles_panitumumab"))

# Allowed units per component (dose units are per administration).
SECOND_LINE_UNITS <- list(
  pathway = "proportion", schedule = "weeks", patient = c("m2", "kg"),
  dose = c("mg/m2", "mg/kg", "mg"), dose_intensity = "proportion",
  price = "NOK per vial", currency = "NOK per EUR", wastage = "proportion",
  resource = c("contacts per cycle", "tests per cycle"), duration = "cycles")

# Model unit costs a resource row may reference.
SECOND_LINE_COST_PARAMETERS <- c("c_other_visit", "c_test_blood")

#' Read the second-line pathway table and check its contract
#'
#' Contract: the required columns exist; every row has a component from
#' SECOND_LINE_UNITS, a unit allowed for it, a source, a DOI or URL and an ISO
#' retrieval date; parameter names are unique; value, lower and upper are
#' finite and ordered (lower <= value <= upper); proportions lie in [0, 1] and
#' every other quantity is positive (the vial-sharing fraction may be 0); every
#' parameter of SECOND_LINE_REQUIRED is present; dose rows name a regimen (a
#' backbone or a biologic), a drug and a positive whole number of
#' administrations per cycle; every dosed drug has at least one price row with
#' a positive vial size; resource rows reference a model unit cost of
#' SECOND_LINE_COST_PARAMETERS.
#'
#' @param path Path of the CSV; default SECOND_LINE_PATHWAY_PATH under the
#'   repository root.
#' @return Data frame of the table with numeric value, lower, upper, vial_mg
#'   and administrations_per_cycle, and attribute `path`.
read_second_line_pathway <- function(path = here::here(SECOND_LINE_PATHWAY_PATH)) {
  if (!file.exists(path)) {
    stop("File not found: ", path, ". The second-line pathway table (issue #209) ",
         "is required for the post-progression treatment cost.", call. = FALSE)
  }
  df_raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                            na.strings = c("", "NA"))
  v_missing <- setdiff(SECOND_LINE_COLUMNS, names(df_raw))
  if (length(v_missing)) {
    stop("second_line_pathway.csv lacks columns: ", paste(v_missing, collapse = ", "), call. = FALSE)
  }
  v_problems <- character(0)
  add <- function(message, rows) {
    if (length(rows)) v_problems <<- c(v_problems, paste0(message, ": ", paste(rows, collapse = ", ")))
  }
  for (column in c("value", "lower", "upper", "vial_mg", "administrations_per_cycle")) {
    df_raw[[column]] <- suppressWarnings(as.numeric(df_raw[[column]]))
  }
  blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))
  add("rows without a parameter name", which(blank(df_raw$parameter)))
  add("duplicated parameters", unique(df_raw$parameter[duplicated(df_raw$parameter)]))
  add("rows without a source", df_raw$parameter[blank(df_raw$source)])
  add("rows without a DOI or URL", df_raw$parameter[blank(df_raw$doi_or_url)])
  v_date <- suppressWarnings(as.Date(as.character(df_raw$retrieved), format = "%Y-%m-%d"))
  add("rows without an ISO retrieval date (YYYY-MM-DD)", df_raw$parameter[is.na(v_date)])
  v_unknown <- !df_raw$component %in% names(SECOND_LINE_UNITS)
  add("unknown components", unique(df_raw$component[v_unknown]))
  v_bad_unit <- !v_unknown & !mapply(function(component, unit) unit %in% SECOND_LINE_UNITS[[component]],
                                     df_raw$component, df_raw$unit)
  add("units not allowed for their component", df_raw$parameter[v_bad_unit])
  v_finite <- is.finite(df_raw$value) & is.finite(df_raw$lower) & is.finite(df_raw$upper)
  add("value, lower or upper not finite", df_raw$parameter[!v_finite])
  add("bounds not ordered (lower <= value <= upper)",
      df_raw$parameter[v_finite & !(df_raw$lower <= df_raw$value & df_raw$value <= df_raw$upper)])
  v_proportion <- df_raw$unit %in% "proportion"
  add("proportions outside [0, 1]",
      df_raw$parameter[v_finite & v_proportion & (df_raw$lower < 0 | df_raw$upper > 1)])
  v_zero_allowed <- df_raw$parameter %in% "vial_sharing" | v_proportion
  add("quantities that must be positive",
      df_raw$parameter[v_finite & !v_zero_allowed & df_raw$lower <= 0])
  for (component in names(SECOND_LINE_REQUIRED)) {
    add(paste0("missing ", component, " parameters"),
        setdiff(SECOND_LINE_REQUIRED[[component]], df_raw$parameter[df_raw$component == component]))
  }
  df_dose <- df_raw[df_raw$component %in% "dose", , drop = FALSE]
  add("dose rows without a backbone or biologic regimen",
      df_dose$parameter[!df_dose$regimen %in% c(SECOND_LINE_BACKBONES, SECOND_LINE_BIOLOGICS[-1])])
  add("dose rows without a drug", df_dose$parameter[blank(df_dose$drug)])
  v_admin <- df_dose$administrations_per_cycle
  add("dose rows without a positive whole number of administrations per cycle",
      df_dose$parameter[!is.finite(v_admin) | v_admin < 1 | v_admin != round(v_admin)])
  df_price <- df_raw[df_raw$component %in% "price", , drop = FALSE]
  add("price rows without a drug", df_price$parameter[blank(df_price$drug)])
  add("price rows without a positive vial size",
      df_price$parameter[!is.finite(df_price$vial_mg) | df_price$vial_mg <= 0])
  add("dosed drugs without a price row", setdiff(unique(df_dose$drug), df_price$drug))
  df_resource <- df_raw[df_raw$component %in% "resource", , drop = FALSE]
  add("resource rows without a model unit cost (cost_parameter)",
      df_resource$parameter[!df_resource$cost_parameter %in% SECOND_LINE_COST_PARAMETERS])
  if (length(v_problems)) {
    stop("second_line_pathway.csv violates its contract:\n",
         paste0("  - ", v_problems, collapse = "\n"), call. = FALSE)
  }
  attr(df_raw, "path") <- path
  df_raw
}

#' Parameters of the pathway that have a range
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @return Character vector of the parameters with lower < upper.
second_line_ranged_parameters <- function(df_pathway) {
  df_pathway$parameter[df_pathway$lower < df_pathway$upper]
}

#' Parameter values of one setting of the pathway
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param setting Named character vector: parameter -> "value", "lower" or
#'   "upper"; parameters not named are at "value".
#' @return Named numeric vector, one entry per parameter of the table.
second_line_values <- function(df_pathway, setting = character(0)) {
  v_unknown <- setdiff(names(setting), df_pathway$parameter)
  if (length(v_unknown)) stop("Unknown second-line parameters: ", paste(v_unknown, collapse = ", "))
  v_out <- stats::setNames(df_pathway$value, df_pathway$parameter)
  for (name in names(setting)) {
    bound <- match.arg(setting[[name]], c("value", "lower", "upper"))
    v_out[[name]] <- df_pathway[[bound]][df_pathway$parameter == name]
  }
  v_out
}

#' Cost of one administered dose of a drug in NOK
#'
#' Without vial sharing the dose is made up from whole vials at the lowest
#' cost (the opened remainder is wasted); with full sharing only the
#' milligrams given are paid, at the lowest price per milligram. A fraction
#' `sharing` between 0 and 1 mixes the two.
#'
#' @param dose_mg Dose in mg.
#' @param vial_mg Vial sizes in mg.
#' @param vial_price Price per vial in NOK, aligned with `vial_mg`.
#' @param sharing Fraction of the wastage avoided by vial sharing (0 to 1).
#' @return Numeric scalar, NOK.
second_line_vial_cost <- function(dose_mg, vial_mg, vial_price, sharing = 0) {
  stopifnot(length(dose_mg) == 1L, is.finite(dose_mg), dose_mg >= 0,
            length(vial_mg) == length(vial_price), length(vial_mg) >= 1L,
            all(vial_mg > 0), all(vial_price >= 0), sharing >= 0, sharing <= 1)
  v_order <- order(vial_mg, decreasing = TRUE)
  v_size <- vial_mg[v_order]
  v_price <- vial_price[v_order]
  # Lowest whole-vial cost covering `remaining` with vials k..n.
  cover <- function(remaining, k) {
    if (remaining <= 1e-9) return(0)
    if (k == length(v_size)) return(ceiling(remaining / v_size[k] - 1e-9) * v_price[k])
    v_counts <- 0:ceiling(remaining / v_size[k] - 1e-9)
    min(vapply(v_counts, function(count)
      count * v_price[k] + cover(remaining - count * v_size[k], k + 1L), numeric(1)))
  }
  whole <- cover(dose_mg, 1L)
  exact <- dose_mg * min(v_price / v_size)
  (1 - sharing) * whole + sharing * exact
}

#' Drug cost per cycle in EUR of one backbone or biologic
#'
#' Sum over the regimen's dose rows of administrations per cycle x the cost
#' of one dose (dose x BSA, weight or flat x relative dose intensity, made up
#' from vials by second_line_vial_cost() with the `vial_sharing` fraction),
#' converted at `nok_per_eur`.
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param regimen A backbone (SECOND_LINE_BACKBONES) or a biologic
#'   ("bevacizumab", "panitumumab").
#' @param values Named numeric vector from second_line_values().
#' @return Numeric scalar, EUR per cycle.
second_line_drug_cost_cycle <- function(df_pathway, regimen, values = second_line_values(df_pathway)) {
  df_dose <- df_pathway[df_pathway$component == "dose", , drop = FALSE]
  df_price <- df_pathway[df_pathway$component == "price", , drop = FALSE]
  df_rows <- df_dose[df_dose$regimen == regimen, , drop = FALSE]
  if (!nrow(df_rows)) stop("No dose rows for second-line regimen '", regimen, "'")
  sum(vapply(seq_len(nrow(df_rows)), function(i) {
    basis <- switch(df_rows$unit[i], "mg/m2" = values[["bsa"]], "mg/kg" = values[["weight"]], "mg" = 1)
    dose_mg <- values[[df_rows$parameter[i]]] * basis * values[["relative_dose_intensity"]]
    df_vials <- df_price[df_price$drug == df_rows$drug[i], , drop = FALSE]
    v_vial_price <- values[df_vials$parameter]
    df_rows$administrations_per_cycle[i] *
      second_line_vial_cost(dose_mg, df_vials$vial_mg, v_vial_price, values[["vial_sharing"]])
  }, numeric(1))) / values[["nok_per_eur"]]
}

#' Second-line regimens with their shares, cost per cycle and cost per course
#'
#' Regimen shares among treated patients: a backbone (FOLFIRI with share
#' `share_folfiri`, Nordic FLIRI otherwise) combined with panitumumab for the
#' RAS/BRAF wild-type patients given anti-EGFR therapy
#' (`share_ras_braf_wt` x `share_antiegfr_if_wt`), with bevacizumab for a share
#' `share_bevacizumab_if_no_antiegfr` of the others, and chemotherapy alone for
#' the rest. Cost per cycle = drugs (doses x relative dose intensity, priced by
#' second_line_vial_cost() and converted at `nok_per_eur`) + contacts x
#' `c_other_visit` + blood tests x `c_test_blood`. Cost per course = cycles x
#' cost per cycle, the cycles depending on the biologic.
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param params List with the model unit costs `c_other_visit` and
#'   `c_test_blood` in EUR (e.g. `l_params_base`).
#' @param values Named numeric vector from second_line_values(); default all
#'   parameters at their value.
#' @return Data frame, one row per backbone x biologic: `regimen`, `backbone`,
#'   `biologic`, `share` (among treated), `drug_cost_cycle`,
#'   `contact_cost_cycle`, `test_cost_cycle`, `cost_cycle`, `cycles`,
#'   `weeks`, `cost_course` (EUR).
second_line_regimens <- function(df_pathway, params, values = second_line_values(df_pathway)) {
  v_cost <- vapply(SECOND_LINE_COST_PARAMETERS, function(name) {
    x <- params[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
      stop("params$", name, " must be a finite nonnegative scalar (the model's unit cost)")
    }
    x
  }, numeric(1))
  drug_cost <- function(regimen) second_line_drug_cost_cycle(df_pathway, regimen, values)
  p_egfr <- values[["share_ras_braf_wt"]] * values[["share_antiegfr_if_wt"]]
  p_bev <- (1 - p_egfr) * values[["share_bevacizumab_if_no_antiegfr"]]
  v_biologic_share <- c(none = 1 - p_egfr - p_bev, bevacizumab = p_bev, panitumumab = p_egfr)
  v_backbone_share <- c(FOLFIRI = values[["share_folfiri"]], FLIRI = 1 - values[["share_folfiri"]])
  df_out <- expand.grid(biologic = SECOND_LINE_BIOLOGICS, backbone = SECOND_LINE_BACKBONES,
                        stringsAsFactors = FALSE)[, c("backbone", "biologic")]
  df_out$regimen <- ifelse(df_out$biologic == "none", df_out$backbone,
                           paste0(df_out$backbone, " + ", df_out$biologic))
  df_out$share <- unname(v_backbone_share[df_out$backbone] * v_biologic_share[df_out$biologic])
  df_out$drug_cost_cycle <- vapply(seq_len(nrow(df_out)), function(i)
    drug_cost(df_out$backbone[i]) +
      if (df_out$biologic[i] == "none") 0 else drug_cost(df_out$biologic[i]), numeric(1))
  df_out$contact_cost_cycle <- unname(values[paste0("contacts_per_cycle_", df_out$backbone)] *
                                        v_cost[["c_other_visit"]])
  df_out$test_cost_cycle <- values[["blood_tests_per_cycle"]] * v_cost[["c_test_blood"]]
  df_out$cost_cycle <- df_out$drug_cost_cycle + df_out$contact_cost_cycle + df_out$test_cost_cycle
  df_out$cycles <- unname(values[paste0("cycles_", df_out$biologic)])
  df_out$weeks <- df_out$cycles * values[["cycle_length_weeks"]]
  df_out$cost_course <- df_out$cycles * df_out$cost_cycle
  df_out[c("regimen", "backbone", "biologic", "share", "drug_cost_cycle", "contact_cost_cycle",
           "test_cost_cycle", "cost_cycle", "cycles", "weeks", "cost_course")]
}

#' Expected second-line cost of one setting of the pathway
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param params List with `c_other_visit` and `c_test_blood` (EUR).
#' @param values Named numeric vector from second_line_values().
#' @return List with `per_progressing` (EUR per patient who progresses:
#'   share treated x mean course cost), `per_treated` (mean course cost among
#'   treated patients), `weeks_treated` (mean weeks on treatment among treated
#'   patients), `per_week_treated` (per_treated / weeks_treated) and
#'   `regimens` (second_line_regimens()).
second_line_cost_at <- function(df_pathway, params, values = second_line_values(df_pathway)) {
  df_regimens <- second_line_regimens(df_pathway, params, values)
  per_treated <- sum(df_regimens$share * df_regimens$cost_course)
  weeks <- sum(df_regimens$share * df_regimens$weeks)
  list(per_progressing = values[["share_treated"]] * per_treated, per_treated = per_treated,
       weeks_treated = weeks, per_week_treated = per_treated / weeks, regimens = df_regimens)
}

#' Expected second-line cost per progressing patient with its low and high
#'
#' The base case puts every parameter at its value. The low (high) scenario
#' puts every ranged parameter at the bound that minimises (maximises) the cost
#' per progressing patient: a coordinate search over the bounds, repeated
#' until no parameter changes. The cost is monotone in every dose, price,
#' size, contact and cycle parameter and linear in each share given the
#' others, so the search ends at a vertex of the parameter box. Like the
#' eligible-population bounds of voi_incidence.csv, low and high are extreme
#' scenarios, not confidence intervals. The model unit costs in `params` are
#' held at their values.
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param params List with `c_other_visit` and `c_test_blood` (EUR).
#' @param max_iter Maximum number of passes of the coordinate search.
#' @return List with `base`, `low`, `high` (each second_line_cost_at() plus
#'   `setting`, the bound chosen per ranged parameter) and `summary`, a data
#'   frame with rows base/low/high and columns `per_progressing`,
#'   `per_treated`, `weeks_treated`, `per_week_treated`.
second_line_cost <- function(df_pathway, params, max_iter = 20L) {
  v_ranged <- second_line_ranged_parameters(df_pathway)
  #' Coordinate search for the extreme in direction `sign` (1 = high).
  extreme <- function(sign) {
    setting <- stats::setNames(rep("value", length(v_ranged)), v_ranged)
    current <- second_line_cost_at(df_pathway, params, second_line_values(df_pathway, setting))
    for (iter in seq_len(max_iter)) {
      changed <- FALSE
      for (name in v_ranged) {
        for (bound in c("lower", "upper")) {
          if (setting[[name]] == bound) next
          candidate <- setting
          candidate[[name]] <- bound
          result <- second_line_cost_at(df_pathway, params, second_line_values(df_pathway, candidate))
          if (sign * (result$per_progressing - current$per_progressing) > 1e-9) {
            setting <- candidate
            current <- result
            changed <- TRUE
          }
        }
      }
      if (!changed) break
    }
    if (changed) stop("second_line_cost(): coordinate search did not settle in ", max_iter, " passes")
    c(current, list(setting = setting))
  }
  base <- c(second_line_cost_at(df_pathway, params), list(setting = character(0)))
  low <- extreme(-1)
  high <- extreme(1)
  v_columns <- c("per_progressing", "per_treated", "weeks_treated", "per_week_treated")
  summary <- do.call(rbind, lapply(list(base = base, low = low, high = high), function(x)
    as.data.frame(x[v_columns])))
  list(base = base, low = low, high = high, summary = summary)
}

#' One-way ranges of the expected second-line cost per progressing patient
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @param params List with `c_other_visit` and `c_test_blood` (EUR).
#' @return Data frame, one row per ranged parameter: `parameter`,
#'   `component`, `lower`, `upper` (the parameter bounds), `cost_lower`,
#'   `cost_upper` (EUR per progressing patient with that parameter at its
#'   bound and every other at its value) and `swing`, sorted by swing.
second_line_one_way <- function(df_pathway, params) {
  v_ranged <- second_line_ranged_parameters(df_pathway)
  df_out <- do.call(rbind, lapply(v_ranged, function(name) {
    row <- df_pathway[df_pathway$parameter == name, ]
    v_cost <- vapply(c("lower", "upper"), function(bound)
      second_line_cost_at(df_pathway, params,
                          second_line_values(df_pathway, stats::setNames(bound, name)))$per_progressing,
      numeric(1))
    data.frame(parameter = name, component = row$component, lower = row$lower, upper = row$upper,
               cost_lower = v_cost[["lower"]], cost_upper = v_cost[["upper"]],
               swing = abs(v_cost[["upper"]] - v_cost[["lower"]]), stringsAsFactors = FALSE)
  }))
  df_out <- df_out[order(-df_out$swing), , drop = FALSE]
  rownames(df_out) <- NULL
  df_out
}

# ===============================================================================
# MODEL INPUTS (issue #209, part 2)
# ===============================================================================

# Table parameters that enter l_params_base, by model field. Shares and cycles
# are sampled in the PSA (parameter_distribution_spec(), group second_line);
# contacts and blood tests are fixed in the PSA and varied in the DSA.
SECOND_LINE_TABLE_FIELDS <- c(
  p_sl_treated = "share_treated",
  p_sl_folfiri = "share_folfiri",
  p_sl_ras_braf_wt = "share_ras_braf_wt",
  p_sl_antiegfr_if_wt = "share_antiegfr_if_wt",
  p_sl_bevacizumab_if_no_antiegfr = "share_bevacizumab_if_no_antiegfr",
  n_sl_cycles_none = "cycles_none",
  n_sl_cycles_bevacizumab = "cycles_bevacizumab",
  n_sl_cycles_panitumumab = "cycles_panitumumab",
  n_sl_contacts_FOLFIRI = "contacts_per_cycle_FOLFIRI",
  n_sl_contacts_FLIRI = "contacts_per_cycle_FLIRI",
  n_sl_blood_tests = "blood_tests_per_cycle")

# Drug cost per cycle (EUR) of each backbone and biologic, by model field.
SECOND_LINE_DRUG_FIELDS <- c(
  c_sl_drug_FOLFIRI = "FOLFIRI", c_sl_drug_FLIRI = "FLIRI",
  c_sl_drug_bevacizumab = "bevacizumab", c_sl_drug_panitumumab = "panitumumab")

#' Second-line model inputs from the pathway table
#'
#' Values are the table's base values; the drug cost per cycle of each
#' backbone and biologic is second_line_drug_cost_cycle() at those values
#' (doses, patient size, relative dose intensity, vial sharing and prices are
#' folded into it, so the model carries four drug costs instead of 30 price
#' and dose rows). Bounds are the table's lower and upper where it gives a
#' range (lower < upper), else NA (the DSA then uses +/- DSA_mult). The
#' standard error of a sampled parameter treats its bounds as an approximate
#' 95% interval, SE = (upper - lower) / (2 x 1.96): an ASSUMPTION, since the
#' bounds are sourced extremes (medians, other trials), not confidence
#' intervals.
#'
#' @param df_pathway Output of read_second_line_pathway().
#' @return List with named numeric vectors `values`, `lower`, `upper` (all
#'   over the fields of SECOND_LINE_TABLE_FIELDS and SECOND_LINE_DRUG_FIELDS)
#'   and `se` (sampled fields only: shares and cycles).
second_line_model_inputs <- function(df_pathway) {
  v_rows <- match(SECOND_LINE_TABLE_FIELDS, df_pathway$parameter)
  if (anyNA(v_rows)) {
    stop("Pathway table lacks: ", paste(SECOND_LINE_TABLE_FIELDS[is.na(v_rows)], collapse = ", "))
  }
  v_fields <- names(SECOND_LINE_TABLE_FIELDS)
  v_values <- stats::setNames(df_pathway$value[v_rows], v_fields)
  v_lower <- stats::setNames(df_pathway$lower[v_rows], v_fields)
  v_upper <- stats::setNames(df_pathway$upper[v_rows], v_fields)
  v_no_range <- !(v_lower < v_upper)
  v_lower[v_no_range] <- NA_real_
  v_upper[v_no_range] <- NA_real_
  values <- second_line_values(df_pathway)
  v_drug <- vapply(SECOND_LINE_DRUG_FIELDS, function(regimen)
    second_line_drug_cost_cycle(df_pathway, regimen, values), numeric(1))
  v_na <- stats::setNames(rep(NA_real_, length(v_drug)), names(v_drug))
  v_sampled <- grep("^(p_sl_|n_sl_cycles_)", v_fields, value = TRUE)
  v_se <- (v_upper[v_sampled] - v_lower[v_sampled]) / (2 * stats::qnorm(0.975))
  if (anyNA(v_se) || any(v_se <= 0)) {
    stop("Every sampled second-line parameter needs a range (lower < upper) in the table.")
  }
  list(values = c(v_values, v_drug), lower = c(v_lower, v_na), upper = c(v_upper, v_na), se = v_se)
}

#' Share of PFS events that are progressions (issue #209)
#'
#' A partitioned survival model sees PFS events, which are progressions and
#' deaths without progression; only progressions start second-line therapy.
#' The share is the count of PFS events classified "progression" over all PFS
#' events ("progression" + "death_within_window") under the #181 rule
#' (derive_pfs_endpoint(): `PFS_rule`), pooled over both arms. Its PSA
#' distribution is Beta(n_progression, n_death): the standard error returned
#' gives exactly these shapes under the method of moments of
#' distribution_parameters(). DSA bounds: the 2.5% quantile of that beta and 1
#' (every PFS event costed).
#'
#' @param data Patient data with `PFS_rule` (the economic cohort,
#'   data_complete).
#' @return List with `value`, `n_progression`, `n_death`, `n_events`, `se`,
#'   `lower`, `upper` and `by_arm` (a progression/death x arm count table
#'   when `Rx` is present, else NULL).
pfs_progression_share <- function(data) {
  if (!"PFS_rule" %in% names(data)) {
    stop("pfs_progression_share() needs PFS_rule (derive_pfs_endpoint()).")
  }
  v_rule <- as.character(data$PFS_rule)
  n_progression <- sum(v_rule == "progression", na.rm = TRUE)
  n_death <- sum(v_rule == "death_within_window", na.rm = TRUE)
  n_events <- n_progression + n_death
  if (n_progression == 0 || n_death == 0) {
    stop("pfs_progression_share() needs at least one progression and one death without ",
         "progression for a proper beta distribution (got ", n_progression, " and ", n_death, ").")
  }
  value <- n_progression / n_events
  by_arm <- NULL
  if ("Rx" %in% names(data)) {
    by_arm <- table(event = factor(v_rule, levels = c("progression", "death_within_window")),
                    arm = data$Rx)
  }
  list(value = value, n_progression = n_progression, n_death = n_death, n_events = n_events,
       se = sqrt(value * (1 - value) / (n_events + 1)),
       lower = stats::qbeta(0.025, n_progression, n_death), upper = 1, by_arm = by_arm)
}
