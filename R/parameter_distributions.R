# ===============================================================================
# NON-SURVIVAL UNCERTAINTY PARAMETERS
# ===============================================================================
# The table below is the single source of truth for PSA, DSA, and EVPPI
# parameter membership. Add a parameter once here; downstream analyses derive
# their parameter vectors from the specification.
#
# Columns (issue #154):
#   psa      TRUE if the parameter is sampled in (or derived within) the PSA.
#   dsa      TRUE if the parameter is varied in the one-way DSA. Published
#            tariffs and list prices are psa = FALSE, dsa = TRUE: their value is
#            known, so sampling them adds decision uncertainty that does not
#            exist and puts "value of research on the FLOX tariff" at the top of
#            the EVPPI ranking, but their level still matters and is tested
#            deterministically.
#   derived  TRUE if the PSA column is a deterministic function of other sampled
#            parameters (apply_derived_psa_parameters()) rather than a draw.
#            Derived parameters are model inputs, so they take part in the DSA,
#            but they are excluded from EVPPI groups because they are exactly
#            collinear with the parameters they are derived from.

# Later treatment sequences (issue #208): the share of patients assigned to an
# arm who start a third sequence and, among them, the share who start a fourth.
# Model inputs are these four shares (later_sequence_resource_use() in
# calculate_outcomes.R reads them by this naming convention); the share
# starting a fourth sequence is p_seq3_<arm> * p_seq4_given_seq3_<arm>, so it
# can never exceed the share starting a third.
LATER_SEQUENCE_SHARE_PARAMETERS <- c("p_seq3_control", "p_seq3_experimental",
                                     "p_seq4_given_seq3_control",
                                     "p_seq4_given_seq3_experimental")

#' Beta shapes of a later-sequence share from the trial counts (issue #208)
#'
#' beta(n_seq3, n_patients - n_seq3) for the third sequence and
#' beta(n_seq4, n_seq3 - n_seq4) for the fourth given the third: by the
#' aggregation property these are the two independent factors of the
#' Dirichlet(n_none, n_seq3_only, n_seq4) posterior of the three outcomes (no
#' third sequence, a third only, a third and a fourth), the same raw-count
#' parameterisation as the joint biomarker masses. A zero count gives a point
#' mass (rbeta() with a zero shape), and beta(0, 0) is set to a point mass at
#' zero.
#'
#' @param parameter One of LATER_SEQUENCE_SHARE_PARAMETERS.
#' @param counts Data frame with columns arm, n_patients, n_seq3, n_seq4
#'   (params$later_sequences$counts).
#' @return Named numeric vector c(shape1, shape2).
later_sequence_beta_shapes <- function(parameter, counts) {
  arm <- sub("^p_seq(3|4_given_seq3)_", "", parameter)
  row <- counts[counts$arm == arm, , drop = FALSE]
  if (nrow(row) != 1L) stop("No later-sequence counts for arm '", arm, "'.")
  v_shapes <- if (grepl("^p_seq3_", parameter)) {
    c(row$n_seq3, row$n_patients - row$n_seq3)
  } else {
    c(row$n_seq4, row$n_seq3 - row$n_seq4)
  }
  if (any(!is.finite(v_shapes)) || any(v_shapes < 0))
    stop("Invalid later-sequence counts for ", parameter, ".")
  if (all(v_shapes == 0)) v_shapes <- c(0, 1)
  c(shape1 = v_shapes[1], shape2 = v_shapes[2])
}

parameter_distribution_spec <- function() {
  v_biomarker_cost_parameters <- unique(unname(biomarker_cost_key()))
  v_prevalence_parameters <- unname(biomarker_prevalence_key())

  # Published list prices and laboratory/procedure tariffs: fixed in the PSA,
  # varied in the one-way DSA only.
  df_fixed_costs <- data.frame(
    parameter = c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                  v_biomarker_cost_parameters),
    distribution = NA_character_,
    cv = NA_real_,
    group = c("drug_costs", "drug_costs", "test_costs", "test_costs",
              rep("test_costs", length(v_biomarker_cost_parameters))),
    psa = FALSE,
    order = c(10:11, 20:21, 22 + seq_along(v_biomarker_cost_parameters)),
    stringsAsFactors = FALSE
  )
  # Visit, baseline, follow-up and end-of-life costs bundle genuine resource-use
  # uncertainty (contact frequency, terminal-care intensity), so they remain
  # probabilistic.
  # CV 0.20 is an ASSUMPTION (issue #56), not an estimate: no source reports
  # the dispersion of these resource-use costs, so a conventional moderate CV
  # is used (gamma 95% interval about 0.65-1.43 times the mean). It is not
  # derived from data or cited literature and is not varied in any scenario.
  df_other_costs <- data.frame(
    parameter = c("c_other_visit", "c_other_baseline", "c_other_follow",
                  "c_other_last"),
    distribution = "gamma",
    cv = 0.20,
    group = "other_costs",
    psa = TRUE,
    order = 30:33,
    stringsAsFactors = FALSE
  )
  # u_p is not sampled directly. Independent beta draws for u_np and u_p put
  # 16% of draws (48% under UTILITY_SOURCE = 0) in the implausible region
  # u_p > u_np. Sampling a non-negative decrement instead makes u_p <= u_np hold
  # by construction, and E[u_decrement] = u_np - u_p leaves the marginal mean of
  # u_p unchanged.
  # CV 0.15 for u_np and u_decrement is an ASSUMPTION (issue #56): the CORRECT
  # utilities come without a usable standard error, so a conventional CV is
  # used (beta 95% interval for u_np = 0.73 about 0.49-0.91). It is not derived
  # from data or cited literature and is not varied in any scenario. In IPD
  # mode (UTILITY_SOURCE = 0, issue #188) patient-clustered bootstrap SEs in
  # params$utility_se replace it (create_parameter_distributions()).
  df_utilities <- data.frame(
    parameter = c("u_np", "u_decrement", "u_p"),
    distribution = c("beta", "gamma", "derived"),
    cv = c(0.15, 0.15, NA_real_),
    group = "utilities",
    psa = TRUE,
    order = 40:42,
    stringsAsFactors = FALSE
  )
  df_prevalence <- data.frame(
    parameter = c(paste0("p_joint_", c("00", "01", "10", "11")), v_prevalence_parameters),
    distribution = c(rep("dirichlet", 3), rep("derived", 1 + length(v_prevalence_parameters))),
    cv = NA_real_,
    group = "prevalence",
    psa = TRUE,
    order = 50 + seq_len(4 + length(v_prevalence_parameters)),
    stringsAsFactors = FALSE
  )

  # Shares of patients starting a third and a fourth treatment sequence, by
  # arm, from the trial's administration records (issue #208). Sampled after
  # every earlier parameter, so those keep their draws.
  df_later_sequences <- data.frame(
    parameter = LATER_SEQUENCE_SHARE_PARAMETERS,
    distribution = "beta_count",
    cv = NA_real_,
    group = "later_sequences",
    psa = TRUE,
    order = 60 + seq_along(LATER_SEQUENCE_SHARE_PARAMETERS),
    stringsAsFactors = FALSE
  )

  # Second-line therapy at incident progression (issue #209). The shares
  # (beta) and the course lengths in cycles (gamma) are sampled with the
  # standard errors that script 05 puts in params$psa_se: from the pathway
  # table's bounds read as an approximate 95% interval (an ASSUMPTION,
  # second_line_model_inputs()), and for the progression share of PFS events
  # from the trial count, Beta(n_progression, n_death). The drug cost per
  # cycle of each backbone and biologic is a list-price product and is fixed
  # in the PSA (#154 rule); contacts and blood tests per cycle are fixed too.
  # Rows come after every earlier row, including the later-sequence shares
  # of issue #208, so the draws of the earlier parameters are unchanged by
  # their addition (generate_psa_samples() draws in spec order).
  v_second_line_sampled <- c("p_pfs_progression", "p_sl_treated", "p_sl_folfiri",
                             "p_sl_ras_braf_wt", "p_sl_antiegfr_if_wt",
                             "p_sl_bevacizumab_if_no_antiegfr", "n_sl_cycles_none",
                             "n_sl_cycles_bevacizumab", "n_sl_cycles_panitumumab")
  v_second_line_fixed <- c("n_sl_contacts_FOLFIRI", "n_sl_contacts_FLIRI", "n_sl_blood_tests",
                           "c_sl_drug_FOLFIRI", "c_sl_drug_FLIRI", "c_sl_drug_bevacizumab",
                           "c_sl_drug_panitumumab")
  df_second_line <- data.frame(
    parameter = c(v_second_line_sampled, v_second_line_fixed),
    distribution = c(ifelse(grepl("^p_", v_second_line_sampled), "beta", "gamma"),
                     rep(NA_character_, length(v_second_line_fixed))),
    cv = NA_real_,
    group = "second_line",
    psa = c(rep(TRUE, length(v_second_line_sampled)), rep(FALSE, length(v_second_line_fixed))),
    order = 70 + seq_len(length(v_second_line_sampled) + length(v_second_line_fixed)),
    stringsAsFactors = FALSE
  )

  df_spec <- rbind(df_fixed_costs, df_other_costs, df_utilities, df_prevalence,
                   df_later_sequences, df_second_line)
  df_spec$derived <- df_spec$distribution %in% "derived"
  # DSA varies marginal prevalences, not individual joint masses. The sampled
  # utility decrement is also excluded because DSA varies state utilities.
  df_spec$dsa <- df_spec$parameter != "u_decrement" & !grepl("^p_joint_", df_spec$parameter)
  df_spec <- df_spec[order(df_spec$order), , drop = FALSE]
  df_spec$order <- NULL
  rownames(df_spec) <- NULL
  if (anyDuplicated(df_spec$parameter)) {
    stop("Parameter distribution specification contains duplicate names.")
  }
  df_spec
}

#' Look up the reporting group of any specified parameter
#'
#' Covers fixed as well as sampled parameters, so deterministic-analysis output
#' can be grouped even for parameters that have no PSA distribution.
#'
#' @param spec Parameter specification.
#' @return Named character vector mapping parameter to group.
parameter_group_lookup <- function(spec = parameter_distribution_spec()) {
  stats::setNames(spec$group, spec$parameter)
}

distribution_parameters <- function(mean_value, distribution, cv) {
  if (identical(distribution, "derived")) {
    return(list(dist = "derived"))
  }
  if (!is.finite(mean_value) || mean_value < 0) {
    stop("Distribution means must be finite and non-negative.")
  }
  if (identical(distribution, "gamma")) {
    shape <- 1 / cv^2
    return(list(dist = "gamma", shape = shape, rate = shape / mean_value))
  }
  if (identical(distribution, "beta")) {
    if (mean_value <= 0 || mean_value >= 1) {
      stop("Beta means must lie strictly between zero and one.")
    }
    variance <- (mean_value * cv)^2
    common <- mean_value * (1 - mean_value) / variance - 1
    if (common <= 0) {
      stop("Beta means and CVs must define positive shape parameters.")
    }
    return(list(
      dist = "beta",
      shape1 = mean_value * common,
      shape2 = (1 - mean_value) * common
    ))
  }
  stop("Unsupported parameter distribution: ", distribution)
}

#' Specification rows that apply to a parameter list
#'
#' A list without any second-line field (test fixtures, lists built before
#' issue #209) costs no second-line therapy (second_line_costed() in
#' calculate_outcomes.R), so the second_line rows are dropped for it. A list
#' with some but not all of them keeps the rows, and the missing-parameter
#' checks of the callers stop.
#'
#' @param spec Parameter specification.
#' @param params Parameter list.
#' @return The specification, without the second_line rows when `params` has
#'   none of their fields.
applicable_parameter_spec <- function(spec, params) {
  v_second_line <- spec$group %in% "second_line"
  if (any(v_second_line) && !any(spec$parameter[v_second_line] %in% names(params))) {
    spec <- spec[!v_second_line, , drop = FALSE]
  }
  spec
}

create_parameter_distributions <- function(params,
                                           spec = parameter_distribution_spec()) {
  spec <- applicable_parameter_spec(spec, params)
  spec <- spec[spec$psa, , drop = FALSE]
  v_missing_parameters <- setdiff(spec$parameter, names(params))
  if (length(v_missing_parameters) > 0) {
    stop("Missing base-case parameters: ",
         paste(v_missing_parameters, collapse = ", "))
  }

  distributions <- lapply(seq_len(nrow(spec)), function(i) {
    if (spec$distribution[i] == "dirichlet") {
      v_counts <- params$joint_counts
      if (is.null(v_counts) || !identical(names(v_counts), c("00", "01", "10", "11")) ||
          any(!is.finite(v_counts)) || any(v_counts < 0) || sum(v_counts) <= 0)
        stop("Joint prevalence uncertainty requires complete-case joint cell counts.")
      return(list(dist = "dirichlet", alpha = v_counts,
                  component = sub("^p_joint_", "", spec$parameter[i])))
    }
    if (spec$distribution[i] == "beta_count") {
      v_counts <- params$later_sequences$counts
      if (is.null(v_counts)) stop("Later-sequence shares require params$later_sequences$counts.")
      v_shapes <- later_sequence_beta_shapes(spec$parameter[i], v_counts)
      return(list(dist = "beta", shape1 = unname(v_shapes[1]), shape2 = unname(v_shapes[2])))
    }
    # In IPD utility mode (UTILITY_SOURCE = 0, issue #188) script 05 supplies
    # patient-clustered bootstrap standard errors in params$utility_se; they
    # replace the assumed CV for the parameters they name. The literature base
    # case has no utility_se and keeps the spec's CV.
    cv <- spec$cv[i]
    se <- params$utility_se[spec$parameter[i]]
    if (length(se) == 1L && !is.na(se)) {
      if (!is.finite(se) || se <= 0) stop("utility_se must be positive and finite.")
      cv <- unname(se / params[[spec$parameter[i]]])
    }
    # Standard errors supplied with the base case (issue #209: the second-line
    # rows, whose spec CV is NA) replace the CV in the same way.
    se <- params$psa_se[spec$parameter[i]]
    if (length(se) == 1L && !is.na(se)) {
      if (!is.finite(se) || se <= 0) stop("psa_se must be positive and finite.")
      cv <- unname(se / params[[spec$parameter[i]]])
    }
    if (!identical(spec$distribution[i], "derived") && is.na(cv)) {
      stop("No dispersion for PSA parameter '", spec$parameter[i],
           "': the spec gives no CV and params$psa_se has no standard error.")
    }
    distribution_parameters(params[[spec$parameter[i]]], spec$distribution[i], cv)
  })
  stats::setNames(distributions, spec$parameter)
}

create_parameter_groups <- function(spec = parameter_distribution_spec()) {
  # EVPPI groups describe what the PSA actually samples. Fixed prices carry no
  # decision uncertainty, and derived columns are collinear with their inputs.
  spec <- spec[spec$psa & !spec$derived, , drop = FALSE]
  grouped <- split(spec$parameter, factor(
    spec$group,
    levels = c("drug_costs", "test_costs", "other_costs", "utilities", "prevalence",
               "later_sequences", "second_line")
  ))
  grouped <- grouped[lengths(grouped) > 0]

  v_cost_group_names <- c("drug_costs", "test_costs", "other_costs")
  v_sampled_cost_groups <- intersect(v_cost_group_names, names(grouped))
  # "all_costs" only earns a row when it spans more than one cost group;
  # otherwise it would duplicate that group in every EVPPI table.
  if (length(v_sampled_cost_groups) > 1) {
    grouped <- append(
      grouped,
      list(all_costs = spec$parameter[spec$group %in% v_cost_group_names]),
      after = match(utils::tail(v_sampled_cost_groups, 1), names(grouped))
    )
  }
  grouped
}

#' Fill in PSA columns that are deterministic functions of sampled parameters
#'
#' Applied by generate_psa_samples() after the sampling loop, so every consumer
#' of the PSA parameter table sees the same derived values.
#'
#' @param samples Data frame of sampled PSA parameters.
#' @return The data frame with derived columns added.
apply_derived_psa_parameters <- function(samples) {
  v_joint_keys <- paste0("p_joint_", c("00", "01", "10"))
  if (all(v_joint_keys %in% names(samples))) {
    samples$p_joint_11 <- pmax(0, 1 - rowSums(samples[v_joint_keys]))
    samples$p_crp <- samples$p_joint_10 + samples$p_joint_11
    samples$p_tmb_braf <- samples$p_joint_01 + samples$p_joint_11
  }
  if (all(c("u_np", "u_decrement") %in% names(samples))) {
    # Non-negative decrement, floored at zero: u_p <= u_np in every draw.
    samples$u_p <- pmax(samples$u_np - samples$u_decrement, 0)
  }
  samples
}

configure_parameter_distributions <- function(params) {
  df_spec <- applicable_parameter_spec(parameter_distribution_spec(), params)
  list(
    distributions = create_parameter_distributions(params, df_spec),
    groups = create_parameter_groups(df_spec),
    spec = df_spec
  )
}

# ===============================================================================
# DETERMINISTIC SENSITIVITY-ANALYSIS RANGES AND STRUCTURAL SCENARIOS
# ===============================================================================
# Single source of truth for the one-way DSA bounds (09_DSA.R, figure4.qmd) and
# for the input-parameter table (table_1.qmd), so the published ranges always
# correspond to an analysis that is actually run (issue #153).

#' One-way DSA bounds: base value +/- mult, capped to the parameter's support
#'
#' The DSA set is spec$dsa, not the PSA parameter set: unit prices that are
#' fixed in the PSA are still varied deterministically (issue #154).
#'
#' Utilities keep their ordering (issue #190). Varying u_np holds the decrement
#' u_np - u_p fixed, so u_p moves with u_np (floored at 0); the linked value is
#' recorded in linked_par, linked_min and linked_max and applied by
#' dsa_endpoint_params(). Varying u_p is capped at the base u_np. Before #190
#' the lower u_np bound (0.584 with the CORRECT utilities) ran with u_p = 0.59,
#' i.e. with utility rising on progression.
#'
#' @param params Base-case parameter list (l_params_base).
#' @param spec Parameter specification from parameter_distribution_spec().
#' Parameters named in `params$dsa_bounds` (issue #209: the second-line
#' shares, cycles and contacts with a sourced range, and the progression
#' share) run at those bounds instead.
#'
#' @param mult Proportional variation (default DSA_mult, +/-20%).
#' @return Data frame with columns pars, min, max, linked_par, linked_min,
#'   linked_max (NA where varying the parameter moves no other parameter).
build_dsa_ranges <- function(params, spec = parameter_distribution_spec(),
                             mult = DSA_mult) {
  spec <- applicable_parameter_spec(spec, params)
  spec <- spec[spec$dsa, , drop = FALSE]
  v_pars <- spec$parameter
  v_missing_parameters <- setdiff(v_pars, names(params))
  if (length(v_missing_parameters) > 0) {
    stop("Missing base-case parameters for the DSA: ",
         paste(v_missing_parameters, collapse = ", "))
  }
  df_ranges <- data.frame(
    pars = v_pars,
    min = unlist(params[v_pars]) * (1 - mult),
    max = unlist(params[v_pars]) * (1 + mult),
    stringsAsFactors = FALSE
  )
  v_cost_groups <- c("drug_costs", "test_costs", "other_costs")
  v_bounded_zero <- spec$group %in% c(v_cost_groups, "prevalence", "later_sequences", "second_line")
  v_bounded_one  <- spec$group %in% c("utilities", "prevalence", "later_sequences") |
    (spec$group %in% "second_line" & grepl("^p_", spec$parameter))
  df_ranges$min[v_bounded_zero] <- pmax(df_ranges$min[v_bounded_zero], 0)
  df_ranges$max[v_bounded_one]  <- pmin(df_ranges$max[v_bounded_one], 1)
  # Sourced bounds (issue #209): a parameter listed in params$dsa_bounds (a
  # data frame with pars, min, max; script 05 fills it from the second-line
  # pathway table and the progression-share beta) runs at those bounds
  # instead of base +/- mult.
  if (!is.null(params$dsa_bounds)) {
    df_bounds <- params$dsa_bounds
    if (!all(c("pars", "min", "max") %in% names(df_bounds)) ||
        any(!is.finite(df_bounds$min) | !is.finite(df_bounds$max))) {
      stop("params$dsa_bounds must be a data frame with finite pars, min and max.")
    }
    v_rows <- match(df_bounds$pars, df_ranges$pars)
    if (anyNA(v_rows)) {
      stop("params$dsa_bounds names parameters outside the DSA: ",
           paste(df_bounds$pars[is.na(v_rows)], collapse = ", "))
    }
    df_ranges$min[v_rows] <- df_bounds$min
    df_ranges$max[v_rows] <- df_bounds$max
  }
  df_ranges$linked_par <- NA_character_
  df_ranges$linked_min <- NA_real_
  df_ranges$linked_max <- NA_real_
  if (all(c("u_np", "u_p") %in% df_ranges$pars)) {
    u_decrement <- params$u_np - params$u_p
    if (u_decrement < 0) stop("Base-case u_p exceeds u_np")
    i_np <- which(df_ranges$pars == "u_np")
    i_p  <- which(df_ranges$pars == "u_p")
    df_ranges$linked_par[i_np] <- "u_p"
    df_ranges$linked_min[i_np] <- max(df_ranges$min[i_np] - u_decrement, 0)
    df_ranges$linked_max[i_np] <- max(df_ranges$max[i_np] - u_decrement, 0)
    df_ranges$max[i_p] <- min(df_ranges$max[i_p], params$u_np)
  }
  if (any(df_ranges$min > df_ranges$max)) {
    stop("Invalid DSA ranges (min > max) for: ",
         paste(df_ranges$pars[df_ranges$min > df_ranges$max], collapse = ", "))
  }
  rownames(df_ranges) <- NULL
  df_ranges
}

#' Parameter list for one endpoint of a one-way DSA run
#'
#' Sets the varied parameter to its bound and any linked parameter to its
#' linked value (issue #190: u_p moves with u_np). u_decrement, which only the
#' PSA reads, is kept equal to u_np - u_p.
#'
#' @param params Base-case parameter list (l_params_base).
#' @param ranges Data frame from build_dsa_ranges().
#' @param i Row of ranges.
#' @param side "min" or "max".
#' @return params with the endpoint applied.
dsa_endpoint_params <- function(params, ranges, i, side = c("min", "max")) {
  side <- match.arg(side)
  params[[ranges$pars[i]]] <- ranges[[side]][i]
  linked_par <- ranges$linked_par[i]
  if (!is.null(linked_par) && !is.na(linked_par)) {
    params[[linked_par]] <- ranges[[paste0("linked_", side)]][i]
  }
  if (all(c("u_np", "u_p", "u_decrement") %in% names(params))) {
    params$u_decrement <- params$u_np - params$u_p
  }
  params
}

#' Deterministic structural scenarios run alongside the one-way DSA
#'
#' Discount rate, model time horizon, post-progression treatment cost and the
#' second treatment sequence. The values here are what 09_DSA.R runs and what
#' table_1.qmd prints as ranges. `sides` names the endpoints that are actually
#' run; a one-sided scenario is one whose other endpoint is the base case.
#'
#' @param base_dr Base-case annual discount rate.
#' @param base_horizon_years Base-case time horizon in years.
#' @param base_pp_cost Base-case post-progression cost per quarter.
#' @return Named list; each element has label, unit, base, min, max, sides.
dsa_structural_scenarios <- function(base_dr, base_horizon_years,
                                     base_pp_cost = 0) {
  list(
    Discount_rate = list(
      label = "Discount rate (costs and QALYs)", unit = "rate",
      base = base_dr, min = 0, max = 0.08, sides = c("min", "max")
    ),
    Time_horizon = list(
      label = "Time horizon", unit = "years",
      base = base_horizon_years, min = 5, max = 20, sides = c("min", "max")
    ),
    # Since issue #209 the base case costs the Norwegian second-line pathway
    # once at incident progression. The former flat EUR 5,000 per quarter of
    # progressed occupancy (issue #154, an illustrative value, not a costed
    # pathway) is kept as a scenario that REPLACES the pathway: adding it would
    # charge second-line therapy twice.
    Post_progression_cost = list(
      label = "Flat post-progression cost per quarter instead of the second-line pathway",
      unit = "EUR", base = base_pp_cost, min = base_pp_cost, max = 5000,
      sides = "max"
    ),
    # No second-line therapy cost (the base case before issue #209): the
    # pathway is switched off and nothing replaces it.
    Second_line_therapy = list(
      label = "Second-line therapy at progression (1 = costed pathway, 0 = not costed)",
      unit = "indicator", base = 1, min = 0, max = 1, sides = "min"
    ),
    # The base case gives every progression-free patient a second eight-cycle
    # sequence at weeks 24-38. The scenario removes it, holding survival at the
    # trial estimate, so it bounds the cost side of that assumption only.
    # Without a second sequence there is no third or fourth (issue #208), so
    # the scenario also removes the later sequences.
    Second_sequence = list(
      label = "Second treatment sequence at weeks 24-38 and the later sequences (1 = given, 0 = omitted)",
      unit = "indicator", base = 1, min = 0, max = 1, sides = "min"
    )
  )
}

#' Endpoints of a structural scenario that 09_DSA.R runs
#'
#' @param scenario One element of dsa_structural_scenarios().
#' @return Character vector of endpoint names ("min", "max").
structural_scenario_sides <- function(scenario) {
  if (is.null(scenario$sides)) c("min", "max") else scenario$sides
}

# ---------------------------------------------------------------------------
# Treatment and monitoring schedules (issues #49, #172)
# ---------------------------------------------------------------------------
# Schedules are 0/1 vectors on the weekly model grid; position i is modeled
# week i - 1. Script 05 and build_horizon_params() both build them here, so a
# horizon shorter than the treatment sequence truncates the schedule instead of
# lengthening the vector, and the recurring monitoring rules work for any
# horizon (seq(13, n, by = 12) fails for n < 13).

# Drug administration positions (modeled weeks 0-38, both sequences).
TREATMENT_SCHEDULE_POSITIONS <- list(
  nivo = c(5, 7, 13, 15, 29, 31, 37, 39),
  FLOX_exp = c(1, 3, 9, 11, 25, 27, 33, 35),
  FLOX_control = c(1, 3, 5, 7, 9, 11, 13, 15, 25, 27, 29, 31, 33, 35, 37, 39)
)

# Monitoring interval of progression-free patients (issue #207; trial team's
# answer MON2 of 15 October 2026, docs/open_questions_2026-10-15-answers.csv:
# "every 8 weeks during treatment and break"): a CT scan at baseline and every
# 8 weeks up to the horizon, including after the last sequence (was every 12
# weeks). Blood tests follow MON1 (blood_schedule()): one per administration
# week on treatment and, off treatment, one with every CT.
PF_CT_INTERVAL_WEEKS <- 8L

#' Positions first, first + by, ... up to n_points (empty when n_points < first)
#'
#' @param first First position.
#' @param by Step between positions.
#' @param n_points Number of grid points.
#' @return Integer vector of positions (empty when n_points < first).
recurring_positions <- function(first, by, n_points) {
  if (n_points < first) integer(0) else seq(first, n_points, by = by)
}

#' Build the six treatment and monitoring schedules on an n-point grid
#'
#' @param n_points Number of model time points (time_horizon + 1).
#' @param nivo,FLOX_exp,FLOX_control Administration positions; positions beyond
#'   the grid are dropped.
#' @return Named list of 0/1 vectors of length n_points: l_nivo, l_FLOX_exp,
#'   l_FLOX_control, l_CT (baseline, then every PF_CT_INTERVAL_WEEKS weeks),
#'   l_blood (baseline, every administration week and every CT;
#'   blood_schedule()) and l_visit (baseline, every administration and every
#'   CT scan; issue #164). The second outpatient contact of a two-day FLOX
#'   administration (issue #207) is not a schedule: calculate_outcomes()
#'   derives it from the FLOX schedule of the regimen.
build_treatment_schedules <- function(n_points,
                                      nivo = TREATMENT_SCHEDULE_POSITIONS$nivo,
                                      FLOX_exp = TREATMENT_SCHEDULE_POSITIONS$FLOX_exp,
                                      FLOX_control = TREATMENT_SCHEDULE_POSITIONS$FLOX_control) {
  if (!is.numeric(n_points) || length(n_points) != 1L || !is.finite(n_points) ||
      n_points < 1 || n_points != round(n_points)) {
    stop("n_points must be a positive whole number. Got: ", paste(n_points, collapse = ", "))
  }
  indicator <- function(positions) {
    v <- rep(0, n_points)
    v[positions[positions <= n_points]] <- 1
    v
  }
  l_nivo <- indicator(nivo)
  l_FLOX_exp <- indicator(FLOX_exp)
  l_FLOX_control <- indicator(FLOX_control)
  l_CT <- indicator(recurring_positions(1, PF_CT_INTERVAL_WEEKS, n_points))
  list(
    l_nivo = l_nivo,
    l_FLOX_exp = l_FLOX_exp,
    l_FLOX_control = l_FLOX_control,
    l_CT = l_CT,
    l_blood = blood_schedule(l_nivo, l_FLOX_exp, l_FLOX_control, l_CT),
    l_visit = visit_schedule(l_nivo, l_FLOX_exp, l_FLOX_control, l_CT)
  )
}

#' Outpatient visit schedule (issue #164)
#'
#' A visit is charged at baseline, at every drug administration and at every
#' CT scan. With the 8-weekly CT (issue #207) the scans at modeled weeks 8, 24
#' and 32 coincide with an administration; the CT rule adds a visit with the
#' week-16 scan of the treatment break and a surveillance visit with each scan
#' after the last administration (week 40 onwards) up to the horizon. The
#' second contact of a two-day FLOX administration is added in
#' calculate_outcomes() (FLOX_CONTACTS_PER_ADMINISTRATION), not here.
#'
#' @param l_nivo,l_FLOX_exp,l_FLOX_control,l_CT 0/1 schedules of equal length.
#' @return 0/1 vector: 1 at baseline, every administration and every CT.
visit_schedule <- function(l_nivo, l_FLOX_exp, l_FLOX_control, l_CT) {
  l_visit <- as.numeric(l_nivo == 1 | l_FLOX_exp == 1 | l_FLOX_control == 1 |
                          l_CT == 1)
  l_visit[1] <- 1
  l_visit
}

#' Blood-test schedule of progression-free patients (issue #207)
#'
#' Trial team's answer MON1 of 15 October 2026 ("blood test every two weeks";
#' the protocol requires local labs within 3 days before each dose): one test
#' per administration week of a sequence, which in both arms falls every two
#' weeks (FLOX or nivolumab). Progression-free time off treatment (the break
#' at weeks 16-22 and after the last sequence) follows the protocol schedule,
#' a blood test every 8 weeks together with the CT; this reading is an
#' interpretation recorded in docs/decisions.md. Both arms have the same
#' administration weeks, so one vector serves both; by construction it
#' coincides with visit_schedule().
#'
#' @param l_nivo,l_FLOX_exp,l_FLOX_control,l_CT 0/1 schedules of equal length.
#' @return 0/1 vector: 1 at baseline, every administration week and every CT.
blood_schedule <- function(l_nivo, l_FLOX_exp, l_FLOX_control, l_CT) {
  l_blood <- as.numeric(l_nivo == 1 | l_FLOX_exp == 1 | l_FLOX_control == 1 |
                          l_CT == 1)
  l_blood[1] <- 1
  l_blood
}

# ---------------------------------------------------------------------------
# Applying a structural scenario to a base-case parameter list (issue #157)
# ---------------------------------------------------------------------------
# Moved here from 09_DSA.R so that the one-way DSA (09_DSA.R, OWSA.qmd,
# Table 1) and the deterministic scenario table (table_s8.qmd) run exactly the
# same scenario code. The discount rate is applied to costs and QALYs together.
# A different time horizon re-predicts the base-case survival curves on the
# new weekly grid from the selected joint models and rebuilds the schedule
# vectors with the same rules as script 05 (drug positions are fixed calendar
# weeks; the 8-weekly CT, with its blood test and visit, recurs to the end of
# the horizon, issue #207). The post-progression cost and second-sequence scenarios test
# the two scope assumptions raised in issue #154.

#' Rebuild the base-case parameter list for a different model horizon
#'
#' @param base_params Base-case parameter list (l_params_base).
#' @param horizon_weeks New horizon in weeks.
#' @param models list(os = <flexsurvreg>, pfs = <flexsurvreg>): the selected
#'   ordering-constrained best-fit models (models$best_fit$os / $pfs).
#' @param strategies_df Strategy table with `id` and `prevalence`.
#' @param data_complete Complete-case data used for population averaging.
#' @return Parameter list with curves and schedules on the new grid.
build_horizon_params <- function(base_params, horizon_weeks, models,
                                 strategies_df, data_complete) {
  if (is.null(models) || is.null(models$os) || is.null(models$pfs)) {
    stop("Fitted best-fit survival models are required for a time-horizon scenario.")
  }
  v_tp <- seq(0, horizon_weeks)
  n <- length(v_tp)
  preds <- generate_population_averaged_predictions(
    models = list(os = models$os, pfs = models$pfs),
    strategies_df = strategies_df,
    data_complete = data_complete,
    time_points = v_tp,
    prevalences = setNames(strategies_df$prevalence, strategies_df$id)
  )

  p <- set_population_predictions(base_params, preds)
  ctrl <- get_control_strategy()
  p$p_os  <- setNames(list(preds[[ctrl]]$os),  paste0(ctrl, "_OS"))
  p$p_pfs <- setNames(list(preds[[ctrl]]$pfs), paste0(ctrl, "_PFS"))
  for (bm in get_biomarkers()) {
    p$p_os[[paste0(bm, "_pos_OS")]]       <- preds[[bm]]$biomarker_positive$os
    p$p_os[[paste0(bm, "_neg_OS")]]       <- preds[[bm]]$biomarker_negative$os
    p$p_os[[paste0(bm, "_weighted_OS")]]  <- preds[[bm]]$os
    p$p_pfs[[paste0(bm, "_pos_PFS")]]     <- preds[[bm]]$biomarker_positive$pfs
    p$p_pfs[[paste0(bm, "_neg_PFS")]]     <- preds[[bm]]$biomarker_negative$pfs
    p$p_pfs[[paste0(bm, "_weighted_PFS")]] <- preds[[bm]]$pfs
  }

  # Keep the base list's administration positions (truncated to the horizon);
  # monitoring recurs to the end of the new horizon.
  schedules <- build_treatment_schedules(
    n,
    nivo = which(base_params$l_nivo == 1),
    FLOX_exp = which(base_params$l_FLOX_exp == 1),
    FLOX_control = which(base_params$l_FLOX_control == 1)
  )
  p[names(schedules)] <- schedules
  p$time_horizon <- horizon_weeks
  p
}

#' Remove the second treatment sequence
#'
#' Zero every administration from the given schedule position onwards and
#' rebuild the visit and blood-test schedules from the remaining
#' administrations and the unchanged CT scans (visit_schedule(),
#' blood_schedule()). The dropped weeks lose their administration visits, the
#' second contacts of their FLOX administrations (derived from the FLOX
#' schedules in calculate_outcomes()) and their on-treatment blood tests; the
#' CT scans with their visit and their blood test (the off-treatment rule,
#' issue #207) are kept. Survival is unchanged, so the scenario bounds the cost
#' of the assumption only. A third or fourth sequence presupposes a second, so
#' the later sequences (issue #208) are removed too: the shares starting a
#' third sequence are set to zero in both arms, which zeroes the fourth.
#'
#' @param p Parameter list.
#' @param first_position First schedule position of the second sequence
#'   (position 25 = modeled week 24).
#' @return p with the second and later sequences removed.
drop_second_sequence <- function(p, first_position = 25L) {
  zap <- function(v) {
    if (length(v) >= first_position) v[seq.int(first_position, length(v))] <- 0
    v
  }
  p$l_nivo         <- zap(p$l_nivo)
  p$l_FLOX_exp     <- zap(p$l_FLOX_exp)
  p$l_FLOX_control <- zap(p$l_FLOX_control)
  p$l_visit <- visit_schedule(p$l_nivo, p$l_FLOX_exp, p$l_FLOX_control, p$l_CT)
  p$l_blood <- blood_schedule(p$l_nivo, p$l_FLOX_exp, p$l_FLOX_control, p$l_CT)
  for (nm in intersect(c("p_seq3_control", "p_seq3_experimental"), names(p))) p[[nm]] <- 0
  p
}

#' Switch off the second-line therapy cost (issue #209)
#'
#' Sets the share of progressing patients treated in second line to 0, which
#' makes second_line_cost_per_progression() zero exactly; a list without
#' second-line fields is returned unchanged.
#'
#' @param p Parameter list.
#' @return p with `p_sl_treated = 0` where the field exists.
drop_second_line_therapy <- function(p) {
  if (!is.null(p$p_sl_treated)) p$p_sl_treated <- 0
  p
}

#' Apply one endpoint of a structural scenario to a parameter list
#'
#' @param scenario_name Name of an entry of dsa_structural_scenarios().
#' @param value The endpoint value (scenario$min or scenario$max).
#' @param base_params Base-case parameter list.
#' @param base_horizon Base-case horizon in weeks (passed back unchanged for
#'   every scenario except Time_horizon).
#' @param models,strategies_df,data_complete Needed for Time_horizon only; see
#'   build_horizon_params().
#' @return list(params, time_horizon) ready for
#'   model_fun(params, time_horizon = time_horizon).
apply_structural_scenario <- function(scenario_name, value, base_params,
                                      base_horizon, models = NULL,
                                      strategies_df = NULL,
                                      data_complete = NULL) {
  p <- base_params
  horizon <- base_horizon
  if (scenario_name == "Discount_rate") {
    p$dr_costs <- value
    p$dr_effects <- value
  } else if (scenario_name == "Time_horizon") {
    horizon <- as.integer(round(value * 52))
    p <- build_horizon_params(base_params, horizon, models, strategies_df,
                              data_complete)
  } else if (scenario_name == "Post_progression_cost") {
    # The flat rate replaces the second-line pathway (issue #209).
    p$c_other_pp <- value
    if (value > 0) p <- drop_second_line_therapy(p)
  } else if (scenario_name == "Second_line_therapy") {
    if (value == 0) p <- drop_second_line_therapy(p)
  } else if (scenario_name == "Second_sequence") {
    if (value == 0) p <- drop_second_sequence(p)
  } else {
    stop("Unknown structural scenario: ", scenario_name)
  }
  list(params = p, time_horizon = horizon)
}

#' Human-readable label for one structural scenario endpoint
#'
#' @param scenario_name Name of an entry of dsa_structural_scenarios().
#' @param value The endpoint value.
#' @return Character label used in Table S8 and the OWSA report.
structural_scenario_label <- function(scenario_name, value) {
  switch(scenario_name,
    Discount_rate = sprintf("Discount rate %g%% (costs and QALYs)", 100 * value),
    Time_horizon = sprintf("Time horizon %g years", value),
    Post_progression_cost = sprintf("Post-progression cost EUR %s per quarter (no second-line pathway)",
                                    format(value, big.mark = ",", scientific = FALSE)),
    Second_line_therapy = if (value == 0) "Second-line therapy not costed" else
      "Second-line therapy costed",
    Second_sequence = if (value == 0) "Second and later treatment sequences omitted" else
      "Second and later treatment sequences given",
    paste(scenario_name, value)
  )
}
