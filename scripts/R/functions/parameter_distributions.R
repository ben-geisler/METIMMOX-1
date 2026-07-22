# ===============================================================================
# NON-SURVIVAL UNCERTAINTY PARAMETERS
# ===============================================================================
# The table below is the single source of truth for PSA, DSA, and EVPPI
# parameter membership. Add a parameter once here; downstream analyses derive
# their parameter vectors from names(param_distributions).

parameter_distribution_spec <- function() {
  biomarker_cost_parameters <- unique(unname(biomarker_cost_key()))
  prevalence_parameters <- unname(biomarker_prevalence_key())
  static_costs <- data.frame(
    parameter = c(
      "c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
      "c_other_visit", "c_other_baseline", "c_other_follow", "c_other_last"
    ),
    distribution = "gamma",
    cv = 0.20,
    group = c(
      "drug_costs", "drug_costs", "test_costs", "test_costs",
      rep("other_costs", 4)
    ),
    order = c(10:11, 20:21, 30:33),
    stringsAsFactors = FALSE
  )
  biomarker_costs <- data.frame(
    parameter = biomarker_cost_parameters,
    distribution = "gamma",
    cv = 0.20,
    group = "test_costs",
    order = 22 + seq_along(biomarker_cost_parameters),
    stringsAsFactors = FALSE
  )
  utilities <- data.frame(
    parameter = c("u_np", "u_p"),
    distribution = "beta",
    cv = 0.15,
    group = "utilities",
    order = 40:41,
    stringsAsFactors = FALSE
  )
  prevalence <- data.frame(
    parameter = prevalence_parameters,
    distribution = "beta",
    cv = 0.15,
    group = "prevalence",
    order = 50 + seq_along(prevalence_parameters),
    stringsAsFactors = FALSE
  )

  spec <- rbind(static_costs, biomarker_costs, utilities, prevalence)
  spec <- spec[order(spec$order), , drop = FALSE]
  spec$order <- NULL
  if (anyDuplicated(spec$parameter)) {
    stop("Parameter distribution specification contains duplicate names.")
  }
  spec
}

distribution_parameters <- function(mean_value, distribution, cv) {
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

create_parameter_distributions <- function(params,
                                           spec = parameter_distribution_spec()) {
  missing_parameters <- setdiff(spec$parameter, names(params))
  if (length(missing_parameters) > 0) {
    stop("Missing base-case parameters: ",
         paste(missing_parameters, collapse = ", "))
  }

  distributions <- lapply(seq_len(nrow(spec)), function(i) {
    distribution_parameters(
      params[[spec$parameter[i]]], spec$distribution[i], spec$cv[i]
    )
  })
  stats::setNames(distributions, spec$parameter)
}

create_parameter_groups <- function(spec = parameter_distribution_spec()) {
  grouped <- split(spec$parameter, factor(
    spec$group,
    levels = c("drug_costs", "test_costs", "other_costs", "utilities", "prevalence")
  ))
  grouped <- grouped[lengths(grouped) > 0]
  append(
    grouped,
    list(all_costs = spec$parameter[spec$group %in% c(
      "drug_costs", "test_costs", "other_costs"
    )]),
    after = match("other_costs", names(grouped))
  )
}

configure_parameter_distributions <- function(params) {
  spec <- parameter_distribution_spec()
  list(
    distributions = create_parameter_distributions(params, spec),
    groups = create_parameter_groups(spec),
    spec = spec
  )
}
