# ===============================================================================
# PREDICTION FUNCTIONS FOR PARAMETRIC SURVIVAL MODELS
# ===============================================================================
# This file contains functions for generating individual and population-level
# survival predictions from fitted parametric models
# ===============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(survival, flexsurv, dplyr, tidyr)
source(here::here("R/prediction_population.R"))

# ===============================================================================
# POPULATION AVERAGING FUNCTIONS
# ===============================================================================
# These functions implement methodologically rigorous population-averaged
# predictions that account for the full distribution of prognostic factors
# (age, sex, biomarkers) in the patient population.
# ===============================================================================

# Helper function to extract survival probabilities from flexsurv predict objects
# Handles the nested structure: pred$.pred[[i]]$.pred_survival
extract_all_survival_probabilities <- function(pred_object) {
  # pred_object is the result of predict(flexsurvreg_model, newdata = data, times = time_points)
  # Structure: tibble with $.pred column containing list of predictions
  # Each element in $.pred has $.pred_survival (vector of survival probs at each time point)

  # Extract all survival probability vectors (one per patient)
  all_surv_probs <- lapply(pred_object$.pred, function(x) x$.pred_survival)

  # Convert to matrix: rows = time points, columns = patients
  surv_matrix <- do.call(cbind, all_surv_probs)

  return(surv_matrix)
}

# Generate a population-average survival curve from a fitted flexsurv model.
# Supports both current tidypredict-style output and the legacy list-of-frames
# structure used by older flexsurv versions.
predict_pop_avg <- function(model, newdata, time_pts) {
  tryCatch({
    pred <- predict(
      model, newdata = newdata, type = "survival", times = time_pts
    )
    if (".pred" %in% names(pred)) {
      surv_probs <- lapply(pred$.pred, function(x) x$.pred_survival)
      surv_matrix <- do.call(cbind, surv_probs)
      return(rowMeans(surv_matrix, na.rm = TRUE))
    }
    if (is.list(pred) && length(pred) > 0 && is.data.frame(pred[[1]])) {
      surv_probs <- sapply(pred, function(x) x$est)
      if (is.matrix(surv_probs)) {
        return(rowMeans(surv_probs, na.rm = TRUE))
      }
      return(surv_probs)
    }
    NULL
  }, error = function(e) NULL)
}

# Convert a survfit object to the minimal data frame used by report plots.
extract_km_data <- function(km_fit) {
  data.frame(
    time = km_fit$time,
    surv = km_fit$surv,
    stringsAsFactors = FALSE
  )
}

# Resolve the joint-model sampling component from the sampling cache
#
# Since issue #156 the cache stores the joint model once under `joint` and lists
# the biomarker strategies that share it under `biomarkers`. Legacy caches held
# one identical component per biomarker key; the first present key is used so
# older objects (and test fixtures) still resolve.
get_joint_sampling_models <- function(sampling_models) {
  if (!is.null(sampling_models$joint)) return(sampling_models$joint)
  biomarkers <- get_biomarkers()
  present <- intersect(biomarkers, names(sampling_models))
  if (length(present) == 0L) {
    stop("sampling_models has no joint-model component; expected `joint` or one of: ",
         paste(biomarkers, collapse = ", "))
  }
  sampling_models[[present[1]]]
}

# Build a flexsurvreg object whose estimates are one sampled coefficient vector
#
# `draw` is on the scale flexsurvreg optimises (log for positive baseline
# parameters, identity for covariate effects), i.e. the scale of `model$opt$par`
# and `model$res.t[, "est"]`. Prediction reads the baseline parameters from
# `res.t` and the covariate effects from `res`, so both tables, the
# coefficient vector and the optimiser output are replaced consistently.
build_sampled_flexsurv_model <- function(model, draw) {
  par_names <- rownames(model$res)
  if (is.null(names(draw))) names(draw) <- par_names
  if (!identical(names(draw), par_names)) {
    draw <- draw[par_names]
  }
  if (anyNA(draw)) stop("Sampled coefficient vector does not cover every model parameter")
  sampled <- model
  sampled$opt$par <- draw[model$optpars]
  sampled$res.t[, "est"] <- draw
  sampled$res[, "est"] <- draw
  for (j in seq_along(model$dlist$pars)) {
    sampled$res[j, "est"] <- model$dlist$inv.transforms[[j]](draw[j])
  }
  sampled$coefficients <- draw
  # The estimate-specific fields below describe the original fit only.
  sampled$res.t[, c("L95%", "U95%", "se")] <- NA_real_
  sampled$res[, c("L95%", "U95%", "se")] <- NA_real_
  sampled$sampled_from_original <- TRUE
  sampled
}

# Return the sampled OS and PFS models for one draw in the legacy sample shape
#
# The result is `list(os = list(model, coefficients, dist), pfs = ..., failed)`,
# which is what model_fun(), the PSA prediction helper and the interaction
# coefficient extraction consumed from the bootstrap cache. A multivariate-normal
# cache (issue #156) stores only the coefficient draws and materialises the
# model object on request; a legacy bootstrap cache returns its stored sample.
sampled_survival_models <- function(component, idx) {
  if (!is.null(component$samples)) {
    return(component$samples[[idx]])
  }
  if (is.null(component$draws)) {
    stop("Sampling component has neither `draws` (issue #156) nor legacy `samples`")
  }
  n <- nrow(component$draws$os)
  if (!is.numeric(idx) || length(idx) != 1L || idx < 1 || idx > n || idx != floor(idx)) {
    stop("idx must be a single integer in 1..", n, "; got ", idx)
  }
  originals <- list(os = component$original_os, pfs = component$original_pfs)
  dists <- list(os = component$dist_os, pfs = component$dist_pfs)
  out <- lapply(c(os = "os", pfs = "pfs"), function(outcome) {
    draw <- component$draws[[outcome]][idx, ]
    list(
      model = build_sampled_flexsurv_model(originals[[outcome]], draw),
      coefficients = draw,
      dist = dists[[outcome]]
    )
  })
  out$failed <- FALSE
  out
}


# Main function for population-averaged predictions across all strategies
# Uses actual patient-level covariate distributions instead of reference patient
generate_population_averaged_predictions <- function(models, strategies_df,
                                                     data_complete, time_points,
                                                     prevalences = NULL,
                                                     quiet = FALSE,
                                                     weights = NULL) {
  # The explicit population is shared by every counterfactual treatment path.
  population <- data_complete
  if (is.null(weights)) weights <- target_population_weights(population, prevalences)
  if (length(weights) != nrow(population)) stop("Population weights have the wrong length.")
  control <- get_control_strategy()
  levels_rx <- levels(population$Rx)
  predict_part <- function(rows, rx) {
    if (!length(rows)) stop("Both biomarker subgroups are required for prediction.")
    newdata <- population[rows, , drop = FALSE]
    newdata$Rx <- factor(rx, levels = levels_rx)
    list(rows = rows,
      os = extract_all_survival_probabilities(predict(models$os, newdata = newdata,
        type = "survival", times = time_points)),
      pfs = extract_all_survival_probabilities(predict(models$pfs, newdata = newdata,
        type = "survival", times = time_points)))
  }
  curves <- setNames(list(predict_part(seq_len(nrow(population)), levels_rx[1])), control)
  for (biomarker in get_biomarkers()) {
    status <- as.numeric(as.character(population[[biomarker]]))
    curves[[biomarker]] <- list(positive = predict_part(which(status == 1), levels_rx[2]),
                                negative = predict_part(which(status == 0), levels_rx[1]))
  }
  if (!quiet) cat("Predicted every strategy in one target population (n =", nrow(population), ").\n")
  standardize_population_curves(curves, population, weights)
}


# Generate the five endpoint-specific curves needed during joint distribution
# selection. Keeping OS and PFS prediction separate lets the selection loop
# compute each candidate distribution once instead of once per pair.
generate_population_averaged_endpoint_curves <- function(model, data_complete,
                                                          time_points) {
  exp_rx <- levels(data_complete$Rx)[2]
  ctrl_rx <- levels(data_complete$Rx)[1]
  biomarkers <- get_biomarkers()

  average_prediction <- function(newdata, rx_level) {
    newdata$Rx <- factor(rx_level, levels = levels(data_complete$Rx))
    prediction <- predict(
      model, newdata = newdata, type = "survival", times = time_points
    )
    rowMeans(extract_all_survival_probabilities(prediction), na.rm = TRUE)
  }

  curves <- setNames(
    list(average_prediction(data_complete, ctrl_rx)),
    get_control_strategy()
  )
  for (biomarker in biomarkers) {
    status <- as.numeric(as.character(data_complete[[biomarker]]))
    positive_data <- data_complete[status == 1, , drop = FALSE]
    negative_data <- data_complete[status == 0, , drop = FALSE]

    if (nrow(positive_data) == 0 || nrow(negative_data) == 0) {
      stop("Both biomarker subgroups are required for ordering selection: ",
           biomarker)
    }

    curves[[paste0(biomarker, "_positive")]] <-
      average_prediction(positive_data, exp_rx)
    curves[[paste0(biomarker, "_negative")]] <-
      average_prediction(negative_data, ctrl_rx)
  }
  curves
}

check_endpoint_curve_ordering <- function(os_curves, pfs_curves,
                                          tolerance = 0) {
  curve_names <- union(names(os_curves), names(pfs_curves))
  details <- lapply(curve_names, function(curve_name) {
    os <- os_curves[[curve_name]]
    pfs <- pfs_curves[[curve_name]]
    if (is.null(os) || is.null(pfs) || length(os) != length(pfs) ||
        any(!is.finite(os)) || any(!is.finite(pfs))) {
      return(data.frame(
        curve = curve_name, n_violations = NA_integer_,
        max_pfs_minus_os = Inf, ordered = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    gaps <- pfs - os
    data.frame(
      curve = curve_name,
      n_violations = sum(gaps > tolerance),
      max_pfs_minus_os = max(gaps),
      ordered = all(gaps <= tolerance),
      stringsAsFactors = FALSE
    )
  })

  details <- do.call(rbind, details)
  list(
    ordered = all(details$ordered),
    n_violations = if (anyNA(details$n_violations)) NA_integer_ else
      sum(details$n_violations),
    max_pfs_minus_os = max(details$max_pfs_minus_os),
    details = details
  )
}

# Check the defining partitioned-survival ordering invariant on population curves.
# The subgroup checks are stricter than checking only prevalence-weighted strategy
# curves: if every subgroup is ordered, each weighted strategy is ordered too.
check_population_survival_ordering <- function(predictions, tolerance = 0) {
  control_strategy <- get_control_strategy()
  curve_sets <- setNames(list(predictions[[control_strategy]]), control_strategy)

  strategy_names <- setdiff(names(predictions), control_strategy)
  for (strategy_name in strategy_names) {
    strategy <- predictions[[strategy_name]]
    curve_sets[[paste0(strategy_name, "_positive")]] <- strategy$biomarker_positive
    curve_sets[[paste0(strategy_name, "_negative")]] <- strategy$biomarker_negative
  }

  check_endpoint_curve_ordering(
    os_curves = lapply(curve_sets, `[[`, "os"),
    pfs_curves = lapply(curve_sets, `[[`, "pfs"),
    tolerance = tolerance
  )
}
