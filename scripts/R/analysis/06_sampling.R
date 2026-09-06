#libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, survival, flexsurv, dplyr, tidyr, mvtnorm)

# Ensure time_points is the same as used in section 3
time_points <- seq(0, time_horizon, by = 1)
time_points_length <- length(time_points)

# ===============================================================================
# SAMPLING CACHE CONFIGURATION
# ===============================================================================

# Define cache directory and file paths
sampling_cache_dir <- cache_dir()
if (!dir.exists(sampling_cache_dir)) {
  dir.create(sampling_cache_dir, recursive = TRUE)
  cat("Created sampling cache directory:", sampling_cache_dir, "\n")
}

# Create cache file path based on n_samples
cache_file <- sampling_cache_path(n_samples)
sampling_seed <- analysis_seed

# ===============================================================================
# MODEL FORMULA SELECTION
# ===============================================================================
# Model formulas are defined in scripts/R/functions/model_configs.R
# This ensures PSA uses the same single joint model as base case.
#
# Economic biomarker strategies are CRP and TMB/BRAF. Both share one joint
# formula, so one bootstrap of the joint model serves every strategy. The
# control arm is NOT a separate age/sex model: it is the same joint fit
# predicted with Rx forced to the control level, exactly as in the base case
# (issue #151). The cache therefore holds one component per biomarker, all
# pointing at the same joint bootstrap; a legacy "control" component, if
# present in an older cache, is ignored.
# ===============================================================================

# Source model configurations if not already loaded
if (!exists("get_model_configs")) {
  source(here::here("scripts/R/functions/model_configs.R"))
}
if (!exists("extract_all_survival_probabilities", mode = "function")) {
  source(here::here("scripts/R/functions/prediction_functions.R"))
}

# Get model configuration
model_config <- get_current_model_config()
cat("Resampling: Using", model_config$label, "\n")
cat("Description:", model_config$description, "\n")

# Create biomarker formula function using central config
create_biomarker_formula <- function(outcome, biomarker) {
  return(get_strategy_formula(biomarker, outcome))
}

# ===============================================================================
# CORRELATED SAMPLING FUNCTION
# ===============================================================================
# MULTIVARIATE-NORMAL COEFFICIENT SAMPLING (issue #156)
# ===============================================================================
# Parameter uncertainty in the survival models is propagated by drawing
# coefficient vectors from the multivariate-normal approximation to the sampling
# distribution of each fitted joint model (mean = maximum-likelihood estimate,
# covariance = inverse observed information, both on flexsurvreg's optimisation
# scale). This replaced the unstratified nonparametric bootstrap, whose
# resamples could contain 0-2 of the 6 control-arm CRP-positive patients and
# produced interaction coefficients between -3.7 and +4.0 and singular fits.
#
# OS and PFS coefficients are drawn independently: the two models are fitted
# separately and no joint covariance between them is estimated. The bootstrap
# preserved the OS-PFS correlation by refitting both to the same resample; that
# property is given up here in exchange for draws that cannot degenerate
# (user decision, issue #156). Within a draw, the same coefficient vector serves
# the control arm (Rx = control) and every biomarker subgroup, so the
# control/biomarker correlation of issue #151 is unchanged.
#
# The cache stores the draws as n_samples x n_parameters matrices and the two
# original fits; sampled_survival_models() (prediction_functions.R) builds a
# flexsurvreg object with the drawn estimates on request.
# ===============================================================================

SAMPLING_METHOD <- "mvn_v1"

sample_survival_coefficients <- function(formula_os, formula_pfs, data,
                                         dist_os = "weibull", dist_pfs = "weibull",
                                         n_samples = n_samples, seed = 123L) {

  n_patients <- nrow(data)
  cat("Generating", n_samples, "multivariate-normal coefficient draws for the",
      "joint OS and PFS models fitted on", n_patients, "patients...\n")
  cat("Using OS distribution:", dist_os, "| PFS distribution:", dist_pfs, "\n")

  original_os <- flexsurvreg(formula_os, data = data, dist = dist_os)
  original_pfs <- flexsurvreg(formula_pfs, data = data, dist = dist_pfs)

  draw_coefficients <- function(model, label) {
    if (is.null(model$cov) || anyNA(model$cov)) {
      stop("Covariance matrix unavailable for the ", label,
           " model (non-converged fit); cannot draw coefficients")
    }
    par_names <- rownames(model$res)
    draws <- matrix(NA_real_, nrow = n_samples, ncol = length(par_names),
                    dimnames = list(NULL, par_names))
    draws[, model$optpars] <- mvtnorm::rmvnorm(n_samples, model$opt$par, model$cov)
    if (length(model$fixedpars) > 0) {
      draws[, model$fixedpars] <- rep(model$res.t[model$fixedpars, "est"],
                                      each = n_samples)
    }
    draws
  }

  # Make the draws independent of prior RNG use and cache branches.
  set.seed(seed)
  start_time <- Sys.time()
  draws <- list(
    os = draw_coefficients(original_os, "OS"),
    pfs = draw_coefficients(original_pfs, "PFS")
  )
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Draw diagnostics: no draw can fail (the normal approximation is always
  # proper), so the reported rate is the share of draws whose treatment-by-
  # biomarker interaction exceeds a ten-fold time ratio in either direction, an
  # interpretive flag for extreme tails rather than a validity criterion.
  interaction_terms <- function(draw_matrix) {
    nm <- colnames(draw_matrix)
    nm[grepl(":", nm, fixed = TRUE) & grepl("Rx", nm, fixed = TRUE)]
  }
  summarise_draws <- function(draw_matrix, outcome) {
    stats <- t(apply(draw_matrix, 2, function(v) c(
      mean = mean(v), sd = stats::sd(v),
      q2.5 = unname(stats::quantile(v, 0.025)),
      q97.5 = unname(stats::quantile(v, 0.975)),
      min = min(v), max = max(v)
    )))
    data.frame(outcome = outcome, coefficient = rownames(stats), stats,
               row.names = NULL, stringsAsFactors = FALSE)
  }
  coefficient_summary <- rbind(summarise_draws(draws$os, "os"),
                               summarise_draws(draws$pfs, "pfs"))
  extreme_bound <- log(10)
  extreme_flags <- lapply(draws, function(m) {
    terms <- interaction_terms(m)
    if (length(terms) == 0) return(rep(FALSE, nrow(m)))
    apply(abs(m[, terms, drop = FALSE]) > extreme_bound, 1, any)
  })
  n_extreme <- sum(extreme_flags$os | extreme_flags$pfs)

  cat(sprintf("Coefficient draws complete: %d draws per outcome in %.1f seconds; 0 failed\n",
              n_samples, total_time))
  cat(sprintf("Extreme-draw rate (any interaction |coefficient| > log(10)): %d/%d (%.2f%%)\n",
              n_extreme, n_samples, 100 * n_extreme / n_samples))
  cat("Interaction coefficient draws (mean, 2.5%, 97.5%):\n")
  ix <- coefficient_summary[grepl(":", coefficient_summary$coefficient, fixed = TRUE), ]
  for (r in seq_len(nrow(ix))) {
    cat(sprintf("  %-4s %-32s %7.3f [%7.3f, %7.3f]\n", ix$outcome[r], ix$coefficient[r],
                ix$mean[r], ix$q2.5[r], ix$q97.5[r]))
  }

  list(
    method = SAMPLING_METHOD,
    draws = draws,
    original_os = original_os,
    original_pfs = original_pfs,
    n_samples = n_samples,
    n_failed = 0L,
    n_extreme = n_extreme,
    extreme_bound = extreme_bound,
    coefficient_summary = coefficient_summary,
    dist_os = dist_os,
    dist_pfs = dist_pfs,
    seed = seed,
    creation_time = Sys.time()
  )
}

# Return TRUE only when the joint component records the expected seed.
sampling_cache_seed_matches <- function(sampling_models, seed = 123L) {
  joint <- tryCatch(get_joint_sampling_models(sampling_models), error = function(e) NULL)
  !is.null(joint) && identical(joint$seed, seed)
}

# ===============================================================================
# DISTRIBUTION RESOLUTION HELPER
# ===============================================================================
# Resolves the AIC-best parametric distributions from the models object.
# Returns a list with $os and $pfs distribution names.
# ===============================================================================

resolve_best_distributions <- function(models) {
  default_dist <- "weibull"

  if (is.null(models)) {
    cat("WARNING: 'models' is NULL. Using default distribution:", default_dist, "\n")
    return(list(os = default_dist, pfs = default_dist))
  }

  os_dist <- models$best_fit$os_distribution
  pfs_dist <- models$best_fit$pfs_distribution

  if (is.null(os_dist)) {
    cat("WARNING: No selected OS distribution found. Using default:", default_dist, "\n")
    os_dist <- default_dist
  }
  if (is.null(pfs_dist)) {
    cat("WARNING: No selected PFS distribution found. Using default:", default_dist, "\n")
    pfs_dist <- default_dist
  }

  cat("Resolved distributions - OS:", os_dist, "| PFS:", pfs_dist, "\n")
  return(list(os = os_dist, pfs = pfs_dist))
}

# ===============================================================================
# CHECK CACHE OR GENERATE SAMPLING MODELS
# ===============================================================================
# The cache is keyed by a fingerprint of everything that determines its content
# (formulas, fitting data, selected distributions, n_samples, seed and sampling
# method; issue #156). A cache whose fingerprint differs is regenerated, so a
# changed formula, data set or distribution can no longer be served from a
# stale file. The file name still carries n_samples for discoverability.

selected_dists <- resolve_best_distributions(models)
biomarkers_to_sample <- get_biomarkers()

# Every economic biomarker strategy must share one joint formula: the cache
# stores that model once, and the interaction-coefficient extraction and the
# control-arm prediction (issue #151) both assume a single joint fit.
formula_keys <- vapply(biomarkers_to_sample, function(biomarker) {
  paste(paste(deparse(create_biomarker_formula("os", biomarker)), collapse = " "),
        paste(deparse(create_biomarker_formula("pfs", biomarker)), collapse = " "),
        sep = " || ")
}, character(1))
if (length(unique(formula_keys)) != 1L) {
  stop("Biomarker strategies use different survival formulas (",
       paste(unique(formula_keys), collapse = "; "),
       "). The shared joint sampling cache requires one formula; extend ",
       "06_sampling.R before adding biomarker-specific formulas.")
}
joint_formulas <- list(
  os = create_biomarker_formula("os", biomarkers_to_sample[1]),
  pfs = create_biomarker_formula("pfs", biomarkers_to_sample[1])
)

expected_sampling <- sampling_cache_fingerprint(
  formulas = joint_formulas,
  data = data_complete,
  distributions = selected_dists,
  n_samples = n_samples,
  seed = sampling_seed,
  method = SAMPLING_METHOD
)

report_fingerprint_differences <- function(cached_inputs, expected_inputs) {
  keys <- union(names(expected_inputs), names(cached_inputs))
  for (key in keys) {
    if (!identical(cached_inputs[[key]], expected_inputs[[key]])) {
      cat("  - ", key, " differs\n", sep = "")
    }
  }
}

sampling_models <- NULL
if (file.exists(cache_file)) {
  cat("\n=== Loading sampling models from cache ===\n")
  cat("Cache file:", cache_file, "\n")
  sampling_models <- readRDS(cache_file)

  if (!is.list(sampling_models) || is.null(sampling_models$joint) ||
      is.null(sampling_models$fingerprint)) {
    cat("Cache predates issue #156 (per-biomarker bootstrap components or no ",
        "fingerprint). Regenerating...\n", sep = "")
    sampling_models <- NULL
  } else if (!identical(sampling_models$fingerprint, expected_sampling$fingerprint)) {
    cat("Cache fingerprint does not match the current inputs. Regenerating...\n")
    report_fingerprint_differences(sampling_models$fingerprint_inputs,
                                   expected_sampling$inputs)
    sampling_models <- NULL
  } else if (!sampling_cache_seed_matches(sampling_models, sampling_seed)) {
    cat("Cache RNG seed is missing or does not match the requested seed (",
        sampling_seed, "). Regenerating...\n", sep = "")
    sampling_models <- NULL
  } else {
    joint_models <- get_joint_sampling_models(sampling_models)
    cat("Cache loaded successfully!\n")
    cat("- Method:", joint_models$method, "\n")
    cat("- Joint-model draws:", joint_models$n_samples, "\n")
    cat("- Biomarkers sharing the joint model:",
        paste(sampling_models$biomarkers, collapse = ", "), "\n")
    cat("- Created:", format(joint_models$creation_time), "\n")
    cat("- Distributions: OS", joint_models$dist_os,
        "| PFS", joint_models$dist_pfs, "\n")
    cat("- RNG seed:", joint_models$seed, "\n")
    cat("- Fingerprint:", sampling_models$fingerprint, "\n")
    rm(joint_models)
  }
} else {
  cat("\n=== No cache found - will generate sampling models ===\n")
}

# ===============================================================================
# SURVIVAL MODEL COEFFICIENT SAMPLING (if not cached)
# ===============================================================================

if (is.null(sampling_models)) {
  cat("\n=== Starting multivariate-normal coefficient sampling ===\n")
  cat("Results will be cached to:", cache_file, "\n\n")

  # One joint model (age, sex, treatment and the biomarker-by-treatment
  # interactions) serves every strategy. No separate control model is fitted:
  # the control arm is this same fit predicted with Rx = control (issue #151).
  joint_component <- sample_survival_coefficients(
    formula_os = joint_formulas$os,
    formula_pfs = joint_formulas$pfs,
    data = data_complete,
    dist_os = selected_dists$os,
    dist_pfs = selected_dists$pfs,
    n_samples = n_samples,
    seed = sampling_seed
  )

  sampling_models <- list(
    joint = joint_component,
    biomarkers = biomarkers_to_sample,
    fingerprint = expected_sampling$fingerprint,
    fingerprint_inputs = expected_sampling$inputs,
    creation_time = Sys.time()
  )
  rm(joint_component)

  cat("\n=== Saving sampling models to cache ===\n")
  tryCatch({
    saveRDS(sampling_models, cache_file)
    cat("Sampling models saved successfully to:", cache_file, "\n")
    cat("File size:", format(file.info(cache_file)$size / 1024^2, digits = 2), "MB\n")
  }, error = function(e) {
    cat("WARNING: Failed to save cache:", conditionMessage(e), "\n")
  })
}
rm(expected_sampling, formula_keys, report_fingerprint_differences)

# Print confirmation of model structure
cat("\n=== Sampling model structure ready ===\n")
cat("- Model:", model_config$label, "\n")
cat("- Biomarker strategies sharing the joint model:",
    paste(sampling_models$biomarkers, collapse = ", "), "\n")
cat("- Formulas defined in: scripts/R/functions/model_configs.R\n")
cat("- Number of coefficient draws:",
    get_joint_sampling_models(sampling_models)$n_samples, "\n")
cat("- Draws: multivariate normal around the fitted joint models; OS and PFS",
    "drawn independently (issue #156)\n")
cat("- Control arm: joint model predicted with Rx = control (issue #151)\n")

# ===============================================================================
# PARAMETER DISTRIBUTIONS FOR UNCERTAINTY ANALYSIS (NON-SURVIVAL PARAMETERS)
# ===============================================================================

# Survival curves are excluded because they come from correlated resampling.
parameter_config <- configure_parameter_distributions(l_params_base)
param_distributions <- parameter_config$distributions
param_groups <- parameter_config$groups
parameter_spec <- parameter_config$spec
rm(parameter_config)

# ===============================================================================
# POPULATION AVERAGING FUNCTION FOR PSA
# ===============================================================================
# With biomarker_name = NULL, predicts the control curve over all patients with
# Rx forced to the control level, i.e. the joint model's standard-of-care
# prediction, mirroring generate_population_averaged_predictions() in the base
# case (issue #151). Otherwise predicts the biomarker-positive/experimental and
# biomarker-negative/control subgroup curves using one shared code path.

generate_psa_population_averaged_predictions <- function(sampling_model_list,
                                                         biomarker_name = NULL,
                                                         outcome = "os",
                                                         sample_idx = 1,
                                                         data_original,
                                                         time_points) {
  sampled_model <- tryCatch(
    sampled_survival_models(sampling_model_list, sample_idx),
    error = function(e) NULL
  )
  if (is.null(sampled_model)) {
    warning("Sampled model ", sample_idx, " is unavailable - returning NULL")
    return(NULL)
  }

  model_obj <- sampled_model[[outcome]]$model
  if (is.null(model_obj)) {
    warning("Model object for ", outcome, " is NULL in sample ", sample_idx)
    return(NULL)
  }

  average_prediction <- function(subgroup_data, rx_level = NULL, label) {
    if (nrow(subgroup_data) == 0) {
      warning("No patients in ", label, " subgroup")
      return(rep(NA_real_, length(time_points)))
    }
    if (!is.null(rx_level)) {
      subgroup_data$Rx <- factor(rx_level, levels = levels(data_original$Rx))
    }
    tryCatch({
      prediction <- predict(
        model_obj, newdata = subgroup_data,
        type = "survival", times = time_points
      )
      rowMeans(extract_all_survival_probabilities(prediction), na.rm = TRUE)
    }, error = function(e) {
      warning("Prediction failed for ", label, " in sim ", sample_idx,
              ": ", conditionMessage(e))
      rep(NA_real_, length(time_points))
    })
  }

  if (is.null(biomarker_name)) {
    return(average_prediction(
      data_original,
      rx_level = levels(data_original$Rx)[1],
      label = "control"
    ))
  }

  status <- as.numeric(as.character(data_original[[biomarker_name]]))
  subgroup_spec <- list(
    positive = list(status = 1, rx = levels(data_original$Rx)[2], suffix = "+"),
    negative = list(status = 0, rx = levels(data_original$Rx)[1], suffix = "-")
  )
  lapply(subgroup_spec, function(subgroup) {
    average_prediction(
      data_original[status == subgroup$status, , drop = FALSE],
      rx_level = subgroup$rx,
      label = paste0(biomarker_name, subgroup$suffix)
    )
  })
}

# Report what the PSA actually samples. Since issue #154 the specification also
# holds parameters that are fixed in the PSA (unit prices) and one that is
# derived rather than drawn (u_p), so the summary is built per group instead of
# assuming a single cost CV and a single utility CV.
sampled_spec <- parameter_spec[parameter_spec$psa & !parameter_spec$derived, ,
                               drop = FALSE]
fixed_spec <- parameter_spec[!parameter_spec$psa, , drop = FALSE]
derived_spec <- parameter_spec[parameter_spec$derived, , drop = FALSE]

cat("\n=== Parameter distributions configured ===\n")
for (group_name in unique(sampled_spec$group)) {
  rows <- sampled_spec[sampled_spec$group == group_name, , drop = FALSE]
  cat(sprintf("- %s: %s (CV = %s): %s\n",
              group_name,
              paste(unique(rows$distribution), collapse = "/"),
              paste(unique(rows$cv), collapse = "/"),
              paste(rows$parameter, collapse = ", ")))
}
if (nrow(derived_spec) > 0) {
  cat("- Derived in the PSA (not drawn):",
      paste(derived_spec$parameter, collapse = ", "), "\n")
}
if (nrow(fixed_spec) > 0) {
  cat("- Fixed in the PSA, varied in the DSA only (issue #154):",
      paste(fixed_spec$parameter, collapse = ", "), "\n")
}
cat("- Survival parameters: multivariate-normal coefficient draws (n =",
    get_joint_sampling_models(sampling_models)$n_samples, ")\n")
cat("- PSA population averaging function ready\n")
cat("\nReady for PSA analysis\n")
