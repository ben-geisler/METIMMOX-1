# One economic target population and common standardisation weights (#166).
# Clinical-only variables, including TLR and scan dates, do not select this cohort.
economic_prediction_population <- function(data, biomarkers = get_biomarkers()) {
  required <- c("Age", "sex", "Rx", biomarkers, "OSwk", "Death", "PFSwk", "Progression")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing prediction-population columns: ", paste(missing, collapse = ", "))
  population <- data[complete.cases(data[, required, drop = FALSE]), , drop = FALSE]
  if (!nrow(population)) stop("The economic prediction population is empty.")
  population
}

joint_biomarker_population <- function(population, biomarkers = get_biomarkers()) {
  status <- as.data.frame(lapply(population[biomarkers], function(x) as.numeric(as.character(x))))
  if (!nrow(status) || anyNA(status) || any(!as.matrix(status) %in% 0:1))
    stop("Prediction population must have complete binary biomarkers.")
  cells <- expand.grid(setNames(rep(list(0:1), length(biomarkers)), biomarkers))
  key <- function(x) apply(x, 1, paste0, collapse = "")
  cells <- cells[order(key(cells)), , drop = FALSE]
  rownames(cells) <- key(cells)
  index <- match(key(status), rownames(cells))
  counts <- setNames(tabulate(index, nrow(cells)), rownames(cells))
  list(cells = cells, index = index, counts = counts, probabilities = counts / sum(counts))
}

# Marginal DSA changes rake the *joint* distribution, retaining its odds ratios.
# The same weights are applied to control and every guided strategy. In the PSA
# the joint masses are sampled together and the marginals are derived from them.
target_population_weights <- function(population, prevalences = NULL,
                                      joint_probabilities = NULL,
                                      biomarkers = get_biomarkers()) {
  joint <- joint_biomarker_population(population, biomarkers)
  mass <- joint$probabilities
  if (!is.null(joint_probabilities)) {
    if (!setequal(names(joint_probabilities), names(mass))) stop("Joint probabilities must name every biomarker cell.")
    mass <- joint_probabilities[names(mass)]
    if (any(!is.finite(mass)) || any(mass < 0) || abs(sum(mass) - 1) > 1e-10)
      stop("Joint probabilities must be non-negative and sum to one.")
    if (any(mass[joint$counts == 0] > 0)) stop("Cannot assign mass to an unobserved biomarker cell.")
  }
  if (!is.null(prevalences)) {
    if (is.null(names(prevalences)) || !all(biomarkers %in% names(prevalences)))
      stop("prevalences must be a named vector containing: ", paste(biomarkers, collapse = ", "))
    targets <- prevalences[biomarkers]
    if (any(!is.finite(targets)) || any(targets < 0 | targets > 1))
      stop("Prevalences must be finite probabilities between zero and one.")
    for (iteration in seq_len(10000L)) {
      if (max(abs(colSums(as.matrix(joint$cells) * mass) - targets)) < 1e-12) break
      for (biomarker in biomarkers) {
        for (value in 0:1) {
          rows <- joint$cells[[biomarker]] == value
          wanted <- if (value == 1) targets[[biomarker]] else 1 - targets[[biomarker]]
          current <- sum(mass[rows])
          if (wanted > 0 && current == 0) stop("Requested prevalences are incompatible with observed joint support.")
          mass[rows] <- if (current == 0) 0 else mass[rows] * wanted / current
        }
      }
    }
    if (max(abs(colSums(as.matrix(joint$cells) * mass) - targets)) >= 1e-10)
      stop("Joint prevalence standardisation did not converge.")
  }
  unname(mass[joint$index] / joint$counts[joint$index])
}

model_population_weights <- function(params, population) {
  biomarkers <- get_biomarkers()
  keys <- paste0("p_joint_", c("00", "01", "10", "11"))
  joint <- if (identical(biomarkers, c("crp", "tmb_braf")) && all(keys %in% names(params)))
    setNames(unlist(params[keys], use.names = FALSE), c("00", "01", "10", "11")) else NULL
  target_population_weights(population,
    setNames(vapply(biomarkers, function(b) params[[paste0("p_", b)]], numeric(1)), biomarkers), joint)
}

weighted_survival_average <- function(curves, weights) {
  if (!is.matrix(curves) || ncol(curves) != length(weights) || any(!is.finite(curves)) ||
      any(!is.finite(weights)) || any(weights < 0)) stop("Invalid population survival predictions or weights.")
  # A zero-mass subgroup is not used in the mixture; retain a finite curve.
  if (sum(weights) == 0) weights <- rep(1, length(weights))
  as.vector(curves %*% (weights / sum(weights)))
}

standardize_population_curves <- function(curves, population, weights) {
  average <- function(part) list(os = weighted_survival_average(part$os, weights[part$rows]),
                                 pfs = weighted_survival_average(part$pfs, weights[part$rows]))
  control <- get_control_strategy()
  out <- setNames(list(c(list(strategy = control), average(curves[[control]]))), control)
  for (biomarker in get_biomarkers()) {
    positive <- average(curves[[biomarker]]$positive)
    negative <- average(curves[[biomarker]]$negative)
    prevalence <- sum(weights[curves[[biomarker]]$positive$rows]) / sum(weights)
    out[[biomarker]] <- list(strategy = biomarker,
      os = prevalence * positive$os + (1 - prevalence) * negative$os,
      pfs = prevalence * positive$pfs + (1 - prevalence) * negative$pfs,
      biomarker_positive = positive, biomarker_negative = negative, prevalence = prevalence)
  }
  attr(out, "population_curves") <- curves
  attr(out, "prediction_population") <- population
  attr(out, "population_weights") <- weights
  out
}

set_population_predictions <- function(params, predictions) {
  control <- get_control_strategy()
  for (outcome in c("os", "pfs")) {
    key <- paste0("p_", outcome)
    suffix <- toupper(outcome)
    params[[key]] <- setNames(list(predictions[[control]][[outcome]]), paste0(control, "_", suffix))
    for (biomarker in get_biomarkers()) {
      params[[key]][[paste0(biomarker, "_pos_", suffix)]] <- predictions[[biomarker]]$biomarker_positive[[outcome]]
      params[[key]][[paste0(biomarker, "_neg_", suffix)]] <- predictions[[biomarker]]$biomarker_negative[[outcome]]
      params[[key]][[paste0(biomarker, "_weighted_", suffix)]] <- predictions[[biomarker]][[outcome]]
      params[[paste0("p_", biomarker)]] <- predictions[[biomarker]]$prevalence
    }
  }
  params$population_curves <- attr(predictions, "population_curves")
  population <- attr(predictions, "prediction_population")
  columns <- c("ID", "Age", "sex", "Rx", get_biomarkers(), "OSwk", "Death", "PFSwk", "Progression")
  params$prediction_population <- population[, intersect(columns, names(population)), drop = FALSE]
  params$population_weights <- attr(predictions, "population_weights")
  params
}
