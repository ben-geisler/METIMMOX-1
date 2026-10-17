# One economic target population and common standardisation weights (#166).
# Clinical-only variables, including TLR and scan dates, do not select this cohort.
economic_prediction_population <- function(data, biomarkers = get_biomarkers()) {
  v_required <- c("Age", "sex", "Rx", biomarkers, "OSwk", "Death", "PFSwk", "Progression")
  v_missing <- setdiff(v_required, names(data))
  if (length(v_missing)) stop("Missing prediction-population columns: ", paste(v_missing, collapse = ", "))
  df_population <- data[complete.cases(data[, v_required, drop = FALSE]), , drop = FALSE]
  if (!nrow(df_population)) stop("The economic prediction population is empty.")
  df_population
}

joint_biomarker_population <- function(population, biomarkers = get_biomarkers()) {
  df_status <- as.data.frame(lapply(population[biomarkers], function(x) as.numeric(as.character(x))))
  if (!nrow(df_status) || anyNA(df_status) || any(!as.matrix(df_status) %in% 0:1))
    stop("Prediction population must have complete binary biomarkers.")
  df_cells <- expand.grid(setNames(rep(list(0:1), length(biomarkers)), biomarkers))
  key <- function(x) apply(x, 1, paste0, collapse = "")
  df_cells <- df_cells[order(key(df_cells)), , drop = FALSE]
  rownames(df_cells) <- key(df_cells)
  v_index <- match(key(df_status), rownames(df_cells))
  v_counts <- setNames(tabulate(v_index, nrow(df_cells)), rownames(df_cells))
  list(cells = df_cells, index = v_index, counts = v_counts, probabilities = v_counts / sum(v_counts))
}

# Marginal DSA changes rake the *joint* distribution, retaining its odds ratios.
# The same weights are applied to control and every guided strategy. In the PSA
# the joint masses are sampled together and the marginals are derived from them.
target_population_weights <- function(population, prevalences = NULL,
                                      joint_probabilities = NULL,
                                      biomarkers = get_biomarkers()) {
  joint <- joint_biomarker_population(population, biomarkers)
  v_mass <- joint$probabilities
  if (!is.null(joint_probabilities)) {
    if (!setequal(names(joint_probabilities), names(v_mass))) stop("Joint probabilities must name every biomarker cell.")
    v_mass <- joint_probabilities[names(v_mass)]
    if (any(!is.finite(v_mass)) || any(v_mass < 0) || abs(sum(v_mass) - 1) > 1e-10)
      stop("Joint probabilities must be non-negative and sum to one.")
    if (any(v_mass[joint$counts == 0] > 0)) stop("Cannot assign mass to an unobserved biomarker cell.")
  }
  if (!is.null(prevalences)) {
    if (is.null(names(prevalences)) || !all(biomarkers %in% names(prevalences)))
      stop("prevalences must be a named vector containing: ", paste(biomarkers, collapse = ", "))
    v_targets <- prevalences[biomarkers]
    if (any(!is.finite(v_targets)) || any(v_targets < 0 | v_targets > 1))
      stop("Prevalences must be finite probabilities between zero and one.")
    for (iteration in seq_len(10000L)) {
      if (max(abs(colSums(as.matrix(joint$cells) * v_mass) - v_targets)) < 1e-12) break
      for (biomarker in biomarkers) {
        for (value in 0:1) {
          v_rows <- joint$cells[[biomarker]] == value
          wanted <- if (value == 1) v_targets[[biomarker]] else 1 - v_targets[[biomarker]]
          current <- sum(v_mass[v_rows])
          if (wanted > 0 && current == 0) stop("Requested prevalences are incompatible with observed joint support.")
          v_mass[v_rows] <- if (current == 0) 0 else v_mass[v_rows] * wanted / current
        }
      }
    }
    if (max(abs(colSums(as.matrix(joint$cells) * v_mass) - v_targets)) >= 1e-10)
      stop("Joint prevalence standardisation did not converge.")
  }
  unname(v_mass[joint$index] / joint$counts[joint$index])
}

model_population_weights <- function(params, population) {
  v_biomarkers <- get_biomarkers()
  v_keys <- paste0("p_joint_", c("00", "01", "10", "11"))
  v_joint <- if (identical(v_biomarkers, c("crp", "tmb_braf")) && all(v_keys %in% names(params)))
    setNames(unlist(params[v_keys], use.names = FALSE), c("00", "01", "10", "11")) else NULL
  target_population_weights(population,
    setNames(vapply(v_biomarkers, function(b) params[[paste0("p_", b)]], numeric(1)), v_biomarkers), v_joint)
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

# ===============================================================================
# JOINT BIOMARKER CELLS (issue #200, paper 3 step 1b)
# ===============================================================================
# Cells are the rows of joint_biomarker_population()$cells, keyed by the
# biomarker values in get_biomarkers() order ("00", "01", "10", "11" for week-4
# CRP x TMB/BRAF, the keys of p_joint_00 ... p_joint_11). A cell curve is the
# weighted average of its patients' curves under one arm; a cell's mass is the
# sum of its patients' weights. Every strategy curve is a mass-weighted
# combination of cell curves (strategy_curves_from_cells()), so the PSA
# predicts each patient once per arm and endpoint.

# Arms in Rx-level order: level 1 of Rx is control, level 2 experimental.
CELL_ARMS <- c("control", "experimental")

#' Column names of the cell x arm curve matrices
#'
#' @param cells Cell keys (e.g. "00").
#' @param arms Arms; default both (CELL_ARMS).
#' @return Character vector `<cell>_<arm>`, cells varying slowest.
cell_curve_columns <- function(cells, arms = CELL_ARMS) {
  as.vector(t(outer(cells, arms, paste, sep = "_")))
}

#' Joint biomarker cell masses of a weighted prediction population
#'
#' @param population Prediction population (complete binary biomarkers).
#' @param weights Patient weights (any non-negative scale); NULL gives equal
#'   weights.
#' @param biomarkers Biomarker columns defining the cells.
#' @return List with `mass` (named by cell key, sums to one), `counts`
#'   (patients per cell), `index` (cell of each patient) and `cells` (the cell
#'   definition data frame).
cell_masses <- function(population, weights = NULL, biomarkers = get_biomarkers()) {
  joint <- joint_biomarker_population(population, biomarkers)
  if (is.null(weights)) weights <- rep(1, nrow(population))
  if (length(weights) != nrow(population) || any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) stop("Invalid population weights for the cell masses.")
  v_mass <- vapply(seq_len(nrow(joint$cells)), function(k) sum(weights[joint$index == k]),
                   numeric(1)) / sum(weights)
  list(mass = setNames(v_mass, rownames(joint$cells)), counts = joint$counts,
       index = joint$index, cells = joint$cells)
}

#' Average patient curves within each joint biomarker cell
#'
#' @param m_curves Matrix of patient curves (time x patient).
#' @param index Cell index of each patient (cell_masses()$index).
#' @param cells Cell keys in index order.
#' @param weights Patient weights; within a cell of zero total weight the
#'   patients are weighted equally, as in weighted_survival_average().
#' @return Matrix time x cell (NA columns for cells without patients).
average_within_cells <- function(m_curves, index, cells, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, ncol(m_curves))
  m_out <- matrix(NA_real_, nrow(m_curves), length(cells), dimnames = list(NULL, cells))
  for (k in seq_along(cells)) {
    v_rows <- which(index == k)
    if (length(v_rows)) m_out[, k] <- weighted_survival_average(m_curves[, v_rows, drop = FALSE],
                                                                 weights[v_rows])
  }
  m_out
}

#' Strategy curves as mass-weighted combinations of cell curves
#'
#' Control is every cell under control; a biomarker-guided strategy is its
#' positive cells under the experimental arm (biomarker-positive subgroup) and
#' its negative cells under control (biomarker-negative subgroup). With masses
#' equal to the summed patient weights this is algebraically the patient-level
#' standardisation of standardize_population_curves(); a subgroup of zero mass
#' falls back to equal patient weights, as weighted_survival_average() does.
#'
#' @param cell_curves List with `os` and `pfs` matrices (time x cell_arm,
#'   columns from cell_curve_columns()).
#' @param cell_weights Output of cell_masses() (`mass`, `counts`, `cells`).
#' @param biomarkers Biomarker ids of the guided strategies.
#' @return List keyed by strategy in the shape of
#'   standardize_population_curves(): control `os`/`pfs`; per biomarker `os`,
#'   `pfs`, `biomarker_positive`, `biomarker_negative` and `prevalence`.
strategy_curves_from_cells <- function(cell_curves, cell_weights, biomarkers = get_biomarkers()) {
  v_mass <- cell_weights$mass
  v_counts <- cell_weights$counts
  combine <- function(v_cells, arm, outcome) {
    v_cells <- v_cells[v_counts[v_cells] > 0]
    if (!length(v_cells)) return(rep(NA_real_, nrow(cell_curves[[outcome]])))
    v_w <- v_mass[v_cells]
    if (sum(v_w) == 0) v_w <- v_counts[v_cells]
    m_part <- cell_curves[[outcome]][, cell_curve_columns(v_cells, arm), drop = FALSE]
    as.vector(m_part %*% (v_w / sum(v_w)))
  }
  v_all <- names(v_mass)
  control <- get_control_strategy()
  out <- setNames(list(list(strategy = control, os = combine(v_all, "control", "os"),
                            pfs = combine(v_all, "control", "pfs"))), control)
  for (biomarker in biomarkers) {
    v_positive <- v_all[cell_weights$cells[[biomarker]] == 1]
    v_negative <- v_all[cell_weights$cells[[biomarker]] == 0]
    positive <- list(os = combine(v_positive, "experimental", "os"),
                     pfs = combine(v_positive, "experimental", "pfs"))
    negative <- list(os = combine(v_negative, "control", "os"),
                     pfs = combine(v_negative, "control", "pfs"))
    prevalence <- sum(v_mass[v_positive])
    out[[biomarker]] <- list(strategy = biomarker,
      os = prevalence * positive$os + (1 - prevalence) * negative$os,
      pfs = prevalence * positive$pfs + (1 - prevalence) * negative$pfs,
      biomarker_positive = positive, biomarker_negative = negative, prevalence = prevalence)
  }
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
  df_population <- attr(predictions, "prediction_population")
  v_columns <- c("ID", "Age", "sex", "Rx", get_biomarkers(), "OSwk", "Death", "PFSwk", "Progression")
  params$prediction_population <- df_population[, intersect(v_columns, names(df_population)), drop = FALSE]
  params$population_weights <- attr(predictions, "population_weights")
  params
}
