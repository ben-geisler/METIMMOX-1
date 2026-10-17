# ===============================================================================
# SURVIVAL COEFFICIENTS AND CELL EFFECTS AS EVPPI PARAMETERS (issue #199)
# ===============================================================================
# extract_interaction_coefficients() (R/evppi_functions.R) returns only the four
# treatment-by-biomarker interactions. Paper 3 needs every coefficient of both
# joint survival models per PSA draw, aligned by the cached model that produced
# the draw (`model_idx`, #156), and the derived treatment effect in each joint
# biomarker cell:
#
#   cell (crp = 0, tmb_braf = 0): Rx
#   cell (1, 0): Rx + Rx:crp
#   cell (0, 1): Rx + Rx:tmb_braf
#   cell (1, 1): Rx + Rx:crp + Rx:tmb_braf
#
# on flexsurv's optimisation scale (log rate for the gamma covariate effects).
# Cell codes list the biomarkers in get_biomarkers() order. The functions live
# in this file rather than in R/evppi_functions.R because that file enters the
# EVPPI and scenario cache fingerprints (calculation_identity()): editing it
# would invalidate both caches without changing any result.
#
# Coefficients are classified structurally, never by position: an interaction
# is a term containing ":" with an Rx part and a biomarker part (either order,
# e.g. "RxExperimental arm:crp1"); the Rx main effect is a single Rx term; a
# prognostic effect is a single biomarker term; everything else (shape, rate,
# Age, sex) is a baseline-survival parameter.
#
# EVPPI groups (survival_evppi_groups()): survival_all (every coefficient of
# both endpoints), cell_effect_<code> (OS and PFS of one cell), rx_main,
# prognostic and baseline_survival. EVPPIs are never summed, differenced or
# expressed as shares; only survival_all is compared with EVPI.
#
# Requires R/model_configs.R (get_biomarkers()) and R/prediction_functions.R
# (get_joint_sampling_models()) for the cache layout; R/evppi_functions.R for
# run_extended_evppi().
# ===============================================================================

#' Endpoints whose coefficients are extracted
SURVIVAL_EVPPI_OUTCOMES <- c("os", "pfs")

#' Classify survival-model coefficient names by role
#'
#' @param coef_names Coefficient names as flexsurvreg emits them.
#' @param outcome Endpoint label ("os" or "pfs").
#' @param biomarkers Biomarker identifiers (get_biomarkers()).
#' @return Data frame with one row per coefficient: `column` (the EVPPI input
#'   name, `coef_<outcome>_<short>`), `outcome`, `coefficient`, `role`
#'   ("baseline", "rx_main", "prognostic", "interaction") and `biomarker`
#'   (NA unless prognostic or interaction).
survival_coefficient_map <- function(coef_names, outcome, biomarkers = get_biomarkers()) {
  if (length(coef_names) == 0 || anyNA(coef_names) || anyDuplicated(coef_names)) {
    stop("Coefficient names must be unique and non-missing")
  }
  # Longest matching biomarker identifier, so that one identifier cannot match
  # another's prefix.
  v_biomarkers <- biomarkers[order(-nchar(biomarkers))]
  term_biomarker <- function(term) {
    v_hit <- v_biomarkers[startsWith(term, v_biomarkers)]
    if (length(v_hit) == 0) NA_character_ else v_hit[1]
  }
  is_rx <- function(term) grepl("^Rx", term)

  df_rows <- lapply(coef_names, function(coef_name) {
    v_terms <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
    v_term_biomarkers <- vapply(v_terms, term_biomarker, character(1))
    v_term_rx <- vapply(v_terms, is_rx, logical(1))
    role <- NA_character_
    biomarker <- NA_character_
    short <- gsub("[^A-Za-z0-9]+", "_", coef_name)
    if (length(v_terms) == 1 && v_term_rx) {
      role <- "rx_main"
      short <- "rx"
    } else if (length(v_terms) == 1 && !is.na(v_term_biomarkers)) {
      role <- "prognostic"
      biomarker <- v_term_biomarkers
      short <- biomarker
    } else if (length(v_terms) == 2 && sum(v_term_rx) == 1 &&
               !is.na(v_term_biomarkers[!v_term_rx])) {
      role <- "interaction"
      biomarker <- unname(v_term_biomarkers[!v_term_rx])
      short <- paste0("rx_x_", biomarker)
    } else if (length(v_terms) == 1) {
      role <- "baseline"
    } else {
      stop("Unrecognised survival coefficient '", coef_name, "' (", outcome, ")")
    }
    data.frame(column = paste0("coef_", outcome, "_", short), outcome = outcome,
               coefficient = coef_name, role = role, biomarker = biomarker,
               stringsAsFactors = FALSE)
  })
  df_map <- do.call(rbind, df_rows)
  if (anyDuplicated(df_map$column)) {
    stop("Survival coefficients map to duplicate columns: ",
         paste(df_map$column[duplicated(df_map$column)], collapse = ", "))
  }
  if (sum(df_map$role == "rx_main") != 1) {
    stop("Expected exactly one treatment main effect in the ", outcome, " model")
  }
  df_map
}

#' Coefficient-draw matrix of one endpoint from a sampling component
#'
#' @param component Joint sampling component (get_joint_sampling_models()):
#'   either `draws$<outcome>` matrices (#156) or legacy `samples` lists.
#' @param outcome Endpoint ("os" or "pfs").
#' @return Numeric matrix, one row per cached draw, named columns.
survival_draw_matrix <- function(component, outcome) {
  if (!is.null(component$draws)) {
    m_draws <- component$draws[[outcome]]
    if (is.null(m_draws)) stop("Sampling component has no ", outcome, " draws")
    m_draws <- as.matrix(m_draws)
  } else if (!is.null(component$samples)) {
    v_names <- names(component$samples[[1]][[outcome]]$coefficients)
    m_draws <- do.call(rbind, lapply(component$samples, function(sample_i) {
      v_coefs <- sample_i[[outcome]]$coefficients
      if (is.null(v_coefs) || !setequal(names(v_coefs), v_names)) {
        return(stats::setNames(rep(NA_real_, length(v_names)), v_names))
      }
      v_coefs[v_names]
    }))
  } else {
    stop("Sampling component has neither `draws` nor legacy `samples`")
  }
  if (is.null(colnames(m_draws))) stop("The ", outcome, " draws have no coefficient names")
  m_draws
}

#' Every survival coefficient of both endpoints per PSA draw
#'
#' @param sampling_models Sampling cache (script 06) or a joint component.
#' @param n_sim Number of cached coefficient draws that `model_idx` indexes
#'   (the cache must hold at least this many).
#' @param model_idx Cached draw used by each PSA row (`psa_params$model_idx`);
#'   defaults to `seq_len(n_sim)`.
#' @param biomarkers Biomarker identifiers.
#' @return Data frame with `length(model_idx)` rows and one column per
#'   coefficient and endpoint (`coef_<outcome>_<short>`, see
#'   survival_coefficient_map()); attributes `coefficient_map` and
#'   `model_idx`. Stops on invalid indices or missing draws.
extract_survival_coefficients <- function(sampling_models, n_sim,
                                          model_idx = seq_len(n_sim),
                                          biomarkers = get_biomarkers()) {
  if (is.null(sampling_models)) stop("sampling_models is NULL")
  component <- if (!is.null(sampling_models$draws) || !is.null(sampling_models$samples)) {
    sampling_models
  } else {
    get_joint_sampling_models(sampling_models)
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || !is.finite(n_sim) || n_sim < 1 ||
      n_sim != floor(n_sim)) {
    stop("n_sim must be one positive integer")
  }
  if (!is.numeric(model_idx) || length(model_idx) == 0 || anyNA(model_idx) ||
      any(model_idx != floor(model_idx)) || any(model_idx < 1 | model_idx > n_sim)) {
    stop("model_idx must be integers in 1..", n_sim)
  }

  v_columns <- list()
  df_maps <- list()
  for (outcome in SURVIVAL_EVPPI_OUTCOMES) {
    m_draws <- survival_draw_matrix(component, outcome)
    if (nrow(m_draws) < n_sim) {
      stop("The ", outcome, " draws hold ", nrow(m_draws), " rows; n_sim is ", n_sim)
    }
    df_map <- survival_coefficient_map(colnames(m_draws), outcome, biomarkers)
    m_selected <- m_draws[model_idx, df_map$coefficient, drop = FALSE]
    if (anyNA(m_selected)) {
      stop("Missing ", outcome, " coefficient draws for ",
           sum(!stats::complete.cases(m_selected)), " PSA rows")
    }
    for (j in seq_len(nrow(df_map))) {
      v_columns[[df_map$column[j]]] <- unname(m_selected[, j])
    }
    df_maps[[outcome]] <- df_map
  }
  df_out <- as.data.frame(v_columns, stringsAsFactors = FALSE)
  df_map_all <- do.call(rbind, df_maps)
  rownames(df_map_all) <- NULL
  attr(df_out, "coefficient_map") <- df_map_all
  attr(df_out, "model_idx") <- as.integer(model_idx)
  df_out
}

#' Add the per-cell treatment effects to extracted survival coefficients
#'
#' @param df_coefs Data frame from extract_survival_coefficients().
#' @param biomarkers Biomarker identifiers; the cell code lists their 0/1
#'   values in this order.
#' @return `df_coefs` with `cell_effect_<code>_<outcome>` columns appended
#'   (Rx main effect plus the interactions of the biomarkers that are positive
#'   in the cell); the `coefficient_map` attribute gains rows with role
#'   "cell_effect" and the cell code in `biomarker`.
add_cell_effects <- function(df_coefs, biomarkers = get_biomarkers()) {
  df_map <- attr(df_coefs, "coefficient_map")
  if (is.null(df_map)) stop("df_coefs must come from extract_survival_coefficients()")
  df_cells <- expand.grid(rep(list(0:1), length(biomarkers)))
  names(df_cells) <- biomarkers
  v_cell_codes <- apply(df_cells, 1, paste, collapse = "")
  df_new_map <- list()
  for (outcome in SURVIVAL_EVPPI_OUTCOMES) {
    df_out_map <- df_map[df_map$outcome == outcome, , drop = FALSE]
    rx_col <- df_out_map$column[df_out_map$role == "rx_main"]
    for (k in seq_len(nrow(df_cells))) {
      v_positive <- biomarkers[unlist(df_cells[k, ]) == 1]
      v_int_cols <- vapply(v_positive, function(b) {
        col <- df_out_map$column[df_out_map$role == "interaction" & df_out_map$biomarker == b]
        if (length(col) != 1) stop("No unique ", outcome, " interaction for biomarker ", b)
        col
      }, character(1))
      v_effect <- df_coefs[[rx_col]]
      for (col in v_int_cols) v_effect <- v_effect + df_coefs[[col]]
      cell_col <- paste0("cell_effect_", v_cell_codes[k], "_", outcome)
      df_coefs[[cell_col]] <- v_effect
      df_new_map[[length(df_new_map) + 1]] <- data.frame(
        column = cell_col, outcome = outcome,
        coefficient = paste(c(df_out_map$coefficient[df_out_map$role == "rx_main"],
                              df_out_map$coefficient[match(v_int_cols, df_out_map$column)]),
                            collapse = " + "),
        role = "cell_effect", biomarker = v_cell_codes[k], stringsAsFactors = FALSE)
    }
  }
  df_map <- rbind(df_map, do.call(rbind, df_new_map))
  rownames(df_map) <- NULL
  attr(df_coefs, "coefficient_map") <- df_map
  df_coefs
}

#' EVPPI parameter groups of the survival coefficients
#'
#' @param df_map Coefficient map (attribute `coefficient_map` after
#'   add_cell_effects()).
#' @return Named list: `survival_all` (every coefficient of both endpoints,
#'   cell effects excluded because they are functions of the coefficients),
#'   `cell_effect_<code>` (OS and PFS), `rx_main`, `prognostic`,
#'   `baseline_survival`.
survival_evppi_groups <- function(df_map) {
  v_coef <- df_map$role != "cell_effect"
  groups <- list(survival_all = df_map$column[v_coef])
  for (code in unique(df_map$biomarker[df_map$role == "cell_effect"])) {
    groups[[paste0("cell_effect_", code)]] <-
      df_map$column[df_map$role == "cell_effect" & df_map$biomarker == code]
  }
  groups$rx_main <- df_map$column[df_map$role == "rx_main"]
  groups$prognostic <- df_map$column[df_map$role == "prognostic"]
  groups$baseline_survival <- df_map$column[df_map$role == "baseline"]
  groups
}

#' PSA parameters with every survival coefficient and cell effect appended
#'
#' Appends the coefficients of the cached model behind each PSA row
#' (`psa_params$model_idx`) and checks that the interaction columns already
#' bound by the PSA (`b_<biomarker>_rx_<outcome>`, when present) are exactly
#' the extracted interactions, which proves the alignment.
#'
#' @param psa_params PSA parameter data frame with `model_idx`.
#' @param sampling_models Sampling cache (script 06).
#' @param biomarkers Biomarker identifiers.
#' @return `psa_params` with the new columns; attributes `coefficient_map`
#'   and `survival_groups` (survival_evppi_groups()).
build_survival_evppi_inputs <- function(psa_params, sampling_models,
                                        biomarkers = get_biomarkers()) {
  if (is.null(psa_params$model_idx)) stop("psa_params has no model_idx column")
  component <- get_joint_sampling_models(sampling_models)
  n_cached <- nrow(survival_draw_matrix(component, SURVIVAL_EVPPI_OUTCOMES[1]))
  df_coefs <- add_cell_effects(
    extract_survival_coefficients(component, n_cached, psa_params$model_idx, biomarkers),
    biomarkers)
  df_map <- attr(df_coefs, "coefficient_map")
  v_clash <- intersect(names(df_coefs), names(psa_params))
  if (length(v_clash)) stop("psa_params already has columns ", paste(v_clash, collapse = ", "))
  for (j in which(df_map$role == "interaction")) {
    legacy_col <- paste0("b_", df_map$biomarker[j], "_rx_", df_map$outcome[j])
    if (legacy_col %in% names(psa_params) &&
        !isTRUE(all(psa_params[[legacy_col]] == df_coefs[[df_map$column[j]]]))) {
      stop("PSA column ", legacy_col, " differs from the extracted ", df_map$column[j],
           ": coefficients are not aligned by model_idx")
    }
  }
  df_out <- cbind(psa_params, df_coefs)
  attr(df_out, "coefficient_map") <- df_map
  attr(df_out, "survival_groups") <- survival_evppi_groups(df_map)
  df_out
}

#' Label of the regression form used for an EVPPI estimate
#'
#' @param method voi::evppi() method.
#' @param n_params Number of parameters in the regression.
#' @return "gam_tensor" (tensor-product spline, up to four parameters),
#'   "gam_additive" (additive splines, more than four) or `method`.
evppi_method_label <- function(method, n_params) {
  if (!identical(method, "gam")) return(method)
  if (n_params <= 4) "gam_tensor" else "gam_additive"
}

#' Extended EVPPI of one scenario: existing groups plus the survival groups
#'
#' @param psa_obj Scenario PSA object.
#' @param psa_params Scenario PSA parameters (with `model_idx`).
#' @param sampling_models Sampling cache.
#' @param wtp Willingness-to-pay threshold.
#' @param base_params Existing single EVPPI parameters (current_evppi_config()$params).
#' @param base_groups Existing groups (current_evppi_config()$groups).
#' @param seed Seed for the standard-error draws.
#' @param B Number of regression-coefficient draws for the standard errors.
#' @param cross_check Named list mapping group names to further voi methods
#'   ("earth", "bart") run as cross-checks of the GAM estimate.
#' @param cross_check_args Named list of extra voi::evppi() arguments per
#'   cross-check method (e.g. list(earth = list(degree = 2))).
#' @param scenario_id,c_drug_nivo Labels copied into the output.
#' @return Data frame with `scenario_id`, `wtp`, `c_drug_nivo`, `group`,
#'   `n_params`, `method`, `evppi`, `evppi_se`, `evpi`, `error`; one row per
#'   group and method. EVPI exactly zero gives every group 0 with method
#'   "zero_evpi_bound" (#168); a small positive EVPI the estimator skips gives
#'   NA with method "not_estimated". Attribute `group_consistency` holds the
#'   check of run_evppi_analysis(), which stops on a violation.
run_extended_evppi <- function(psa_obj, psa_params, sampling_models, wtp,
                               base_params, base_groups, seed, B = 1000,
                               cross_check = list(), cross_check_args = list(),
                               scenario_id = NA_character_, c_drug_nivo = NA_real_) {
  df_inputs <- build_survival_evppi_inputs(psa_params, sampling_models)
  survival_groups <- attr(df_inputs, "survival_groups")
  v_survival_params <- unique(unlist(survival_groups))
  groups <- c(base_groups, survival_groups)
  v_params <- unique(c(base_params, v_survival_params))

  outcomes <- psa_outcome_matrices(psa_obj)
  evpi <- calculate_evpi_from_nmb(outcomes$effect * wtp - outcomes$cost)
  make_rows <- function(group, n_params, method, evppi, evppi_se, error = "") {
    data.frame(scenario_id = scenario_id, wtp = wtp, c_drug_nivo = c_drug_nivo,
               group = group, n_params = as.integer(n_params), method = method,
               evppi = evppi, evppi_se = evppi_se, evpi = evpi, error = error,
               stringsAsFactors = FALSE)
  }

  if (!is.finite(evpi) || evpi <= 0.01) {
    zero <- isTRUE(evpi == 0)
    df_out <- do.call(rbind, lapply(names(groups), function(g) make_rows(
      g, length(groups[[g]]), if (zero) "zero_evpi_bound" else "not_estimated",
      if (zero) 0 else NA_real_, if (zero) 0 else NA_real_,
      if (zero) "" else "EVPI below the estimator's threshold; not estimated")))
    attr(df_out, "group_consistency") <- NULL
    return(df_out)
  }

  df_results <- run_evppi_analysis(psa_obj, df_inputs, wtp = wtp, evppi_params = v_params,
                                   param_groups = groups, seed = seed, B = B)
  df_consistency <- attr(df_results, "group_consistency")
  df_group_rows <- df_results[grepl("^\\[GROUP\\] ", df_results$parameter), , drop = FALSE]
  df_out <- do.call(rbind, lapply(names(groups), function(g) {
    df_row <- df_group_rows[df_group_rows$parameter == paste0("[GROUP] ", g), , drop = FALSE]
    if (nrow(df_row) == 0) {
      return(make_rows(g, length(groups[[g]]), "not_estimated", NA_real_, NA_real_,
                       "Group not estimated (fewer than two members with variation)"))
    }
    make_rows(g, df_row$n_params, evppi_method_label(df_row$method, df_row$n_params),
              df_row$evppi, df_row$evppi_se, df_row$error)
  }))

  for (g in names(cross_check)) {
    for (method in cross_check[[g]]) {
      result <- do.call(calculate_evppi_regression, c(list(
        psa_obj = psa_obj, psa_params = df_inputs, param_names = groups[[g]],
        wtp = wtp, method = method, se = TRUE, B = B, seed = seed),
        cross_check_args[[method]]))
      df_out <- rbind(df_out, make_rows(
        g, length(groups[[g]]), method, result$evppi, result$evppi_se,
        if (is.null(result$error)) "" else result$error))
    }
  }
  attr(df_out, "group_consistency") <- df_consistency
  df_out
}

#' Is the union of the survival coefficients at least each nested group?
#'
#' Every survival group (cell effects, rx_main, prognostic, baseline_survival,
#' the interaction groups) is a function of `survival_all`, so its EVPPI cannot
#' exceed the union's beyond Monte Carlo error. The additive-spline regression
#' of the union can understate a joint EVPPI; this table shows where it does.
#' No difference is reported.
#'
#' @param df_evppi Output of run_extended_evppi() for one scenario.
#' @param v_nested Groups that are functions of `survival_all`.
#' @param method Method of the union row to compare (default "gam_additive").
#' @param tol_se Number of combined standard errors tolerated.
#' @return Data frame with `group`, `group_evppi`, `group_se`, `union_evppi`,
#'   `union_se`, `union_method`, `union_at_least_group` (logical).
check_survival_union <- function(df_evppi, v_nested, method = "gam_additive", tol_se = 2) {
  df_union <- df_evppi[df_evppi$group == "survival_all" & df_evppi$method == method, , drop = FALSE]
  if (nrow(df_union) != 1) stop("No unique survival_all row for method ", method)
  df_gam <- df_evppi[df_evppi$group %in% v_nested & grepl("^gam", df_evppi$method), , drop = FALSE]
  data.frame(
    group = df_gam$group, group_evppi = df_gam$evppi, group_se = df_gam$evppi_se,
    union_evppi = df_union$evppi, union_se = df_union$evppi_se, union_method = method,
    union_at_least_group = df_union$evppi + tol_se * sqrt(df_union$evppi_se^2 + df_gam$evppi_se^2) >=
      df_gam$evppi,
    stringsAsFactors = FALSE)
}
