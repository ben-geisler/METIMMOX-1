# ===============================================================================
# VALUE OF INFORMATION ON THE NIVOLUMAB PRICE AXIS (issue #199, paper 3 step 1a)
# ===============================================================================
# Script 12 runs one PSA per nivolumab price (list price and biosimilar price),
# both with `analysis_seed`. The nivolumab price is a fixed input (#154) that
# enters only the drug cost of progression-free patients on the nivolumab
# schedule, so in every draw i and for every strategy s
#
#   C_is(p) = C_is(p0) + (p - p0) * N_is,
#
# where N_is is the discounted number of nivolumab administrations weighted by
# progression-free occupancy, and the QALYs, model indices and every sampled
# non-price parameter are identical at both prices. check_price_pairing()
# verifies this per draw and derives N_is from the two evaluated prices; the
# remaining functions evaluate EVPI on a price grid from the paired draws.
#
# EVPI is mean_i max_s NMB_is - max_s mean_i NMB_is (calculate_evpi_from_nmb()
# in R/evppi_functions.R). The QALY-scale EVPI, the same quantity computed on
# QALYs, is the limit of EVPI / WTP as the threshold grows and does not depend
# on the price.
#
# Requires R/evppi_functions.R (calculate_evpi_from_nmb(),
# calculate_population_evppi(), psa_outcome_matrices()) and, for the defaults,
# R/model_configs.R (get_control_strategy()).
# ===============================================================================

#' Directory for the value-of-information outputs of paper 3
#'
#' @return Path of `data/output/voi` under the repository root (not created).
voi_output_dir <- function() here::here("data", "output", "voi")

#' One PSA pair per nivolumab price from the scenario cache
#'
#' Scenarios that share a price reuse one PSA in script 12; the function keeps
#' that PSA once per price and stops if two scenarios at the same price carry
#' different PSA outcomes or parameters.
#'
#' @param scenario_cache Scenario cache as returned by load_scenario_cache().
#' @return List with one element per price, ordered from the highest price
#'   down, each `list(price, scenario_ids, psa_obj, psa_params)`.
scenario_price_pairs <- function(scenario_cache) {
  df_scenarios <- scenario_cache$scenarios
  if (is.null(df_scenarios) || !all(c("scenario_id", "c_drug_nivo") %in% names(df_scenarios))) {
    stop("scenario_cache$scenarios must have scenario_id and c_drug_nivo columns")
  }
  v_prices <- sort(unique(df_scenarios$c_drug_nivo), decreasing = TRUE)
  lapply(v_prices, function(price) {
    v_ids <- df_scenarios$scenario_id[df_scenarios$c_drug_nivo == price]
    first <- scenario_cache$all_scenario_results[[v_ids[1]]]
    if (is.null(first)) stop("Scenario cache has no result for scenario ", v_ids[1])
    for (id in v_ids[-1]) {
      other <- scenario_cache$all_scenario_results[[id]]
      if (is.null(other) ||
          !identical(other$psa_obj$cost, first$psa_obj$cost) ||
          !identical(psa_outcome_matrices(other$psa_obj)$effect,
                     psa_outcome_matrices(first$psa_obj)$effect) ||
          !identical(other$psa_params, first$psa_params)) {
        stop("Scenarios ", v_ids[1], " and ", id, " share the nivolumab price ",
             price, " but carry different PSA results")
      }
    }
    list(price = price, scenario_ids = v_ids,
         psa_obj = first$psa_obj, psa_params = first$psa_params)
  })
}

#' Largest absolute difference between two numeric matrices or data frames
#'
#' @param x,y Numeric matrices or data frames of the same shape.
#' @return Maximum absolute elementwise difference (0 for identical inputs,
#'   Inf when the shapes differ or a missing value is unmatched).
max_abs_difference <- function(x, y) {
  m_x <- as.matrix(x)
  m_y <- as.matrix(y)
  if (!identical(dim(m_x), dim(m_y))) return(Inf)
  if (!identical(is.na(m_x), is.na(m_y))) return(Inf)
  if (length(m_x) == 0 || all(is.na(m_x))) return(0)
  max(abs(m_x - m_y), na.rm = TRUE)
}

#' Verify the per-draw price pairing of two PSA pairs and derive N_i
#'
#' Checks, for two PSA pairs that differ only in the nivolumab price, that the
#' strategies, draw identifiers, `model_idx`, QALYs, life-years (when stored
#' as `ly_disc`/`ly_undisc`, issue #197) and every sampled non-price parameter
#' are identical per
#' draw, derives the per-draw administrations
#' `N_is = (C_is(p1) - C_is(p0)) / (p1 - p0)`, and checks that they are finite,
#' non-negative, zero for the control strategy, at most `max_administrations`
#' and, when `df_threshold` is supplied, that their probabilistic means
#' (incremental versus control) equal its `discounted_administrations` column.
#'
#' @param pair_ref Reference pair `list(psa_obj, psa_params)` at price `p0`.
#' @param pair_alt Second pair `list(psa_obj, psa_params)` at price `p1`.
#' @param price_ref Reference nivolumab price p0 (the list price).
#' @param price_alt Second nivolumab price p1 (must differ from p0).
#' @param control_strategy Strategy that receives no nivolumab (its N_i must be
#'   zero); NULL skips that check.
#' @param max_administrations Upper bound on N_i (the number of scheduled
#'   administrations, `sum(l_params_base$l_nivo)`); NULL skips the check.
#' @param df_threshold Optional data frame read from
#'   outputs/tables/cea_threshold_price.csv; its probabilistic rows are
#'   compared with the mean incremental administrations.
#' @param price_params Names of PSA parameter columns that are allowed to
#'   differ between the prices (default: `c_drug_nivo`).
#' @param tol Absolute tolerance for identity and for the bounds.
#' @param tol_threshold Relative tolerance for the comparison with
#'   `df_threshold` (the CSV stores 15 significant digits).
#' @return Object of class `voi_price_pairing`: a list with `price_ref`,
#'   `price_alt`, `strategies`, `n_sim`, `sim`, `model_idx`, `cost_ref`
#'   (matrix of costs at p0), `effect` (QALY matrix), `life_years` (discounted
#'   life-year matrix or
#'   NULL), `administrations` (matrix of N_i), `checks` (data frame with
#'   `check`, `passed`, `detail`) and `administrations_vs_table` (data frame
#'   or NULL). Stops listing every failed check.
check_price_pairing <- function(pair_ref, pair_alt, price_ref, price_alt,
                                control_strategy = NULL,
                                max_administrations = NULL,
                                df_threshold = NULL,
                                price_params = "c_drug_nivo",
                                tol = 1e-9, tol_threshold = 1e-9) {
  if (!is.numeric(price_ref) || !is.numeric(price_alt) || length(price_ref) != 1 ||
      length(price_alt) != 1 || !is.finite(price_ref) || !is.finite(price_alt) ||
      price_ref == price_alt) {
    stop("price_ref and price_alt must be two different finite numbers")
  }
  m_out_ref <- psa_outcome_matrices(pair_ref$psa_obj)
  m_out_alt <- psa_outcome_matrices(pair_alt$psa_obj)
  df_params_ref <- as.data.frame(pair_ref$psa_params)
  df_params_alt <- as.data.frame(pair_alt$psa_params)

  df_checks <- data.frame(check = character(), passed = logical(), detail = character(),
                          stringsAsFactors = FALSE)
  add_check <- function(check, passed, detail) {
    df_checks[nrow(df_checks) + 1, ] <<- list(check, isTRUE(passed), detail)
  }

  v_strategies <- colnames(m_out_ref$cost)
  add_check("Strategies identical",
            identical(v_strategies, colnames(m_out_alt$cost)) &&
              identical(colnames(m_out_ref$effect), colnames(m_out_alt$effect)) &&
              identical(v_strategies, colnames(m_out_ref$effect)),
            paste(v_strategies, collapse = ", "))

  n_sim <- nrow(m_out_ref$cost)
  same_rows <- nrow(m_out_alt$cost) == n_sim && nrow(m_out_ref$effect) == n_sim &&
    nrow(m_out_alt$effect) == n_sim && nrow(df_params_ref) == n_sim &&
    nrow(df_params_alt) == n_sim
  add_check("Number of draws identical", same_rows,
            sprintf("%d and %d draws", n_sim, nrow(m_out_alt$cost)))

  if (!all(df_checks$passed)) {
    stop("Price pairing failed:\n", paste0("  - ", df_checks$check[!df_checks$passed], ": ",
                                           df_checks$detail[!df_checks$passed], collapse = "\n"))
  }

  first_mismatch <- function(x, y) {
    v_bad <- which(!(x == y) | xor(is.na(x), is.na(y)))
    if (length(v_bad) == 0) "none" else sprintf("%d draws differ (first: draw %d)",
                                                length(v_bad), v_bad[1])
  }
  add_check("Draw identifiers (sim) identical",
            identical(as.numeric(df_params_ref$sim), as.numeric(df_params_alt$sim)) &&
              !is.null(df_params_ref$sim),
            if (is.null(df_params_ref$sim) || is.null(df_params_alt$sim)) "sim column missing"
            else first_mismatch(df_params_ref$sim, df_params_alt$sim))
  v_model_ref <- pair_ref$psa_obj$model_idx
  v_model_alt <- pair_alt$psa_obj$model_idx
  model_ok <- !is.null(v_model_ref) && !is.null(v_model_alt) &&
    identical(as.integer(v_model_ref), as.integer(v_model_alt)) &&
    identical(as.integer(df_params_ref$model_idx), as.integer(v_model_ref)) &&
    identical(as.integer(df_params_alt$model_idx), as.integer(v_model_alt))
  add_check("model_idx identical", model_ok,
            if (is.null(v_model_ref) || is.null(v_model_alt)) "model_idx missing"
            else first_mismatch(v_model_ref, v_model_alt))

  qaly_diff <- max_abs_difference(m_out_ref$effect, m_out_alt$effect)
  add_check("QALYs identical per draw", qaly_diff <= tol,
            sprintf("max |difference| %.3g", qaly_diff))

  m_life_years <- NULL
  for (key in c("ly_disc", "ly_undisc")) {
    m_ly_ref <- pair_ref$psa_obj[[key]]
    m_ly_alt <- pair_alt$psa_obj[[key]]
    label <- paste0("Life-years identical per draw (", key, ")")
    if (is.null(m_ly_ref) && is.null(m_ly_alt)) {
      add_check(label, TRUE, "not stored in these PSA objects")
    } else {
      ly_diff <- if (is.null(m_ly_ref) || is.null(m_ly_alt)) Inf
                 else max_abs_difference(m_ly_ref, m_ly_alt)
      add_check(label, ly_diff <= tol, sprintf("max |difference| %.3g", ly_diff))
      if (key == "ly_disc") m_life_years <- as.matrix(m_ly_ref)
    }
  }

  v_param_cols <- setdiff(names(df_params_ref), c("sim", "model_idx", price_params))
  v_param_cols_alt <- setdiff(names(df_params_alt), c("sim", "model_idx", price_params))
  if (!setequal(v_param_cols, v_param_cols_alt)) {
    add_check("Non-price parameters identical per draw", FALSE,
              paste("column sets differ:",
                    paste(union(setdiff(v_param_cols, v_param_cols_alt),
                                setdiff(v_param_cols_alt, v_param_cols)), collapse = ", ")))
  } else {
    v_param_diff <- vapply(v_param_cols, function(col) {
      max_abs_difference(df_params_ref[[col]], df_params_alt[[col]])
    }, numeric(1))
    v_bad_cols <- v_param_cols[v_param_diff > tol]
    add_check("Non-price parameters identical per draw", length(v_bad_cols) == 0,
              if (length(v_bad_cols) == 0) {
                sprintf("%d columns, max |difference| %.3g", length(v_param_cols),
                        if (length(v_param_diff)) max(v_param_diff) else 0)
              } else {
                paste("differ:", paste(sprintf("%s (%.3g)", v_bad_cols,
                                               v_param_diff[v_bad_cols]), collapse = ", "))
              })
  }

  m_cost_ref <- m_out_ref$cost
  m_cost_alt <- m_out_alt$cost[, v_strategies, drop = FALSE]
  m_admin <- (m_cost_alt - m_cost_ref) / (price_alt - price_ref)

  n_nonfinite <- sum(!is.finite(m_admin))
  add_check("Administrations N_i finite", n_nonfinite == 0,
            sprintf("%d non-finite values", n_nonfinite))
  min_admin <- suppressWarnings(min(m_admin, na.rm = TRUE))
  add_check("Administrations N_i non-negative", is.finite(min_admin) && min_admin >= -tol,
            sprintf("minimum %.6g", min_admin))
  if (!is.null(control_strategy)) {
    if (!control_strategy %in% v_strategies) {
      add_check("Control strategy receives no nivolumab", FALSE,
                paste("control strategy", control_strategy, "not in the PSA"))
    } else {
      max_control <- max(abs(m_admin[, control_strategy]))
      add_check("Control strategy receives no nivolumab", max_control <= tol,
                sprintf("max |N_i| %.3g for %s", max_control, control_strategy))
    }
  }
  if (!is.null(max_administrations)) {
    max_admin <- suppressWarnings(max(m_admin, na.rm = TRUE))
    add_check("Administrations N_i at most the scheduled number",
              is.finite(max_admin) && max_admin <= max_administrations + tol,
              sprintf("maximum %.6g of %g scheduled", max_admin, max_administrations))
  }
  # Reconstruction at the second price from the reference costs and N_i.
  recon_diff <- max_abs_difference(m_cost_ref + (price_alt - price_ref) * m_admin, m_cost_alt)
  add_check("C_i(p1) = C_i(p0) + (p1 - p0) N_i reproduced", recon_diff <= tol * max(1, abs(price_alt)),
            sprintf("max |difference| %.3g", recon_diff))

  df_admin_table <- NULL
  if (!is.null(df_threshold)) {
    df_prob <- df_threshold[df_threshold$analysis == "probabilistic", , drop = FALSE]
    if (nrow(df_prob) == 0 || is.null(control_strategy) ||
        !all(df_prob$strategy %in% v_strategies)) {
      add_check("Mean administrations match cea_threshold_price.csv", FALSE,
                "no probabilistic rows, unknown strategies or no control strategy")
    } else {
      v_mean_inc <- vapply(df_prob$strategy, function(s) {
        mean(m_admin[, s] - m_admin[, control_strategy])
      }, numeric(1))
      df_admin_table <- data.frame(
        strategy = df_prob$strategy,
        mean_incremental_administrations = unname(v_mean_inc),
        table_discounted_administrations = df_prob$discounted_administrations,
        stringsAsFactors = FALSE)
      df_admin_table$relative_difference <- abs(v_mean_inc - df_prob$discounted_administrations) /
        pmax(1, abs(df_prob$discounted_administrations))
      add_check("Mean administrations match cea_threshold_price.csv",
                all(df_admin_table$relative_difference <= tol_threshold),
                paste(sprintf("%s %.10g vs %.10g", df_admin_table$strategy,
                              df_admin_table$mean_incremental_administrations,
                              df_admin_table$table_discounted_administrations), collapse = "; "))
    }
  }

  if (!all(df_checks$passed)) {
    stop("Price pairing failed:\n", paste0("  - ", df_checks$check[!df_checks$passed], ": ",
                                           df_checks$detail[!df_checks$passed], collapse = "\n"))
  }

  structure(list(
    price_ref = price_ref, price_alt = price_alt, strategies = v_strategies,
    n_sim = n_sim, sim = df_params_ref$sim, model_idx = v_model_ref,
    cost_ref = m_cost_ref, effect = m_out_ref$effect, life_years = m_life_years,
    administrations = m_admin, checks = df_checks,
    administrations_vs_table = df_admin_table
  ), class = "voi_price_pairing")
}

#' Price pairing of the scenario cache (list price versus biosimilar price)
#'
#' @param scenario_cache Scenario cache as returned by load_scenario_cache().
#' @param ... Further arguments passed to check_price_pairing().
#' @return `voi_price_pairing` object with the highest cached price as the
#'   reference and the lowest as the second price.
scenario_price_pairing <- function(scenario_cache, ...) {
  pairs <- scenario_price_pairs(scenario_cache)
  if (length(pairs) != 2L) {
    stop("Expected exactly two nivolumab prices in the scenario cache, found ", length(pairs))
  }
  check_price_pairing(pairs[[1]], pairs[[2]], price_ref = pairs[[1]]$price,
                      price_alt = pairs[[2]]$price, ...)
}

#' Per-draw costs at any nivolumab price
#'
#' @param pairing `voi_price_pairing` object from check_price_pairing().
#' @param price Nivolumab price (finite, non-negative).
#' @return Matrix of costs, one row per draw and one column per strategy:
#'   `C_i(p) = C_i(p0) + (p - p0) N_i`.
costs_at_price <- function(pairing, price) {
  if (!inherits(pairing, "voi_price_pairing")) stop("pairing must come from check_price_pairing()")
  if (!is.numeric(price) || length(price) != 1 || !is.finite(price) || price < 0) {
    stop("price must be one finite non-negative number")
  }
  pairing$cost_ref + (price - pairing$price_ref) * pairing$administrations
}

#' Per-patient EVPI and decision probabilities at one price and threshold
#'
#' @param pairing `voi_price_pairing` object.
#' @param price Nivolumab price.
#' @param wtp Willingness-to-pay threshold (EUR per QALY).
#' @return List with `evpi`, `prob_optimal` (named share of draws in which each
#'   strategy has the highest NMB; ties go to the first strategy),
#'   `expected_nmb` (named) and `optimal_strategy` (highest expected NMB).
evpi_at_price <- function(pairing, price, wtp) {
  if (!is.numeric(wtp) || length(wtp) != 1 || !is.finite(wtp) || wtp < 0) {
    stop("wtp must be one finite non-negative number")
  }
  m_nmb <- pairing$effect * wtp - costs_at_price(pairing, price)
  v_best <- max.col(m_nmb, ties.method = "first")
  v_prob <- tabulate(v_best, nbins = ncol(m_nmb)) / nrow(m_nmb)
  v_expected <- colMeans(m_nmb)
  list(
    evpi = calculate_evpi_from_nmb(m_nmb),
    prob_optimal = stats::setNames(v_prob, pairing$strategies),
    expected_nmb = stats::setNames(v_expected, pairing$strategies),
    optimal_strategy = pairing$strategies[which.max(v_expected)]
  )
}

#' QALY-scale EVPI (the infinite-threshold limit, independent of the price)
#'
#' @param pairing `voi_price_pairing` object (or any list with a QALY matrix
#'   `effect`).
#' @return mean_i max_s QALY_is - max_s mean_i QALY_is.
qaly_scale_evpi <- function(pairing) {
  calculate_evpi_from_nmb(as.matrix(pairing$effect))
}

#' Probabilistic adoption threshold price from cea_threshold_price.csv
#'
#' The adoption threshold of paper 2 is the nivolumab price at which a guided
#' strategy's probabilistic incremental NMB versus standard of care is zero at
#' the base-case threshold. When several strategies have a finite threshold
#' price, the highest is returned (the first adoption as the price falls).
#'
#' @param df_threshold Data frame read from outputs/tables/cea_threshold_price.csv.
#' @return List with `price` and `strategy`.
adoption_threshold_price <- function(df_threshold) {
  df_prob <- df_threshold[df_threshold$analysis == "probabilistic" &
                            is.finite(df_threshold$threshold_price), , drop = FALSE]
  if (nrow(df_prob) == 0) stop("cea_threshold_price.csv has no finite probabilistic threshold price")
  idx <- which.max(df_prob$threshold_price)
  list(price = df_prob$threshold_price[idx], strategy = df_prob$strategy[idx])
}

#' Step of the nivolumab price axis (EUR per administration)
#'
#' Used by steps 1a and 1b in script 14 and value_of_research.qmd, which read one
#' saved step 1a result. EUR 50 since issue #210 (base price EUR 1,392; EUR 250
#' with the former EUR 13,923).
VOI_PRICE_STEP <- 50

#' Nivolumab price grid for the price axis
#'
#' @param list_price Upper end of the grid (the list price).
#' @param marked_prices Further prices that must be grid points (biosimilar
#'   price, adoption threshold price).
#' @param step Grid step from zero.
#' @return Sorted unique numeric vector from 0 to `list_price`.
price_axis_grid <- function(list_price, marked_prices = numeric(0), step = VOI_PRICE_STEP) {
  v_marked <- as.numeric(marked_prices)
  if (!is.finite(list_price) || list_price <= 0 || any(!is.finite(v_marked)) ||
      any(v_marked < 0 | v_marked > list_price)) {
    stop("Grid prices must be finite and within 0..list_price")
  }
  sort(unique(c(seq(0, list_price, by = step), list_price, v_marked)))
}

#' Per-patient and population EVPI on a price grid
#'
#' @param pairing `voi_price_pairing` object.
#' @param v_prices Nivolumab prices (price_axis_grid()).
#' @param v_wtp Willingness-to-pay thresholds.
#' @param population List with `annual_incidence`, `research_horizon` and
#'   `discount_rate` (current_population_inputs() reads the globals).
#' @param marked_prices Named numeric vector of prices to label in the
#'   `price_marker` column (e.g. c(list = 1392, biosimilar = 464)).
#' @return Data frame with one row per price and threshold: `c_drug_nivo`,
#'   `wtp`, `price_marker`, `evpi_per_patient`, `evpi_population_millions`,
#'   `evpi_qaly_equivalent` (EVPI / WTP), `evpi_qaly_scale` and
#'   `evpi_qaly_scale_population` (QALYs; price- and threshold-independent),
#'   `optimal_strategy`, `prob_optimal_<strategy>` per strategy, the
#'   population inputs and `n_sim`.
run_price_axis_evpi <- function(pairing, v_prices, v_wtp, population,
                                marked_prices = numeric(0)) {
  for (field in c("annual_incidence", "research_horizon", "discount_rate")) {
    if (is.null(population[[field]])) stop("population$", field, " is required")
  }
  scale_population <- function(x) calculate_population_evppi(
    x, annual_incidence = population$annual_incidence,
    research_horizon = population$research_horizon,
    discount_rate = population$discount_rate)
  evpi_qaly <- qaly_scale_evpi(pairing)
  label_price <- function(price) {
    v_hits <- names(marked_prices)[abs(marked_prices - price) < 1e-9]
    if (length(v_hits) == 0) "" else paste(v_hits, collapse = "; ")
  }
  rows <- list()
  for (wtp in v_wtp) {
    for (price in v_prices) {
      result <- evpi_at_price(pairing, price, wtp)
      df_row <- data.frame(
        c_drug_nivo = price, wtp = wtp, price_marker = label_price(price),
        evpi_per_patient = result$evpi,
        evpi_population_millions = scale_population(result$evpi),
        evpi_qaly_equivalent = if (wtp > 0) result$evpi / wtp else NA_real_,
        evpi_qaly_scale = evpi_qaly,
        evpi_qaly_scale_population = scale_population(evpi_qaly) * 1e6,
        optimal_strategy = result$optimal_strategy,
        stringsAsFactors = FALSE)
      for (s in pairing$strategies) {
        df_row[[paste0("prob_optimal_", s)]] <- result$prob_optimal[[s]]
      }
      df_row$annual_incidence <- population$annual_incidence
      df_row$research_horizon_years <- population$research_horizon
      df_row$discount_rate_research <- population$discount_rate
      df_row$n_sim <- pairing$n_sim
      rows[[length(rows) + 1]] <- df_row
    }
  }
  df_out <- do.call(rbind, rows)
  rownames(df_out) <- NULL
  df_out
}

# ===============================================================================
# STEP 1A DRIVER: pairing, price axis, reproduction checks, extended EVPPI
# ===============================================================================
# Needs R/voi_survival_parameters.R (run_extended_evppi()) and R/cache_paths.R
# (cache_fingerprint()). The result is saved under data/output/voi with a
# fingerprint of every input, so the PDF and GFM passes of the technical report
# compute it once; the pipeline caches are only read.

#' EVPI of every scenario reproduced from the price axis
#'
#' @param df_axis Output of run_price_axis_evpi().
#' @param scenario_cache Scenario cache.
#' @param df_table_s7 Optional data frame read from outputs/tables/table_s7.csv
#'   (rows with `group == "evpi"`).
#' @param tol Absolute tolerance (EUR per patient).
#' @return Data frame with `scenario_id`, `wtp`, `c_drug_nivo`, `evpi_axis`,
#'   `evpi_scenario_psa`, `evpi_table_s7` and `reproduced`.
reproduce_scenario_evpi <- function(df_axis, scenario_cache, df_table_s7 = NULL, tol = 1e-6) {
  df_rows <- lapply(seq_len(nrow(scenario_cache$scenarios)), function(k) {
    info <- scenario_cache$scenarios[k, ]
    result <- scenario_cache$all_scenario_results[[info$scenario_id]]
    outcomes <- psa_outcome_matrices(result$psa_obj)
    evpi_psa <- calculate_evpi_from_nmb(outcomes$effect * info$wtp - outcomes$cost)
    v_hit <- abs(df_axis$c_drug_nivo - info$c_drug_nivo) < 1e-9 & df_axis$wtp == info$wtp
    evpi_axis <- if (sum(v_hit) == 1) df_axis$evpi_per_patient[v_hit] else NA_real_
    evpi_s7 <- NA_real_
    if (!is.null(df_table_s7)) {
      v_s7 <- df_table_s7$scenario_id == info$scenario_id & df_table_s7$group == "evpi"
      if (sum(v_s7) == 1) evpi_s7 <- df_table_s7$evpi[v_s7]
    }
    reproduced <- is.finite(evpi_axis) && abs(evpi_axis - evpi_psa) <= tol &&
      (is.null(df_table_s7) || (is.finite(evpi_s7) && abs(evpi_axis - evpi_s7) <= tol))
    data.frame(scenario_id = info$scenario_id, wtp = info$wtp, c_drug_nivo = info$c_drug_nivo,
               evpi_axis = evpi_axis, evpi_scenario_psa = evpi_psa, evpi_table_s7 = evpi_s7,
               reproduced = reproduced, stringsAsFactors = FALSE)
  })
  do.call(rbind, df_rows)
}

#' Identity of the step 1a inputs (for the saved results)
#'
#' @param scenario_cache,sampling_models,df_threshold,df_table_s7 Inputs of
#'   run_voi_step1a().
#' @param base_config,settings,population Inputs of run_voi_step1a().
#' @return List with `fingerprint` and `inputs`: `upstream` holds the
#'   sampling-cache fingerprint and the scenario cache's fingerprint, result
#'   hash and PSA pair generations (voi_upstream_fingerprints()), and
#'   `implementation` the `voi` and `evppi` categories of
#'   calculation_identity() (#205).
voi_step1a_identity <- function(scenario_cache, sampling_models, df_threshold, df_table_s7,
                                base_config, settings, population) {
  inputs <- list(
    schema = "voi_step1a_v2",
    upstream = voi_upstream_fingerprints(sampling_models, scenario_cache = scenario_cache),
    threshold_table = df_threshold, table_s7 = df_table_s7,
    base_config = base_config, settings = settings, population = population,
    implementation = list(voi = calculation_identity("voi"), evppi = calculation_identity("evppi")))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Run paper 3 step 1a on the scenario cache
#'
#' @param scenario_cache Validated scenario cache (load_scenario_cache()).
#' @param sampling_models Sampling cache (script 06).
#' @param df_threshold Data frame of outputs/tables/cea_threshold_price.csv.
#' @param df_table_s7 Data frame of outputs/tables/table_s7.csv, or NULL.
#' @param base_config Existing EVPPI configuration (current_evppi_config()).
#' @param population Population inputs (current_population_inputs()).
#' @param settings List with `wtp` (thresholds of the price axis), `step`
#'   (grid step), `seed`, `B`, `control_strategy`, `max_administrations`,
#'   `cross_check` and `cross_check_args` (see run_extended_evppi()).
#' @param compare_published Compare N_i with cea_threshold_price.csv and EVPI
#'   with table_s7.csv (TRUE for reported results).
#' @return List with `pairing_checks`, `administrations_vs_table`,
#'   `administrations_summary`, `adoption_threshold`, `marked_prices`,
#'   `price_axis`, `evpi_reproduction`, `evppi`, `union_checks`,
#'   `group_consistency`, `population` and `prices`. Stops on a failed
#'   pairing check, a failed reproduction or a group-consistency violation.
run_voi_step1a <- function(scenario_cache, sampling_models, df_threshold, df_table_s7,
                           base_config, population, settings, compare_published = TRUE) {
  pairing <- scenario_price_pairing(
    scenario_cache, control_strategy = settings$control_strategy,
    max_administrations = settings$max_administrations,
    df_threshold = if (compare_published) df_threshold else NULL)
  threshold <- adoption_threshold_price(df_threshold)
  v_marked <- c(list_price = pairing$price_ref, biosimilar = pairing$price_alt,
                adoption_threshold = threshold$price)
  v_grid <- price_axis_grid(pairing$price_ref, v_marked, step = settings$step)
  df_axis <- run_price_axis_evpi(pairing, v_grid, settings$wtp, population, v_marked)

  df_repro <- reproduce_scenario_evpi(df_axis, scenario_cache,
                                      if (compare_published) df_table_s7 else NULL)
  if (!all(df_repro$reproduced)) {
    df_bad <- df_repro[!df_repro$reproduced, , drop = FALSE]
    stop("Price-axis EVPI does not reproduce the scenario EVPI:\n",
         paste(sprintf("  %s: axis %.8f, scenario PSA %.8f, Table S7 %.8f",
                       df_bad$scenario_id, df_bad$evpi_axis, df_bad$evpi_scenario_psa,
                       df_bad$evpi_table_s7), collapse = "\n"))
  }

  df_map <- attr(build_survival_evppi_inputs(
    scenario_cache$all_scenario_results[[1]]$psa_params, sampling_models), "coefficient_map")
  v_nested <- c(grep("^cell_effect_", names(survival_evppi_groups(df_map)), value = TRUE),
                "rx_main", "prognostic", "baseline_survival",
                grep("^interaction_", names(base_config$groups), value = TRUE))
  evppi_rows <- list()
  union_rows <- list()
  consistency <- list()
  for (k in seq_len(nrow(scenario_cache$scenarios))) {
    info <- scenario_cache$scenarios[k, ]
    result <- scenario_cache$all_scenario_results[[info$scenario_id]]
    df_evppi <- run_extended_evppi(
      result$psa_obj, result$psa_params, sampling_models, wtp = info$wtp,
      base_params = base_config$params, base_groups = base_config$groups,
      seed = settings$seed, B = settings$B, cross_check = settings$cross_check,
      cross_check_args = settings$cross_check_args,
      scenario_id = info$scenario_id, c_drug_nivo = info$c_drug_nivo)
    consistency[[info$scenario_id]] <- attr(df_evppi, "group_consistency")
    if (any(df_evppi$method == "gam_additive" & df_evppi$group == "survival_all")) {
      union_rows[[info$scenario_id]] <- cbind(scenario_id = info$scenario_id,
                                              check_survival_union(df_evppi, v_nested))
    }
    evppi_rows[[info$scenario_id]] <- df_evppi
  }
  df_evppi_all <- do.call(rbind, evppi_rows)
  rownames(df_evppi_all) <- NULL
  df_union_all <- if (length(union_rows)) do.call(rbind, union_rows) else NULL
  if (!is.null(df_union_all)) rownames(df_union_all) <- NULL

  df_admin_summary <- do.call(rbind, lapply(pairing$strategies, function(s) {
    v_n <- pairing$administrations[, s]
    data.frame(strategy = s, mean = mean(v_n), min = min(v_n),
               p2_5 = unname(stats::quantile(v_n, 0.025, type = 7)),
               p97_5 = unname(stats::quantile(v_n, 0.975, type = 7)), max = max(v_n),
               stringsAsFactors = FALSE)
  }))

  list(pairing_checks = pairing$checks,
       administrations_vs_table = pairing$administrations_vs_table,
       administrations_summary = df_admin_summary,
       adoption_threshold = threshold, marked_prices = v_marked,
       price_axis = df_axis, evpi_reproduction = df_repro,
       evppi = df_evppi_all, union_checks = df_union_all,
       group_consistency = consistency, population = population,
       prices = c(reference = pairing$price_ref, alternative = pairing$price_alt))
}

#' Step 1a results, computed once per input identity
#'
#' Reads `cache_file` when its stored identity equals the current one and
#' otherwise runs run_voi_step1a() (always comparing with the published
#' tables) and saves the result. Writes `voi_price_axis.csv` and
#' `evppi_extended.csv` next to it in either case.
#'
#' @param scenario_cache,sampling_models,df_threshold,df_table_s7 See
#'   run_voi_step1a().
#' @param base_config,population,settings See run_voi_step1a().
#' @param cache_file Path of the saved results (under data/output/voi).
#' @return The list returned by run_voi_step1a(), with `identity` and
#'   `computed` (TRUE when computed in this call).
load_or_run_voi_step1a <- function(scenario_cache, sampling_models, df_threshold, df_table_s7,
                                   base_config, population, settings,
                                   cache_file = file.path(voi_output_dir(),
                                                          "voi_step1a_results.rds")) {
  identity <- voi_step1a_identity(scenario_cache, sampling_models, df_threshold, df_table_s7,
                                  base_config, settings, population)
  results <- NULL
  if (file.exists(cache_file)) {
    saved <- readRDS(cache_file)
    # Names the differing input (e.g. upstream$scenario$fingerprint) before recomputing (#205).
    if (check_voi_identity(saved$identity$fingerprint, saved$identity$inputs, identity,
                           paste0("step 1a results (", basename(cache_file), ")"),
                           "Recomputing step 1a.", required = FALSE)) results <- saved
  }
  if (is.null(results)) {
    results <- run_voi_step1a(scenario_cache, sampling_models, df_threshold, df_table_s7,
                              base_config, population, settings, compare_published = TRUE)
    results$identity <- identity
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(results, cache_file)
    results$computed <- TRUE
  } else {
    results$computed <- FALSE
  }
  write_voi_step1a_csvs(results, dirname(cache_file))
  results
}

#' Write the step 1a CSV outputs
#'
#' @param results List returned by run_voi_step1a().
#' @param directory Output directory (created if needed).
#' @return Paths of `voi_price_axis.csv` and `evppi_extended.csv`, invisibly.
write_voi_step1a_csvs <- function(results, directory = voi_output_dir()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  v_paths <- file.path(directory, c("voi_price_axis.csv", "evppi_extended.csv"))
  utils::write.csv(results$price_axis, v_paths[1], row.names = FALSE)
  utils::write.csv(results$evppi, v_paths[2], row.names = FALSE)
  invisible(v_paths)
}
