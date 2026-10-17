# ===============================================================================
# CELL-LEVEL POLICY SPACE AND ITS PRICE-AXIS EVPI (issue #200, paper 3 step 1b)
# ===============================================================================
# The PSA stores, per draw, the outcomes of the four joint biomarker cells
# (week-4 CRP x TMB/BRAF; keys "00", "01", "10", "11" with the CRP digit first)
# under FLOX and under FLOX + nivolumab, and the draw's joint cell weights
# (psa_obj$cells, built by build_psa_cell_component() in R/psa_functions.R).
# A policy assigns nivolumab to a subset of the cells; with 2 biomarkers there
# are 2^4 = 16 policies. Per draw i and policy k
#
#   Y_ik = sum_c w_ic * Y_ic[a_kc] + test cost_k        (Y = cost)
#   Y_ik = sum_c w_ic * Y_ic[a_kc]                      (QALYs, life-years,
#                                                         administrations)
#
# where a_kc is the arm policy k gives cell c. A policy pays the diagnostic
# test of every biomarker it conditions on (c_test_CRP for CRP, c_test_NGS for
# TMB/BRAF, both when both); treating nobody or everyone needs no test.
#
# The three paper 2 strategies are policies "0000" (standard of care), "0011"
# (CRP-guided) and "0101" (TMB/BRAF-guided). model_fun() applies the ordering
# clamp (PFS capped at OS) and pmax() on the death transitions to each
# strategy subgroup curve, which is a weighted average of cell curves, while
# the cell outcomes clamp each cell's own curves. The two therefore agree
# exactly (to floating point) whenever no cell curve of the draw crosses, and
# differ by the clamp contribution otherwise; strategy_reproduction_gaps()
# reports that gap per draw. Life-years are the area under OS alone, which is
# linear in the cell curves and never clamped, so their gap is zero (to
# floating point) in every draw, crossing or not (issue #197).
#
# Cost is linear in the nivolumab price with slope equal to the discounted
# administrations, so policy_price_pairing() turns the policy space into the
# `voi_price_pairing` object of R/voi_price_axis.R, whose EVPI functions work
# for any number of policies.
#
# Requires R/model_configs.R, R/evppi_functions.R (calculate_evpi_from_nmb(),
# calculate_population_evppi()) and R/voi_price_axis.R.
# ===============================================================================

#' Does a cell assignment condition on a biomarker?
#'
#' A policy conditions on biomarker b when two cells that differ only in b
#' receive different arms; only then is the test of b needed.
#'
#' @param v_assign 0/1 assignment of nivolumab per cell (cell order of
#'   `df_cells`).
#' @param df_cells Cell definition (one column per biomarker, values 0/1).
#' @param biomarker Biomarker column.
#' @return Logical scalar.
policy_conditions_on <- function(v_assign, df_cells, biomarker) {
  v_other <- setdiff(names(df_cells), biomarker)
  v_key <- if (length(v_other)) apply(df_cells[v_other], 1, paste0, collapse = "") else
    rep("", nrow(df_cells))
  any(vapply(split(v_assign, v_key), function(x) length(unique(x)) > 1, logical(1)))
}

#' The policy space: every assignment of nivolumab to the joint cells
#'
#' @param df_cells Cell definition data frame (rows named by cell key, one
#'   0/1 column per biomarker), e.g. psa_obj$cells$cell_definition.
#' @param biomarkers Biomarker columns, in cell-key order.
#' @return Data frame with one row per policy: `policy` (id `nivo_<pattern>`),
#'   `pattern` (one digit per cell, 1 = nivolumab), `label` (the six named
#'   policies; otherwise "Cells <keys>"), `named` (logical), `paper2_strategy`
#'   (paper 2 strategy id or NA), `n_cells_treated`, one logical
#'   `tests_<biomarker>` column per biomarker, and the assignment as a list
#'   column `assign`.
policy_space_definition <- function(df_cells, biomarkers = get_biomarkers()) {
  v_cells <- rownames(df_cells)
  n_cells <- length(v_cells)
  m_assign <- as.matrix(expand.grid(rep(list(0:1), n_cells)))[, n_cells:1, drop = FALSE]
  m_assign <- m_assign[order(apply(m_assign, 1, paste0, collapse = "")), , drop = FALSE]
  v_pattern <- apply(m_assign, 1, paste0, collapse = "")
  df_metadata <- get_strategy_metadata(c(get_control_strategy(), biomarkers))
  v_rules <- list(
    nobody = rep(0L, n_cells),
    everyone = rep(1L, n_cells),
    either_positive = as.integer(rowSums(as.matrix(df_cells[biomarkers])) > 0),
    both_positive = as.integer(rowSums(as.matrix(df_cells[biomarkers])) == length(biomarkers)))
  v_named_patterns <- c(
    setNames(paste0(v_rules$nobody, collapse = ""), df_metadata$report_label[1]),
    setNames(paste0(v_rules$everyone, collapse = ""), "Nivolumab for everyone"),
    setNames(vapply(biomarkers, function(b) paste0(df_cells[[b]], collapse = ""), character(1)),
             df_metadata$report_label[match(biomarkers, df_metadata$id)]),
    setNames(paste0(v_rules$either_positive, collapse = ""), "Either biomarker positive"),
    setNames(paste0(v_rules$both_positive, collapse = ""), "Both biomarkers positive"))
  v_paper2 <- setNames(c(get_control_strategy(), biomarkers),
                       c(v_named_patterns[[1]], v_named_patterns[2 + seq_along(biomarkers)]))
  df_out <- data.frame(policy = paste0("nivo_", v_pattern), pattern = v_pattern,
                       stringsAsFactors = FALSE)
  df_out$label <- ifelse(v_pattern %in% v_named_patterns,
    names(v_named_patterns)[match(v_pattern, v_named_patterns)],
    vapply(seq_len(nrow(m_assign)), function(k)
      paste("Cells", paste(v_cells[m_assign[k, ] == 1], collapse = ", ")), character(1)))
  df_out$named <- v_pattern %in% v_named_patterns
  df_out$paper2_strategy <- unname(v_paper2[v_pattern])
  df_out$n_cells_treated <- rowSums(m_assign)
  for (biomarker in biomarkers) {
    df_out[[paste0("tests_", biomarker)]] <- vapply(seq_len(nrow(m_assign)), function(k)
      policy_conditions_on(m_assign[k, ], df_cells, biomarker), logical(1))
  }
  df_out$assign <- lapply(seq_len(nrow(m_assign)), function(k) unname(m_assign[k, ]))
  rownames(df_out) <- NULL
  df_out
}

#' Diagnostic-test cost of every policy
#'
#' @param df_policies Output of policy_space_definition().
#' @param params Parameter list holding the scalar test prices
#'   (biomarker_cost_key(): c_test_CRP, c_test_NGS).
#' @param biomarkers Biomarkers.
#' @return Numeric vector, one test cost per policy.
policy_test_costs <- function(df_policies, params, biomarkers = get_biomarkers()) {
  v_prices <- vapply(biomarkers, function(b) {
    price <- params[[biomarker_cost_key(b)[[1]]]]
    if (!is.numeric(price) || length(price) != 1L || !is.finite(price) || price < 0)
      stop("No finite diagnostic-test price for ", b, " (", biomarker_cost_key(b)[[1]], ")")
    price
  }, numeric(1))
  m_tests <- as.matrix(df_policies[paste0("tests_", biomarkers)])
  as.vector(m_tests %*% v_prices)
}

#' Outcomes of every policy for every draw
#'
#' @param a_values Cell outcomes, draw x cell x arm x outcome (or one draw,
#'   cell x arm x outcome).
#' @param m_weights Joint cell weights, draw x cell (or a vector for one draw).
#' @param v_test_costs Test cost per policy (policy_test_costs()).
#' @param df_policies Output of policy_space_definition().
#' @return Named list over the outcomes (cost, qalys, ly_disc, ly_undisc,
#'   nivo_admin) of draw x policy matrices; test costs are in `cost` only.
#'   A cell of zero weight contributes nothing, also when its outcomes are NA
#'   (a cell without patients).
policy_outcome_matrices <- function(a_values, m_weights, v_test_costs, df_policies) {
  if (length(dim(a_values)) == 3L) {
    a_values <- array(a_values, dim = c(1L, dim(a_values)),
                      dimnames = c(list(NULL), dimnames(a_values)))
  }
  m_weights <- matrix(m_weights, ncol = dim(a_values)[2])
  v_outcomes <- dimnames(a_values)[[4]]
  n_draws <- dim(a_values)[1]
  if (nrow(m_weights) != n_draws) stop("Cell weights and cell outcomes differ in draws.")
  if (length(v_test_costs) != nrow(df_policies)) stop("One test cost per policy is required.")
  out <- setNames(lapply(v_outcomes, function(outcome) {
    m_out <- vapply(seq_len(nrow(df_policies)), function(k) {
      v_assign <- df_policies$assign[[k]]
      v_total <- numeric(n_draws)
      for (c_idx in seq_along(v_assign)) {
        v_w <- m_weights[, c_idx]
        v_y <- a_values[, c_idx, v_assign[c_idx] + 1L, outcome]
        v_total <- v_total + ifelse(v_w > 0, v_w * v_y, 0)
      }
      v_total
    }, numeric(n_draws))
    m_out <- matrix(m_out, nrow = n_draws, dimnames = list(NULL, df_policies$policy))
    if (outcome == "cost") m_out <- sweep(m_out, 2, v_test_costs, "+")
    m_out
  }), v_outcomes)
  out
}

#' Cost, QALYs, life-years and administrations of the policies for one draw
#'
#' @param cell_values Cell outcomes of one draw (cell x arm x outcome array,
#'   as in model_fun(return_cells = TRUE)).
#' @param weights Joint cell weights of the draw (named by cell key).
#' @param params Parameter list with the test prices (c_test_CRP, c_test_NGS).
#' @param df_policies Policy space; default all 16 assignments.
#' @return Data frame with one row per policy: `policy`, `pattern`, `label`,
#'   `test_cost`, `cost`, `qalys`, `ly_disc`, `ly_undisc`, `nivo_admin`.
policy_outcomes <- function(cell_values, weights, params, df_policies = NULL) {
  if (is.null(df_policies)) {
    v_cells <- dimnames(cell_values)[[1]]
    df_cells <- expand.grid(setNames(rep(list(0:1), length(get_biomarkers())), get_biomarkers()))
    rownames(df_cells) <- apply(df_cells, 1, paste0, collapse = "")
    df_policies <- policy_space_definition(df_cells[v_cells, , drop = FALSE])
  }
  v_test_costs <- policy_test_costs(df_policies, params)
  l_out <- policy_outcome_matrices(cell_values, weights[dimnames(cell_values)[[1]]],
                                   v_test_costs, df_policies)
  df_out <- df_policies[c("policy", "pattern", "label")]
  df_out$test_cost <- v_test_costs
  for (outcome in names(l_out)) df_out[[outcome]] <- as.vector(l_out[[outcome]][1, ])
  df_out
}

#' Policy space of a PSA object with cell-level outcomes
#'
#' @param psa_obj PSA object with the `cells` component (psa_pair_v2).
#' @param params Parameter list with the test prices (l_params_base).
#' @param exclude Policy ids to leave out (e.g. "nivo_1111" for the 15-policy
#'   sensitivity analysis).
#' @return List with `policies` (policy_space_definition() rows, plus
#'   `test_cost`), `outcomes` (policy_outcome_matrices()), `price_ref` (the
#'   nivolumab price of the cell costs), `n_sim` and `cells` (the component).
policy_space_from_psa <- function(psa_obj, params, exclude = character(0)) {
  cells <- psa_obj$cells
  if (is.null(cells) || !identical(cells$schema, "psa_cells_v1")) {
    stop("psa_obj has no cell-level component; rerun analysis/10_PSA.R (issue #200).")
  }
  df_policies <- policy_space_definition(cells$cell_definition, cells$biomarkers)
  df_policies <- df_policies[!df_policies$policy %in% exclude, , drop = FALSE]
  v_test_costs <- policy_test_costs(df_policies, params, cells$biomarkers)
  df_policies$test_cost <- v_test_costs
  v_price <- unique(cells$nivolumab_price)
  list(policies = df_policies,
       outcomes = policy_outcome_matrices(cells$values, cells$weights, v_test_costs, df_policies),
       price_ref = if (length(v_price) == 1L) v_price else cells$nivolumab_price,
       n_sim = nrow(cells$weights), cells = cells)
}

#' Price pairing object of a policy space
#'
#' @param policy_space Output of policy_space_from_psa().
#' @param keep Policy ids to keep; NULL keeps all.
#' @return Object of class `voi_price_pairing` (R/voi_price_axis.R) with the
#'   policies as strategies: `cost_ref` at `price_ref`, `effect` (QALYs),
#'   `life_years` (discounted), `administrations`.
policy_price_pairing <- function(policy_space, keep = NULL) {
  if (length(policy_space$price_ref) != 1L) {
    stop("The cell costs were evaluated at more than one nivolumab price.")
  }
  v_keep <- if (is.null(keep)) policy_space$policies$policy else keep
  if (!all(v_keep %in% policy_space$policies$policy)) stop("Unknown policies: ",
    paste(setdiff(v_keep, policy_space$policies$policy), collapse = ", "))
  pick <- function(outcome) policy_space$outcomes[[outcome]][, v_keep, drop = FALSE]
  structure(list(
    price_ref = policy_space$price_ref, price_alt = NA_real_, strategies = v_keep,
    n_sim = policy_space$n_sim, cost_ref = pick("cost"), effect = pick("qalys"),
    life_years = pick("ly_disc"), administrations = pick("nivo_admin"),
    checks = NULL, administrations_vs_table = NULL), class = "voi_price_pairing")
}

#' Gap between the cell-based policies and the PSA strategy outcomes
#'
#' @param policy_space Output of policy_space_from_psa() for `psa_obj`.
#' @param psa_obj The PSA object (cost, effectiveness, ly_disc, ly_undisc).
#' @return Data frame with one row per draw and paper 2 strategy: `draw`,
#'   `strategy`, `policy`, the gaps (policy minus PSA strategy) `gap_cost`,
#'   `gap_qalys`, `gap_ly_disc`, `gap_ly_undisc`, and `cell_violation`
#'   (TRUE when any cell curve of the draw crossed before the cap).
strategy_reproduction_gaps <- function(policy_space, psa_obj) {
  df_map <- policy_space$policies[!is.na(policy_space$policies$paper2_strategy), , drop = FALSE]
  outcomes <- psa_outcome_matrices(psa_obj)
  v_violation <- apply(policy_space$cells$violations, 1, function(x) any(x > 0, na.rm = TRUE))
  do.call(rbind, lapply(seq_len(nrow(df_map)), function(k) {
    strategy <- df_map$paper2_strategy[k]
    policy <- df_map$policy[k]
    data.frame(
      draw = seq_len(policy_space$n_sim), strategy = strategy, policy = policy,
      gap_cost = policy_space$outcomes$cost[, policy] - outcomes$cost[, strategy],
      gap_qalys = policy_space$outcomes$qalys[, policy] - outcomes$effect[, strategy],
      gap_ly_disc = policy_space$outcomes$ly_disc[, policy] - psa_obj$ly_disc[[strategy]],
      gap_ly_undisc = policy_space$outcomes$ly_undisc[, policy] - psa_obj$ly_undisc[[strategy]],
      cell_violation = v_violation, stringsAsFactors = FALSE)
  }))
}

#' Summary of the reproduction gaps per strategy
#'
#' @param df_gaps Output of strategy_reproduction_gaps().
#' @param tol_cost,tol_effect Absolute tolerances (EUR; QALYs and life-years)
#'   for "reproduced" in draws without a cell crossing (floating point).
#' @return Data frame per strategy: maximum absolute gaps, draws with any
#'   cell crossing, maximum gaps in draws without a crossing, and
#'   `reproduced_without_crossing` (all such gaps within `tol` relative).
summarise_reproduction_gaps <- function(df_gaps, tol_cost = 1e-6, tol_effect = 1e-10) {
  do.call(rbind, lapply(split(df_gaps, df_gaps$strategy), function(d) {
    v_clean <- !d$cell_violation
    max_clean <- function(x) if (any(v_clean)) max(abs(x[v_clean])) else 0
    data.frame(strategy = d$strategy[1], policy = d$policy[1], n_draws = nrow(d),
      n_draws_cell_crossing = sum(d$cell_violation),
      max_abs_gap_cost = max(abs(d$gap_cost)), max_abs_gap_qalys = max(abs(d$gap_qalys)),
      max_abs_gap_ly_disc = max(abs(d$gap_ly_disc)),
      max_abs_gap_cost_no_crossing = max_clean(d$gap_cost),
      max_abs_gap_qalys_no_crossing = max_clean(d$gap_qalys),
      mean_gap_cost = mean(d$gap_cost), mean_gap_qalys = mean(d$gap_qalys),
      reproduced_without_crossing = max_clean(d$gap_cost) <= tol_cost &&
        max_clean(d$gap_qalys) <= tol_effect && max_clean(d$gap_ly_disc) <= tol_effect &&
        max_clean(d$gap_ly_undisc) <= tol_effect,
      stringsAsFactors = FALSE)
  }))
}

#' Summary of the cell-level outcomes (voi_cells.csv)
#'
#' @param cells The `cells` component of a PSA object.
#' @return Data frame with one row per cell, arm and outcome: the biomarker
#'   values, mean cell weight, mean, SD and 95% percentile interval (type 7) of
#'   the outcome, and the mean and interval of the increment experimental
#'   minus control in that cell.
summarise_cell_outcomes <- function(cells) {
  rows <- list()
  for (cell in cells$cells) for (arm in cells$arms) for (outcome in cells$outcomes) {
    v_y <- cells$values[, cell, arm, outcome]
    v_inc <- cells$values[, cell, "experimental", outcome] - cells$values[, cell, "control", outcome]
    df_row <- data.frame(cell = cell, stringsAsFactors = FALSE)
    for (b in cells$biomarkers) df_row[[b]] <- cells$cell_definition[cell, b]
    df_row$arm <- arm
    df_row$outcome <- outcome
    df_row$mean_weight <- mean(cells$weights[, cell])
    df_row$mean <- mean(v_y)
    df_row$sd <- stats::sd(v_y)
    df_row$p2_5 <- unname(stats::quantile(v_y, 0.025, type = 7, na.rm = TRUE))
    df_row$p97_5 <- unname(stats::quantile(v_y, 0.975, type = 7, na.rm = TRUE))
    df_row$mean_increment <- mean(v_inc)
    df_row$increment_p2_5 <- unname(stats::quantile(v_inc, 0.025, type = 7, na.rm = TRUE))
    df_row$increment_p97_5 <- unname(stats::quantile(v_inc, 0.975, type = 7, na.rm = TRUE))
    df_row$n_draws_crossing <- sum(cells$violations[, cell, arm] > 0, na.rm = TRUE)
    df_row$n_sim <- nrow(cells$weights)
    rows[[length(rows) + 1]] <- df_row
  }
  df_out <- do.call(rbind, rows)
  rownames(df_out) <- NULL
  df_out
}

#' Summary of the policies (voi_policies.csv)
#'
#' @param policy_space Output of policy_space_from_psa().
#' @param wtp Willingness-to-pay threshold for the NMB columns.
#' @param prices Named nivolumab prices at which to report the expected NMB
#'   and the probability of being optimal.
#' @param df_gap_summary Optional summarise_reproduction_gaps() output, joined
#'   on the paper 2 strategy.
#' @return Data frame with one row per policy.
summarise_policies <- function(policy_space, wtp, prices, df_gap_summary = NULL) {
  df_out <- policy_space$policies
  df_out$assign <- NULL
  for (outcome in names(policy_space$outcomes)) {
    df_out[[paste0("mean_", outcome)]] <- colMeans(policy_space$outcomes[[outcome]])
  }
  pairing <- policy_price_pairing(policy_space)
  for (k in seq_along(prices)) {
    result <- evpi_at_price(pairing, prices[[k]], wtp)
    df_out[[paste0("expected_nmb_", names(prices)[k])]] <- unname(result$expected_nmb[df_out$policy])
    df_out[[paste0("prob_optimal_", names(prices)[k])]] <- unname(result$prob_optimal[df_out$policy])
  }
  df_out$wtp <- wtp
  df_out$price_ref <- policy_space$price_ref
  df_out$n_sim <- policy_space$n_sim
  if (!is.null(df_gap_summary)) {
    v_cols <- setdiff(names(df_gap_summary), c("strategy", "policy"))
    v_match <- match(df_out$policy, df_gap_summary$policy)
    for (col in v_cols) df_out[[paste0("reproduction_", col)]] <- df_gap_summary[[col]][v_match]
  }
  df_out
}

#' Run paper 3 step 1b on a PSA object with cell-level outcomes
#'
#' @param psa_obj PSA object of the price at which the cell costs were
#'   evaluated (validated pair, psa_pair_v2).
#' @param params Parameter list with the test prices (l_params_base).
#' @param population Population inputs (current_population_inputs()).
#' @param settings List with `wtp` (thresholds of the price axis), `wtp_base`
#'   (threshold of voi_policies.csv), `step` (price grid step),
#'   `marked_prices` (named prices that must be grid points, e.g.
#'   c(biosimilar = 464)), and `everyone` (id of the treat-all policy,
#'   excluded in the 15-policy space).
#' @return List with `cells` (voi_cells.csv), `policies` (voi_policies.csv),
#'   `evpi_price` (voi_evpi_price.csv: 16 and 15 policies), `gaps` (per draw)
#'   and `gap_summary`.
run_voi_step1b <- function(psa_obj, params, population, settings) {
  policy_space <- policy_space_from_psa(psa_obj, params)
  df_gaps <- strategy_reproduction_gaps(policy_space, psa_obj)
  df_gap_summary <- summarise_reproduction_gaps(df_gaps)
  v_marked <- c(list_price = policy_space$price_ref, settings$marked_prices)
  v_grid <- price_axis_grid(policy_space$price_ref, v_marked, step = settings$step)
  run_space <- function(keep, label) {
    df_axis <- run_price_axis_evpi(policy_price_pairing(policy_space, keep), v_grid,
                                   settings$wtp, population, v_marked)
    cbind(policy_space = label, df_axis, stringsAsFactors = FALSE)
  }
  v_all <- policy_space$policies$policy
  df_16 <- run_space(v_all, sprintf("%d policies", length(v_all)))
  df_15 <- run_space(setdiff(v_all, settings$everyone),
                     sprintf("%d policies (everyone excluded)", length(v_all) - 1L))
  # The 15-policy rows lack the prob_optimal column of the excluded policy.
  for (col in setdiff(names(df_16), names(df_15))) df_15[[col]] <- NA_real_
  list(cells = summarise_cell_outcomes(policy_space$cells),
       policies = summarise_policies(policy_space, settings$wtp_base, v_marked, df_gap_summary),
       evpi_price = rbind(df_16, df_15[names(df_16)]),
       gaps = df_gaps, gap_summary = df_gap_summary)
}

#' Write the step 1b CSV outputs
#'
#' @param results List returned by run_voi_step1b().
#' @param directory Output directory (created if needed); default
#'   data/output/voi.
#' @return Paths of `voi_cells.csv`, `voi_policies.csv` and
#'   `voi_evpi_price.csv`, invisibly.
write_voi_step1b_csvs <- function(results, directory = voi_output_dir()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  v_paths <- file.path(directory, c("voi_cells.csv", "voi_policies.csv", "voi_evpi_price.csv"))
  utils::write.csv(results$cells, v_paths[1], row.names = FALSE)
  utils::write.csv(results$policies, v_paths[2], row.names = FALSE)
  utils::write.csv(results$evpi_price, v_paths[3], row.names = FALSE)
  invisible(v_paths)
}
