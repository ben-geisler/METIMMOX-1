# ===============================================================================
# PSA ON THE DRAWS OF AN EVIDENCE STATE (paper 3, E2; issue #203)
# ===============================================================================
# The published-trials state E2 is the closed-form (conjugate) update of the
# joint survival-coefficient distribution; build_evidence_state() regenerates
# its draws in the sampling-cache layout. Ruling (9) of 11 October 2026: before
# any E2 number is reported, the net benefit of E2 comes from a PSA on those
# draws, not from resampling the rows of the E0 PSA. run_evidence_psa() runs
# the production PSA (build_psa_obj(), the code of analysis/10_PSA.R) on the
# state's draws, with the economic-parameter draws, seed, base-case
# parameters, strategies and time grid of E0, so the two states differ only in
# the survival evidence. Draw i of a state regenerated with the base seed is
# draw i of E0 moved to the posterior, so the survival draws are paired too.
#
# The PSA is a paper 3 cache (psa_evidence_<label>_<utility label>.rds), never
# the paper 2 PSA path. Its identity (evidence_psa_identity()) records the PSA
# inputs with the state's fingerprint as the sampling fingerprint, the
# upstream fingerprints (base sampling cache, E0 PSA generation, evidence
# state) and the `voi` implementation; the loader names the differing input.
# The importance resample of the E0 PSA stays as the reported check
# (evidence_state_comparison()).
#
# Requires R/model_fun.R, R/psa_functions.R (build_psa_obj()), R/cache_paths.R
# and R/cache_provenance.R (psa_cache_fingerprint(), bind_psa_pair(),
# validate_psa_pair(), voi_upstream_fingerprints()), R/voi_cache.R
# (save_voi_results(), load_voi_results()), R/voi_price_axis.R
# (costs_at_price()) and R/evppi_functions.R (calculate_evpi_from_nmb()).
# ===============================================================================

#' Schema of an evidence-state PSA cache
EVIDENCE_PSA_SCHEMA <- "voi_evidence_psa_v1"

#' Path of the PSA cache of an evidence state
#'
#' @param label Evidence-state label (letters, digits, "_" and "-").
#' @param util_label Utility-source label (resolve_util_label()).
#' @param directory Cache directory.
#' @return File path "psa_evidence_<label>_<util_label>.rds"; never the paper 2
#'   PSA path.
evidence_psa_cache_path <- function(label, util_label = resolve_util_label(), directory = cache_dir()) {
  if (!grepl("^[A-Za-z0-9_-]+$", label)) stop("Evidence-state labels may contain letters, digits, '_' and '-' only")
  file.path(directory, paste0("psa_evidence_", label, "_", util_label, ".rds"))
}

#' Evaluate code with another survival-draw object as `sampling_models`
#'
#' model_fun() reads the survival draws of a PSA from the free variable
#' `sampling_models` in its own environment. This sets it to `state` for the
#' evaluation of `code` and restores the previous value (or removes the
#' binding) afterwards, also on error.
#'
#' @param state Object in the sampling-cache layout (`joint` component).
#' @param code Expression, evaluated lazily after the swap.
#' @param envir Environment holding the binding (default model_fun()'s).
#' @return The value of `code`.
with_sampling_models <- function(state, code, envir = environment(model_fun)) {
  if (is.null(state$joint)) stop("The evidence state has no joint sampling component")
  had <- exists("sampling_models", envir = envir, inherits = FALSE)
  old <- if (had) get("sampling_models", envir = envir, inherits = FALSE) else NULL
  assign("sampling_models", state, envir = envir)
  on.exit(if (had) assign("sampling_models", old, envir = envir) else
    rm(list = "sampling_models", envir = envir), add = TRUE)
  code
}

#' Check that two PSAs share their economic-parameter draws
#'
#' @param base_params Retained parameter rows of the E0 PSA (`sim` column).
#' @param state_params Retained parameter rows of the evidence-state PSA.
#' @param parameters Economic parameter columns (names(param_distributions)).
#' @return Number of draws (`sim`) both PSAs retain; stops if any shared draw
#'   differs in an economic parameter.
check_common_economic_draws <- function(base_params, state_params, parameters) {
  v_missing <- setdiff(parameters, intersect(names(base_params), names(state_params)))
  if (length(v_missing)) stop("Economic parameters missing from a PSA: ", paste(v_missing, collapse = ", "))
  v_sim <- intersect(base_params$sim, state_params$sim)
  if (!length(v_sim)) stop("The two PSAs share no draw")
  df_a <- base_params[match(v_sim, base_params$sim), parameters, drop = FALSE]
  df_b <- state_params[match(v_sim, state_params$sim), parameters, drop = FALSE]
  v_differ <- parameters[!vapply(parameters, function(p) identical(df_a[[p]], df_b[[p]]), logical(1))]
  if (length(v_differ)) {
    stop("The evidence-state PSA does not reuse the E0 economic draws (", paste(v_differ, collapse = ", "),
         "); same seed, n_sim and parameter distributions are required", call. = FALSE)
  }
  length(v_sim)
}

#' Run the production PSA on the draws of an evidence state
#'
#' build_psa_obj() (life-years and cell outcomes included) with the state's
#' survival draws as `sampling_models` and everything else as in the E0 PSA;
#' the result is stamped with psa_cache_fingerprint() on the state's
#' fingerprint and bound (bind_psa_pair()).
#'
#' @param state Evidence state (build_evidence_state()$state: `joint`,
#'   `fingerprint`, `evidence_state`).
#' @param l_params_base,param_distributions,strategies,time_horizon,cl,n_sim,seed
#'   Inputs of the E0 PSA (analysis/10_PSA.R).
#' @param base_psa_params Retained parameter rows of the E0 PSA, or NULL to
#'   skip the shared-draw check.
#' @param psa_fingerprint psa_cache_fingerprint() of the state (default
#'   computed here).
#' @param fallback_threshold Maximum initial fallback rate.
#' @param envir Environment of the `sampling_models` binding.
#' @return List with `psa_obj` and `psa_params` (a bound psa_pair_v2) and
#'   `diagnostics` (label, state fingerprint, requested and retained draws,
#'   fallback and dropped counts, draws shared with E0, minutes).
run_evidence_psa <- function(state, l_params_base, param_distributions, strategies, time_horizon, cl, n_sim,
                             seed, base_psa_params = NULL, psa_fingerprint = NULL, fallback_threshold = 0.02,
                             envir = environment(model_fun)) {
  if (!identical(state$joint$method, "mvn_joint_v2") || is.null(state$fingerprint)) {
    stop("An evidence-state PSA needs an mvn_joint_v2 state with a fingerprint")
  }
  if (is.null(psa_fingerprint)) {
    psa_fingerprint <- psa_cache_fingerprint(state$fingerprint, l_params_base, param_distributions, strategies,
                                             n_sim, seed, time_horizon, cl)
  }
  started <- Sys.time()
  build <- with_sampling_models(state, build_psa_obj(
    l_params_base = l_params_base, param_distributions = param_distributions, strategies = strategies,
    time_horizon = time_horizon, cl = cl, n_sim = n_sim, seed = seed, fallback_threshold = fallback_threshold,
    extended_outcomes = TRUE), envir = envir)
  psa_obj <- build$psa_obj
  psa_obj$fingerprint <- psa_fingerprint$fingerprint
  psa_obj$fingerprint_inputs <- psa_fingerprint$inputs
  psa_obj$sampling_method <- state$joint$method
  psa_obj$evidence_state <- state$evidence_state
  pair <- bind_psa_pair(psa_obj, build$psa_params)
  n_common <- if (is.null(base_psa_params)) NA_integer_ else
    check_common_economic_draws(base_psa_params, pair$psa_params, names(param_distributions))
  list(psa_obj = pair$psa_obj, psa_params = pair$psa_params,
       diagnostics = list(evidence_state = state$evidence_state, state_fingerprint = state$fingerprint,
                          requested_n_sim = n_sim, n_sim = pair$psa_obj$n_sim,
                          fallback_count = build$psa_results$fallback_count,
                          dropped_count = build$psa_results$dropped_count,
                          draws_shared_with_e0 = n_common,
                          minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))))
}

#' Identity of an evidence-state PSA cache
#'
#' @param state Evidence state (`fingerprint`, `evidence_state`).
#' @param psa_inputs psa_cache_fingerprint()$inputs of the state (the PSA
#'   inputs, with the state's fingerprint as the sampling fingerprint and the
#'   `psa` implementation).
#' @param base_sampling_models Base (E0) sampling cache, or NULL.
#' @param base_psa_obj Bound E0 PSA object (its generation fixes the shared
#'   economic draws), or NULL.
#' @return List with `fingerprint` and `inputs` (`schema`, `evidence_state`,
#'   `psa`, `upstream`, `implementation`).
evidence_psa_identity <- function(state, psa_inputs, base_sampling_models = NULL, base_psa_obj = NULL) {
  v_label <- state$evidence_state
  if (!is.character(v_label) || length(v_label) != 1L || !nzchar(v_label)) stop("The state has no label")
  inputs <- list(schema = EVIDENCE_PSA_SCHEMA, evidence_state = v_label, psa = psa_inputs,
                 upstream = voi_upstream_fingerprints(sampling_models = base_sampling_models, psa_obj = base_psa_obj,
                                                      evidence_states = stats::setNames(list(state), v_label)),
                 implementation = list(voi = calculation_identity("voi")))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Load an evidence-state PSA and validate it
#'
#' @param path Cache path (evidence_psa_cache_path()).
#' @param expected Current identity (evidence_psa_identity()).
#' @param required TRUE stops on an absent or stale cache, naming the
#'   differing input; FALSE returns NULL (the producer then reruns it).
#' @return List with `psa_obj`, `psa_params`, `diagnostics`, or NULL.
load_evidence_psa <- function(path, expected, required = FALSE) {
  results <- load_voi_results(path, expected, "evidence-state PSA",
                              "Rerun analysis/14_VOI.R with the published-trials state (E2).", required)
  if (is.null(results)) return(NULL)
  validate_psa_pair(results$psa_obj, results$psa_params)
  results
}

#' The PSA of an evidence state: from its cache when current, otherwise run and saved
#'
#' @param state Evidence state (build_evidence_state()$state).
#' @param l_params_base,param_distributions,strategies,time_horizon,cl,n_sim,seed
#'   Inputs of the E0 PSA.
#' @param base_sampling_models Base (E0) sampling cache.
#' @param base_psa_obj,base_psa_params Bound E0 PSA pair.
#' @param path Cache path (default evidence_psa_cache_path(state label)).
#' @param envir Environment of the `sampling_models` binding.
#' @return run_evidence_psa() result plus `identity`, `path` and
#'   `from_cache`.
load_or_run_evidence_psa <- function(state, l_params_base, param_distributions, strategies, time_horizon, cl,
                                     n_sim, seed, base_sampling_models, base_psa_obj, base_psa_params,
                                     path = evidence_psa_cache_path(state$evidence_state),
                                     envir = environment(model_fun)) {
  base_psa <- normalizePath(c(psa_obj_path(directory = dirname(path)), psa_params_path(directory = dirname(path))),
                            mustWork = FALSE)
  if (normalizePath(path, mustWork = FALSE) %in% base_psa) stop("Refusing to use the E0 PSA path for an evidence state")
  psa_fp <- psa_cache_fingerprint(state$fingerprint, l_params_base, param_distributions, strategies, n_sim, seed,
                                  time_horizon, cl)
  identity <- evidence_psa_identity(state, psa_fp$inputs, base_sampling_models, base_psa_obj)
  cached <- load_evidence_psa(path, identity, required = FALSE)
  if (!is.null(cached)) return(c(cached, list(identity = identity, path = path, from_cache = TRUE)))
  result <- run_evidence_psa(state, l_params_base, param_distributions, strategies, time_horizon, cl, n_sim, seed,
                             base_psa_params = base_psa_params, psa_fingerprint = psa_fp, envir = envir)
  save_voi_results(result, identity, path)
  c(result, list(identity = identity, path = path, from_cache = FALSE))
}

#' Compare an evidence state's PSA with the importance resample of the E0 PSA
#'
#' The check of ruling (9) of 11 October 2026: per price, the EVPI and the
#' mean incremental NMB of every policy versus the first (standard of care)
#' from the resampled E0 rows and from the state's own PSA. They differ by
#' Monte Carlo error and by the linearisation of the conjugate update.
#'
#' @param pairing_base voi_price_pairing of the E0 policies.
#' @param rows Resampled E0 rows (evidence_state_rows()$index).
#' @param pairing_state voi_price_pairing of the state's PSA (same policies).
#' @param prices Nivolumab prices.
#' @param wtp Willingness-to-pay threshold.
#' @return Data frame per price: `evpi_resampled`, `evpi_psa`,
#'   `max_abs_inmb_difference` (largest difference of a policy's mean
#'   incremental NMB), `policy_max_difference`.
evidence_state_comparison <- function(pairing_base, rows, pairing_state, prices, wtp) {
  if (!identical(pairing_base$strategies, pairing_state$strategies)) stop("The two PSAs have different policies")
  rows_df <- lapply(prices, function(p) {
    m_base <- as.matrix(pairing_base$effect * wtp - costs_at_price(pairing_base, p))[rows, , drop = FALSE]
    m_state <- as.matrix(pairing_state$effect * wtp - costs_at_price(pairing_state, p))
    v_inmb <- function(m) colMeans(m - m[, 1])
    v_diff <- abs(v_inmb(m_base) - v_inmb(m_state))
    data.frame(price = p, wtp = wtp, evpi_resampled = calculate_evpi_from_nmb(m_base),
               evpi_psa = calculate_evpi_from_nmb(m_state), max_abs_inmb_difference = max(v_diff),
               policy_max_difference = names(v_diff)[which.max(v_diff)], stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows_df)
  rownames(df) <- NULL
  df
}
