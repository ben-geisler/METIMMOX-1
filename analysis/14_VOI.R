# ===============================================================================
# 14 VALUE OF INFORMATION FOR PAPER 3 (issues #199, #200, #202, #203)
# ===============================================================================
# Runs in the pipeline process after 12 (make.R), or alone from a session in
# which analysis/02-06 have been sourced (globals, l_params_base,
# sampling_models). It only READS the paper 2 caches and stops with a pipeline
# instruction when one is missing or stale; it never regenerates the sampling,
# PSA, EVPPI or scenario caches. It writes paper 3 caches of its own (step 1a,
# the E2 evidence state and its PSA, steps 5-6), each next to the base caches
# under its own name and fingerprint.
#
#   step 1a  price pairing, price-axis EVPI, extended EVPPI
#            (load_or_run_voi_step1a(); data/output/voi/voi_step1a_results.rds,
#            shared with reports/technical/value_of_research.qmd)
#   step 1b  cell-level outcomes, 16-policy space, EVPI on the price axis
#            (run_voi_step1b(); voi_cells.csv, voi_policies.csv, voi_evpi_price.csv)
#   E0, E2   METIMMOX alone; METIMMOX plus one result per trial and endpoint
#            of data/external/voi_external_evidence.csv chosen by the named
#            rule of select_evidence_rows() (pMMR/MSS over ITT, latest cut,
#            95% interval; ruling (2) of 11 October 2026), built with
#            build_evidence_state() (tau elicitation pending; saved as
#            sampling_models_n<n>_evidence_E2.rds) and evaluated with a PSA on
#            its own draws with the economic draws and seed of E0
#            (load_or_run_evidence_psa(); psa_evidence_E2_<label>.rds, about 20
#            minutes when not current; ruling (9) of 11 October 2026); the
#            importance resample of the E0 PSA (evidence_state_rows()) is
#            reported as the check (evidence_state_comparison())
#   step 5   candidate designs S1, S2, S3 (METIMMOX-3 as in protocol v1.1 of 1
#            October 2026, #213: non-randomised, stratified by TMB/BRAF), S3R
#            (its randomised counterpart), S4,
#            the power-based S6 and the S4 + Sx combinations
#   step 6   EVSI by design, n and price; maximum justifiable research cost,
#            ENBS, optimal n, research threshold price (evsi_by_design.csv,
#            voi_designs.csv, voi_research_value.csv, enbs_by_design.csv,
#            voi_threshold_price.csv; voi_power_sizes.csv with the
#            engine-based and Schoenfeld power-based n under the conventional
#            and the decision-relevant hypothesis; voi_decision_effects.csv;
#            ENBS for the oncology cost band and, in the appendix, Speich 2018)
#
# Steps 5-6 are cached in data/tidy/evsi_results_research_<label>.rds with a
# fingerprint over the PSA generation, the sampling cache, the evidence file,
# the settings, every named design/cost parameter and the code
# (voi_research_fingerprint()); a rerun with unchanged inputs only rewrites
# the CSVs. Settings can be overridden with options(metimmox.voi14_settings =
# list(...)).
#
# Grid (ruling (10) of 11 October 2026): the default run is the REDUCED grid
# of voi14_grid_settings(): the marked prices only (list, biosimilar,
# adoption threshold), both evidence states, Norway only and n in
# VOI14_N_GRID_REDUCED (40-1280). The full grid (price axis, VOI_N_GRID,
# Norway, Nordic, Europe) is options(metimmox.voi14_settings = list(grid =
# "full")). Steps 1a and 1b always run first.
#
# METIMMOX-only switch (user, 12 October 2026): evidence_states = "E0" skips
# the published-trials state E2 (neither built nor evaluated, no E2 PSA, not
# in the fingerprint). From the command line:
#   Rscript make.R --from 14 --no-vignettes --no-reports --no-tests --voi-states E0
# (make.R passes it as METIMMOX_VOI14_STATES); in a session:
#   options(metimmox.voi14_settings = list(evidence_states = "E0"))
#
# Runtime: not yet measured on live caches; step 1a takes about 50 minutes
# when its result is not saved, the E2 PSA about 20 minutes when not current.
# voi14_runtime_estimate() times one information call and one EVSI
# evaluation and prints the extrapolated runtime of steps 5-6 before they run.
# ===============================================================================

if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(here, dampack, flexsurv, survival, mvtnorm, voi, mgcv, Matrix)

for (f in c("model_fun", "calculate_outcomes", "prediction_functions", "prediction_population",
            "parameter_distributions", "psa_functions", "evppi_functions", "scenario_analysis", "cache_paths",
            "cache_provenance", "cea_helpers", "joint_survival_sampling", "pfs_os_violation_diagnostics",
            "voi_price_axis", "voi_survival_parameters", "voi_policies", "design_information",
            "evidence_update", "evsi_functions", "voi_cache", "voi_designs", "voi_research_value",
            "voi_evidence_psa")) {
  source(here::here("R", paste0(f, ".R")))
}

v_required <- c("l_params_base", "sampling_models", "analysis_seed", "WTP", "n_sim", "strategies",
                "param_distributions", "annual_incidence_norway", "annual_incidence_scenarios",
                "research_horizon_years", "discount_rate_research", "time_horizon", "cl")
v_absent <- v_required[!vapply(v_required, exists, logical(1))]
if (length(v_absent)) {
  stop("analysis/14_VOI.R needs a session in which analysis/02-06 have been sourced; missing: ",
       paste(v_absent, collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------------------------
# Validated inputs (stop on missing or stale caches)
# ------------------------------------------------------------------------------
psa_obj <- load_current_psa_cache()
psa_params <- load_psa_params_cache(required = TRUE, seed = analysis_seed, n_sim = psa_obj$n_sim)
validate_psa_pair(psa_obj, psa_params)
if (is.null(psa_obj$cells) || !identical(psa_obj$cells$schema, "psa_cells_v1")) {
  stop("The PSA cache has no cell-level component; rerun analysis/10_PSA.R (issues #197, #200).", call. = FALSE)
}
scenario_cache <- load_scenario_cache()
if (!identical(scenario_cache$failed_draw_policy, psa_failed_draw_policy())) {
  stop("Scenario cache uses another failed-draw policy; rerun analysis/12_scenario_EVPPIs.R.", call. = FALSE)
}
threshold_file <- here::here("outputs", "tables", "cea_threshold_price.csv")
s7_file <- here::here("outputs", "tables", "table_s7.csv")
for (file in c(threshold_file, s7_file)) {
  if (!file.exists(file)) {
    stop("Missing ", file, "; render cea_numbers.qmd and table_s7.qmd first.", call. = FALSE)
  }
}
df_threshold <- utils::read.csv(threshold_file)
df_table_s7 <- utils::read.csv(s7_file)
component <- get_joint_sampling_models(sampling_models)
if (!identical(component$method, "mvn_joint_v2")) {
  stop("The sampling cache is not mvn_joint_v2; rerun analysis/06_sampling.R.", call. = FALSE)
}
if (!"model_idx" %in% names(psa_params)) stop("psa_params has no model_idx; rerun analysis/10_PSA.R.", call. = FALSE)
reference_population <- l_params_base$prediction_population
if (is.null(reference_population)) stop("l_params_base has no prediction_population; rerun analysis/05.", call. = FALSE)

# ------------------------------------------------------------------------------
# Step 1a: price axis and extended EVPPI (same settings as value_of_research.qmd,
# so both read one saved result)
# ------------------------------------------------------------------------------
cat("\n=== 14 VoI: step 1a (price axis, extended EVPPI) ===\n")
v_scenario_wtp <- sort(unique(scenario_cache$scenarios$wtp))
settings_1a <- list(
  wtp = v_scenario_wtp, step = VOI_PRICE_STEP, seed = analysis_seed, B = 1000L,
  control_strategy = get_control_strategy(),
  max_administrations = sum(l_params_base$l_nivo),
  cross_check = list(survival_all = c("earth", "bart")),
  cross_check_args = list(earth = list(degree = 2)))
voi_step1a <- load_or_run_voi_step1a(scenario_cache, sampling_models, df_threshold, df_table_s7,
                                     base_config = current_evppi_config(),
                                     population = current_population_inputs(), settings = settings_1a)
price_list <- unname(voi_step1a$prices["reference"])
price_biosimilar <- unname(voi_step1a$prices["alternative"])
price_adoption <- voi_step1a$adoption_threshold$price

# ------------------------------------------------------------------------------
# Step 1b: cells and the 16-policy space
# ------------------------------------------------------------------------------
cat("\n=== 14 VoI: step 1b (cells, 16 policies) ===\n")
settings_1b <- list(wtp = v_scenario_wtp, wtp_base = WTP, step = VOI_PRICE_STEP,
                    marked_prices = c(biosimilar = price_biosimilar, adoption_threshold = price_adoption),
                    everyone = "nivo_1111")
voi_step1b <- run_voi_step1b(psa_obj, l_params_base, current_population_inputs(), settings_1b)
write_voi_step1b_csvs(voi_step1b)
policy_space <- policy_space_from_psa(psa_obj, l_params_base)
pairing16 <- policy_price_pairing(policy_space)

# ------------------------------------------------------------------------------
# Steps 5-6 settings
# ------------------------------------------------------------------------------
v_override <- getOption("metimmox.voi14_settings")
voi14_grid <- voi14_grid_settings(
  if (is.null(v_override$grid)) "reduced" else v_override$grid,
  marked_prices = c(list = price_list, biosimilar = price_biosimilar, adoption_threshold = price_adoption),
  price_axis = price_axis_grid(price_list, c(price_biosimilar, price_adoption), step = 200),
  evidence_states = voi14_requested_states())
voi14_settings <- c(voi14_grid, list(
  wtp = WTP,
  horizons = sort(unique(c(research_horizon_years, 5, 15))),
  primary_horizon = research_horizon_years,
  rate = discount_rate_research,
  methods = EVSI_METHODS,
  primary_method = "gam_regression",
  seed = analysis_seed,
  B = 1000L,
  tau_m2 = TAU_M2_DEFAULT,
  tau_soc = TAU_SOC_DEFAULT,
  include_power = TRUE,
  cost_bands = trial_cost_bands(),
  reference_rate = accrual_per_year(METIMMOX_ACCRUAL_PATIENTS, METIMMOX_ACCRUAL_START, METIMMOX_ACCRUAL_END),
  e2_tau_shared = E2_TAU_SHARED_DEFAULT,
  e2_tau_trial = E2_TAU_TRIAL_DEFAULT,
  e2_maintenance_landmark_weeks = E2_MAINTENANCE_LANDMARK_WEEKS,
  e2_min_ess = E2_MIN_ESS))
if (length(v_override)) voi14_settings[names(v_override)] <- v_override
voi14_settings$evidence_states <- voi14_evidence_states(voi14_settings$evidence_states)
run_e2 <- "E2" %in% voi14_settings$evidence_states
if (!run_e2) {
  # METIMMOX alone: the E2 settings do not enter the run or its fingerprint.
  voi14_settings <- voi14_settings[!startsWith(names(voi14_settings), "e2_")]
}
eligible <- c(norway = annual_incidence_norway, annual_incidence_scenarios)
if (!all(c("nordic", "europe") %in% names(eligible))) {
  stop("annual_incidence_scenarios must name 'nordic' and 'europe' (analysis/02).", call. = FALSE)
}
eligible <- voi14_eligible(eligible, voi14_settings$accrual)
cat("\n=== 14 VoI: ", voi14_settings$grid, " grid: ", length(voi14_settings$prices), " price(s), ",
    length(voi14_settings$n_grid), " sample sizes, accrual ", paste(voi14_settings$accrual, collapse = ", "),
    ", evidence states ", paste(voi14_settings$evidence_states, collapse = ", "),
    if (!run_e2) " (METIMMOX alone)", " ===\n", sep = "")

evidence_raw <- NULL
evidence_selection <- NULL
evidence_rows <- NULL
if (run_e2) {
  evidence_raw <- read_external_evidence()
  evidence_selection <- select_evidence_rows(evidence_raw, subgroups = TRUE,
                                             maintenance_landmark_weeks = voi14_settings$e2_maintenance_landmark_weeks)
  evidence_rows <- evidence_selection$selected
}

research_identity <- voi_research_fingerprint(list(
  upstream = voi_upstream_fingerprints(sampling_models, psa_obj, scenario_cache),
  model_idx = psa_params$model_idx,
  evidence_md5 = if (run_e2) attr(evidence_raw, "md5"), evidence_rows = evidence_rows,
  settings = voi14_settings, eligible = eligible, specs = candidate_designs(),
  parameters = voi_research_parameters(),
  test_prices = list(crp = l_params_base$c_test_CRP, ngs = l_params_base$c_test_NGS),
  reference_population = reference_population_hash(reference_population)))
research_cache <- evsi_cache_path(paste0("research_", resolve_util_label()))
# A stale cache names the differing inputs (e.g. upstream$psa$generation_id) and is recomputed (#205).
voi_research <- load_voi_results(research_cache, research_identity, "steps 5-6 results",
                                 "Recomputing steps 5-6.", required = FALSE)
if (!is.null(voi_research)) {
  cat("\n=== 14 VoI: steps 5-6 read from", basename(research_cache), "(inputs unchanged) ===\n")
}

if (is.null(voi_research)) {
  # Evidence states: E0 is the PSA itself; E2 (unless switched off) has its own PSA.
  theta_psa <- joint_draw_matrix(component)[psa_params$model_idx, , drop = FALSE]
  states <- list(E0 = list(component = component, rows = NULL))
  if (run_e2) {
    cat("\n=== 14 VoI: evidence state E2 (", nrow(evidence_rows), " external results; tau elicitation pending) ===\n",
        sep = "")
    e2 <- build_evidence_state(component, reference_population, evidence = evidence_rows,
                               tau_shared = voi14_settings$e2_tau_shared, tau_trial = voi14_settings$e2_tau_trial,
                               label = "E2", base_fingerprint = sampling_models$fingerprint,
                               evidence_md5 = attr(evidence_raw, "md5"), seed = analysis_seed)
    save_evidence_state(e2$state)
    cat("E2: conjugate update in", e2$update$iterations, "linearisations; importance ESS",
        round(e2$check$ess), "of", e2$check$n_draws, "draws\n")
    # E2's net benefit from a PSA on its own draws, with the economic draws and
    # seed of E0 (ruling (9) of 11 October 2026); reused when current, about
    # 20 minutes otherwise.
    cat("\n=== 14 VoI: PSA on the E2 draws ===\n")
    e2_psa <- load_or_run_evidence_psa(e2$state, l_params_base, param_distributions, strategies, time_horizon, cl,
                                       n_sim, analysis_seed, base_sampling_models = sampling_models,
                                       base_psa_obj = psa_obj, base_psa_params = psa_params)
    cat("E2 PSA:", if (e2_psa$from_cache) "read from" else "run and saved to", basename(e2_psa$path), "-",
        e2_psa$psa_obj$n_sim, "draws;", e2_psa$diagnostics$draws_shared_with_e0,
        "economic draws shared with E0\n")
    pairing_e2 <- policy_price_pairing(policy_space_from_psa(e2_psa$psa_obj, l_params_base))
    theta_e2 <- joint_draw_matrix(e2$state$joint)[e2_psa$psa_params$model_idx, , drop = FALSE]
    states$E2 <- list(component = e2$state$joint, pairing = pairing_e2, theta = theta_e2)
    # The check: the importance resample of the E0 PSA against E2's own PSA.
    e2_rows <- tryCatch(evidence_state_rows(e2$check$log_weights, psa_params$model_idx, analysis_seed,
                                            voi14_settings$e2_min_ess),
                        error = function(e) {
                          message("E2 check not computed: ", conditionMessage(e))
                          NULL
                        })
    e2_check <- NULL
    if (!is.null(e2_rows)) {
      e2_check <- evidence_state_comparison(pairing16, e2_rows$index, pairing_e2, voi14_settings$prices,
                                            voi14_settings$wtp)
      cat("E2 check (resampled E0 PSA, ESS", round(e2_rows$ess), "of", nrow(psa_params), "rows, against E2's PSA):\n")
      print(e2_check, row.names = FALSE)
    }
  }
  voi14_runtime <- voi14_runtime_estimate(pairing16, theta_psa, candidate_designs(), component, reference_population,
                                          params = l_params_base, settings = voi14_settings,
                                          n_states = length(states))
  cat("\n=== 14 VoI: steps 5-6 (EVSI by design, n and price; research value) ===\n")
  voi_research <- run_voi_research_value(pairing16, theta_psa, states, candidate_designs(), reference_population,
                                         params = l_params_base, cells = psa_obj$cells, eligible = eligible,
                                         settings = voi14_settings)
  voi_research$runtime_estimate <- voi14_runtime
  voi_research$evidence_selection <- evidence_selection
  voi_research$evidence_states <- if (run_e2) list(
    E2 = list(observations = e2$observations, check = e2$check$table, ess_draws = e2$check$ess,
              ess_psa_rows = if (is.null(e2_rows)) NA_real_ else e2_rows$ess, iterations = e2$update$iterations,
              fingerprint = e2$state$fingerprint, tau_shared = voi14_settings$e2_tau_shared,
              tau_trial = voi14_settings$e2_tau_trial,
              psa = list(generation_id = e2_psa$psa_obj$pair_metadata$generation_id,
                         fingerprint = e2_psa$identity$fingerprint, from_cache = e2_psa$from_cache,
                         diagnostics = e2_psa$diagnostics),
              resample_check = e2_check))
  save_voi_results(voi_research, research_identity, research_cache)
}
v_csv <- write_voi_research_csvs(voi_research)

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
cat("\n=== 14 VoI: summary (WTP EUR", format(WTP, big.mark = ","), "per QALY; evidence states",
    paste(voi14_settings$evidence_states, collapse = ", "), ") ===\n")
df_evpi_list <- voi_research$evpi[abs(voi_research$evpi$price - price_list) < 1e-9, ]
cat("Per-patient EVPI over the 16 policies at the list price (EUR ", format(price_list, big.mark = ","), "): ",
    paste(sprintf("%s %.2f", df_evpi_list$evidence_state, df_evpi_list$evpi), collapse = "; "), "\n", sep = "")
df_thr <- voi_research$threshold[voi_research$threshold$population == "norway" &
                                   voi_research$threshold$horizon == voi14_settings$primary_horizon, ]
for (i in seq_len(nrow(df_thr))) {
  r <- df_thr[i, ]
  cat(sprintf("Research threshold price, %s, nivolumab %s, %s cost band, %s end: %s\n", r$evidence_state,
              r$nivolumab, r$cost_band, r$cost_end, if (is.na(r$threshold_price)) "none on the price grid" else
                sprintf("EUR %s (%s, n* = %s)", format(round(r$threshold_price), big.mark = ","), r$best_design,
                        format(r$n_star, big.mark = ","))))
}
df_s6 <- voi_research$designs[startsWith(voi_research$designs$design, "S6_") &
                                voi_research$designs$accrual_scenario == "norway", ]
cat("Power-based n, rate ratio 0.5 (Norway accrual; engine, Schoenfeld in brackets):",
    paste(sprintf("%s %s (%s)", df_s6$design, df_s6$n_power, df_s6$n_power_schoenfeld), collapse = ", "), "\n")
cat("Outputs:", paste(basename(v_csv), collapse = ", "), "in", dirname(v_csv[1]), "\n")
