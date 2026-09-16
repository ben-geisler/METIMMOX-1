# Adversarial review probes. Run from the repository root with Rscript.
# Reads local trial data/caches; writes no data or caches and prints aggregates only.
# These probes describe current behavior, including bugs, rather than asserting
# that the current behavior is correct. No full PSA or report render is run.

quiet_source <- function(path) {
  invisible(capture.output(suppressMessages(suppressWarnings(source(path)))) )
}
show <- function(label, x) {
  cat("\n### ", label, "\n", sep = "")
  print(x, row.names = FALSE)
}
attempt <- function(expr) {
  tryCatch(force(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
}

for (path in c("analysis/02_setup_and_global_variables.R",
               "analysis/03_biomarker_strategies.R",
               "analysis/04_parametric_survival_analysis.R",
               "analysis/05_basecase_input_parameters.R",
               "R/model_fun.R", "R/calculate_outcomes.R",
               "R/psa_functions.R", "R/cea_helpers.R")) quiet_source(path)

base <- model_fun(l_params_base)
show("Current base case", base)
show("Prevalence populations", data.frame(
  biomarker = get_biomarkers(),
  canonical = vapply(get_biomarkers(), function(b) l_params_base[[paste0("p_", b)]], numeric(1)),
  complete_case = vapply(get_biomarkers(), function(b) mean(data_complete[[b]] == 1), numeric(1))
))

# Same treatment, same schedules, no screening: increments must be zero for a
# comparison of policies applied to the same population.
null_params <- l_params_base
null_params$c_drug_nivo <- 0
null_params$l_FLOX_exp <- null_params$l_FLOX_control
null_params$c_test_CRP <- null_params$c_test_NGS <- 0
null_params <- sync_biomarker_test_costs(null_params)
for (b in get_biomarkers()) {
  pos <- data_complete[data_complete[[b]] == 1, , drop = FALSE]
  pos$Rx <- factor(levels(data_complete$Rx)[1], levels = levels(data_complete$Rx))
  null_params$p_os[[paste0(b, "_pos_OS")]] <- predict_pop_avg(models$best_fit$os, pos, time_points)
  null_params$p_pfs[[paste0(b, "_pos_PFS")]] <- predict_pop_avg(models$best_fit$pfs, pos, time_points)
}
null_result <- model_fun(null_params)
null_result$increment_cost <- null_result$Cost - null_result$Cost[1]
null_result$increment_qaly <- null_result$Effect - null_result$Effect[1]
show("Identical-treatment comparison at canonical prevalences", null_result)
matched_params <- null_params
for (b in get_biomarkers()) matched_params[[paste0("p_", b)]] <- mean(data_complete[[b]] == 1)
matched_result <- model_fun(matched_params)
show("Identical-treatment comparison with matched population weights", data.frame(
  strategy = matched_result$Strategy,
  increment_cost = matched_result$Cost - matched_result$Cost[1],
  increment_qaly = matched_result$Effect - matched_result$Effect[1]
))
null_params$p_tmb_braf <- null_params$p_tmb_braf * 1.2
varied_null <- model_fun(null_params)
show("Identical treatment after a supported +20% TMB/BRAF prevalence variation", data.frame(
  strategy = varied_null$Strategy,
  increment_qaly = varied_null$Effect - varied_null$Effect[1]
))

# Occupancy grid integration and short-horizon support.
unit_params <- l_params_base
unit_params$u_np <- unit_params$u_p <- 1
unit_params$dr_costs <- unit_params$dr_effects <- 0
show("No mortality and unit utilities", do.call(rbind, lapply(c(1L, 2L, 5L, 12L, 13L, 53L, 521L), function(n) {
  result <- attempt(calculate_outcomes(unit_params, rep(1, n), rep(0, n), rep(0, n),
                                       "standard", NULL, rep(1, n), rep(1, n), 1 / 52)$qalys_total)
  data.frame(grid_points = n, expected_years = (n - 1) / 52, actual = as.character(result))
})))

# The wrapper should honor the parameters it receives.
monthly_params <- l_params_base
monthly_params$cl <- 1 / 12
show("run_basecase ignores params$cl", data.frame(
  strategy = base$Strategy,
  wrapper = run_basecase(monthly_params, verbose = FALSE)$base_results$Effect,
  explicit_cl = model_fun(monthly_params, cl = monthly_params$cl)$Effect
))
short_params <- build_horizon_params(l_params_base, 260,
  list(os = models$best_fit$os, pfs = models$best_fit$pfs), strategies_df, data_complete)
show("run_basecase on a valid five-year parameter list", attempt(run_basecase(short_params, verbose = FALSE)))

show("PSA with no survival draw index returns deterministic results", list(
  equal_to_base = isTRUE(all.equal(model_fun(l_params_base, determpsa = "psa"), base)),
  fallback_used = attr(model_fun(l_params_base, determpsa = "psa"), "fallback_used")
))

# Non-finite model outcomes bypass the initial-pass fallback threshold.
real_model_fun <- model_fun
model_fun <- function(params, time_horizon, cl, determpsa, return_traces, sim_idx) {
  data.frame(Strategy = get_strategies(), Cost = if (sim_idx <= 3) NA_real_ else 100,
             Effect = 1)
}
nonfinite_result <- suppressWarnings(run_psa_analysis(
  psa_params = data.frame(sim = seq_len(100)), l_params_base = l_params_base,
  param_distributions = list(), strategies = get_strategies(),
  time_horizon = 520, cl = 1 / 52, n_sim = 100, fallback_threshold = 0.02
))
model_fun <- real_model_fun
show("Three invalid draws out of 100 under a 2% limit", nonfinite_result[c(
  "fallback_count", "fallback_rate", "dropped_count", "n_sim")])

# Range checks on the survival curves are absent.
bad_curve <- l_params_base
bad_curve$p_os$control_OS[2] <- 1.1
bad_result <- model_fun(bad_curve, return_traces = TRUE)
show("Invalid survival probability accepted", data.frame(
  week = 1,
  dead_occupancy = bad_result$traces$control$p_d[2],
  total_cost = bad_result$results$Cost[1],
  total_qaly = bad_result$results$Effect[1]
))

# Changing a canonical scalar test price does not change its duplicate lookup.
test_price <- l_params_base
test_price$c_test_NGS <- test_price$c_test_NGS + 1000
unsynced <- model_fun(test_price)
synced <- model_fun(sync_biomarker_test_costs(test_price))
show("Direct test-price change versus explicit lookup synchronization", data.frame(
  strategy = base$Strategy, unsynced_cost_change = unsynced$Cost - base$Cost,
  synced_cost_change = synced$Cost - base$Cost
))

# The expectation hash has no dependency on the model implementation.
sampling_models <- readRDS(sampling_cache_path())
parameter_config <- configure_parameter_distributions(l_params_base)
param_distributions <- parameter_config$distributions
hash_now <- function() psa_cache_fingerprint(sampling_models$fingerprint, l_params_base,
  param_distributions, strategies, n_sim, analysis_seed, time_horizon, cl)
before_hash <- hash_now()
real_outcomes <- calculate_outcomes
calculate_outcomes <- function(...) {
  result <- real_outcomes(...)
  result$qalys_total <- result$qalys_total * 2
  result
}
after_hash <- hash_now()
changed_model <- model_fun(l_params_base)
calculate_outcomes <- real_outcomes
show("Changed model logic has the same PSA fingerprint", list(
  identical_fingerprint = identical(before_hash$fingerprint, after_hash$fingerprint),
  control_qaly_before = base$Effect[1], control_qaly_after = changed_model$Effect[1]
))

# The two PSA files are not bound to each other by a draw-content fingerprint.
psa <- readRDS(psa_obj_path())
psa_pars <- readRDS(psa_params_path())
show("PSA parameter-cache metadata", list(
  seed = attr(psa_pars, "seed"), fingerprint = attr(psa_pars, "fingerprint"),
  model_indices_match = identical(psa$model_idx, psa_pars$model_idx)
))
tmp_cache_dir <- tempfile("metimmox-review-cache-")
dir.create(tmp_cache_dir)
reordered_pars <- psa_pars[rev(seq_len(nrow(psa_pars))), , drop = FALSE]
saveRDS(reordered_pars, psa_params_path(directory = tmp_cache_dir))
accepted <- attempt(load_psa_params_cache(directory = tmp_cache_dir, required = TRUE,
  verbose = FALSE, seed = analysis_seed, n_sim = psa$n_sim))
show("A reversed PSA parameter table passes the report loader", is.data.frame(accepted))
unlink(psa_params_path(directory = tmp_cache_dir))
unlink(tmp_cache_dir)

# The EVPPI cache stamps only the PSA fingerprint, although WTP changes EVPI.
quiet_source("R/evppi_functions.R")
show("EVPI changes with WTP while the PSA fingerprint is unchanged", do.call(rbind,
  lapply(c(51000, 100000, 150000, 500000), function(w) {
    data.frame(wtp = w, evpi = calculate_evpi_from_nmb(as.matrix(psa$effect) * w - as.matrix(psa$cost)))
  })))

# Constructed interactions are constant conditional on BOTH their parents.
quiet_source("R/dag_helpers.R")
quiet_source("R/assoc_tests.R")
quiet_source("R/tlr_landmark.R")
quiet_source("R/dag_association_tests.R")
coverage <- dag_ci_coverage(dag)
show("Deterministic-child CI statements classified as empirical tests", coverage[
  grepl("TLR _", coverage$statement) & grepl("Tx", coverage$statement),
  c("statement", "definitional")])
toy <- expand.grid(T_num = 0:1, crp = 0:1, replicate = seq_len(100))
toy$TxCRP <- toy$T_num * toy$crp
toy$tlr <- as.numeric(toy$replicate <= ifelse(toy$TxCRP == 1, 90, 10))
toy_test <- run_firth_logistic_test(toy, tlr ~ T_num + crp + TxCRP, "TxCRP",
  "TLR _||_ TxCRP | {CRP, T}", "Firth logistic")
show("Interaction test rejects a CI that holds by construction", data.frame(
  conditional_unique_TxCRP = max(tapply(toy$TxCRP, interaction(toy$T_num, toy$crp), function(x) length(unique(x)))),
  p_value = toy_test$p_value
))

# Production PSA subgroup calls use 'data', whereas the control and the base
# case use 'data_complete'. Endpoint missingness can change the prediction set.
parsed_sampling <- parse("analysis/06_sampling.R")
for (expr in parsed_sampling) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("generate_psa_population_averaged_predictions"))) eval(expr)
}
point_component <- list(samples = list(list(os = list(model = models$best_fit$os),
  pfs = list(model = models$best_fit$pfs), failed = FALSE)))
incomplete <- data_complete
incomplete$OSwk[which(incomplete$crp == 1)[1]] <- NA
filtered <- incomplete[complete.cases(incomplete[, c("Age", "sex", "Rx", "crp", "tmb_braf", "OSwk", "Death", "PFSwk", "Progression")]), ]
curves_full <- generate_psa_population_averaged_predictions(point_component, "crp", "os", 1, incomplete, time_points)
curves_filtered <- generate_psa_population_averaged_predictions(point_component, "crp", "os", 1, filtered, time_points)
show("Endpoint-missing patient is included only in PSA subgroup prediction", list(
  max_positive_curve_difference = max(abs(curves_full$positive - curves_filtered$positive))
))

# Static parse of every live R script and executable R chunks in reports and
# vignettes. Purl extracts code only; it does not execute report chunks.
r_files <- unlist(lapply(c("R", "analysis", "tests", "publish"), function(d) {
  list.files(d, "\\.R$", full.names = TRUE, recursive = TRUE)
}))
qmd_files <- unlist(lapply(c("reports", "outputs/vignettes"), function(d) {
  list.files(d, "\\.qmd$", full.names = TRUE, recursive = TRUE)
}))
parse_errors <- character()
for (f in r_files) {
  result <- attempt(parse(f))
  if (is.character(result)) parse_errors <- c(parse_errors, paste(f, result))
}
for (f in qmd_files) {
  scratch_r <- tempfile(fileext = ".R")
  result <- attempt({
    suppressWarnings(knitr::purl(f, output = scratch_r, documentation = 0, quiet = TRUE))
    parse(scratch_r)
  })
  unlink(scratch_r)
  if (is.character(result)) parse_errors <- c(parse_errors, paste(f, result))
}
show("Static syntax coverage", list(R_files = length(r_files), Qmd_files = length(qmd_files),
                                    errors = parse_errors))

cat("\nReview probes completed. This probe script did not change production code or caches.\n")
