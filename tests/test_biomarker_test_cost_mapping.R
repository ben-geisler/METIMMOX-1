# Scalar diagnostic prices must affect every public calculation path (#173).
source("R/model_configs.R")
source("R/calculate_outcomes.R")
source("R/model_fun.R")
source("R/cea_helpers.R")

params <- list(
  time_horizon = 520, cl = 1 / 52,  # model_fun() reads its grid from params (#172)
  dr_costs = 0.04, dr_effects = 0.04, u_np = 1, u_p = 1,
  c_drug_nivo = 0, c_drug_FLOX = 0, c_test_CT = 0, c_test_blood = 16,
  c_test_CRP = 17, c_test_NGS = 2500,
  # Deliberately stale legacy lookup must have no effect.
  c_test_biomarker = list(crp = 999, tmb_braf = 999),
  c_other_visit = 0, c_other_baseline = 0, c_other_follow = 0, c_other_last = 0,
  p_crp = 0.3, p_tmb_braf = 0.6)
for (nm in c("l_nivo", "l_FLOX_exp", "l_FLOX_control", "l_CT", "l_visit"))
  params[[nm]] <- rep(0, 521)
params$l_blood <- c(1, rep(0, 520))
for (endpoint in c("OS", "PFS")) {
  v_keys <- paste0(c("control", "crp_pos", "crp_neg", "tmb_braf_pos", "tmb_braf_neg"), "_", endpoint)
  params[[paste0("p_", tolower(endpoint))]] <- setNames(rep(list(rep(1, 521)), 5), v_keys)
}
df_base <- model_fun(params)
stopifnot(identical(df_base$Cost, c(16, 33, 2516)))
for (biomarker in get_biomarkers()) {
  changed <- params
  changed[[biomarker_cost_key(biomarker)[[1L]]]] <-
    changed[[biomarker_cost_key(biomarker)[[1L]]]] + 1000
  v_expected <- 1000 * as.numeric(df_base$Strategy == biomarker)
  # One-time diagnostic charge at t=0: discount factor is exactly 1.
  stopifnot(max(abs(model_fun(changed)$Cost - df_base$Cost - v_expected)) < 1e-10)
  stopifnot(max(abs(run_basecase(changed, verbose = FALSE)$base_results$Cost - df_base$Cost - v_expected)) < 1e-10)
}
without_lookup <- params
without_lookup$c_test_biomarker <- NULL
stopifnot(identical(model_fun(without_lookup), df_base))

# Direct calculate_outcomes() callers (including enriched analysis) use scalars.
#' Run calculate_outcomes() on a 14-point all-progression-free trace for one
#' biomarker with constant discount weight `dw`; returns its outcome list.
outcomes <- function(p, biomarker, dw = 1) calculate_outcomes(
  p, rep(1, 14), rep(0, 14), rep(0, 14), "standard", biomarker,
  rep(dw, 14), rep(1, 14), 1 / 52)
params$c_test_CRP <- 0
stopifnot(outcomes(params, "crp")$costs_total == 16,
          outcomes(params, NULL)$costs_total == 16)
changed <- params
changed$c_test_NGS <- changed$c_test_NGS + 1000
stopifnot(abs(outcomes(changed, "tmb_braf", 0.9)$costs_total -
                outcomes(params, "tmb_braf", 0.9)$costs_total - 900) < 1e-10)
stopifnot(inherits(tryCatch(outcomes(params, "unknown"), error = identity), "error"))

# PSA economic rows reach the real model without a caller-side sync.
source("R/psa_functions.R")
actual_model <- model_fun
#' Wrapper that forwards only `params` to the real model_fun(), so the PSA loop
#' runs the deterministic model on each economic draw.
model_fun <- function(params, ...) actual_model(params)
df_draws <- data.frame(sim = 1, c_test_CRP = 31, c_test_NGS = 4100)
invisible(capture.output(result <- run_psa_analysis(
  df_draws, params, list(c_test_CRP = list(), c_test_NGS = list()),
  get_strategies(), 520, 1 / 52, 1)))
stopifnot(identical(as.numeric(result$cost[1, ]), c(16, 47, 4116)))
cat("Diagnostic-price regression tests passed.\n")
