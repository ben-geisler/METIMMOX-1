# White-box regression tests for WB-004 / BB-SW / BB-EQ2.

source("scripts/R/functions/model_configs.R")
source("scripts/R/functions/calculate_outcomes.R")

make_params <- function() {
  list(
    u_np = 1,
    u_p = 1,
    c_drug_nivo = 0,
    c_drug_FLOX = 0,
    c_test_CT = 0,
    c_test_blood = 0,
    c_test_CRP = 17,
    c_test_NGS = 2500,
    c_test_biomarker = list(crp = 17, tmb_braf = 2500),
    c_other_visit = 0,
    c_other_baseline = 0,
    c_other_follow = 0,
    c_other_last = 0,
    l_nivo = 0,
    l_FLOX_exp = 0,
    l_FLOX_control = 0,
    l_CT = 0,
    l_blood = 0,
    l_visit = 0
  )
}

run_outcomes <- function(params, biomarker) {
  n_cycles <- 14
  calculate_outcomes(
    params = params,
    p_pf = rep(1, n_cycles),
    p_p = rep(0, n_cycles),
    p_d = rep(0, n_cycles),
    treatment_type = "control",
    biomarker = biomarker,
    v_dw_c = rep(1, n_cycles),
    v_dw_e = rep(1, n_cycles),
    cl = 1
  )
}

params <- make_params()

# The diagnostic charge comes from the mapping, not legacy scalar fields.
params$c_test_CRP <- 999
params$c_test_NGS <- 999
stopifnot(run_outcomes(params, "crp")$costs_total == 17)
stopifnot(run_outcomes(params, "tmb_braf")$costs_total == 2500)

# Swapping the data mapping swaps costs without changing calculation logic.
params$c_test_biomarker <- rev(params$c_test_biomarker)
names(params$c_test_biomarker) <- c("crp", "tmb_braf")
stopifnot(run_outcomes(params, "crp")$costs_total == 2500)
stopifnot(run_outcomes(params, "tmb_braf")$costs_total == 17)

# CRP diagnostic cost can be zeroed independently of routine blood monitoring.
params <- make_params()
params$c_test_biomarker$crp <- 0
params$c_test_blood <- 16
params$l_blood <- c(1, rep(0, 13))
stopifnot(run_outcomes(params, "crp")$costs_total == 16)
stopifnot(run_outcomes(params, NULL)$costs_total == 16)

# Unknown biomarkers fail explicitly rather than silently omitting the charge.
unknown_error <- tryCatch(
  {
    run_outcomes(params, "unknown")
    NULL
  },
  error = identity
)
stopifnot(inherits(unknown_error, "error"))
stopifnot(grepl("No diagnostic-test cost configured", conditionMessage(unknown_error)))

# PSA draws update both sampled scalars and the lookup consumed by the model.
source("scripts/R/functions/psa_functions.R")
captured_psa_params <- NULL
model_fun <- function(params, ...) {
  captured_psa_params <<- params
  data.frame(Strategy = "control", Cost = 0, Effect = 0)
}
psa_draws <- data.frame(sim = 1, c_test_CRP = 31, c_test_NGS = 4100)
invisible(capture.output(run_psa_analysis(
  psa_params = psa_draws,
  l_params_base = make_params(),
  param_distributions = list(c_test_CRP = list(), c_test_NGS = list()),
  strategies = "control",
  time_horizon = 1,
  cl = 1,
  n_sim = 1
)))
stopifnot(captured_psa_params$c_test_biomarker$crp == 31)
stopifnot(captured_psa_params$c_test_biomarker$tmb_braf == 4100)

cat("Biomarker diagnostic-test cost mapping tests passed.\n")
