# ===============================================================================
# TEST: model_fun() input contract (issues #172, #49, #60)
# ===============================================================================
# 1. Time settings (#172): model_fun() takes time_horizon and cl from the
#    parameter list; an explicit argument must agree with it, and a missing
#    setting must be supplied explicitly. A valid five-year list therefore runs
#    through model_fun(params) and run_basecase(params).
# 2. Schedules (#172, #49): build_treatment_schedules() builds every schedule on
#    any horizon, truncating positions beyond it, and reproduces the script-05
#    schedules at 521 points; build_horizon_params() uses it for short horizons;
#    model_fun() rejects schedules whose length differs from the horizon.
# 3. Survival curves (#60): curves must be finite, within [0, 1], start at 1 and
#    be nonincreasing, in every mode and in partitioned_survival_states(); an
#    invalid PSA prediction is treated as a failed draw (fallback_used).
# 4. Scalars (#60): utilities, discount rates and prevalences must be finite
#    numeric scalars.
#
# Synthetic data only; no trial data or cache is read.
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_model_input_contract.R
# ===============================================================================

suppressPackageStartupMessages(library(dampack))
source("R/model_configs.R")
source("R/prediction_population.R")
source("R/prediction_functions.R")
source("R/model_fun.R")
source("R/calculate_outcomes.R")
source("R/parameter_distributions.R")
source("R/cea_helpers.R")

expect_error_matching <- function(expr, pattern) {
  err <- tryCatch({ force(expr); NULL }, error = identity)
  if (is.null(err)) stop("Expected an error matching '", pattern, "', got none.")
  if (!grepl(pattern, conditionMessage(err))) {
    stop("Expected an error matching '", pattern, "', got: ", conditionMessage(err))
  }
  invisible(err)
}

curve_keys <- function() {
  ctrl <- get_control_strategy()
  os <- paste0(ctrl, "_OS"); pfs <- paste0(ctrl, "_PFS")
  for (bm in get_biomarkers()) {
    os <- c(os, paste0(bm, c("_pos_OS", "_neg_OS")))
    pfs <- c(pfs, paste0(bm, c("_pos_PFS", "_neg_PFS")))
  }
  list(os = os, pfs = pfs)
}

make_params <- function(horizon, cl = 1 / 52) {
  t <- seq(0, horizon)
  keys <- curve_keys()
  p <- list(cl = cl, time_horizon = horizon, dr_costs = 0.04, dr_effects = 0.04,
            u_np = 0.73, u_p = 0.59, c_drug_nivo = 1000, c_drug_FLOX = 100,
            c_test_CT = 10, c_test_blood = 5, c_test_CRP = 16, c_test_NGS = 2518,
            c_other_visit = 10, c_other_baseline = 100, c_other_follow = 20,
            c_other_pp = 0, c_other_last = 1000, p_crp = 0.34, p_tmb_braf = 0.44)
  p$p_os <- setNames(lapply(seq_along(keys$os), function(i) exp(-0.004 * i * t)), keys$os)
  p$p_pfs <- setNames(lapply(seq_along(keys$pfs), function(i) exp(-0.010 * i * t)), keys$pfs)
  c(p, build_treatment_schedules(horizon + 1))
}

# ---------------------------------------------------------------------------
# 1. Time settings come from the parameter list (#172)
# ---------------------------------------------------------------------------
p260 <- make_params(260)
bare <- model_fun(p260)
stopifnot(nrow(bare) == 3L, all(is.finite(bare$Cost)), all(is.finite(bare$Effect)),
          identical(bare, model_fun(p260, time_horizon = 260, cl = 1 / 52)))
basecase <- run_basecase(p260, verbose = FALSE)
stopifnot(identical(basecase$base_results, bare))

expect_error_matching(model_fun(p260, time_horizon = 520), "conflicts with params\\$time_horizon")
expect_error_matching(model_fun(p260, cl = 1 / 12), "conflicts with params\\$cl")

no_time <- p260; no_time$time_horizon <- NULL; no_time$cl <- NULL
expect_error_matching(model_fun(no_time), "needs 'time_horizon'")
expect_error_matching(model_fun(no_time, time_horizon = 260), "needs 'cl'")
stopifnot(identical(model_fun(no_time, time_horizon = 260, cl = 1 / 52), bare))

# A monthly cycle length in params is honoured without an explicit argument.
monthly <- p260; monthly$cl <- 1 / 12
stopifnot(identical(model_fun(monthly), model_fun(monthly, cl = 1 / 12)),
          all(model_fun(monthly)$Effect > bare$Effect))

bad_horizon <- p260; bad_horizon$time_horizon <- 260.5
expect_error_matching(model_fun(bad_horizon), "time_horizon must be a non-negative whole number")
bad_cl <- p260; bad_cl$cl <- -1
expect_error_matching(model_fun(bad_cl), "cl must be a positive finite")

# ---------------------------------------------------------------------------
# 2. Schedules on every horizon (#172, #49)
# ---------------------------------------------------------------------------
schedule_names <- c("l_nivo", "l_FLOX_exp", "l_FLOX_control", "l_CT", "l_blood", "l_visit")
for (n in 1:60) {
  s <- build_treatment_schedules(n)
  stopifnot(identical(names(s), schedule_names),
            all(vapply(s, length, integer(1)) == n),
            all(unlist(s) %in% c(0, 1)), s$l_visit[1] == 1, s$l_CT[1] == 1, s$l_blood[1] == 1)
}
expect_error_matching(build_treatment_schedules(0), "positive whole number")

# At the base-case 521 points the builder reproduces the script-05 construction.
legacy <- local({
  n <- 521
  l_nivo <- rep(0, n); l_nivo[c(5, 7, 13, 15, 29, 31, 37, 39)] <- 1
  l_FLOX_exp <- rep(0, n); l_FLOX_exp[c(1, 3, 9, 11, 25, 27, 33, 35)] <- 1
  l_FLOX_control <- rep(0, n)
  l_FLOX_control[c(1, 3, 5, 7, 9, 11, 13, 15, 25, 27, 29, 31, 33, 35, 37, 39)] <- 1
  l_CT <- rep(0, n); l_CT[1] <- 1; l_CT[seq(13, n, by = 12)] <- 1
  l_blood <- rep(0, n); l_blood[1] <- 1; l_blood[seq(5, n, by = 4)] <- 1
  l_visit <- rep(0, n); l_visit[1] <- 1
  l_visit[which(l_nivo == 1 | l_FLOX_exp == 1 | l_FLOX_control == 1)] <- 1
  l_visit[l_CT == 1] <- 1  # surveillance visit with every CT (issue #164)
  list(l_nivo = l_nivo, l_FLOX_exp = l_FLOX_exp, l_FLOX_control = l_FLOX_control,
       l_CT = l_CT, l_blood = l_blood, l_visit = l_visit)
})
stopifnot(identical(build_treatment_schedules(521), legacy))
# Issue #164: the CT rule adds visits only after the last administration
# (position 39); every in-treatment CT already falls on an administration.
s521 <- build_treatment_schedules(521)
added <- which(s521$l_visit == 1 & !(s521$l_nivo == 1 | s521$l_FLOX_exp == 1 |
                                        s521$l_FLOX_control == 1))
stopifnot(identical(added, seq(49L, 521L, by = 12L)),
          all(s521$l_visit[s521$l_CT == 1] == 1))

# A 30-point horizon keeps 30 points (script 05 used to lengthen it to 39).
stopifnot(all(vapply(build_treatment_schedules(30), length, integer(1)) == 30L))

# The model runs on horizons shorter than every schedule rule.
for (h in c(0, 1, 3, 4, 11, 12, 30, 38)) {
  r <- model_fun(make_params(h))
  stopifnot(all(is.finite(r$Cost)), all(is.finite(r$Effect)))
}

# Schedules of the wrong length or non-binary values are rejected.
long_schedule <- p260; long_schedule$l_nivo <- c(long_schedule$l_nivo, 0)
expect_error_matching(model_fun(long_schedule), "Schedule 'l_nivo' has length 262")
short_schedule <- p260; short_schedule$l_CT <- short_schedule$l_CT[1:39]
expect_error_matching(model_fun(short_schedule), "Schedule 'l_CT' has length 39")
fractional <- p260; fractional$l_blood[2] <- 0.5
expect_error_matching(model_fun(fractional), "Schedule 'l_blood' must contain only 0 and 1")

# build_horizon_params() on short horizons (mock survival models, as in
# test_canonical_prevalence.R).
predict.mock_survival_model <- function(object, newdata, type, times, ...) {
  rate <- (0.01 + newdata$Age / 3000 + as.numeric(as.character(newdata$crp)) / 100 +
             as.numeric(as.character(newdata$tmb_braf)) / 50) * object$rate
  structure(list(.pred = lapply(rate, function(r)
    data.frame(.pred_survival = exp(-r * times)))), class = "data.frame",
    row.names = seq_along(rate))
}
mock_data <- data.frame(ID = 1:8, Age = seq(45, 80, 5), sex = factor(rep(0:1, 4)),
  Rx = factor(rep(c("control", "experimental"), 4)),
  crp = c(0, 0, 0, 0, 1, 1, 1, 1), tmb_braf = c(0, 0, 0, 1, 0, 1, 1, 1),
  OSwk = 100, Death = 1, PFSwk = 50, Progression = 1)
mock_complete <- economic_prediction_population(mock_data)
mock_models <- list(os = structure(list(rate = 1), class = "mock_survival_model"),
                    pfs = structure(list(rate = 2), class = "mock_survival_model"))
mock_strategies <- data.frame(id = get_strategies(), prevalence = c(1, 0.5, 0.5))
mock_base <- set_population_predictions(make_params(52), generate_population_averaged_predictions(
  mock_models, mock_strategies, mock_complete, 0:52, quiet = TRUE))
for (h in c(2, 4, 12, 30, 104)) {
  hp <- build_horizon_params(mock_base, h, mock_models, mock_strategies, mock_complete)
  stopifnot(hp$time_horizon == h,
            all(vapply(hp[schedule_names], length, integer(1)) == h + 1),
            sum(hp$l_nivo) == sum(c(5, 7, 13, 15, 29, 31, 37, 39) <= h + 1))
  r <- model_fun(hp)
  stopifnot(all(is.finite(r$Cost)), all(is.finite(r$Effect)))
}

# ---------------------------------------------------------------------------
# 3. Survival-curve contract (#60)
# ---------------------------------------------------------------------------
ctrl_os <- paste0(get_control_strategy(), "_OS")
above_one <- p260; above_one$p_os[[ctrl_os]][2] <- 1.1
expect_error_matching(model_fun(above_one), "leaves \\[0, 1\\]")
expect_error_matching(model_fun(above_one, determpsa = "curves"), "leaves \\[0, 1\\]")

late_start <- p260
late_start$p_os <- lapply(late_start$p_os, function(x) 0.8 * x)
late_start$p_pfs <- lapply(late_start$p_pfs, function(x) 0.8 * x)
expect_error_matching(model_fun(late_start), "starts at 0.8")

rising <- p260; rising$p_os[[ctrl_os]][3:4] <- c(0.5, 0.9)
expect_error_matching(model_fun(rising), "increases at")

missing_value <- p260; missing_value$p_os[[ctrl_os]][10] <- NA
expect_error_matching(model_fun(missing_value), "non-finite")

negative <- p260; negative$p_pfs[[1]][261] <- -0.01
expect_error_matching(model_fun(negative), "leaves \\[0, 1\\]")

# partitioned_survival_states() checks the same contract for direct callers
# (the enriched-population path), labelled or not.
expect_error_matching(partitioned_survival_states(c(0.8, 0.7), c(0.8, 0.6)), "starts at 0.8")
expect_error_matching(partitioned_survival_states(c(1, 1.1), c(1, 0.5), "enriched"),
                      "Survival curve 'enriched OS' leaves")
states <- partitioned_survival_states(c(1, 0.9, 0.8), c(1, 0.7, 0.5))
stopifnot(isTRUE(all.equal(states$p_pf + states$p_p + states$p_d, rep(1, 3))))
# Numerical noise within the tolerance is accepted.
stopifnot(is.null(survival_curve_problem(c(1 - 1e-13, 0.5, 0.5 + 1e-13, 0))))

# Invalid PSA predictions are failed draws, not accepted curves: model_fun()
# flags fallback_used and keeps the base-case curves, so the PSA loop replaces
# the draw. The sampled-prediction machinery is mocked.
n_samples <- 1L
sampling_models <- list(joint = list())
get_joint_sampling_models <- function(x) x$joint
sampled_survival_models <- function(component, idx) list(failed = FALSE)
model_population_weights <- function(params, population) NULL
psa_mock_curve <- NULL
generate_psa_population_averaged_predictions <- function(sampling_model_list, biomarker_name = NULL,
                                                         outcome, sample_idx, data_original,
                                                         weights, time_points) {
  curve <- psa_mock_curve(outcome, time_points)
  if (is.null(biomarker_name)) curve else list(positive = curve, negative = curve)
}
psa_params <- p260; psa_params$prediction_population <- mock_complete

psa_mock_curve <- function(outcome, t) exp(-(if (outcome == "os") 0.003 else 0.006) * t)
valid_draw <- model_fun(psa_params, determpsa = "psa", sim_idx = 1)
stopifnot(!isTRUE(attr(valid_draw, "fallback_used")),
          !isTRUE(all.equal(valid_draw$Effect, bare$Effect)))

psa_mock_curve <- function(outcome, t) {
  curve <- exp(-(if (outcome == "os") 0.003 else 0.006) * t)
  if (outcome == "os") curve[2] <- 1.1
  curve
}
invalid_draw <- suppressWarnings(model_fun(psa_params, determpsa = "psa", sim_idx = 1))
stopifnot(isTRUE(attr(invalid_draw, "fallback_used")),
          isTRUE(all.equal(as.data.frame(invalid_draw), as.data.frame(bare),
                           check.attributes = FALSE)))

# ---------------------------------------------------------------------------
# 4. Scalar parameters (#60)
# ---------------------------------------------------------------------------
for (case in list(list("u_np", NA_real_), list("u_p", Inf), list("u_np", c(0.7, 0.8)),
                  list("u_p", "0.5"), list("dr_costs", NaN), list("dr_effects", Inf),
                  list("dr_costs", c(0.04, 0.04)), list("p_crp", NA_real_))) {
  p <- p260; p[[case[[1]]]] <- case[[2]]
  expect_error_matching(model_fun(p), paste0("'", case[[1]], "' must be a finite numeric scalar"))
}

cat("PASS: time settings from params, schedules on every horizon, survival-curve",
    "and scalar contracts (issues #172, #49, #60).\n")
