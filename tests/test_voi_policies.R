# ===============================================================================
# TEST: cell-level PSA outcomes, life-years and the 16-policy space
#       (issues #200 and #197, paper 3 step 1b)
# ===============================================================================
# Synthetic (mock exponential survival models, 12 patients in all four joint
# biomarker cells; no trial data or cache):
#   1. Cell-to-strategy identity: the strategy curves built from the cell curves
#      (strategy_curves_from_cells()) equal the patient-level standardisation of
#      generate_population_averaged_predictions() and the former PSA helper of
#      script 06, for joint weights that differ from the observed cell shares,
#      for weights that vary within a cell and for a biomarker subgroup of
#      zero mass.
#   2. model_fun() in PSA mode (cell path) reproduces model_fun() on the
#      patient-level curves (curves mode), costs, QALYs and life-years.
#   3. Administration slope: the discounted nivolumab administrations equal
#      the derivative of cost in the nivolumab price, per cell and per strategy,
#      including the third and fourth treatment sequences (issue #208), whose
#      expected cost per patient is the same in every cell of an arm.
#   4. Life-years: trapezoidal area under OS; with both utilities at 1 the
#      discounted QALYs equal the discounted life-years; at most the horizon.
#   5. Policy space: 16 assignments, the six named labels, the tests each
#      policy needs and their costs (treat-all and nobody need no test); the
#      three paper 2 strategies reproduce the strategy outcomes exactly when no
#      cell curve crosses, and the clamp gap is reported (and non-zero) when
#      one does. The fixture costs second-line therapy at incident
#      progression (issue #209): its cell and arm component equals the hand
#      calculation and the cell combination reproduces the strategies'.
#   6. End to end: build_psa_obj() stores life-years and the cell component,
#      bind_psa_pair()/validate_psa_pair() accept it, the policy space
#      reproduces the PSA, the 3-policy EVPI equals the PSA EVPI, the price
#      pairing of two PSAs holds for life-years, the policy costs at a second
#      price equal a PSA run at that price, and the step 1b driver writes its
#      three CSVs.
# Live (skip until the regenerated caches exist):
#   7. The three paper 2 policies reproduce the current PSA per draw (gap zero
#      in draws without a cell crossing).
#   8. From the list-price scenario PSA's cells, EVPI over the three paper 2
#      strategies at the biosimilar price and EUR 51,000 equals the biosimilar
#      scenario PSA's EVPI and Table S7 (EUR 576.72 at the time of writing).
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_policies.R
# ===============================================================================

suppressPackageStartupMessages(library(dampack))
source("R/model_configs.R")
source("R/prediction_functions.R")
source("R/model_fun.R")
source("R/calculate_outcomes.R")
source("R/parameter_distributions.R")
source("R/psa_functions.R")
source("R/cache_paths.R")
source("R/cache_provenance.R")
source("R/evppi_functions.R")
source("R/voi_price_axis.R")
source("R/voi_policies.R")
source("R/later_sequences.R")

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Mock predict() method: exponential survival with a per-patient rate that
#' rises with age, CRP and TMB/BRAF; under the experimental arm the rate is
#' multiplied by `rx_mult[cell]` (cell key "<crp><tmb_braf>"). Returns
#' flexsurv-style nested `.pred` survival frames.
predict.mock_survival_model <- function(object, newdata, type, times, ...) {
  v_crp <- as.numeric(as.character(newdata$crp))
  v_tmb <- as.numeric(as.character(newdata$tmb_braf))
  v_rate <- (0.004 + newdata$Age / 20000 + v_crp / 500 + v_tmb / 400) * object$rate
  v_exp <- newdata$Rx == levels(newdata$Rx)[2]
  v_rate[v_exp] <- v_rate[v_exp] * object$rx_mult[paste0(v_crp, v_tmb)[v_exp]]
  structure(list(.pred = lapply(v_rate, function(r)
    data.frame(.pred_survival = exp(-r * times)))), class = "data.frame",
    row.names = seq_along(v_rate))
}

#' Mock OS/PFS model pair; `pfs_mult` scales the PFS rate in the experimental
#' arm per cell (a value below 1 / 2 makes that cell's PFS cross its OS).
mock_models <- function(os_rate = 1, rx_mult = c("00" = 1.1, "01" = 0.9, "10" = 0.7, "11" = 0.8),
                        pfs_mult = rx_mult) {
  list(os = structure(list(rate = os_rate, rx_mult = rx_mult), class = "mock_survival_model"),
       pfs = structure(list(rate = 2 * os_rate, rx_mult = pfs_mult), class = "mock_survival_model"))
}

horizon <- 104
time_points <- seq(0, horizon)
df_population <- data.frame(
  ID = 1:12, Age = c(45, 52, 60, 71, 48, 66, 74, 55, 69, 50, 63, 77),
  sex = factor(c(0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0)),
  Rx = factor(rep(c("control", "experimental"), 6), levels = c("control", "experimental")),
  crp = c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
  tmb_braf = c(0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1),
  OSwk = 100, Death = 1, PFSwk = 50, Progression = 1)
df_population <- economic_prediction_population(df_population)
joint <- joint_biomarker_population(df_population)
stopifnot(identical(rownames(joint$cells), c("00", "01", "10", "11")),
          identical(unname(joint$counts), c(4L, 3L, 2L, 3L)))

# Second-line inputs of the fixture (issue #209), roughly the base case.
SECOND_LINE_FIXTURE <- list(
  p_pfs_progression = 0.98, p_sl_treated = 0.8, p_sl_folfiri = 0.5, p_sl_ras_braf_wt = 0.362,
  p_sl_antiegfr_if_wt = 0.5, p_sl_bevacizumab_if_no_antiegfr = 0.25, n_sl_cycles_none = 9.8,
  n_sl_cycles_bevacizumab = 12, n_sl_cycles_panitumumab = 10.8, n_sl_contacts_FOLFIRI = 2,
  n_sl_contacts_FLIRI = 2, n_sl_blood_tests = 1, c_sl_drug_FOLFIRI = 822, c_sl_drug_FLIRI = 623,
  c_sl_drug_bevacizumab = 887, c_sl_drug_panitumumab = 1729)
SECOND_LINE_FIXTURE_SE <- c(p_pfs_progression = 0.0196, p_sl_treated = 0.09, p_sl_folfiri = 0.15,
  p_sl_ras_braf_wt = 0.03, p_sl_antiegfr_if_wt = 0.14, p_sl_bevacizumab_if_no_antiegfr = 0.13,
  n_sl_cycles_none = 0.6, n_sl_cycles_bevacizumab = 1, n_sl_cycles_panitumumab = 0.6)

#' Base parameter list on the 104-week grid with the production schedules.
make_params <- function(models = mock_models(), joint_mass = c(0.3, 0.2, 0.15, 0.35)) {
  p <- list(time_horizon = horizon, cl = 1 / 52, dr_costs = 0.04, dr_effects = 0.04,
            u_np = 0.73, u_p = 0.59, u_decrement = 0.14, c_drug_nivo = 13923, c_drug_FLOX = 100,
            c_test_CT = 150, c_test_blood = 20, c_test_CRP = 16, c_test_NGS = 2518,
            c_other_visit = 213.46, c_other_baseline = 213.46, c_other_follow = 213.46,
            c_other_pp = 0, c_other_last = 5000, p_crp = 0.5, p_tmb_braf = 0.5)
  p <- c(p, build_treatment_schedules(horizon + 1))
  # Second-line therapy at incident progression (issue #209), so every
  # identity below also covers the second-line cost.
  p <- c(p, SECOND_LINE_FIXTURE)
  p$psa_se <- SECOND_LINE_FIXTURE_SE
  p$joint_counts <- joint$counts
  for (k in seq_along(joint_mass)) p[[paste0("p_joint_", rownames(joint$cells)[k])]] <- joint_mass[k]
  p$p_crp <- joint_mass[3] + joint_mass[4]
  p$p_tmb_braf <- joint_mass[2] + joint_mass[4]
  # Third and fourth treatment sequences (issue #208): synthetic starters in
  # both arms, all within the 104-week grid except one fourth sequence that is
  # truncated at the horizon.
  later <- later_sequence_inputs(df_population, data.frame(
    ID = df_population$ID, Seq3 = c(1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0),
    Seq3wk = c(40, 50, NA, 45.6, NA, 60, 55, 70, NA, NA, 48, NA),
    Seq4 = c(0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0),
    Seq4wk = c(NA, 70, NA, NA, NA, 95, NA, NA, NA, NA, 72, NA)))
  p$later_sequences <- later$later_sequences
  p[names(later$shares)] <- later$shares
  v_weights <- model_population_weights(p, df_population)
  preds <- generate_population_averaged_predictions(models, NULL, df_population, time_points,
                                                    quiet = TRUE, weights = v_weights)
  set_population_predictions(p, preds)
}

# ---------------------------------------------------------------------------
# 1. Cell-to-strategy identity
# ---------------------------------------------------------------------------
#' Largest absolute difference between the strategy curves of two curve sets.
curve_set_difference <- function(a, b) {
  v_diff <- c(abs(a$control$os - b$control$os), abs(a$control$pfs - b$control$pfs))
  for (bm in get_biomarkers()) {
    for (side in c("biomarker_positive", "biomarker_negative")) for (outcome in c("os", "pfs"))
      v_diff <- c(v_diff, abs(a[[bm]][[side]][[outcome]] - b[[bm]][[side]][[outcome]]))
    v_diff <- c(v_diff, abs(a[[bm]]$os - b[[bm]]$os), abs(a[[bm]]$prevalence - b[[bm]]$prevalence))
  }
  max(v_diff)
}
models <- mock_models()
params <- make_params(models)
v_weights <- params$population_weights
preds_patient <- generate_population_averaged_predictions(models, NULL, df_population, time_points,
                                                          quiet = TRUE, weights = v_weights)
cell_curves <- predict_cell_curves(models$os, models$pfs, df_population, time_points, v_weights)
stopifnot(identical(colnames(cell_curves$os), c("00_control", "00_experimental", "01_control",
  "01_experimental", "10_control", "10_experimental", "11_control", "11_experimental")),
  max(abs(cell_curves$mass - c(0.3, 0.2, 0.15, 0.35))) < 1e-12)
preds_cells <- strategy_curves_from_cells(cell_curves, cell_masses(df_population, v_weights))
diff_identity <- curve_set_difference(preds_patient, preds_cells)
stopifnot(diff_identity < 1e-14)

# Weights that vary within a cell, and a TMB/BRAF-positive subgroup of zero mass.
v_uneven <- seq(1, 2.1, by = 0.1)
v_zero_tmb <- ifelse(df_population$tmb_braf == 1, 0, 1)
for (v_w in list(v_uneven / sum(v_uneven), v_zero_tmb / sum(v_zero_tmb))) {
  a <- generate_population_averaged_predictions(models, NULL, df_population, time_points,
                                                quiet = TRUE, weights = v_w)
  b <- strategy_curves_from_cells(predict_cell_curves(models$os, models$pfs, df_population,
                                                      time_points, v_w),
                                  cell_masses(df_population, v_w))
  stopifnot(curve_set_difference(a, b) < 1e-14)
}

# The former PSA helper of script 06 (subgroup predictions) gives the same curves.
for (expr in parse("analysis/06_sampling.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("generate_psa_population_averaged_predictions"))) eval(expr)
}
sampling_fixture <- list(samples = list(list(os = list(model = models$os),
                                             pfs = list(model = models$pfs), failed = FALSE)))
diff_old_psa <- max(sapply(c("os", "pfs"), function(outcome) {
  v_control <- generate_psa_population_averaged_predictions(sampling_fixture, NULL, outcome, 1,
                                                            df_population, time_points, v_weights)
  v_d <- abs(v_control - preds_cells$control[[outcome]])
  for (bm in get_biomarkers()) {
    sub <- generate_psa_population_averaged_predictions(sampling_fixture, bm, outcome, 1,
                                                        df_population, time_points, v_weights)
    v_d <- c(v_d, abs(sub$positive - preds_cells[[bm]]$biomarker_positive[[outcome]]),
             abs(sub$negative - preds_cells[[bm]]$biomarker_negative[[outcome]]))
  }
  max(v_d)
}))
stopifnot(diff_old_psa < 1e-14)
cat(sprintf("PASS: strategy curves are joint-weight combinations of the 8 cell curves (max |diff| %.1e; former PSA helper %.1e)\n",
            diff_identity, diff_old_psa))

# ---------------------------------------------------------------------------
# 2. model_fun() PSA (cell path) = model_fun() on the patient-level curves
# ---------------------------------------------------------------------------
n_samples <- 1L
sampling_models <- list(joint = sampling_fixture)
df_curves <- model_fun(params, determpsa = "curves")
df_psa <- model_fun(params, determpsa = "psa", sim_idx = 1, return_cells = TRUE)
stopifnot(!isTRUE(attr(df_psa, "fallback_used")),
          max(abs(df_psa$Cost - df_curves$Cost)) < 1e-8,
          max(abs(df_psa$Effect - df_curves$Effect)) < 1e-12,
          max(abs(as.matrix(df_psa[, c("LY_disc", "LY_undisc")]) -
                  as.matrix(df_curves[, c("LY_disc", "LY_undisc")]))) < 1e-12,
          is.null(attr(df_psa, "life_years")),
          is.null(attr(df_curves, "cells")), !is.null(attr(df_psa, "cells")))
expect_error_matching(model_fun(params, return_cells = TRUE), "only available with determpsa = 'psa'")
cat("PASS: model_fun() PSA mode on the cell path reproduces the patient-level curves\n")

# ---------------------------------------------------------------------------
# 3. Administration slope
# ---------------------------------------------------------------------------
cell_results <- attr(df_psa, "cells")
params_up <- params
params_up$c_drug_nivo <- params$c_drug_nivo + 1000
cell_results_up <- attr(model_fun(params_up, determpsa = "psa", sim_idx = 1, return_cells = TRUE),
                        "cells")
m_slope <- (cell_results_up$values[, , "cost"] - cell_results$values[, , "cost"]) / 1000
v_dw <- discount_weights(params, horizon + 1)
later_exp <- later_sequence_resource_use(params, "experimental", horizon + 1)
later_ctl <- later_sequence_resource_use(params, "control", horizon + 1)
n_later_nivo <- sum(later_exp$nivolumab * v_dw$cost)
stopifnot(max(abs(m_slope - cell_results$values[, , "nivo_admin"])) < 1e-9,
          all(cell_results$values[, "control", "nivo_admin"] == 0),
          all(cell_results$values[, "experimental", "nivo_admin"] > n_later_nivo),
          n_later_nivo > 0,
          all(cell_results$values[, "experimental", "nivo_admin"] <=
                sum(params$l_nivo) + 2 * length(params$later_sequences$regimen$experimental$nivolumab)))
states <- partitioned_survival_states(cell_curves$os[, "11_experimental"],
                                      cell_curves$pfs[, "11_experimental"])
outcome <- calculate_outcomes(params, states$p_pf, states$p_p, states$p_d, "experimental", NULL,
                              v_dw$cost, v_dw$effect, 1 / 52)
stopifnot(abs(outcome$nivo_administrations_total -
                sum(params$l_nivo * states$p_pf * v_dw$cost) - n_later_nivo) < 1e-12)
# The later sequences are expected costs per patient assigned to the arm,
# whatever the health state: removing them lowers every cell of an arm by the
# same amount and leaves QALYs and life-years unchanged.
params_none <- params
params_none$later_sequences <- NULL
cell_none <- attr(model_fun(params_none, determpsa = "psa", sim_idx = 1, return_cells = TRUE), "cells")
cost_later <- function(u) sum((u$FLOX * params$c_drug_FLOX + u$nivolumab * params$c_drug_nivo +
                                 u$contacts * params$c_other_visit + u$blood * params$c_test_blood) * v_dw$cost)
m_cost_gap <- cell_results$values[, , "cost"] - cell_none$values[, , "cost"]
stopifnot(max(abs(m_cost_gap[, "experimental"] - cost_later(later_exp))) < 1e-8,
          max(abs(m_cost_gap[, "control"] - cost_later(later_ctl))) < 1e-8,
          cost_later(later_ctl) > 0,
          identical(cell_results$values[, , c("qalys", "ly_disc", "ly_undisc")],
                    cell_none$values[, , c("qalys", "ly_disc", "ly_undisc")]),
          max(abs(cell_results$values[, "experimental", "nivo_admin"] -
                    cell_none$values[, "experimental", "nivo_admin"] - n_later_nivo)) < 1e-12)
# Sequence shares: p_seq4 = p_seq3 x p_seq4_given_seq3 places that mass at the
# fourth-sequence start weeks; a start beyond the grid drops its doses.
v_start4 <- params$later_sequences$start_weeks
v_start4 <- v_start4[v_start4$arm == "experimental" & v_start4$sequence == 4L, ]
share4 <- params$p_seq3_experimental * params$p_seq4_given_seq3_experimental
v_nivo_weeks <- params$later_sequences$regimen$experimental$nivolumab
n_in_grid <- sum(sapply(seq_len(nrow(v_start4)), function(i)
  v_start4$n[i] * sum(v_start4$week[i] + v_nivo_weeks <= horizon))) / sum(v_start4$n)
stopifnot(share4 > 0, n_in_grid < length(v_nivo_weeks))
cat("PASS: discounted nivolumab administrations are the cost slope in the price, later sequences included\n")

# ---------------------------------------------------------------------------
# 4. Life-years
# ---------------------------------------------------------------------------
params_u1 <- params
params_u1$u_np <- params_u1$u_p <- 1
df_u1 <- model_fun(params_u1)
# Life-years are result columns of model_fun() (issue #197), not an attribute.
# The deterministic fixture has no OS/PFS crossing (det mode would stop), so
# the QALYs at utility 1 equal the discounted life-years.
df_ly <- df_u1
v_control_os <- params$p_os$control_OS
v_trapezoid <- sum(v_control_os * (1 / 52) * interval_weights(horizon + 1))
stopifnot(max(abs(df_u1$Effect - df_ly$LY_disc)) < 1e-12,
          all(df_ly$LY_undisc >= df_ly$LY_disc), all(df_ly$LY_undisc <= horizon / 52),
          abs(df_ly$LY_undisc[df_ly$Strategy == "control"] - v_trapezoid) < 1e-12,
          identical(model_fun(params)$Strategy, get_strategies()),
          is.null(attr(model_fun(params), "life_years")))
cell_u1 <- attr(model_fun(params_u1, determpsa = "psa", sim_idx = 1, return_cells = TRUE), "cells")
stopifnot(max(abs(cell_u1$values[, , "qalys"] - cell_u1$values[, , "ly_disc"])) < 1e-12)
cat("PASS: life-years are the trapezoidal area under OS and equal QALYs at utility 1\n")

# ---------------------------------------------------------------------------
# 5. Policy space, test costs, reproduction of the three strategies
# ---------------------------------------------------------------------------
df_policies <- policy_space_definition(joint$cells)
stopifnot(nrow(df_policies) == 16L, !anyDuplicated(df_policies$pattern),
          identical(df_policies$pattern[c(1, 16)], c("0000", "1111")))
v_named <- setNames(df_policies$label, df_policies$pattern)[c("0000", "1111", "0011", "0101", "0111", "0001")]
stopifnot(identical(unname(v_named), c("Standard of Care", "Nivolumab for everyone", "CRP-guided",
                                       "TMB/BRAF-guided", "Either biomarker positive",
                                       "Both biomarkers positive")),
          sum(df_policies$named) == 6L,
          identical(df_policies$paper2_strategy[match(c("0000", "0011", "0101"), df_policies$pattern)],
                    c("control", "crp", "tmb_braf")))
v_tests <- policy_test_costs(df_policies, params)
names(v_tests) <- df_policies$pattern
stopifnot(v_tests[["0000"]] == 0, v_tests[["1111"]] == 0,
          v_tests[["0011"]] == 16, v_tests[["1100"]] == 16,
          v_tests[["0101"]] == 2518, v_tests[["1010"]] == 2518,
          v_tests[["0111"]] == 16 + 2518, v_tests[["0001"]] == 16 + 2518,
          v_tests[["0110"]] == 16 + 2518, v_tests[["0010"]] == 16 + 2518)

df_policy <- policy_outcomes(cell_results$values, cell_results$mass, params)
v_rows <- match(c("nivo_0000", "nivo_0011", "nivo_0101"), df_policy$policy)
v_strat_rows <- match(c("control", "crp", "tmb_braf"), df_psa$Strategy)
df_ly_psa <- df_psa
gap_cost <- max(abs(df_policy$cost[v_rows] - df_psa$Cost[v_strat_rows]))
gap_qaly <- max(abs(df_policy$qalys[v_rows] - df_psa$Effect[v_strat_rows]))
stopifnot(all(cell_results$violations == 0), gap_cost < 1e-8, gap_qaly < 1e-12,
          max(abs(df_policy$ly_disc[v_rows] - df_ly_psa$LY_disc[v_strat_rows])) < 1e-12,
          max(abs(df_policy$ly_undisc[v_rows] - df_ly_psa$LY_undisc[v_strat_rows])) < 1e-12,
          df_policy$test_cost[df_policy$pattern == "1111"] == 0,
          df_policy$nivo_admin[df_policy$pattern == "0000"] == 0)
# Treat-all = every cell under the experimental arm, no test.
v_mass <- cell_results$mass
stopifnot(abs(df_policy$cost[df_policy$pattern == "1111"] -
                sum(v_mass * cell_results$values[, "experimental", "cost"])) < 1e-8)
cat(sprintf("PASS: 16 policies, test costs, paper 2 strategies reproduced (max gap EUR %.1e, %.1e QALYs)\n",
            gap_cost, gap_qaly))

# Second-line therapy (issue #209) is part of every cell and strategy cost
# above: per cell and arm it is the discounted PFS decrement x progression
# share x cost per progressing patient, the cell combination reproduces the
# strategies' second-line cost, and switching it off leaves QALYs unchanged.
df_psa_off <- model_fun(drop_second_line_therapy(params), determpsa = "psa", sim_idx = 1,
                        return_cells = TRUE)
cell_off <- attr(df_psa_off, "cells")
m_sl_cell <- cell_results$values[, , "cost"] - cell_off$values[, , "cost"]
c_progression <- second_line_cost_per_progression(params)
m_sl_hand <- m_sl_cell
for (cell in rownames(m_sl_hand)) for (arm in colnames(m_sl_hand)) {
  column <- cell_curve_columns(cell, arm)
  v_pf <- pmin(cell_curves$pfs[, column], cell_curves$os[, column])
  m_sl_hand[cell, arm] <- sum(c(0, -diff(v_pf)) * v_dw$cost) * params$p_pfs_progression * c_progression
}
df_policy_off <- policy_outcomes(cell_off$values, cell_off$mass, params)
v_sl_strategy <- df_psa$Cost[v_strat_rows] - df_psa_off$Cost[v_strat_rows]
stopifnot(c_progression > 0, all(m_sl_cell > 0), max(abs(m_sl_cell - m_sl_hand)) < 1e-8,
          max(abs((df_policy$cost[v_rows] - df_policy_off$cost[v_rows]) - v_sl_strategy)) < 1e-8,
          identical(cell_results$values[, , "qalys"], cell_off$values[, , "qalys"]),
          identical(df_psa$Effect, df_psa_off$Effect))
cat(sprintf("PASS: second-line cost per cell and arm, summed to the strategies (EUR %.0f to %.0f per cell)\n",
            min(m_sl_cell), max(m_sl_cell)))

# A cell whose experimental PFS crosses its OS: cell-level clamping differs
# from the subgroup-level clamp of model_fun(), and the gap is reported.
models_cross <- mock_models(pfs_mult = c("00" = 1.1, "01" = 0.9, "10" = 0.7, "11" = 0.2))
params_cross <- make_params(models_cross)
sampling_models <- list(joint = list(samples = list(list(os = list(model = models_cross$os),
  pfs = list(model = models_cross$pfs), failed = FALSE))))
df_cross <- suppressWarnings(model_fun(params_cross, determpsa = "psa", sim_idx = 1,
                                       return_cells = TRUE))
cell_cross <- attr(df_cross, "cells")
df_policy_cross <- policy_outcomes(cell_cross$values, cell_cross$mass, params_cross)
gap_cross <- df_policy_cross$qalys[match("nivo_0011", df_policy_cross$policy)] -
  df_cross$Effect[df_cross$Strategy == "crp"]
stopifnot(cell_cross$violations["11", "experimental"] > 0,
          sum(cell_cross$violations) == cell_cross$violations["11", "experimental"],
          abs(gap_cross) > 1e-6,
          abs(df_policy_cross$qalys[1] - df_cross$Effect[df_cross$Strategy == "control"]) < 1e-12)
# Life-years come from OS alone (issue #197), which is linear in the cell curves
# and never clamped, so the weighted cell life-years reproduce every strategy
# exactly even with the crossing cell (the cost/QALY clamp gap above remains).
v_rows_cross <- match(c("nivo_0000", "nivo_0011", "nivo_0101"), df_policy_cross$policy)
v_strat_cross <- match(c("control", "crp", "tmb_braf"), df_cross$Strategy)
gap_ly_cross <- max(abs(df_policy_cross$ly_disc[v_rows_cross] - df_cross$LY_disc[v_strat_cross]),
                    abs(df_policy_cross$ly_undisc[v_rows_cross] - df_cross$LY_undisc[v_strat_cross]))
stopifnot(gap_ly_cross < 1e-9)
cat(sprintf("PASS: a crossing cell curve is flagged and its clamp gap reported (%.2e QALYs; life-year gap %.1e)\n",
            gap_cross, gap_ly_cross))

# ---------------------------------------------------------------------------
# 6. End to end through build_psa_obj()
# ---------------------------------------------------------------------------
n_draw <- 8L
set.seed(200)
draw_models <- lapply(seq_len(n_draw), function(i) mock_models(
  os_rate = runif(1, 0.8, 1.2),
  rx_mult = setNames(runif(4, 0.6, 1.2), c("00", "01", "10", "11"))))
n_samples <- n_draw
sampling_models <- list(joint = list(samples = lapply(draw_models, function(m)
  list(os = list(model = m$os), pfs = list(model = m$pfs), failed = FALSE))))
config <- configure_parameter_distributions(params)
#' Build and bind a synthetic extended PSA at nivolumab price `price`.
build_pair <- function(price) {
  p <- params
  p$c_drug_nivo <- price
  build <- NULL
  capture.output(build <- suppressWarnings(build_psa_obj(p, config$distributions,
    get_strategies(), horizon, 1 / 52, n_draw, seed = 7L, currency = "EUR",
    extended_outcomes = TRUE)))
  build$psa_obj$fingerprint <- paste0("synthetic-", price)
  bind_psa_pair(build$psa_obj, build$psa_params)
}
pair_list <- build_pair(13923)
pair_bio <- build_pair(4641)
po <- pair_list$psa_obj
validate_psa_pair(po, pair_list$psa_params)
stopifnot(identical(dim(po$cells$values), c(n_draw, 4L, 2L, 5L)),
          identical(dim(po$ly_disc), c(n_draw, 3L)),
          max(abs(rowSums(po$cells$weights) - 1)) < 1e-12,
          max(abs(as.matrix(po$cells$weights) -
                  as.matrix(pair_list$psa_params[paste0("p_joint_", c("00", "01", "10", "11"))]))) < 1e-12,
          all(po$cells$nivolumab_price == 13923))
policy_space <- policy_space_from_psa(po, params)
df_gaps <- strategy_reproduction_gaps(policy_space, po)
df_gap_summary <- summarise_reproduction_gaps(df_gaps)
stopifnot(nrow(df_gaps) == 3L * n_draw, !any(df_gaps$cell_violation),
          all(df_gap_summary$reproduced_without_crossing),
          max(abs(df_gaps$gap_cost)) < 1e-8, max(abs(df_gaps$gap_qalys)) < 1e-12)
pairing3 <- policy_price_pairing(policy_space, c("nivo_0000", "nivo_0011", "nivo_0101"))
outcomes_po <- psa_outcome_matrices(po)
for (wtp in c(51000, 150000)) {
  stopifnot(abs(evpi_at_price(pairing3, 13923, wtp)$evpi -
                  calculate_evpi_from_nmb(outcomes_po$effect * wtp - outcomes_po$cost)) < 1e-6)
}
# The step 1a pairing of two PSA prices now also checks the life-years.
pairing_scen <- check_price_pairing(pair_list, pair_bio, 13923, 4641, control_strategy = "control")
stopifnot(all(pairing_scen$checks$passed),
          any(grepl("ly_undisc", pairing_scen$checks$check)))
# Policy costs at the biosimilar price from the list-price cells equal the
# policy costs of the PSA run at the biosimilar price.
policy_space_bio <- policy_space_from_psa(pair_bio$psa_obj, params)
pairing16 <- policy_price_pairing(policy_space)
stopifnot(max(abs(costs_at_price(pairing16, 4641) - policy_space_bio$outcomes$cost)) < 1e-6,
          max(abs(policy_space$outcomes$qalys - policy_space_bio$outcomes$qalys)) == 0)
expect_error_matching(policy_space_from_psa(list(cells = NULL), params), "no cell-level component")
cat("PASS: extended PSA object binds, reproduces the strategies, and pairs across prices\n")

settings <- list(wtp = c(51000, 100000), wtp_base = 51000, step = 2000,
                 marked_prices = c(biosimilar = 4641), everyone = "nivo_1111")
population <- list(annual_incidence = 519, research_horizon = 10, discount_rate = 0.04)
results <- run_voi_step1b(po, params, population, settings)
v_prob_cols <- grep("^prob_optimal_", names(results$evpi_price), value = TRUE)
df_16 <- results$evpi_price[results$evpi_price$policy_space == "16 policies", ]
df_15 <- results$evpi_price[grepl("everyone excluded", results$evpi_price$policy_space), ]
stopifnot(length(v_prob_cols) == 16L, nrow(df_16) == nrow(df_15),
          max(abs(rowSums(df_16[v_prob_cols]) - 1)) < 1e-12,
          all(is.na(df_15$prob_optimal_nivo_1111)),
          max(abs(rowSums(df_15[setdiff(v_prob_cols, "prob_optimal_nivo_1111")]) - 1)) < 1e-12,
          all(df_16$evpi_per_patient >= -1e-9), all(c(0, 4641, 13923) %in% df_16$c_drug_nivo),
          nrow(results$policies) == 16L, nrow(results$cells) == 4L * 2L * 5L,
          all(results$policies$reproduction_reproduced_without_crossing[
            !is.na(results$policies$paper2_strategy)]))
out_dir <- file.path(tempdir(), "voi_step1b_test")
v_paths <- write_voi_step1b_csvs(results, out_dir)
stopifnot(all(file.exists(v_paths)),
          identical(basename(v_paths), c("voi_cells.csv", "voi_policies.csv", "voi_evpi_price.csv")),
          nrow(utils::read.csv(v_paths[3])) == nrow(results$evpi_price))
unlink(out_dir, recursive = TRUE)
cat("PASS: step 1b driver: 16- and 15-policy EVPI on the price grid, three CSVs written\n")

# ---------------------------------------------------------------------------
# 7./8. Live checks (after the regeneration of the PSA and scenario caches)
# ---------------------------------------------------------------------------
#' Read and validate a PSA pair; NULL with a SKIP line when absent or legacy.
read_live_pair <- function(obj_file, params_file) {
  if (!file.exists(obj_file) || !file.exists(params_file)) {
    cat("SKIP: live PSA reproduction (", basename(obj_file), " absent)\n", sep = "")
    return(NULL)
  }
  po_live <- readRDS(obj_file)
  pp_live <- readRDS(params_file)
  ok <- tryCatch({ validate_psa_pair(po_live, pp_live); TRUE }, error = function(e) {
    cat("SKIP: live PSA reproduction (", conditionMessage(e), ")\n", sep = "")
    FALSE
  })
  if (ok) list(psa_obj = po_live, psa_params = pp_live) else NULL
}
live <- read_live_pair(psa_obj_path("correct"), psa_params_path("correct"))
if (!is.null(live)) {
  live_params <- list(c_test_CRP = live$psa_obj$cells$test_costs[["crp"]],
                      c_test_NGS = live$psa_obj$cells$test_costs[["tmb_braf"]])
  space_live <- policy_space_from_psa(live$psa_obj, live_params)
  gaps_live <- summarise_reproduction_gaps(strategy_reproduction_gaps(space_live, live$psa_obj))
  print(gaps_live)
  stopifnot(all(gaps_live$reproduced_without_crossing))
  cat("PASS: live PSA: the three paper 2 policies reproduce the strategies in every draw without a cell crossing\n")
}

scenario_file <- scenario_evppi_path("correct")
s7_file <- here::here("outputs", "tables", "table_s7.csv")
if (!file.exists(scenario_file) || !file.exists(s7_file)) {
  cat("SKIP: live biosimilar EVPI (scenario cache or Table S7 absent)\n")
} else {
  scenario_cache <- readRDS(scenario_file)
  price_pairs <- tryCatch(scenario_price_pairs(scenario_cache), error = function(e) NULL)
  has_cells <- !is.null(price_pairs) && length(price_pairs) == 2L &&
    all(vapply(price_pairs, function(x) !is.null(x$psa_obj$cells), logical(1)))
  if (!has_cells) {
    cat("SKIP: live biosimilar EVPI (scenario cache predates the cell outcomes of #200;",
        "rerun analysis/12_scenario_EVPPIs.R)\n")
  } else {
    for (x in price_pairs) validate_psa_pair(x$psa_obj, x$psa_params)
    po_ref <- price_pairs[[1]]$psa_obj
    price_alt <- price_pairs[[2]]$price
    live_params <- list(c_test_CRP = po_ref$cells$test_costs[["crp"]],
                        c_test_NGS = po_ref$cells$test_costs[["tmb_braf"]])
    space_ref <- policy_space_from_psa(po_ref, live_params)
    pairing_live <- policy_price_pairing(space_ref, c("nivo_0000", "nivo_0011", "nivo_0101"))
    evpi_cells <- evpi_at_price(pairing_live, price_alt, 51000)$evpi
    out_alt <- psa_outcome_matrices(price_pairs[[2]]$psa_obj)
    evpi_psa <- calculate_evpi_from_nmb(out_alt$effect * 51000 - out_alt$cost)
    df_s7 <- utils::read.csv(s7_file)
    evpi_s7 <- df_s7$evpi[df_s7$scenario_id == "s1_decreased_nivo" & df_s7$group == "evpi"]
    gaps_ref <- summarise_reproduction_gaps(strategy_reproduction_gaps(space_ref, po_ref))
    cat(sprintf("Live EVPI at EUR %s and EUR 51,000: cells %.6f, scenario PSA %.6f, Table S7 %.2f (576.72 before the regeneration); draws with a cell crossing: %s\n",
                format(price_alt, big.mark = ","), evpi_cells, evpi_psa, evpi_s7,
                paste(gaps_ref$n_draws_cell_crossing, collapse = "/")))
    tol_live <- if (all(gaps_ref$n_draws_cell_crossing == 0)) 1e-6 else
      max(abs(c(gaps_ref$max_abs_gap_cost, 51000 * gaps_ref$max_abs_gap_qalys))) + 1e-6
    stopifnot(length(evpi_s7) == 1, abs(evpi_cells - evpi_psa) <= tol_live,
              abs(evpi_psa - evpi_s7) < 0.005)
    cat("PASS: live biosimilar EVPI from the list-price cells equals the scenario PSA and Table S7\n")
  }
}

cat("PASS: test_voi_policies.R\n")
