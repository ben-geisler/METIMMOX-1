# ===============================================================================
# TEST: second-line therapy in the model (issue #209, part 2)
# ===============================================================================
# Synthetic (no trial data or cache):
#   1. The closed form second_line_cost_per_progression() on the model inputs
#      of second_line_model_inputs() equals second_line_cost_at() on the
#      pathway table, at the base values and with shares, cycles and contacts
#      moved; contacts and blood tests follow the model unit costs.
#   2. Accrual at incident progression: the cost of each interval is the PFS
#      decrement x the progression share x the cost per progressing patient,
#      charged at the grid point closing the interval and discounted there;
#      undiscounted it sums to share x cost x (1 - final PFS); it does not
#      depend on post-progression survival.
#   3. Zero second-line cost reproduces the former results exactly: a list
#      with p_sl_treated = 0 gives the same outcomes as a list without any
#      second-line field (the pre-#209 model); a partial list stops.
#   4. pfs_progression_share() is the count ratio of "progression" over
#      "progression" + "death_within_window" PFS events, its standard error
#      gives Beta(n_progression, n_death) in distribution_parameters(), and
#      it stops without a death or without a progression.
#   5. Structural scenarios: Second_line_therapy = 0 removes the cost; the
#      flat Post_progression_cost rate replaces it.
# Live (skips without data/tidy/METIMMOX.rds or the pathway table):
#   6. Script 05's p_pfs_progression equals the count ratio computed here from
#      data_complete$PFS_rule; its second-line cost equals second_line_cost()
#      on the sourced table; switching the pathway off reproduces the model
#      without second-line fields exactly; the base-case increase in cost is
#      the discounted second-line cost.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_second_line_model.R
# ===============================================================================

source("R/model_configs.R")
source("R/calculate_outcomes.R")
source("R/second_line_costs.R")
source("R/parameter_distributions.R")

#' Stop unless two numbers agree within a relative tolerance.
expect_near <- function(x, y, tol = 1e-10) {
  if (!isTRUE(all(abs(x - y) <= tol * pmax(1, abs(y))))) {
    stop("Expected ", paste(signif(y, 10), collapse = ", "), ", got ", paste(signif(x, 10), collapse = ", "))
  }
  invisible(TRUE)
}

# --- 1. Closed form against the pathway table ---------------------------------
#' One row of a synthetic pathway table.
row <- function(component, parameter, value, lower = value, upper = value, unit,
                regimen = NA, drug = NA, vial_mg = NA, administrations_per_cycle = NA,
                cost_parameter = NA) {
  data.frame(component = component, parameter = parameter, regimen = regimen, drug = drug,
             value = value, lower = lower, upper = upper, unit = unit, vial_mg = vial_mg,
             administrations_per_cycle = administrations_per_cycle,
             cost_parameter = cost_parameter, source = "Synthetic", doi_or_url = "https://example.org",
             retrieved = "2026-10-17", note = "", stringsAsFactors = FALSE)
}
df_table <- rbind(
  row("pathway", "share_treated", 0.8, 0.6, 0.9, unit = "proportion"),
  row("pathway", "share_folfiri", 0.4, 0.2, 0.8, unit = "proportion"),
  row("pathway", "share_ras_braf_wt", 0.35, 0.3, 0.5, unit = "proportion"),
  row("pathway", "share_antiegfr_if_wt", 0.5, 0.25, 0.75, unit = "proportion"),
  row("pathway", "share_bevacizumab_if_no_antiegfr", 0.3, 0, 0.5, unit = "proportion"),
  row("schedule", "cycle_length_weeks", 2, unit = "weeks"),
  row("patient", "bsa", 1.9, 1.8, 2, unit = "m2"),
  row("patient", "weight", 70, unit = "kg"),
  row("dose", "dose_FOLFIRI_irinotecan", 180, unit = "mg/m2", regimen = "FOLFIRI",
      drug = "irinotecan", administrations_per_cycle = 1),
  row("dose", "dose_FLIRI_irinotecan", 180, unit = "mg/m2", regimen = "FLIRI",
      drug = "irinotecan", administrations_per_cycle = 1),
  row("dose", "dose_FLIRI_folinate", 60, unit = "mg", regimen = "FLIRI",
      drug = "folinate", administrations_per_cycle = 2),
  row("dose", "dose_bevacizumab", 5, unit = "mg/kg", regimen = "bevacizumab",
      drug = "bevacizumab", administrations_per_cycle = 1),
  row("dose", "dose_panitumumab", 6, unit = "mg/kg", regimen = "panitumumab",
      drug = "panitumumab", administrations_per_cycle = 1),
  row("dose_intensity", "relative_dose_intensity", 1, 0.9, 1, unit = "proportion"),
  row("price", "price_irinotecan_100mg", 1500, unit = "NOK per vial", drug = "irinotecan", vial_mg = 100),
  row("price", "price_irinotecan_40mg", 650, unit = "NOK per vial", drug = "irinotecan", vial_mg = 40),
  row("price", "price_folinate_100mg", 330, unit = "NOK per vial", drug = "folinate", vial_mg = 100),
  row("price", "price_bevacizumab_400mg", 9700, unit = "NOK per vial", drug = "bevacizumab", vial_mg = 400),
  row("price", "price_panitumumab_100mg", 3800, unit = "NOK per vial", drug = "panitumumab", vial_mg = 100),
  row("currency", "nok_per_eur", 10.9, unit = "NOK per EUR"),
  row("wastage", "vial_sharing", 0, 0, 1, unit = "proportion"),
  row("resource", "contacts_per_cycle_FOLFIRI", 2, 1, 2, unit = "contacts per cycle",
      cost_parameter = "c_other_visit"),
  row("resource", "contacts_per_cycle_FLIRI", 2, unit = "contacts per cycle",
      cost_parameter = "c_other_visit"),
  row("resource", "blood_tests_per_cycle", 1, unit = "tests per cycle", cost_parameter = "c_test_blood"),
  row("duration", "cycles_none", 9.8, 8, 10.4, unit = "cycles"),
  row("duration", "cycles_bevacizumab", 12, 9.1, 12.9, unit = "cycles"),
  row("duration", "cycles_panitumumab", 10.8, 9, 11.4, unit = "cycles"))
path <- tempfile(fileext = ".csv")
utils::write.csv(df_table, path, row.names = FALSE)
df_pathway <- read_second_line_pathway(path)
inputs <- second_line_model_inputs(df_pathway)
unit_costs <- list(c_other_visit = 213.46, c_test_blood = 16)

#' Parameter list with the second-line model inputs, a progression share and
#' the model unit costs.
sl_params <- function(values = inputs$values, share = 0.9, costs = unit_costs) {
  c(as.list(values), list(p_pfs_progression = share), costs)
}
stopifnot(setequal(c(names(inputs$values), "p_pfs_progression"), SECOND_LINE_MODEL_PARAMETERS),
          setequal(names(inputs$se), c(grep("^p_sl_", SECOND_LINE_MODEL_PARAMETERS, value = TRUE),
                                       grep("^n_sl_cycles_", SECOND_LINE_MODEL_PARAMETERS, value = TRUE))),
          all(is.na(inputs$lower[c("n_sl_contacts_FLIRI", "n_sl_blood_tests", "c_sl_drug_FLIRI")])),
          inputs$lower[["n_sl_contacts_FOLFIRI"]] == 1)
expect_near(inputs$se[["p_sl_treated"]], (0.9 - 0.6) / (2 * qnorm(0.975)))
p_base <- sl_params()
expect_near(second_line_cost_per_progression(p_base),
            second_line_cost_at(df_pathway, unit_costs)$per_progressing)
# Shares, cycles and contacts moved: the closed form follows the table.
v_setting <- c(share_treated = "upper", share_folfiri = "lower", share_ras_braf_wt = "upper",
               share_antiegfr_if_wt = "lower", share_bevacizumab_if_no_antiegfr = "upper",
               cycles_none = "lower", cycles_bevacizumab = "upper", cycles_panitumumab = "lower",
               contacts_per_cycle_FOLFIRI = "lower")
v_moved <- second_line_values(df_pathway, v_setting)
v_model_moved <- inputs$values
for (field in names(SECOND_LINE_TABLE_FIELDS)) {
  v_model_moved[[field]] <- v_moved[[SECOND_LINE_TABLE_FIELDS[[field]]]]
}
expect_near(second_line_cost_per_progression(sl_params(v_model_moved)),
            second_line_cost_at(df_pathway, unit_costs, v_moved)$per_progressing)
# Contacts at the model's visit cost: the PSA's sampled c_other_visit moves it.
costs_up <- list(c_other_visit = 300, c_test_blood = 16)
expect_near(second_line_cost_per_progression(sl_params(costs = costs_up)),
            second_line_cost_at(df_pathway, costs_up)$per_progressing)
cat("PASS: closed-form cost per progressing patient equals second_line_cost_at() on the table\n")

# --- 2. Accrual at incident progression and discounting -----------------------
zero_params <- list(u_np = 0.7, u_p = 0.5, c_drug_nivo = 0, c_drug_FLOX = 0, c_test_CT = 0,
                    c_test_blood = 16, c_test_CRP = 0, c_test_NGS = 0, c_other_visit = 213.46,
                    c_other_baseline = 0, c_other_follow = 0, c_other_pp = 0, c_other_last = 0,
                    l_nivo = 0, l_FLOX_exp = 0, l_FLOX_control = 0, l_CT = 0, l_blood = 0, l_visit = 0)
p_sl <- c(zero_params[setdiff(names(zero_params), names(unit_costs))], sl_params())
v_pfs <- c(1, 0.9, 0.9, 0.6, 0.35, 0.2)
v_os <- c(1, 1, 0.95, 0.8, 0.6, 0.5)
states <- partitioned_survival_states(v_os, v_pfs)
v_dw <- 1 / 1.04^(seq_along(v_pfs) - 1)  # yearly grid: discounting is visible
c_progression <- second_line_cost_per_progression(p_sl)
out <- calculate_outcomes(p_sl, states$p_pf, states$p_p, states$p_d, "standard", NULL,
                          v_dw, v_dw, 1)
v_expected <- c(0, 0.1, 0, 0.3, 0.25, 0.15) * 0.9 * c_progression
expect_near(out$second_line_costs, v_expected)
expect_near(out$costs_total, sum(v_expected * v_dw))
expect_near(sum(out$second_line_costs), 0.9 * c_progression * (1 - tail(v_pfs, 1)))
stopifnot(sum(v_expected * v_dw) < sum(v_expected))
# Longer post-progression survival (higher OS, same PFS) leaves it unchanged.
states_long <- partitioned_survival_states(pmax(v_os, c(1, 1, 1, 0.95, 0.9, 0.9)), v_pfs)
out_long <- calculate_outcomes(p_sl, states_long$p_pf, states_long$p_p, states_long$p_d,
                               "standard", NULL, v_dw, v_dw, 1)
expect_near(out_long$costs_total, out$costs_total)
# Flat PFS: no progression, no cost; a one-point grid spans no time.
stopifnot(all(calculate_outcomes(p_sl, rep(1, 4), rep(0, 4), rep(0, 4), "standard", NULL,
                                 rep(1, 4), rep(1, 4), 1)$second_line_costs == 0),
          calculate_outcomes(p_sl, 1, 0, 0, "standard", NULL, 1, 1, 1)$second_line_costs == 0)
cat("PASS: second-line cost charged once at incident progression and discounted there\n")

# --- 3. Zero second-line cost reproduces the former model exactly -------------
p_old <- zero_params
p_zero <- p_sl
p_zero$p_sl_treated <- 0
out_old <- calculate_outcomes(p_old, states$p_pf, states$p_p, states$p_d, "experimental", "crp",
                              v_dw, v_dw, 1)
out_zero <- calculate_outcomes(p_zero, states$p_pf, states$p_p, states$p_d, "experimental", "crp",
                               v_dw, v_dw, 1)
stopifnot(identical(out_zero[setdiff(names(out_zero), "second_line_costs")],
                    out_old[setdiff(names(out_old), "second_line_costs")]),
          all(out_zero$second_line_costs == 0), all(out_old$second_line_costs == 0),
          !second_line_costed(p_old), second_line_costed(p_sl))
p_partial <- p_old
p_partial$p_sl_treated <- 0.8
stopifnot(grepl("Second-line parameters incomplete",
                tryCatch(second_line_costed(p_partial), error = conditionMessage)))
p_bad <- p_sl
p_bad$p_sl_folfiri <- 1.2
stopifnot(grepl("must be in \\[0, 1\\]",
                tryCatch(second_line_cost_per_progression(p_bad), error = conditionMessage)))
cat("PASS: zero second-line cost reproduces the pre-#209 outcomes exactly\n")

# --- 4. Progression share of PFS events ----------------------------------------
v_levels <- c("progression", "death_within_window", "death_censored_at_assessment",
              "censored_alive", "anchor_missing")
df_rule <- data.frame(
  PFS_rule = factor(c(rep("progression", 7), rep("death_within_window", 2),
                      rep("death_censored_at_assessment", 3), "censored_alive"), levels = v_levels),
  Rx = factor(c(rep("A", 4), rep("B", 3), "A", "B", "A", "B", "B", "A")))
share <- pfs_progression_share(df_rule)
stopifnot(share$n_progression == 7, share$n_death == 2, share$n_events == 9,
          share$upper == 1, share$lower < share$value,
          identical(as.vector(share$by_arm), c(4L, 1L, 3L, 1L)))
expect_near(share$value, 7 / 9)
beta <- distribution_parameters(share$value, "beta", share$se / share$value)
expect_near(c(beta$shape1, beta$shape2), c(7, 2))
expect_near(share$lower, qbeta(0.025, 7, 2))
stopifnot(grepl("at least one progression and one death",
                tryCatch(pfs_progression_share(df_rule[df_rule$PFS_rule != "death_within_window", ]),
                         error = conditionMessage)),
          grepl("needs PFS_rule", tryCatch(pfs_progression_share(df_rule["Rx"]), error = conditionMessage)))
cat("PASS: progression share is the count ratio of PFS events with Beta(n_progression, n_death)\n")

# --- 5. Structural scenarios ---------------------------------------------------
sc_off <- apply_structural_scenario("Second_line_therapy", 0, p_sl, 5)$params
sc_flat <- apply_structural_scenario("Post_progression_cost", 5000, p_sl, 5)$params
stopifnot(second_line_cost_per_progression(sc_off) == 0, sc_flat$c_other_pp == 5000,
          second_line_cost_per_progression(sc_flat) == 0,
          identical(apply_structural_scenario("Second_line_therapy", 1, p_sl, 5)$params, p_sl))
out_flat <- calculate_outcomes(sc_flat, states$p_pf, states$p_p, states$p_d, "standard", NULL,
                               v_dw, v_dw, 1)
stopifnot(all(out_flat$second_line_costs == 0), sum(out_flat$post_progression_costs) > 0)
cat("PASS: no-second-line and flat-rate scenarios switch the pathway off\n")

# --- 6. Live -------------------------------------------------------------------
if (!file.exists("data/tidy/METIMMOX.rds") || !file.exists(SECOND_LINE_PATHWAY_PATH)) {
  cat("SKIP: live checks need data/tidy/METIMMOX.rds and", SECOND_LINE_PATHWAY_PATH, "\n")
} else {
  suppressPackageStartupMessages({
    library(dplyr); library(flexsurv); library(survival)
  })
  options(metimmox.sampling_allow_regenerate = FALSE)
  invisible(capture.output(suppressMessages(suppressWarnings({
    for (script in c("02_setup_and_global_variables", "03_biomarker_strategies",
                     "04_parametric_survival_analysis", "05_basecase_input_parameters")) {
      source(file.path("analysis", paste0(script, ".R")))
    }
  }))))
  source("R/model_fun.R"); source("R/calculate_outcomes.R"); source("R/prediction_functions.R")
  v_rule <- as.character(data_complete$PFS_rule)
  n_prog <- sum(v_rule == "progression")
  n_events <- n_prog + sum(v_rule == "death_within_window")
  expect_near(l_params_base$p_pfs_progression, n_prog / n_events)
  expect_near(second_line_cost_per_progression(l_params_base),
              second_line_cost(read_second_line_pathway(), l_params_base)$summary["base", "per_progressing"])
  df_base <- model_fun(l_params_base)
  df_off <- model_fun(drop_second_line_therapy(l_params_base))
  p_without <- l_params_base[setdiff(names(l_params_base), SECOND_LINE_MODEL_PARAMETERS)]
  df_without <- model_fun(p_without)
  stopifnot(identical(df_off, df_without),
            identical(df_base[c("Effect", "LY_disc", "LY_undisc")], df_off[c("Effect", "LY_disc", "LY_undisc")]),
            all(df_base$Cost > df_off$Cost))
  # Control: the cost increase is the discounted second-line cost of its PFS curve.
  v_dw_live <- discount_weights(l_params_base, time_horizon + 1)$cost
  v_pfs_control <- l_params_base$p_pfs[[paste0(get_control_strategy(), "_PFS")]]
  expected_control <- sum(c(0, -diff(v_pfs_control)) * v_dw_live) *
    l_params_base$p_pfs_progression * second_line_cost_per_progression(l_params_base)
  expect_near(df_base$Cost[df_base$Strategy == get_control_strategy()] -
                df_off$Cost[df_off$Strategy == get_control_strategy()], expected_control, tol = 1e-9)
  cat(sprintf("PASS (live): progression share %d/%d, EUR %.0f per progressing patient, control +EUR %.0f\n",
              n_prog, n_events, second_line_cost_per_progression(l_params_base), expected_control))
}
