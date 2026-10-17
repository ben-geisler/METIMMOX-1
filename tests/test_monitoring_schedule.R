# ===============================================================================
# TEST: monitoring and contact schedule from the trial team's answers (issue #207)
# ===============================================================================
# Answers MON1-MON3 of 15 October 2026:
# 1. CT every 8 weeks from baseline while progression-free, to the horizon
#    (PF_CT_INTERVAL_WEEKS), on the 521-point grid and on every short horizon.
# 2. Blood test at every administration week on treatment and, off treatment
#    (break, after the last sequence), with every CT; identical in both arms.
# 3. Outpatient contacts: a FLOX administration is two contacts, a nivolumab
#    administration one; visits with every CT are kept (#164). Counted per arm
#    through calculate_outcomes()$visit_contacts and the visit cost.
# 4. Second_sequence scenario: dropping weeks 24-38 removes their
#    administrations, their blood tests and both FLOX contacts, and keeps the
#    CT scans with their visit and blood test.
# 5. Progressed state: follow-up visit and CT every 2.5 months (4.8 per
#    patient-year); the post-progression treatment cost stays quarterly.
# 6. Every schedule stays a 0/1 vector of time_horizon + 1 points.
#
# Synthetic inputs only; no trial data or cache is read.
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_monitoring_schedule.R
# ===============================================================================

source("R/model_configs.R")
source("R/parameter_distributions.R")
source("R/calculate_outcomes.R")

#' Assert that two numeric values or vectors agree within 1e-10.
equal <- function(x, y) stopifnot(length(x) == length(y), max(abs(x - y)) < 1e-10)

# ---------------------------------------------------------------------------
# 1-2. Positions on the base-case grid
# ---------------------------------------------------------------------------
stopifnot(PF_CT_INTERVAL_WEEKS == 8L, FLOX_CONTACTS_PER_ADMINISTRATION == 2L,
          PROGRESSED_MONITORING_INTERVAL_MONTHS == 2.5)
equal(PROGRESSED_MONITORING_PER_YEAR, 4.8)

s <- build_treatment_schedules(521)
v_admin <- which(s$l_nivo == 1 | s$l_FLOX_exp == 1 | s$l_FLOX_control == 1)
stopifnot(identical(v_admin, c(seq(1L, 15L, by = 2L), seq(25L, 39L, by = 2L))),
          # Both arms have an administration every two weeks of a sequence.
          identical(which(s$l_nivo == 1 | s$l_FLOX_exp == 1), which(s$l_FLOX_control == 1)),
          identical(which(s$l_CT == 1), seq(1L, 521L, by = 8L)),
          identical(which(s$l_blood == 1), sort(union(v_admin, which(s$l_CT == 1)))),
          # Off treatment the blood test is drawn with the CT only.
          identical(which(s$l_blood == 1 & !seq_len(521) %in% v_admin),
                    c(17L, seq(41L, 521L, by = 8L))),
          identical(s$l_blood, s$l_visit))

# Counts over the first 52 weeks (positions 1-53) and the horizon.
s53 <- build_treatment_schedules(53)
stopifnot(sum(s53$l_CT) == 7, sum(s53$l_blood) == 19, sum(s53$l_visit) == 19,
          sum(s53$l_visit + s53$l_FLOX_exp) == 27, sum(s53$l_visit + s53$l_FLOX_control) == 35,
          sum(s$l_CT) == 66, sum(s$l_blood) == 78, sum(s$l_visit) == 78)

# Every horizon: 0/1 vectors of exactly n points; CT every 8 weeks; blood and
# visits on administrations and CTs only.
for (n in 1:60) {
  sn <- build_treatment_schedules(n)
  stopifnot(all(vapply(sn, length, integer(1)) == n), all(unlist(sn) %in% c(0, 1)),
            identical(which(sn$l_CT == 1), seq(1L, n, by = 8L)),
            identical(sn$l_blood, sn$l_visit),
            identical(sn$l_blood, s$l_blood[seq_len(n)]))
}

# ---------------------------------------------------------------------------
# 3. Contacts per arm in calculate_outcomes()
# ---------------------------------------------------------------------------
params <- c(s, list(u_np = 1, u_p = 1, c_drug_nivo = 0, c_drug_FLOX = 0,
  c_test_CT = 0, c_test_blood = 0, c_test_CRP = 0, c_test_NGS = 0,
  c_other_visit = 1, c_other_baseline = 0, c_other_follow = 0,
  c_other_pp = 0, c_other_last = 0))
#' calculate_outcomes() of a cohort that stays progression-free, undiscounted.
run_pf <- function(p, arm, n = length(p$l_visit)) {
  calculate_outcomes(p, rep(1, n), rep(0, n), rep(0, n), arm, NULL,
                     rep(1, n), rep(1, n), 1 / 52)
}
out_exp <- run_pf(params, "experimental")
out_ctrl <- run_pf(params, "standard")
equal(out_exp$visit_contacts, s$l_visit + s$l_FLOX_exp)
equal(out_ctrl$visit_contacts, s$l_visit + s$l_FLOX_control)
# FLOX days are two contacts, nivolumab administrations one.
stopifnot(all(out_exp$visit_contacts[s$l_FLOX_exp == 1] == 2),
          all(out_exp$visit_contacts[s$l_nivo == 1] == 1),
          all(out_ctrl$visit_contacts[s$l_FLOX_control == 1] == 2),
          all(out_exp$visit_contacts[s$l_CT == 1 & !seq_len(521) %in% v_admin] == 1))
equal(sum(out_exp$visit_costs), 86)
equal(sum(out_ctrl$visit_costs), 94)
# Blood tests are the same in both arms: one per l_blood point.
params_blood <- params; params_blood$c_other_visit <- 0; params_blood$c_test_blood <- 1
equal(sum(run_pf(params_blood, "experimental")$test_costs), 78)
equal(sum(run_pf(params_blood, "standard")$test_costs), 78)

# ---------------------------------------------------------------------------
# 4. Second_sequence scenario
# ---------------------------------------------------------------------------
dropped <- apply_structural_scenario("Second_sequence", 0, params, 520)$params
v_late <- 25:521
stopifnot(all(dropped$l_nivo[v_late] == 0), all(dropped$l_FLOX_exp[v_late] == 0),
          all(dropped$l_FLOX_control[v_late] == 0),
          identical(dropped$l_CT, s$l_CT),
          # Blood tests and visits from week 24 only with the CT.
          identical(which(dropped$l_blood[v_late] == 1) + 24L, seq(25L, 521L, by = 8L)),
          identical(which(dropped$l_visit[v_late] == 1) + 24L, seq(25L, 521L, by = 8L)),
          identical(dropped$l_blood[1:24], s$l_blood[1:24]),
          identical(dropped$l_visit[1:24], s$l_visit[1:24]),
          all(unlist(dropped[names(s)]) %in% c(0, 1)),
          all(vapply(dropped[names(s)], length, integer(1)) == 521))
d_exp <- run_pf(dropped, "experimental")
d_ctrl <- run_pf(dropped, "standard")
# Weeks 24-38 held 8 administrations in each arm, 2 of them CT weeks (24, 32):
# 6 visits lost in each arm plus the second FLOX contacts (4 experimental, 8 control).
equal(sum(out_exp$visit_costs) - sum(d_exp$visit_costs), 6 + 4)
equal(sum(out_ctrl$visit_costs) - sum(d_ctrl$visit_costs), 6 + 8)
stopifnot(all(d_exp$visit_contacts[v_late] <= 1), all(d_ctrl$visit_contacts[v_late] <= 1))
equal(sum(run_pf(params_blood, "experimental")$test_costs) -
        sum(run_pf(apply_structural_scenario("Second_sequence", 0, params_blood, 520)$params,
                   "experimental")$test_costs), 6)
# The base-case endpoint of the scenario leaves the schedules unchanged.
stopifnot(identical(apply_structural_scenario("Second_sequence", 1, params, 520)$params, params))

# ---------------------------------------------------------------------------
# 5. Progressed-state rates
# ---------------------------------------------------------------------------
p_prog <- params
p_prog$c_other_visit <- 0
p_prog$c_other_follow <- 100
p_prog$c_test_CT <- 10
p_prog$c_other_pp <- 1000
n <- 53  # one year
out_prog <- calculate_outcomes(p_prog, rep(0, n), rep(1, n), rep(0, n), "standard", NULL,
                               rep(1, n), rep(1, n), 1 / 52)
equal(sum(out_prog$follow_up_costs), 100 * 4.8)
equal(sum(out_prog$progressed_imaging_costs), 10 * 4.8)
equal(sum(out_prog$post_progression_costs), 1000 * 4)

cat("PASS: 8-weekly CT, blood tests per administration and with every CT, two contacts per FLOX administration, Second_sequence scenario and 2.5-monthly progressed-state monitoring (issue #207).\n")
