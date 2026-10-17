# ===============================================================================
# TEST: candidate study designs of paper 3 (issue #203, plan step 5)
# ===============================================================================
# Synthetic joint gamma component (simulated patients without identifiers; no
# trial data, no caches):
#
#   1. Design table: S1, S2, S3, S3R and S4 enrol the right cells, randomise
#      at the right week, pay the right screening tests (S1 CRP + NGS, S2 CRP,
#      S3 and S3R NGS, S4 none), and identify the targets of the information
#      engine (S1 {Rx, Rx:crp, Rx:tmb_braf}; S2 {Rx:tmb_braf, cell_effect_10};
#      S3R {Rx:crp, cell_effect_01}; S3 {control log rate, CRP, experimental
#      log rate of the TMB/BRAF-positive cells, Rx:crp}; S4 the
#      experimental-arm log rates); design_table() writes one row per design
#      plus the S6 power-based rows.
#   2. Accrual and reporting-year arithmetic: METIMMOX's rate from its dates,
#      the Nordic/European rates scaled by the eligible population, the
#      enrolment fraction (enrolled mass, times P(PFS > 4) for randomisation at
#      week 4), and reporting year = L/52 + n / rate + minimum follow-up.
#   3. Schoenfeld sizes against a hand calculation: 4 z^2 / log(RR)^2 events
#      for a main effect, four times that for an interaction; the heuristic
#      variance of balanced event tables; power_based_n() returns the
#      smallest n meeting the variance target.
#   4. Heterogeneity (Woodbury): tau = 0 leaves the information unchanged;
#      a finite tau equals (P^-1 + tau^2 D D')^-1 for non-singular P; tau = Inf
#      removes the information in the biased directions.
#   5. Combined-design information: information_target_covariance()
#      reproduces design_information()'s V_d(n) (recoloured and independent);
#      the S4 + S1 information is the sum of both; on the shared targets the
#      combined covariance is no larger (Loewner order) than S1's; the S4
#      single arm aliases rate with Rx and identifies its targets only.
#   6. Power from the engine (ruling (4) of 11 October 2026): with gamma
#      shape 1 (exponential PFS) the engine's power-based n
#      (power_based_n_information()) and Schoenfeld's (power_based_n()) agree
#      within 10% for S3R's conditional effect and both S1 interactions; the
#      engine's n is the smallest meeting the variance target; the primary n
#      is the engine's, the largest over S1's interactions.
#   7. Decision-relevant effect: the CRP interaction shift at which the
#      CRP-guided policy breaks even is more beneficial at a higher nivolumab
#      price; the incremental NMB is zero at the root; cost is linear in the
#      price; no root at a prohibitive price (NA with the reason); the design
#      layer gives both hypotheses by both methods, n only for a benefit.
#   8. METIMMOX-3 (S3, protocol v1.1, #213; stratified, non-randomised):
#      experimental data in the TMB/BRAF-positive cells and control data in
#      the negative cells only; information rank 7 of 9 per endpoint, Rx:crp
#      identified and the treatment effect in the positive cells not; with
#      tau_soc = 0 the start-up and standard-of-care periods sum to the whole
#      control follow-up; tau_soc discounts the start-up and the standard of
#      care together (one shared bias) and removes information on the
#      control log rate; tau_m3 = 0 leaves the experimental information
#      unchanged, a larger tau_m3 removes information on the CRP interaction
#      and never on the control cells, tau_m3 = Inf leaves Rx:crp
#      unidentified; the experimental share is the positive cells' mass; the
#      planned n reports at the protocol's data cut-off (planned last
#      inclusion plus the minimum follow-up); the size is described as
#      descriptive.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_designs.R
# ===============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(flexsurv)
})
for (f in c("model_configs", "evppi_functions", "joint_survival_sampling", "pfs_os_violation_diagnostics",
            "voi_survival_parameters", "design_information", "evsi_functions", "voi_price_axis",
            "voi_designs", "voi_research_value")) {
  source(here::here("R", paste0(f, ".R")))
}

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Synthetic joint gamma component in the sampling-cache layout (no identifiers).
make_component <- function(n = 300, seed = 11, n_draws = 1000, rho = 0.4) {
  set.seed(seed)
  df <- data.frame(Age = stats::rnorm(n, 62, 9), sex = factor(stats::rbinom(n, 1, 0.45)),
                   Rx = factor(stats::rbinom(n, 1, 0.5), labels = c("Control", "Experimental arm")),
                   crp = factor(stats::rbinom(n, 1, 0.35)), tmb_braf = factor(stats::rbinom(n, 1, 0.45)))
  v_eta <- log(0.02) + 0.01 * (df$Age - 62) - 0.2 * (df$Rx != "Control") -
    0.3 * (df$crp == "1") * (df$Rx != "Control")
  v_t <- stats::rgamma(n, 1.6, exp(v_eta)); v_c <- stats::runif(n, 30, 150)
  df$OSwk <- pmin(v_t, v_c); df$Death <- as.numeric(v_t <= v_c)
  v_t2 <- pmin(v_t, stats::rgamma(n, 1.2, exp(v_eta + 0.8)))
  df$PFSwk <- pmin(v_t2, v_c); df$Prog <- as.numeric(v_t2 <= v_c)
  rhs <- "Age + sex + Rx + crp * Rx + tmb_braf * Rx"
  fit_os <- flexsurvreg(stats::as.formula(paste("Surv(OSwk, Death) ~", rhs)), data = df, dist = "gamma")
  fit_pfs <- flexsurvreg(stats::as.formula(paste("Surv(PFSwk, Prog) ~", rhs)), data = df, dist = "gamma")
  m_ro <- covariance_root(fit_os$cov); m_rp <- covariance_root(fit_pfs$cov)
  m_joint <- rbind(cbind(fit_os$cov, rho * m_ro %*% m_rp), cbind(rho * m_rp %*% m_ro, fit_pfs$cov))
  v_names <- c(paste0("os::", rownames(fit_os$res)), paste0("pfs::", rownames(fit_pfs$res)))
  dimnames(m_joint) <- list(v_names, v_names)
  list(component = list(method = "mvn_joint_v2", joint_covariance = list(covariance = m_joint),
                        draws = draw_joint_survival_coefficients(fit_os, fit_pfs, m_joint, n_draws, seed),
                        original_os = fit_os, original_pfs = fit_pfs, n_samples = n_draws, n_failed = 0L,
                        dist_os = "gamma", dist_pfs = "gamma", seed = seed),
       reference = df[, c("Age", "sex", "Rx", "crp", "tmb_braf")])
}

synthetic <- make_component()
component <- synthetic$component
df_ref <- synthetic$reference
draw <- split_joint_draw(component_mean_draw(component), component)
specs <- candidate_designs()
params <- list(c_test_CRP = 16, c_test_NGS = 2518)

# ---------------------------------------------------------------------------
# 1. Design table
# ---------------------------------------------------------------------------
stopifnot(identical(names(specs), c("S1", "S2", "S3", "S3R", "S4")),
          identical(specs$S1$enrolled_cells, c("00", "01", "10", "11")),
          identical(specs$S2$enrolled_cells, c("10", "11")),     # CRP digit first: CRP-normal
          identical(specs$S3R$enrolled_cells, c("01", "11")),    # TMB/BRAF-positive
          identical(specs$S3$enrolled_cells, c("00", "01", "10", "11")),
          identical(specs$S3$experimental_cells, c("01", "11")),
          identical(specs$S3R$experimental_cells, specs$S3R$enrolled_cells),
          specs$S1$randomisation_week == 4, specs$S2$randomisation_week == 4,
          specs$S3R$randomisation_week == 0, specs$S3$randomisation_week == 0, specs$S4$p_experimental == 1,
          specs$S3$type == "stratified", is.na(specs$S3$p_experimental), is.na(specs$S3$power_target),
          specs$S3$n_proposal == 77, specs$S3R$n_proposal == 77, specs$S4$n_fixed == 80,
          abs(specs$S3$min_follow_up_years - 4) < 1e-12, specs$S3$startup_weeks == 4,
          grepl("^descriptive", specs$S3$sample_size_basis), grepl("^fixed", specs$S4$sample_size_basis),
          stratifying_biomarker(specs$S3) == "tmb_braf")
stopifnot(design_screening_cost(specs$S1, params) == 16 + 2518,
          design_screening_cost(specs$S2, params) == 16,
          design_screening_cost(specs$S3, params) == 2518,
          design_screening_cost(specs$S3R, params) == 2518,
          design_screening_cost(specs$S4, params) == 0)
expect_error_matching(design_screening_cost(specs$S2, list(c_test_NGS = 1)), "c_test_CRP")

built <- build_evsi_designs(specs, registered_per_year = 200, n_grid = c(80, 400), draw = draw,
                            component = component, reference_population = df_ref, include_power = FALSE)
v_targets <- lapply(built$designs, function(d) rownames(d$targets))
stopifnot(identical(v_targets$S1, c(paste0("os::", c("rx", "rx_x_crp", "rx_x_tmb_braf")),
                                    paste0("pfs::", c("rx", "rx_x_crp", "rx_x_tmb_braf")))),
          setequal(v_targets$S2, c("os::rx_x_tmb_braf", "os::cell_effect_10", "pfs::rx_x_tmb_braf",
                                   "pfs::cell_effect_10")),
          setequal(v_targets$S3R, c("os::rx_x_crp", "os::cell_effect_01", "pfs::rx_x_crp", "pfs::cell_effect_01")),
          identical(v_targets$S3, c(paste0("os::", c("ctl_base", "ctl_crp", "exp_tmb_braf", "rx_x_crp")),
                                    paste0("pfs::", c("ctl_base", "ctl_crp", "exp_tmb_braf", "rx_x_crp")))),
          setequal(v_targets$S4, c(paste0("os::", c("exp_base", "exp_crp", "exp_tmb_braf")),
                                   paste0("pfs::", c("exp_base", "exp_crp", "exp_tmb_braf")))),
          identical(sort(names(built$designs)),
                    sort(c("S1", "S2", "S3", "S3R", "S4", "S4+S1", "S4+S2", "S4+S3", "S4+S3R"))),
          length(v_targets$`S4+S1`) == 12, length(v_targets$`S4+S2`) == 10, length(v_targets$`S4+S3`) == 10,
          identical(names(built$designs$S1$covariances), c("80", "400")),
          identical(names(built$designs$S3R$covariances), c("77", "80", "400")),
          identical(names(built$designs$S3$covariances), c("77", "80", "400")))
# Enriched designs: the information matrices themselves are singular (rank 7 of 9).
stopifnot(built$information$S3R[["80"]]$endpoints$pfs$rank == 7,
          built$information$S1[["80"]]$endpoints$pfs$rank == 9)

df_accrual <- accrual_scenarios(c(norway = 500, nordic = 2000, europe = 45000), reference_rate = 20)
df_table <- design_table(specs, df_accrual, params, df_ref, NULL, component, draw,
                         power = list(norway = list(S1 = list(n = 300, events = 290),
                                                    S2 = list(n = 70, events = 69), S3R = list(n = 90, events = 88))))
row_s3 <- df_table[df_table$design == "S3" & df_table$accrual_scenario == "europe", ]
stopifnot(nrow(df_table) == 3 * 5 + 3, all(c("S6_S1", "S6_S2", "S6_S3R") %in% df_table$design),
          !"S6_S3" %in% df_table$design,
          all(df_table$screening_cost_per_registered[df_table$design == "S1"] == 2534),
          df_table$n_power[df_table$design == "S6_S1"] == 300,
          grepl("rx_x_crp", df_table$identified_targets[df_table$design == "S1"][1]),
          all(is.na(df_table$randomised_per_year[df_table$design %in% c("S3", "S4")])),
          # S3 keeps its own single-centre rate in every accrual scenario.
          abs(row_s3$registered_per_year - specs$S3$accrual_per_year) < 1e-12,
          abs(row_s3$min_follow_up_years - 4) < 1e-12, row_s3$tau_soc == TAU_SOC_DEFAULT,
          row_s3$tau_m3 == TAU_M3_DEFAULT, grepl("^descriptive", row_s3$sample_size_basis),
          all(is.na(df_table$tau_m3[df_table$type != "stratified"])), row_s3$experimental_cells == "01;11",
          abs(row_s3$p_experimental - mean(design_cell_codes(df_ref) %in% c("01", "11"))) < 1e-12,
          abs(row_s3$reporting_year_proposal - stratified_reporting_year(77, specs$S3)) < 1e-12,
          grepl("ctl_base", row_s3$identified_targets))
cat("1. designs: cells, randomisation week, screening tests, identified targets, design table\n")

# ---------------------------------------------------------------------------
# 2. Accrual and reporting-year arithmetic
# ---------------------------------------------------------------------------
rate_metimmox <- accrual_per_year(80, as.Date("2018-01-01"), as.Date("2022-01-01"))
stopifnot(abs(rate_metimmox - 80 / (1461 / 365.25)) < 1e-12)
expect_error_matching(accrual_per_year(80, as.Date("2020-01-01"), as.Date("2019-01-01")), "must follow")
stopifnot(identical(df_accrual$accrual_scenario, c("norway", "nordic", "europe")),
          max(abs(df_accrual$registered_per_year - 20 * c(1, 4, 90))) < 1e-12,
          max(abs(df_accrual$scale - c(1, 4, 90))) < 1e-12)
default_rate <- accrual_scenarios(c(norway = 500))$registered_per_year
stopifnot(abs(default_rate - accrual_per_year(METIMMOX_ACCRUAL_PATIENTS, METIMMOX_ACCRUAL_START,
                                              METIMMOX_ACCRUAL_END)) < 1e-12)
# Enrolment fractions by hand.
v_cell <- design_cell_codes(df_ref)
enrol_s3 <- design_enrolment_fraction(specs$S3R, df_ref, NULL, component, draw)
stopifnot(abs(enrol_s3$fraction - mean(v_cell %in% c("01", "11"))) < 1e-12)
v_sel <- stats::pgamma(4, shape = exp(draw$pfs[["shape"]]),
                       rate = exp(gamma_linear_predictor(component$original_pfs, draw$pfs, df_ref, "Control")),
                       lower.tail = FALSE)
enrol_s1 <- design_enrolment_fraction(specs$S1, df_ref, NULL, component, draw)
enrol_s2 <- design_enrolment_fraction(specs$S2, df_ref, NULL, component, draw)
stopifnot(abs(enrol_s1$fraction - mean(v_sel)) < 1e-12,
          abs(enrol_s2$fraction - mean(v_sel * (v_cell %in% c("10", "11")))) < 1e-12,
          abs(sum(enrol_s2$cell_mix_randomisation) - 1) < 1e-12)
# Reporting year: L / 52 + n / rate + F.
stopifnot(abs(design_reporting_year(100, 50, 2, 4) - (4 / 52 + 2 + 2)) < 1e-12,
          abs(design_reporting_year(c(50, 100), 25, 2) - c(4, 6)) < 1e-12)
meta_s1 <- built$meta[built$meta$design == "S1", ]
stopifnot(max(abs(meta_s1$randomised_per_year - 200 * enrol_s1$fraction)) < 1e-12,
          max(abs(meta_s1$reporting_year - (4 / 52 + meta_s1$n / (200 * enrol_s1$fraction) + 2))) < 1e-12,
          max(abs(meta_s1$n_registered - meta_s1$n / enrol_s1$fraction)) < 1e-9,
          all(built$meta$reporting_year[built$meta$design == "S4+S1"] >= s4_reporting_year()))
stopifnot(abs(s4_reporting_year(as.Date("2026-01-01"), as.Date("2028-01-01"), 0.5) - (730 / 365.25 + 0.5)) < 1e-12)
cat("2. accrual rates, enrolment fractions and reporting years agree with hand calculations\n")

# ---------------------------------------------------------------------------
# 3. Schoenfeld sizes
# ---------------------------------------------------------------------------
z <- stats::qnorm(0.975) + stats::qnorm(0.8)
d_main <- 4 * z^2 / log(0.5)^2
stopifnot(abs(schoenfeld_events(0.5, 0.8, 0.05, "main") - d_main) < 1e-10,
          abs(schoenfeld_events(0.5, 0.8, 0.05, "interaction") - 4 * d_main) < 1e-10,
          abs(d_main - 65.345) < 1e-3)
# Balanced event tables: the heuristic variance at the Schoenfeld total equals the target.
df_bal_int <- data.frame(arm = rep(c("control", "experimental"), each = 4), cell = rep(c("00", "01", "10", "11"), 2),
                         events = 4 * d_main / 8)
df_bal_main <- data.frame(arm = c("control", "experimental"), cell = "11", events = d_main / 2)
stopifnot(abs(power_target_variance(df_bal_int, specs$S1) - power_variance_target()) < 1e-12,
          abs(power_target_variance(df_bal_main, specs$S3R) - power_variance_target()) < 1e-12)
expect_error_matching(power_target_variance(df_bal_main, specs$S4), "no power target")
pw_s3 <- power_based_n(specs$S3R, 100, draw, component, df_ref)
v_at <- function(n) power_target_variance(design_expected_events(as_trial_design(specs$S3R, 100), n, draw, component,
                                                                 df_ref), specs$S3R)
stopifnot(is.finite(pw_s3$n), pw_s3$variance <= pw_s3$target_variance, v_at(pw_s3$n - 1) > pw_s3$target_variance,
          pw_s3$events >= d_main * 0.999)
cat("3. Schoenfeld events (", round(d_main, 2), " main, ", round(4 * d_main, 2),
    " interaction) and the smallest n meeting the target (S3R: ", pw_s3$n, ")\n", sep = "")

# ---------------------------------------------------------------------------
# 4. Heterogeneity-adjusted information
# ---------------------------------------------------------------------------
v_names <- rownames(component$original_pfs$res)
m_p <- solve(component$original_pfs$cov); dimnames(m_p) <- list(v_names, v_names)
m_d <- treatment_effect_directions(v_names, "pfs")
stopifnot(ncol(m_d) == 3, all(colSums(m_d) == 1))
stopifnot(identical(heterogeneity_adjusted_information(m_p, m_d, 0), m_p))
m_closed <- solve(solve(m_p) + 0.3^2 * m_d %*% t(m_d))
stopifnot(max(abs(heterogeneity_adjusted_information(m_p, m_d, 0.3) - m_closed)) < 1e-8 * max(abs(m_p)))
m_inf <- heterogeneity_adjusted_information(m_p, m_d, Inf)
stopifnot(max(abs(t(m_d) %*% m_inf %*% m_d)) < 1e-8 * max(abs(m_p)),
          max(abs(heterogeneity_adjusted_information(m_p, m_d, 1e4) - m_inf)) < 1e-4 * max(abs(m_p)))
expect_error_matching(heterogeneity_adjusted_information(m_p, m_d, -1), "non-negative")
cat("4. heterogeneity adjustment: tau = 0 identity, Woodbury closed form, tau = Inf removes the directions\n")

# ---------------------------------------------------------------------------
# 5. Combined-design information
# ---------------------------------------------------------------------------
info_s1 <- built$information$S1[["400"]]
l_info_s1 <- lapply(info_s1$endpoints, function(e) e$information)
stopifnot(max(abs(information_target_covariance(l_info_s1, info_s1$targets, component) - info_s1$covariance)) < 1e-10)
info_s1_ind <- design_information(as_trial_design(specs$S1, 200 * enrol_s1$fraction), 400, draw, component, df_ref,
                                  cross_endpoint = "independent",
                                  step = design_node_step(104 + 400 / (200 * enrol_s1$fraction) * 52))
stopifnot(max(abs(information_target_covariance(lapply(info_s1_ind$endpoints, function(e) e$information),
                                                info_s1_ind$targets, component, "independent") -
                    info_s1_ind$covariance)) < 1e-10)
info_s4 <- built$info_s4
comb <- combined_design_information(info_s4, info_s1, s4_targets(component, df_ref), component, "S4+S1")
for (outcome in c("os", "pfs")) {
  stopifnot(max(abs(comb$information[[outcome]] - info_s4$endpoints[[outcome]]$information -
                      info_s1$endpoints[[outcome]]$information)) < 1e-10 * max(abs(comb$information[[outcome]])))
}
v_shared <- rownames(info_s1$targets)
m_gap <- info_s1$covariance[v_shared, v_shared] - comb$covariance[v_shared, v_shared]
stopifnot(min(eigen((m_gap + t(m_gap)) / 2, symmetric = TRUE, only.values = TRUE)$values) > -1e-10,
          max(abs(comb$covariance - built$designs$`S4+S1`$covariances[["400"]])) < 1e-12)
# The single arm: rate and Rx are aliased, its own targets are identified.
pinv_s4 <- information_pseudo_inverse(info_s4$endpoints$pfs$information)
v_alias <- stats::setNames(numeric(length(v_names)), v_names)
v_alias[c("rate", "RxExperimental arm")] <- c(1, -1)
stopifnot(pinv_s4$rank < length(v_names),
          max(abs(info_s4$endpoints$pfs$information %*% v_alias)) < 1e-6 * max(abs(info_s4$endpoints$pfs$information)),
          info_s4$n_late < info_s4$n, info_s4$n_late > 0)
expect_error_matching(information_target_covariance(lapply(info_s4$endpoints, function(e) e$information),
                                                    info_s1$targets, component, label = "S4"),
                      "does not identify")
# A larger tau_m2 removes information in the late period.
info_s4_hi <- single_arm_information(80, draw, component, df_ref, tau_m2 = 1)
v_rx <- stats::setNames(numeric(length(v_names)), v_names); v_rx[c("rate", "RxExperimental arm")] <- 1
stopifnot(sum(v_rx * (info_s4$endpoints$pfs$information %*% v_rx)) >
            sum(v_rx * (info_s4_hi$endpoints$pfs$information %*% v_rx)))
cat("5. combined information sums, reproduces V_d(n), shrinks the shared covariance; S4 aliases rate and Rx\n")


# ---------------------------------------------------------------------------
# 6. Power from the engine (ruling (4) of 11 October 2026)
# ---------------------------------------------------------------------------
# Exponential case (gamma shape 1 on PFS): the engine's standard error and
# Schoenfeld's events give power-based sizes within 10% of each other.
draw_exp <- draw
# A smaller reference population (all four cells present) keeps the searches fast.
df_small <- df_ref[seq_len(80), ]
stopifnot(length(unique(design_cell_codes(df_small))) == 4)
draw_exp$pfs[["shape"]] <- 0
l_ps <- list()
for (id in c("S3R", "S1")) {
  spec <- specs[[id]]
  df_ps <- design_power_sizes(spec, 100, draw_exp, component, df_small)
  l_ps[[id]] <- df_ps
  for (b in unique(df_ps$biomarker)) {
    n_eng <- df_ps$n[df_ps$biomarker == b & df_ps$method == "engine"]
    n_sch <- df_ps$n[df_ps$biomarker == b & df_ps$method == "schoenfeld"]
    if (!is.finite(n_eng) || !is.finite(n_sch) || abs(n_eng / n_sch - 1) > 0.10) {
      stop(id, " ", b, ": engine n ", n_eng, " and Schoenfeld n ", n_sch, " differ by more than 10%")
    }
    cat(sprintf("6. %s (%s, %s): engine n %d, Schoenfeld n %d (ratio %.3f)\n", id, spec$power_target, b,
                as.integer(n_eng), as.integer(n_sch), n_eng / n_sch))
  }
  stopifnot(setequal(df_ps$method, c("engine", "schoenfeld")), all(df_ps$hypothesis == "conventional"),
            all(abs(df_ps$effect - log(POWER_RATE_RATIO)) < 1e-12))
}
# The engine's n is the smallest meeting the variance target.
tg_s3 <- design_power_targets(specs$S3R, component, df_small, NULL, draw_exp)[[1]]
stopifnot(identical(names(design_power_targets(specs$S1, component, df_small, NULL, draw)), c("crp", "tmb_braf")),
          tg_s3$biomarker == "tmb_braf", tg_s3$policy == "0101",
          abs(tg_s3$target[1, "pfs::RxExperimental arm:tmb_braf1"] - 1) < 1e-12,
          abs(tg_s3$target[1, "pfs::RxExperimental arm"] - 1) < 1e-12,
          guided_policy_pattern("crp") == "0011")
trial_s3 <- as_trial_design(specs$S3R, 100)
pw_eng <- power_based_n_information(trial_s3, tg_s3$target, log(0.5), draw = draw_exp, component = component,
                                    reference_population = df_small)
v_var <- function(n) diag(design_information(trial_s3, n, draw_exp, component, df_small, targets = tg_s3$target,
                                             step = design_node_step(104 + n / 100 * 52))$covariance)
stopifnot(pw_eng$variance <= pw_eng$target_variance, v_var(pw_eng$n - 1) > pw_eng$target_variance,
          abs(pw_eng$se - sqrt(pw_eng$variance)) < 1e-12)
expect_error_matching(power_based_n_information(trial_s3, tg_s3$target, 0, draw = draw_exp, component = component,
                                                reference_population = df_small), "non-zero")
pp <- primary_power_size(l_ps$S1)
stopifnot(pp$method == "engine", pp$n == pp$n_engine, is.finite(pp$n_schoenfeld))
cat("6. engine-based power: smallest n meeting the variance target; primary n is the engine's, max over S1's interactions\n")

# ---------------------------------------------------------------------------
# 7. Decision-relevant effect
# ---------------------------------------------------------------------------
for (f in c("prediction_functions", "prediction_population", "calculate_outcomes", "parameter_distributions",
            "voi_policies")) {
  source(here::here("R", paste0(f, ".R")))
}
horizon_test <- 104
params_test <- c(list(time_horizon = horizon_test, cl = 1 / 52, dr_costs = 0.04, dr_effects = 0.04, u_np = 0.73,
                      u_p = 0.59, u_decrement = 0.14, c_drug_nivo = 13923, c_drug_FLOX = 100, c_test_CT = 150,
                      c_test_blood = 20, c_test_CRP = 16, c_test_NGS = 2518, c_other_visit = 213.46,
                      c_other_baseline = 213.46, c_other_follow = 213.46, c_other_pp = 0, c_other_last = 5000),
                 build_treatment_schedules(horizon_test + 1))
df_pop <- df_small
v_shift <- c(os = "RxExperimental arm:crp1", pfs = "RxExperimental arm:crp1")
v_prices <- c(2000, 4000)
df_dre <- decision_relevant_effect("0011", v_prices, 1e5, component, params_test, df_pop, v_shift)
stopifnot(nrow(df_dre) == 2, all(df_dre$method == "uniroot"),
          df_dre$effect[2] < df_dre$effect[1],          # a larger price needs a larger benefit
          all(abs(df_dre$effect - (df_dre$effect_at_mean + df_dre$shift)) < 1e-12))
# At the root the incremental NMB is zero (re-evaluated directly).
par_root <- component_mean_draw(component)
for (o in c("os", "pfs")) par_root[[o]][v_shift[[o]]] <- par_root[[o]][v_shift[[o]]] + df_dre$shift[1]
inmb_root <- policy_incremental_nmb(guided_policy_cell_outcomes(par_root, component, params_test, df_pop), "0011",
                                    v_prices[1], 1e5, params_test)
stopifnot(abs(inmb_root) < 50)
# The price slope: cost linear in the price with slope the administrations.
co0 <- guided_policy_cell_outcomes(component_mean_draw(component), component, params_test, df_pop)
params_hi <- params_test; params_hi$c_drug_nivo <- 20000
co_hi <- guided_policy_cell_outcomes(component_mean_draw(component), component, params_hi, df_pop)
stopifnot(abs(policy_incremental_nmb(co0, "0011", 20000, 1e5, params_test) -
                policy_incremental_nmb(co_hi, "0011", 20000, 1e5, params_hi)) < 1e-6)
# No root at a prohibitive price: NA with the reason.
df_none <- decision_relevant_effect("nivo_0011", 1e6, 1e5, component, params_test, df_pop, v_shift)
stopifnot(is.na(df_none$effect), df_none$method == "no_root", grepl("not cost-effective", df_none$reason))
cat(sprintf("7. decision-relevant CRP effect: rate ratio %.3f at EUR %s and %.3f at EUR %s; none at EUR 1,000,000\n",
            df_dre$rate_ratio[1], format(v_prices[1]), df_dre$rate_ratio[2], format(v_prices[2])))
# Through the design layer: both hypotheses and both methods for S2 (CRP-enriched).
df_dec <- design_decision_effects(specs["S2"], v_prices, 1e5, component, params_test, df_pop)
df_ps2 <- design_power_sizes(specs$S2, 100, draw, component, df_small, decision_effects = df_dec)
d_dec <- df_ps2[df_ps2$hypothesis == "decision_relevant", ]
stopifnot(nrow(df_dec) == 2, all(df_dec$design == "S2"), all(df_dec$biomarker == "crp"),
          max(abs(df_dec$shift - df_dre$shift)) < 1e-3,
          setequal(paste(df_ps2$hypothesis, df_ps2$method),
                   c("conventional engine", "conventional schoenfeld", "decision_relevant engine",
                     "decision_relevant schoenfeld")),
          setequal(d_dec$price, v_prices),
          # n only for a benefit; a smaller benefit than rate ratio 0.5 needs more patients
          all(is.finite(d_dec$n) == (d_dec$effect < 0)),
          all(grepl("not a benefit", d_dec$reason[d_dec$effect >= 0])),
          any(is.finite(d_dec$n)),
          all(d_dec$n[is.finite(d_dec$n) & d_dec$effect > log(0.5)] >
                max(df_ps2$n[df_ps2$hypothesis == "conventional"])))
cat("7. design layer: S2 power-based n under the conventional and the decision-relevant hypothesis, both methods\n")

# ---------------------------------------------------------------------------
# 8. METIMMOX-3 as planned (S3: stratified, non-randomised)
# ---------------------------------------------------------------------------
s3 <- specs$S3
alloc_s3 <- design_allocation(s3, df_ref, NULL, component, draw)
stopifnot(abs(alloc_s3$p_experimental - mean(v_cell %in% c("01", "11"))) < 1e-12,
          identical(names(alloc_s3$experimental_mix), c("01", "11")),
          abs(sum(alloc_s3$experimental_mix) - 1) < 1e-12,
          alloc_s3$enrolment$fraction == 1)
# Reporting year: at the planned n, the planned last inclusion plus the minimum follow-up, i.e. the
# protocol's data cut-off.
stopifnot(abs(stratified_reporting_year(77, s3) -
                (as.numeric(difftime(METIMMOX3_ACCRUAL_END, VOI_DECISION_DATE, units = "days")) / 365.25 +
                   METIMMOX3_MIN_FOLLOW_UP_YEARS)) < 1e-12,
          abs(stratified_reporting_year(77, s3) -
                as.numeric(difftime(METIMMOX3_DATA_CUTOFF, VOI_DECISION_DATE, units = "days")) / 365.25) < 1e-12,
          abs(s3$accrual_per_year - accrual_per_year(77, METIMMOX3_ACCRUAL_START, METIMMOX3_ACCRUAL_END)) < 1e-12,
          all(diff(stratified_reporting_year(c(40, 77, 200), s3)) > 0))
expect_error_matching(stratified_reporting_year(77, specs$S1), "stratified designs")
info_s3 <- built$information$S3[["77"]]
# Arms by cell: experimental only in the TMB/BRAF-positive cells, control only in the negative ones.
df_ev <- info_s3$endpoints$pfs$expected_events
stopifnot(setequal(paste(df_ev$arm, df_ev$cell), c("experimental 01", "experimental 11", "control 00", "control 10")),
          all(df_ev$period[df_ev$arm == "experimental"] == "whole"),
          setequal(df_ev$period[df_ev$arm == "control"], c("startup", "soc")),
          abs(info_s3$n_experimental + info_s3$n_control - 77) < 1e-12,
          info_s3$endpoints$pfs$n_at_risk_soc < info_s3$n_control)
# Rank 7 of 9: Rx, the TMB/BRAF main effect and Rx:tmb_braf enter only as one sum.
pinv_s3 <- information_pseudo_inverse(info_s3$endpoints$pfs$information)
v_alias3 <- stats::setNames(numeric(length(v_names)), v_names)
v_alias3[c("RxExperimental arm", "tmb_braf1")] <- c(1, -1)
stopifnot(pinv_s3$rank == 7, information_pseudo_inverse(info_s3$endpoints$os$information)$rank == 7,
          max(abs(info_s3$endpoints$pfs$information %*% v_alias3)) <
            1e-6 * max(abs(info_s3$endpoints$pfs$information)))
# The CRP interaction is identified; the treatment effect in the TMB/BRAF-positive cells is not.
stopifnot(all(is.finite(diag(info_s3$covariance))), all(diag(info_s3$covariance) > 0))
expect_error_matching(information_target_covariance(lapply(info_s3$endpoints, function(e) e$information),
                                                    built$designs$S3R$targets, component, label = "S3"),
                      "does not identify: .*cell_effect_01")
# With tau_soc = 0 the start-up and standard-of-care periods partition the control follow-up.
info_s3_0 <- stratified_design_information(s3, 77, draw, component, df_ref, tau_soc = 0)
d_ctl_full <- single_arm_trial_design("control_full", c("00", "10"), info_s3_0$min_follow_up_weeks,
                                      info_s3_0$accrual_weeks, arm = "control")
for (outcome in c("os", "pfs")) {
  e_full <- endpoint_design_information(d_ctl_full, info_s3_0$n_control, draw, component, df_ref, NULL, outcome,
                                        design_node_step(info_s3_0$min_follow_up_weeks + info_s3_0$accrual_weeks))
  m_parts <- info_s3_0$endpoints[[outcome]]$control_startup + info_s3_0$endpoints[[outcome]]$control_soc
  stopifnot(max(abs(m_parts - e_full$information)) < 1e-3 * max(abs(e_full$information)))
}
# A larger tau_soc removes information on the control log rate, never on the experimental cells.
v_rate <- stats::setNames(numeric(length(v_names)), v_names); v_rate["rate"] <- 1
q_rate <- function(info) sum(v_rate * (info$endpoints$pfs$information %*% v_rate))
# (tau_soc = Inf removes the control log rate entirely: with the start-up under the same bias nothing
# identifies the control baseline any more.)
info_s3_inf <- stratified_design_information(s3, 77, draw, component, df_ref, tau_soc = TAU_SOC_RANGE[["high"]])
expect_error_matching(stratified_design_information(s3, 77, draw, component, df_ref, tau_soc = Inf),
                      "does not identify: .*ctl_base")
stopifnot(q_rate(info_s3_0) > q_rate(info_s3) && q_rate(info_s3) > q_rate(info_s3_inf),
          max(abs(info_s3_0$endpoints$pfs$experimental - info_s3_inf$endpoints$pfs$experimental)) == 0,
          max(abs(info_s3$covariance - built$designs$S3$covariances[["77"]])) < 1e-12)
# tau_soc acts on the start-up and the standard of care together (one shared bias of the control regimen).
for (outcome in c("os", "pfs")) {
  e <- info_s3$endpoints[[outcome]]
  m_shared <- heterogeneity_adjusted_information(e$control_startup + e$control_soc,
                                                 control_regimen_directions(v_names), TAU_SOC_DEFAULT)
  stopifnot(max(abs(e$control - m_shared)) < 1e-10 * max(abs(m_shared)),
            max(abs(e$information - e$experimental - e$control)) < 1e-10 * max(abs(e$information)))
}
stopifnot(sum(v_rate * (info_s3$endpoints$pfs$control %*% v_rate)) <
            sum(v_rate * ((info_s3$endpoints$pfs$control_startup +
                             heterogeneity_adjusted_information(info_s3$endpoints$pfs$control_soc,
                                                                control_regimen_directions(v_names),
                                                                TAU_SOC_DEFAULT)) %*% v_rate)))
# tau_m3: the regimen heterogeneity of the experimental cells (tislelizumab, bevacizumab).
info_m3_0 <- stratified_design_information(s3, 77, draw, component, df_ref, tau_m3 = 0)
info_m3_hi <- stratified_design_information(s3, 77, draw, component, df_ref, tau_m3 = TAU_M3_RANGE[["high"]])
v_int <- stats::setNames(numeric(length(v_names)), v_names); v_int[v_names[grepl(":", v_names) & grepl("crp", v_names)]] <- 1
q_int <- function(info) sum(v_int * (info$endpoints$pfs$information %*% v_int))
stopifnot(max(abs(info_m3_0$endpoints$pfs$experimental - info_m3_0$endpoints$pfs$experimental_unadjusted)) == 0,
          q_int(info_m3_0) > q_int(info_s3) && q_int(info_s3) > q_int(info_m3_hi),
          max(abs(info_m3_0$endpoints$pfs$control - info_m3_hi$endpoints$pfs$control)) == 0,
          info_s3$tau_m3 == TAU_M3_DEFAULT,
          all(diag(info_m3_hi$covariance) > diag(info_s3$covariance) - 1e-12),
          information_pseudo_inverse(info_m3_hi$endpoints$pfs$information)$rank == 7)
expect_error_matching(stratified_design_information(s3, 77, draw, component, df_ref, tau_m3 = Inf),
                      "does not identify: .*rx_x_crp")
# Accrual-independent designs and the cost table's experimental share.
stopifnot(identical(accrual_independent_designs(specs), c("S4", "S3", "S4+S3")))
df_cost_s3 <- design_cost_table(specs, params, c(S3 = 5, S3R = 6), "norway", c("S3", "S4+S3", "S3R"),
                                p_experimental = c(S3 = alloc_s3$p_experimental))
stopifnot(identical(df_cost_s3$cost_design, c("S3", "S3", "S3R")),
          all(df_cost_s3$p_experimental[1:2] == alloc_s3$p_experimental), df_cost_s3$p_experimental[3] == 0.5)
expect_error_matching(design_cost_table(specs, params, c(S3 = 5), "norway", "S3"), "No experimental share")
cat(sprintf("8. S3 (protocol v1.1): %.1f experimental and %.1f control patients at n = 77, rank 7 of 9, reports in %.2f years\n",
            info_s3$n_experimental, info_s3$n_control, stratified_reporting_year(77, s3)))

cat("\nAll design tests passed.\n")
