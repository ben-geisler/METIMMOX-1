# ===============================================================================
# TEST: research value, ENBS and the research threshold price (issue #203,
#       paper 3 plan step 6)
# ===============================================================================
# Synthetic inputs only (no trial data, no caches):
#
#   1. population_discount_factor(): start year 1 reproduces the annuity of
#      calculate_population_evppi(); a later start drops the early years; a
#      start after the horizon gives zero. design_value_start_year() rounds
#      the reporting year up and starts the year after.
#   2. research_value_table(): population value = EVSI x eligible x factor;
#      a design reporting after the horizon has zero value; the EVPI rows use
#      start year 1.
#   3. Trial-cost bands: the Norwegian METIMMOX band (NOK 50 million; #213,
#      proposed main-text default), the oncology band (Moore et al. 2020)
#      and the Speich 2018 range in 2023 euros against a hand conversion. ENBS and n* on a synthetic EVSI
#      curve (Speich band) against a brute-force calculation written out in
#      the test (low, mid and high trial cost; nivolumab donated and charged); the share of the cost range with
#      ENBS > 0; n* is NA when no n pays; S4 is costed at S4_RESEARCH_COST.
#   4. Research threshold price on a synthetic price axis against a
#      brute-force search; NA when no price pays.
#   5. drop_dominated_policies() leaves EVPI and EVSI unchanged; the
#      sampling-importance resample has ESS n under equal weights.
#   6. run_evsi_grid() and sequential_increments() on scalar synthetic
#      designs: zero EVSI where EVPI is zero, "reports_after_horizon" rows,
#      increments equal EVSI(S4 + S1) - EVSI(S4).
#   7. select_external_evidence() keeps one ITT row per trial and endpoint
#      (synthetic table; the real table when data/external is present).
#   8. CSV writers: the five files are written and read back with their rows.
#   9. Live (skips without the trial data and unless the cache directory
#      holds a valid psa_pair_v2 with the cell component, #205): the policy
#      NMB at the list price has the reference policy first.
#  10. End to end: run_voi_research_value() on a synthetic component, every
#      candidate (S3 as planned and S3R included), combined and incremental
#      designs, E0 and E2; the accrual-independent designs (S4, S3) are
#      identical across accrual scenarios.
#  11. Script 14 grids, the METIMMOX-only switch (evidence states: E0
#      always, E2 optional, from the setting or METIMMOX_VOI14_STATES) and
#      the runtime estimate.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_research_value.R
# ===============================================================================

suppressPackageStartupMessages(library(survival))
for (f in c("model_configs", "evppi_functions", "joint_survival_sampling", "pfs_os_violation_diagnostics",
            "voi_survival_parameters", "design_information", "evsi_functions", "voi_price_axis",
            "evidence_update", "voi_designs", "voi_research_value")) {
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

# ---------------------------------------------------------------------------
# 1. Discount factor and start year
# ---------------------------------------------------------------------------
for (rate in c(0.035, 0.04)) for (h in c(5, 10, 15)) {
  evpi <- 123.4
  stopifnot(abs(evpi * 519 * population_discount_factor(1, h, rate) / 1e6 -
                  calculate_population_evppi(evpi, 519, h, rate)) < 1e-12)
}
stopifnot(abs(population_discount_factor(4, 10, 0.04) - sum(1.04^-(4:10))) < 1e-14,
          population_discount_factor(11, 10, 0.04) == 0, population_discount_factor(10, 10, 0) == 1,
          identical(population_discount_factor(c(1, 11), 10, 0.04) > 0, c(TRUE, FALSE)))
expect_error_matching(population_discount_factor(1.5, 10, 0.04), "whole years")
stopifnot(identical(design_value_start_year(c(0, 2.3, 3, 9.0000001)), c(1, 4, 4, 11)))
cat("1. start year 1 reproduces calculate_population_evppi(); later starts drop years; after the horizon 0\n")

# ---------------------------------------------------------------------------
# 2. Research value table
# ---------------------------------------------------------------------------
v_n <- c(40, 80, 160, 320, 640, 1280)
evsi_curve <- function(n, price) pmax(0, 8000 - price) / 8000 * 1000 * (1 - exp(-n / 200))
v_prices <- c(0, 2000, 4000, 6000, 8000, 10000)
rate_rand <- 100
meta <- data.frame(design = "S1", base_design = "S1", role = "standalone", n = v_n,
                   n_role = c("grid", "grid", "power", "grid", "grid", "grid"),
                   reporting_year = design_reporting_year(v_n, rate_rand), enrolment_fraction = 0.9,
                   randomised_per_year = rate_rand, registered_per_year = rate_rand / 0.9,
                   n_registered = v_n / 0.9, evaluated = TRUE, accrual_scenario = "norway", stringsAsFactors = FALSE)
meta <- rbind(meta, data.frame(design = "S4", base_design = "S4", role = "standalone", n = 80, n_role = "fixed",
                               reporting_year = 1.8, enrolment_fraction = 1, randomised_per_year = NA,
                               registered_per_year = NA, n_registered = 80, evaluated = TRUE,
                               accrual_scenario = "norway", stringsAsFactors = FALSE))
df_evsi <- do.call(rbind, lapply(v_prices, function(p) rbind(
  data.frame(design = "S1", n = v_n, n_role = meta$n_role[1:6], price = p, evidence_state = "E0",
             accrual_scenario = "norway", evsi = evsi_curve(v_n, p), mcse = 1, method = "gam_regression",
             evpi = 1500, error = "", stringsAsFactors = FALSE),
  data.frame(design = "S4", n = 80, n_role = "fixed", price = p, evidence_state = "E0", accrual_scenario = "norway",
             evsi = evsi_curve(30, p), mcse = 1, method = "gam_regression", evpi = 1500, error = "",
             stringsAsFactors = FALSE))))
df_evpi <- data.frame(price = v_prices, evidence_state = "E0", evpi = 1500)
populations <- population_scenarios(c(norway = 500))
stopifnot(identical(populations$accrual_scenario, "norway"))
df_value <- research_value_table(df_evsi, meta, df_evpi, populations, horizons = c(5, 10), rate = 0.04)
r <- df_value[df_value$design == "S1" & df_value$price == 2000 & df_value$horizon == 10 & df_value$n == 160, ]
stopifnot(nrow(r) == 1, r$value_start_year == design_value_start_year(160 / 100 + 2),
          abs(r$max_justifiable_cost - evsi_curve(160, 2000) * 500 * sum(1.04^-(5:10))) < 1e-8)
r_late <- df_value[df_value$design == "S1" & df_value$n == 1280 & df_value$horizon == 10, ]
stopifnot(all(r_late$reporting_year > 10), all(r_late$max_justifiable_cost == 0), all(r_late$discount_factor == 0))
r_evpi <- df_value[df_value$design == "EVPI" & df_value$horizon == 10, ]
stopifnot(nrow(r_evpi) == length(v_prices), all(r_evpi$value_start_year == 1),
          max(abs(r_evpi$max_justifiable_cost - 1500 * 500 * population_discount_factor(1, 10, 0.04))) < 1e-8)
cat("2. population value = EVSI x eligible x factor; reporting after the horizon gives 0; EVPI from year 1\n")

# ---------------------------------------------------------------------------
# 3. ENBS and n*
# ---------------------------------------------------------------------------
design_costs <- rbind(
  data.frame(design = "S1", accrual_scenario = "norway", cost_design = "S1", screening_cost = 2534,
             p_experimental = 0.5, administrations = 7.5, fixed_research_cost = NA_real_),
  data.frame(design = "S4", accrual_scenario = "norway", cost_design = "S4", screening_cost = 0,
             p_experimental = 1, administrations = 0, fixed_research_cost = S4_RESEARCH_COST))
# Cost bands: the Norwegian METIMMOX band (#213, proposed main-text default:
# NOK 50 million per 80 enrolled at face value to per 68 analysed inflated
# 2018-2023, at 10.90 NOK/EUR); USD bands inflated with the CPI-U to 2023,
# converted at the ECB 2023 rate; the Speich 2018 range is the band used for
# the brute force.
bands <- trial_cost_bands()
to_eur <- function(x, year) x * 304.702 / c("2017" = 245.120, "2019" = 255.657)[[as.character(year)]] / 1.0813
stopifnot(identical(names(bands), c("norway", "oncology", "speich")), TRIAL_COST_BAND_MAIN == "norway",
          abs(bands$norway$per_patient[["low"]] - 50e6 / 80 / 10.90) < 1e-8,
          abs(bands$norway$per_patient[["high"]] - 50e6 / 68 * 130.4 / 108.4 / 10.90) < 1e-8,
          abs(bands$norway$per_patient_point - 50e6 / 68 / 10.90) < 1e-8,
          bands$norway$per_patient[["low"]] < bands$norway$per_patient_point,
          bands$norway$per_patient_point < bands$norway$per_patient[["high"]],
          all(bands$norway$fixed == 0), bands$norway$price_year == 2023,
          grepl("NOK 50 million", bands$norway$source),
          abs(bands$oncology$per_patient[["low"]] - to_eur(80880, 2019)) < 1e-8,
          abs(bands$oncology$per_patient[["high"]] - to_eur(155714, 2019)) < 1e-8,
          all(bands$oncology$fixed == 0), bands$oncology$price_year == 2023,
          abs(bands$speich$per_patient[["low"]] - to_eur(43, 2017)) < 1e-10,
          abs(bands$speich$per_patient[["high"]] - to_eur(103254, 2017)) < 1e-8,
          abs(bands$speich$fixed[["high"]] - to_eur(2e5, 2017)) < 1e-8,
          grepl("10.1136/bmjopen-2020-038863", bands$oncology$source))
df_onc <- trial_cost_parameters(c(0, 0.5, 1), bands$oncology)
stopifnot(max(abs(df_onc$per_patient - bands$oncology$per_patient[["low"]] *
                    (bands$oncology$per_patient[["high"]] / bands$oncology$per_patient[["low"]])^c(0, 0.5, 1))) < 1e-8,
          all(df_onc$fixed == 0))
# The default band is the main-text (Norwegian) band.
df_nor <- trial_cost_parameters(c(0, 1))
stopifnot(max(abs(df_nor$per_patient - bands$norway$per_patient[c("low", "high")])) < 1e-8)
expect_error_matching(usd_to_eur2023(1, 2005), "No CPI-U value")
band_s <- bands$speich
enbs <- enbs_analysis(df_value, design_costs, u_grid = seq(0, 1, length.out = 21), band = band_s, cost_band = "speich")
enbs_onc <- enbs_analysis(df_value, design_costs, u_grid = seq(0, 1, length.out = 21))
# The default band is the main-text (Norwegian) band.
stopifnot(all(enbs$curves$cost_band == "speich"), all(enbs_onc$optimal$cost_band == "norway"),
          all(abs(enbs_onc$optimal$per_patient_cost[enbs_onc$optimal$u == 0] - bands$norway$per_patient[["low"]]) < 1e-8))
df_cost <- trial_cost_parameters(c(0, 0.5, 1), band_s)
stopifnot(abs(df_cost$per_patient[1] - to_eur(43, 2017)) < 1e-12,
          abs(df_cost$per_patient[3] - to_eur(103254, 2017)) < 1e-9,
          abs(df_cost$per_patient[2] - to_eur(sqrt(43 * 103254), 2017)) < 1e-9,
          df_cost$fixed[1] == 0, abs(df_cost$fixed[3] - to_eur(2e5, 2017)) < 1e-9)
for (price in c(0, 4000)) for (niv in c("donated", "charged")) for (j in 1:3) {
  v_mjc <- evsi_curve(v_n, price) * 500 * population_discount_factor(design_value_start_year(v_n / 100 + 2), 10, 0.04)
  v_cost <- df_cost$fixed[j] + df_cost$per_patient[j] * v_n + 2534 * v_n / 0.9 +
    (if (niv == "charged") v_n * 0.5 * 7.5 * price else 0)
  v_enbs <- v_mjc - v_cost
  n_star <- if (max(v_enbs) > 0) v_n[which.max(v_enbs)] else NA
  row <- enbs$optimal[enbs$optimal$design == "S1" & enbs$optimal$price == price & enbs$optimal$horizon == 10 &
                        enbs$optimal$nivolumab == niv & abs(enbs$optimal$u - df_cost$u[j]) < 1e-12, ]
  stopifnot(nrow(row) == 1, identical(row$n_star, as.numeric(n_star)),
            abs(row$enbs_star - max(max(v_enbs), 0)) < 1e-6)
  curve <- enbs$curves[enbs$curves$design == "S1" & enbs$curves$price == price & enbs$curves$horizon == 10 &
                         enbs$curves$nivolumab == niv & enbs$curves$cost_point == c("low", "mid", "high")[j], ]
  stopifnot(max(abs(curve$enbs - v_enbs)) < 1e-6, sum(curve$is_n_star) == as.integer(!is.na(n_star)))
}
# Share of the cost range with ENBS > 0 at n*(u), brute force.
v_u <- seq(0, 1, length.out = 21)
df_u <- trial_cost_parameters(v_u, band_s)
v_mjc0 <- evsi_curve(v_n, 0) * 500 * population_discount_factor(design_value_start_year(v_n / 100 + 2), 10, 0.04)
share0 <- mean(vapply(seq_along(v_u), function(j)
  max(v_mjc0 - df_u$fixed[j] - df_u$per_patient[j] * v_n - 2534 * v_n / 0.9) > 0, logical(1)))
row0 <- enbs$optimal[enbs$optimal$design == "S1" & enbs$optimal$price == 0 & enbs$optimal$horizon == 10 &
                       enbs$optimal$nivolumab == "donated", ]
stopifnot(all(abs(row0$share_cost_range_positive - share0) < 1e-12), share0 > 0, share0 < 1,
          all(row0$n_power == 160))
# No n pays at a price with zero EVSI; S4 costs nothing.
row_none <- enbs$optimal[enbs$optimal$design == "S1" & enbs$optimal$price == 10000, ]
row_s4 <- enbs$optimal[enbs$optimal$design == "S4" & enbs$optimal$price == 0 & enbs$optimal$horizon == 10, ]
stopifnot(all(is.na(row_none$n_star)), all(row_none$enbs_star == 0), all(row_none$share_cost_range_positive == 0),
          all(row_s4$trial_cost_star == 0), all(row_s4$n_star == 80))
cat("3. ENBS curves, n* and the share of the cost range with ENBS > 0 match a brute-force calculation\n")

# ---------------------------------------------------------------------------
# 4. Research threshold price
# ---------------------------------------------------------------------------
df_thr <- research_threshold_price(enbs$optimal)
for (niv in c("donated", "charged")) for (end in c("low", "high")) {
  u <- c(low = 0, high = 1)[[end]]
  cp <- trial_cost_parameters(u, band_s)
  v_pay <- vapply(v_prices, function(price) {
    v_mjc <- evsi_curve(v_n, price) * 500 *
      population_discount_factor(design_value_start_year(v_n / 100 + 2), 10, 0.04)
    v_c <- cp$fixed + cp$per_patient * v_n + 2534 * v_n / 0.9 + (if (niv == "charged") v_n * 3.75 * price else 0)
    mjc_s4 <- evsi_curve(30, price) * 500 * population_discount_factor(design_value_start_year(1.8), 10, 0.04)
    max(v_mjc - v_c) > 0 || mjc_s4 > 0
  }, logical(1))
  expected <- if (any(v_pay)) max(v_prices[v_pay]) else NA
  row <- df_thr[df_thr$nivolumab == niv & df_thr$cost_end == end & df_thr$horizon == 10, ]
  stopifnot(nrow(row) == 1, identical(row$threshold_price, as.numeric(expected)))
}
# Without the free pooling design the high end has no paying price on this axis.
df_thr_s1 <- research_threshold_price(enbs$optimal[enbs$optimal$design == "S1", ])
r_hi <- df_thr_s1[df_thr_s1$cost_end == "high" & df_thr_s1$nivolumab == "charged" & df_thr_s1$horizon == 5, ]
stopifnot(nrow(r_hi) == 1, is.na(r_hi$threshold_price) || r_hi$threshold_price < 8000,
          all(df_thr$prices_evaluated == length(v_prices)))
cat("4. research threshold price equals a brute-force search over the price grid\n")

# ---------------------------------------------------------------------------
# 5. Dominated policies and the importance resample
# ---------------------------------------------------------------------------
set.seed(3)
n_sim <- 3000
v_theta <- stats::rnorm(n_sim, -0.1, 0.3)
m_nmb <- cbind(nivo_0000 = 0, nivo_0011 = 1000 * (v_theta + stats::rnorm(n_sim, 0, 0.2)),
               nivo_0001 = 1000 * (v_theta + stats::rnorm(n_sim, 0, 0.2)))
m_nmb <- cbind(m_nmb, nivo_0111 = m_nmb[, "nivo_0011"] - 50, nivo_1111 = pmin(m_nmb[, 2], m_nmb[, 3]) - 1)
m_pruned <- drop_dominated_policies(m_nmb)
stopifnot(identical(colnames(m_pruned), c("nivo_0000", "nivo_0011", "nivo_0001")),
          setequal(attr(m_pruned, "dropped"), c("nivo_0111", "nivo_1111")),
          abs(calculate_evpi_from_nmb(m_pruned) - calculate_evpi_from_nmb(m_nmb)) < 1e-10)
scalar_design <- function(name, n_grid, sigma = 1) {
  evsi_design(name, matrix(1, 1, 1, dimnames = list("os::x", "os::x")),
              stats::setNames(lapply(n_grid, function(n) matrix(sigma^2 / n)), n_grid))
}
m_theta <- cbind("os::x" = v_theta)
e_full <- evsi_by_design(m_nmb, m_theta, scalar_design("d", 10), 10, methods = "gaussian_approximation", B = 50)
e_pruned <- evsi_by_design(m_pruned, m_theta, scalar_design("d", 10), 10, methods = "gaussian_approximation", B = 50)
stopifnot(abs(e_full$evsi - e_pruned$evsi) < 1e-8)
rs <- evidence_resample_index(rep(0, 500), seed = 1)
stopifnot(abs(rs$ess - 500) < 1e-8, length(rs$index) == 500)
rs2 <- evidence_resample_index(c(0, rep(-50, 99)), seed = 1)
stopifnot(rs2$ess < 1.01, all(rs2$index == 1))
expect_error_matching(evidence_resample_index(c(0, -Inf)), "finite")
cat("5. dominated policies dropped without changing EVPI or EVSI; resample ESS correct\n")

# ---------------------------------------------------------------------------
# 6. EVSI grid and sequential increments
# ---------------------------------------------------------------------------
m_cost_ref <- matrix(rep(c(0, 500, 500), each = n_sim), n_sim, dimnames = list(NULL, colnames(m_pruned)))
m_admin <- matrix(rep(c(0, 0.1, 0.1), each = n_sim), n_sim, dimnames = list(NULL, colnames(m_pruned)))
pairing <- structure(list(price_ref = 1000, price_alt = NA_real_, strategies = colnames(m_pruned), n_sim = n_sim,
                          cost_ref = m_cost_ref, effect = (m_pruned + m_cost_ref) / 51000,
                          administrations = m_admin), class = "voi_price_pairing")
l_nmb <- policy_nmb_at_prices(pairing, c(0, 1000, 1e6), 51000)
stopifnot(identical(names(l_nmb), c("0", "1000", "1e+06")) || identical(names(l_nmb), c("0", "1000", "1000000")))
stopifnot(max(abs(l_nmb[["1000"]] - m_pruned)) < 1e-6)
designs <- list(S4 = scalar_design("S4", 80, 4), `S4+S1` = scalar_design("S4+S1", c(40, 160), 1))
meta6 <- data.frame(design = c("S4", "S4+S1", "S4+S1", "S4+S1"), n = c(80, 40, 160, 4000),
                    n_role = c("fixed", "grid", "grid", "grid"), evaluated = c(TRUE, TRUE, TRUE, FALSE),
                    stringsAsFactors = FALSE)
df_grid <- run_evsi_grid(l_nmb, m_theta, designs, meta6, "E0", "norway", methods = "gam_regression", B = 50)
stopifnot(all(df_grid$method[df_grid$n == 4000] == "reports_after_horizon"), all(is.na(df_grid$evsi[df_grid$n == 4000])),
          all(df_grid$evsi[as.character(df_grid$price) == "1e+06" & df_grid$n != 4000] == 0),
          all(df_grid$method[df_grid$price == 1e6 & df_grid$n != 4000] == "zero_evpi_bound"),
          all(is.finite(df_grid$evsi[df_grid$price == 1000 & df_grid$n != 4000])))
df_inc <- sequential_increments(df_grid)
for (n in c(40, 160)) {
  e_c <- df_grid$evsi[df_grid$design == "S4+S1" & df_grid$n == n & df_grid$price == 1000]
  e_4 <- df_grid$evsi[df_grid$design == "S4" & df_grid$price == 1000]
  stopifnot(abs(df_inc$evsi[df_inc$n == n & df_inc$price == 1000] - (e_c - e_4)) < 1e-12)
}
stopifnot(all(df_inc$design == "S4->S1"), all(design_role(c("S1", "S4+S1", "S4->S1")) ==
                                                c("standalone", "combined", "increment")))
cat("6. EVSI grid: zero bound, reports-after-horizon rows; increments EVSI(S4+S1) - EVSI(S4)\n")

# ---------------------------------------------------------------------------
# 7. External-evidence selection
# ---------------------------------------------------------------------------
df_lik <- data.frame(label = paste("row", 1:6), trial = c("A", "A", "A", "B", "B", "C"),
                     population = c("ITT", "ITT", "pMMR subgroup", "ITT, BRAF wild-type", "ITT", "ITT"),
                     endpoint = c("pfs", "pfs", "pfs", "os", "pfs", "os"), log_hr = -0.1, se = 0.2, z = 1.96,
                     p_experimental = c(2 / 3, 2 / 3, 2 / 3, 0.5, 0.5, NA),
                     follow_up_weeks = c(50, 100, 100, 80, 80, 60), maintenance = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE),
                     source_row = 1:6, stringsAsFactors = FALSE)
sel <- select_external_evidence(df_lik)
stopifnot(identical(sel$label, c("row 2", "row 4", "row 5")), all(sel$landmark_weeks == c(0, 16, 16)))
sel_nm <- select_external_evidence(df_lik, maintenance_landmark_weeks = NA)
stopifnot(identical(sel_nm$label, "row 2"))
if (file.exists(here::here(EXTERNAL_EVIDENCE_PATH))) {
  sel_live <- suppressMessages(select_external_evidence(external_evidence_likelihood(read_external_evidence())))
  stopifnot(nrow(sel_live) > 0, !anyDuplicated(paste(sel_live$trial, sel_live$endpoint)),
            all(grepl("^ITT", sel_live$population)))
  cat("7. evidence selection: one ITT row per trial and endpoint (", nrow(sel_live), " rows of the table)\n", sep = "")
} else {
  cat("7. evidence selection: one ITT row per trial and endpoint (synthetic; SKIP: evidence table absent)\n")
}

# ---------------------------------------------------------------------------
# 8. CSV writers
# ---------------------------------------------------------------------------
out_dir <- file.path(tempdir(), "voi_research_test")
results <- list(evsi = df_grid, designs = data.frame(design = c("S1", "S4"), n_fixed = c(NA, 80)),
                research_value = df_value, enbs = enbs$curves, threshold = df_thr)
v_paths <- write_voi_research_csvs(results, out_dir)
stopifnot(identical(unname(basename(v_paths)), c("evsi_by_design.csv", "voi_designs.csv", "voi_research_value.csv",
                                                 "enbs_by_design.csv", "voi_threshold_price.csv")),
          all(file.exists(v_paths)),
          nrow(utils::read.csv(v_paths[["research_value"]])) == nrow(df_value),
          nrow(utils::read.csv(v_paths[["enbs"]])) == nrow(enbs$curves),
          all(c("design", "n", "price", "evidence_state", "evsi", "mcse", "method") %in%
                names(utils::read.csv(v_paths[["evsi"]]))))
expect_error_matching(write_voi_research_csvs(results[1:3], out_dir), "results lacks")
unlink(out_dir, recursive = TRUE)
cat("8. the five CSV writers write and read back\n")

# ---------------------------------------------------------------------------
# 9. Live check (after the regeneration with the cell component)
# ---------------------------------------------------------------------------
source(here::here("R", "cache_paths.R"))
#' The bound PSA object of the cache directory (option metimmox.cache_dir,
#' set by tests/run_all.R), or NULL with a SKIP line when the trial data or the
#' PSA pair are absent or the pair is not a valid psa_pair_v2 with the cell
#' component (#205: a stale pair must skip, not run).
read_live_psa_obj <- function() {
  skip <- function(why) {
    cat("SKIP: live policy NMB (", why, ")\n", sep = "")
    NULL
  }
  if (!file.exists(here::here("data", "tidy", "METIMMOX.rds"))) return(skip("confidential trial data absent"))
  v_files <- c(psa_obj_path("correct"), psa_params_path("correct"))
  if (!all(file.exists(v_files))) {
    return(skip(paste0(paste(basename(v_files), collapse = " or "), " absent in ", cache_dir())))
  }
  po_live <- readRDS(v_files[1])
  problem <- tryCatch({ validate_psa_pair(po_live, readRDS(v_files[2])); NULL },
                      error = function(e) conditionMessage(e))
  if (!is.null(problem)) return(skip(problem))
  po_live
}
psa_live <- read_live_psa_obj()
if (!is.null(psa_live)) {
  source(here::here("R", "voi_policies.R"))
  space <- policy_space_from_psa(psa_live, list(c_test_CRP = psa_live$cells$test_costs[["crp"]],
                                                c_test_NGS = psa_live$cells$test_costs[["tmb_braf"]]))
  l_live <- policy_nmb_at_prices(policy_price_pairing(space), space$price_ref, 51000)
  stopifnot(colnames(l_live[[1]])[1] == "nivo_0000", nrow(l_live[[1]]) == space$n_sim)
  cat("9. live: policy NMB at the list price keeps", ncol(l_live[[1]]), "undominated policies\n")
}

# ---------------------------------------------------------------------------
# 10. End to end: run_voi_research_value() on a synthetic component and PSA
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(flexsurv))
set.seed(11)
n_pat <- 300
df_pat <- data.frame(Age = stats::rnorm(n_pat, 62, 9), sex = factor(stats::rbinom(n_pat, 1, 0.45)),
                     Rx = factor(stats::rbinom(n_pat, 1, 0.5), labels = c("Control", "Experimental arm")),
                     crp = factor(stats::rbinom(n_pat, 1, 0.35)), tmb_braf = factor(stats::rbinom(n_pat, 1, 0.45)))
v_eta <- log(0.02) + 0.01 * (df_pat$Age - 62) - 0.2 * (df_pat$Rx != "Control") -
  0.3 * (df_pat$crp == "1") * (df_pat$Rx != "Control")
v_t <- stats::rgamma(n_pat, 1.6, exp(v_eta)); v_cens <- stats::runif(n_pat, 30, 150)
df_pat$OSwk <- pmin(v_t, v_cens); df_pat$Death <- as.numeric(v_t <= v_cens)
v_t2 <- pmin(v_t, stats::rgamma(n_pat, 1.2, exp(v_eta + 0.8)))
df_pat$PFSwk <- pmin(v_t2, v_cens); df_pat$Prog <- as.numeric(v_t2 <= v_cens)
fit_os <- flexsurvreg(Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx, data = df_pat, dist = "gamma")
fit_pfs <- flexsurvreg(Surv(PFSwk, Prog) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx, data = df_pat, dist = "gamma")
m_joint <- as.matrix(Matrix::bdiag(fit_os$cov, fit_pfs$cov))
v_joint <- c(paste0("os::", rownames(fit_os$res)), paste0("pfs::", rownames(fit_pfs$res)))
dimnames(m_joint) <- list(v_joint, v_joint)
n_d <- 2000
component <- list(method = "mvn_joint_v2", joint_covariance = list(covariance = m_joint),
                  draws = draw_joint_survival_coefficients(fit_os, fit_pfs, m_joint, n_d, 11L),
                  original_os = fit_os, original_pfs = fit_pfs, n_samples = n_d, dist_os = "gamma",
                  dist_pfs = "gamma", seed = 11L)
m_draws <- joint_draw_matrix(component)
model_idx <- sample.int(n_d)
theta_psa <- m_draws[model_idx, , drop = FALSE]
# Three policies; the CRP-guided one gains with a more negative CRP-positive cell effect.
v_cell <- theta_psa[, "pfs::RxExperimental arm"] + theta_psa[, "pfs::RxExperimental arm:crp1"]
v_qaly_gain <- -0.4 * v_cell + stats::rnorm(n_d, 0, 0.02)
m_eff <- cbind(nivo_0000 = 1, nivo_0011 = 1 + v_qaly_gain, nivo_1111 = 1 + 0.5 * v_qaly_gain - 0.05)
m_admin <- cbind(nivo_0000 = 0, nivo_0011 = 2.5, nivo_1111 = 7)
m_admin <- matrix(rep(m_admin, each = n_d), n_d, dimnames = list(NULL, colnames(m_eff)))
m_cost <- 13923 * m_admin + cbind(0, 16, 0)[rep(1, n_d), ]
pairing_e2e <- structure(list(price_ref = 13923, price_alt = NA_real_, strategies = colnames(m_eff), n_sim = n_d,
                              cost_ref = m_cost, effect = m_eff, administrations = m_admin),
                         class = "voi_price_pairing")
cells <- list(values = array(7, dim = c(n_d, 4, 2, 1),
                             dimnames = list(NULL, c("00", "01", "10", "11"), c("control", "experimental"),
                                             "nivo_admin")))
rs <- evidence_state_rows(-0.5 * ((theta_psa[, "pfs::RxExperimental arm"] + 0.2) / 0.5)^2, seq_len(n_d), 11L,
                          min_ess = 100)
expect_error_matching(evidence_state_rows(rep(0, 10), 11, 1), "outside")
settings_e2e <- list(wtp = 51000, prices = c(0, 4000, 13923), n_grid = c(80, 240), horizons = c(10, 5),
                     primary_horizon = 10, rate = 0.04, methods = "gaussian_approximation",
                     primary_method = "gaussian_approximation", seed = 11L, B = 50, tau_m2 = TAU_M2_DEFAULT,
                     include_power = FALSE, accrual = c("norway", "europe"), reference_rate = 20)
res <- run_voi_research_value(pairing_e2e, theta_psa,
                              list(E0 = list(component = component, rows = NULL),
                                   E2 = list(component = component, rows = rs$index)),
                              candidate_designs(), df_pat[, c("Age", "sex", "Rx", "crp", "tmb_braf")],
                              params = list(c_test_CRP = 16, c_test_NGS = 2518), cells = cells,
                              eligible = c(norway = 500, europe = 45000), settings = settings_e2e, verbose = FALSE)
stopifnot(all(c("S1", "S2", "S3", "S3R", "S4", "S4+S1", "S4->S1", "S4->S3", "S4->S3R") %in% res$evsi$design),
          all(c("E0", "E2") %in% res$evsi$evidence_state),
          all(res$evsi$evsi[res$evsi$price == 13923] == 0 | is.na(res$evsi$evsi[res$evsi$price == 13923])),
          any(res$evsi$evsi[res$evsi$price == 0] > 0, na.rm = TRUE),
          all(res$evpi$evpi[res$evpi$price == 0] > 0),
          # S4 and the single-centre S3 are computed once per state and copied to the other accrual scenario.
          identical(res$evsi$evsi[res$evsi$design == "S4" & res$evsi$accrual_scenario == "norway"],
                    res$evsi$evsi[res$evsi$design == "S4" & res$evsi$accrual_scenario == "europe"]),
          identical(res$evsi$evsi[res$evsi$design == "S3" & res$evsi$accrual_scenario == "norway"],
                    res$evsi$evsi[res$evsi$design == "S3" & res$evsi$accrual_scenario == "europe"]),
          all(res$research_value$max_justifiable_cost >= 0 | !is.finite(res$research_value$evsi_per_patient) |
                res$research_value$evsi_per_patient < 0),
          all(res$enbs$horizon == 10), nrow(res$threshold) > 0,
          all(c("S1", "S2", "S3", "S3R", "S4") %in% res$designs$design),
          all(res$designs$administrations_per_experimental[res$designs$type %in% c("rct", "stratified")] == 7),
          # S3 pays per-patient costs for every enrolled patient and the drug for its experimental share.
          all(res$design_costs$p_experimental[res$design_costs$cost_design == "S3"] ==
                res$designs$p_experimental[res$designs$design == "S3"][1]))
# S3 reports at the same time in both accrual scenarios (single centre at its own rate).
ry3 <- res$meta[res$meta$design == "S3" & res$meta$n == 77 & res$meta$evidence_state == "E0", ]
stopifnot(nrow(ry3) == 2, length(unique(ry3$reporting_year)) == 1)
# Europe accrues faster than Norway, so a design reports earlier there.
ry <- res$meta[res$meta$design == "S1" & res$meta$n == 240 & res$meta$evidence_state == "E0", ]
stopifnot(ry$reporting_year[ry$accrual_scenario == "europe"] < ry$reporting_year[ry$accrual_scenario == "norway"])
out_dir <- file.path(tempdir(), "voi_research_e2e")
v_paths <- write_voi_research_csvs(res, out_dir)
stopifnot(all(file.exists(v_paths)))
unlink(out_dir, recursive = TRUE)
cat("10. end to end: EVSI of S1, S2, S3, S3R, S4, combined and incremental designs under E0 and E2, research value, ENBS,\n",
    "    threshold price and the design table on a synthetic component\n", sep = "")

# A state with its own PSA (pairing and theta; R/voi_evidence_psa.R): given E0's
# own PSA it reproduces E0 exactly; rows and an own pairing are exclusive.
settings_own <- modifyList(settings_e2e, list(prices = 0, n_grid = 80, accrual = "norway", reference_rate = 500))
specs_own <- candidate_designs()[c("S2", "S4")]
run_own <- function(states) {
  run_voi_research_value(pairing_e2e, theta_psa, states, specs_own, df_pat[, c("Age", "sex", "Rx", "crp", "tmb_braf")],
                         params = list(c_test_CRP = 16, c_test_NGS = 2518), cells = cells,
                         eligible = c(norway = 500), settings = settings_own, verbose = FALSE)
}
res_own <- run_own(list(E0 = list(component = component, rows = NULL),
                        E2 = list(component = component, pairing = pairing_e2e, theta = theta_psa)))
ev0 <- res_own$evsi[res_own$evsi$evidence_state == "E0", c("design", "n", "price", "evsi", "method")]
ev2 <- res_own$evsi[res_own$evsi$evidence_state == "E2", c("design", "n", "price", "evsi", "method")]
rownames(ev0) <- rownames(ev2) <- NULL
stopifnot(nrow(ev0) > 0, identical(ev0, ev2), any(ev0$evsi > 0, na.rm = TRUE),
          identical(res_own$evpi$evpi[1], res_own$evpi$evpi[2]))
expect_error_matching(run_own(list(E0 = list(component = component, rows = NULL),
                                   E2 = list(component = component, rows = 1:10, pairing = pairing_e2e,
                                             theta = theta_psa))), "either rows or its own pairing")
expect_error_matching(run_own(list(E0 = list(component = component, rows = NULL),
                                   E2 = list(component = component, pairing = pairing_e2e,
                                             theta = theta_psa[1:10, ]))), "one row per PSA draw")
cat("10. a state with its own PSA reproduces E0 when given E0's PSA; rows and an own pairing are exclusive\n")

# Every cost band enters the ENBS and the threshold price; no power sizes without include_power.
stopifnot(setequal(unique(res$enbs$cost_band), c("norway", "oncology", "speich")),
          setequal(unique(res$threshold$cost_band), c("norway", "oncology", "speich")),
          is.null(res$power_sizes), is.null(res$decision_effects), identical(names(res$cost_bands), names(bands)))

# ---------------------------------------------------------------------------
# 11. Script 14 grid (ruling (10) of 11 October 2026) and runtime estimate
# ---------------------------------------------------------------------------
marked <- c(list = 13923, biosimilar = 4641, adoption_threshold = 2500)
grid_red <- voi14_grid_settings("reduced", marked, seq(0, 13923, by = 2000))
grid_full <- voi14_grid_settings("full", marked, seq(0, 13923, by = 2000))
stopifnot(identical(grid_red$prices, c(2500, 4641, 13923)), identical(grid_red$n_grid, c(40, 80, 160, 320, 640, 1280)),
          identical(grid_red$accrual, "norway"), identical(grid_red$decision_prices, grid_red$prices),
          identical(grid_full$n_grid, VOI_N_GRID), identical(grid_full$accrual, c("norway", "nordic", "europe")),
          all(marked %in% grid_full$prices), length(grid_full$prices) > length(grid_red$prices))
expect_error_matching(voi14_grid_settings("medium", marked, 0), "should be one of")
# METIMMOX-only switch (12 October 2026): E0 always, E2 optional, canonical order.
stopifnot(identical(grid_red$evidence_states, c("E0", "E2")), identical(grid_full$evidence_states, c("E0", "E2")),
          identical(voi14_grid_settings("reduced", marked, 0, evidence_states = "E0")$evidence_states, "E0"),
          identical(voi14_evidence_states(c("E2", "E0")), c("E0", "E2")),
          identical(voi14_evidence_states(NULL), c("E0", "E2")),
          identical(voi14_requested_states(""), c("E0", "E2")),
          identical(voi14_requested_states("E0"), "E0"),
          identical(voi14_requested_states(" E0 , E2 "), c("E0", "E2")))
expect_error_matching(voi14_evidence_states("E2"), "must include E0")
expect_error_matching(voi14_evidence_states(c("E0", "E1")), "Unknown evidence state")
expect_error_matching(voi14_requested_states("E1"), "Unknown evidence state")
v_elig <- c(norway = 519, norway_low = 400, norway_high = 650, nordic = 2000, europe = 45000)
stopifnot(identical(names(voi14_eligible(v_elig, "norway")), c("norway", "norway_low", "norway_high")),
          identical(voi14_eligible(v_elig, c("norway", "nordic", "europe")), v_elig))
rt <- voi14_runtime_estimate(pairing_e2e, theta_psa, candidate_designs(), component,
                             df_pat[, c("Age", "sex", "Rx", "crp", "tmb_braf")], params = NULL,
                             settings = c(settings_e2e[setdiff(names(settings_e2e), "accrual")], grid_red["accrual"],
                                          list(grid = "reduced", decision_prices = 4000)),
                             verbose = FALSE)
stopifnot(identical(rt$component, c("design information", "power-based sizes", "EVSI evaluations",
                                    "decision-relevant effects")),
          all(rt$seconds_each[1:3] > 0), rt$count[1] == 2 * 1 * 4 * 2, rt$count[2] == 0, rt$count[4] == 0,
          rt$count[3] == 2 * (1 * 2 * 4 * 2 + 1) * 3, is.finite(attr(rt, "total_minutes")))
cat(sprintf("11. script 14 grids (reduced: %d prices, n %s, Norway; full: %d prices); runtime estimate %.2f min\n",
            length(grid_red$prices), paste(range(grid_red$n_grid), collapse = "-"), length(grid_full$prices),
            attr(rt, "total_minutes")))

cat("\nAll research-value tests passed.\n")
