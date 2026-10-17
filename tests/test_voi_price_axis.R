# ===============================================================================
# TEST: per-draw nivolumab price pairing and price-axis EVPI (issue #199)
# ===============================================================================
# Paper 3 evaluates EVPI at any nivolumab price from the two scenario PSAs
# (list and biosimilar price), which must agree per draw in everything except
# the cost, and differ in cost by the price difference times the per-draw
# administrations N_i: C_i(p) = C_i(p0) + (p - p0) N_i.
#
#   1. Synthetic pair with known N_i: the pairing passes, N_i is recovered and
#      costs at a third price are reproduced.
#   2. Every kind of mismatch stops with a message naming the check: QALYs,
#      model_idx, a non-price parameter, negative N_i, nivolumab cost in the
#      control strategy, and the published mean administrations.
#   3. EVPI is 0 at a price where one strategy has the highest NMB in every
#      draw, equals a hand-computed value on a 3-draw example, and the
#      QALY-scale EVPI equals its hand-computed value.
#   4. The price grid contains the marked prices, sorted and unique.
#   5. Live (when the scenario cache is present and current): the pairing holds
#      on the scenario cache, the mean administrations match
#      outputs/tables/cea_threshold_price.csv, and the EVPI at the biosimilar
#      price and EUR 51,000 equals Table S7 within 1e-6. Skips when the cache
#      is absent or predates the cache fingerprints of #167.
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_voi_price_axis.R
# ===============================================================================

source(here::here("R", "model_configs.R"))
source(here::here("R", "evppi_functions.R"))
source(here::here("R", "voi_price_axis.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error))) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Synthetic PSA pair at one price from costs excluding nivolumab, QALYs and N.
make_pair <- function(m_cost_other, m_effect, m_admin, price, model_idx, df_params) {
  m_cost <- m_cost_other + price * m_admin
  psa_obj <- list(strategies = colnames(m_cost), n_sim = nrow(m_cost),
                  cost = as.data.frame(m_cost), effectiveness = as.data.frame(m_effect),
                  model_idx = model_idx)
  list(psa_obj = psa_obj, psa_params = df_params)
}

# ---------------------------------------------------------------------------
# 1. Known N_i reproduces C(p)
# ---------------------------------------------------------------------------
set.seed(1)
n_sim <- 50L
v_strategies <- c("control", "crp", "tmb_braf")
m_effect <- cbind(control = runif(n_sim, 1.3, 1.5), crp = runif(n_sim, 1.3, 1.6),
                  tmb_braf = runif(n_sim, 1.2, 1.6))
m_cost_other <- cbind(control = runif(n_sim, 2e4, 3e4), crp = runif(n_sim, 2e4, 3e4),
                      tmb_braf = runif(n_sim, 2e4, 3e4))
m_admin <- cbind(control = 0, crp = runif(n_sim, 1, 3.5), tmb_braf = runif(n_sim, 1, 4.5))
v_model_idx <- sample.int(200L, n_sim)
df_params <- data.frame(sim = seq_len(n_sim), u_np = runif(n_sim), c_other_last = rgamma(n_sim, 25),
                        model_idx = v_model_idx)
price_list <- 13923
price_bio <- 4641
pair_list <- make_pair(m_cost_other, m_effect, m_admin, price_list, v_model_idx, df_params)
pair_bio <- make_pair(m_cost_other, m_effect, m_admin, price_bio, v_model_idx, df_params)

pairing <- check_price_pairing(pair_list, pair_bio, price_list, price_bio,
                               control_strategy = "control", max_administrations = 8)
stopifnot(inherits(pairing, "voi_price_pairing"), all(pairing$checks$passed),
          max(abs(pairing$administrations - m_admin)) < 1e-9,
          max(abs(costs_at_price(pairing, 1000) - (m_cost_other + 1000 * m_admin))) < 1e-7,
          max(abs(costs_at_price(pairing, 0) - m_cost_other)) < 1e-7)
expect_error_matching(costs_at_price(pairing, -1), "non-negative")

# The comparison with the published mean administrations passes on agreement.
df_threshold <- data.frame(
  analysis = c("deterministic", "probabilistic", "probabilistic"),
  strategy = c("crp", "crp", "tmb_braf"),
  discounted_administrations = c(99, mean(m_admin[, "crp"]), mean(m_admin[, "tmb_braf"])),
  threshold_price = c(1000, 800, NA))
pairing_table <- check_price_pairing(pair_list, pair_bio, price_list, price_bio,
                                     control_strategy = "control", df_threshold = df_threshold)
stopifnot(nrow(pairing_table$administrations_vs_table) == 2,
          identical(adoption_threshold_price(df_threshold),
                    list(price = 800, strategy = "crp")))
cat("PASS: known administrations reproduce costs at any price\n")

# ---------------------------------------------------------------------------
# 2. Mismatches are detected
# ---------------------------------------------------------------------------
pair_bad <- pair_bio
pair_bad$psa_obj$effectiveness$crp[7] <- pair_bad$psa_obj$effectiveness$crp[7] + 1e-4
expect_error_matching(check_price_pairing(pair_list, pair_bad, price_list, price_bio),
                      "QALYs identical per draw")

pair_bad <- pair_bio
pair_bad$psa_obj$model_idx[3] <- pair_bad$psa_obj$model_idx[3] + 1L
pair_bad$psa_params$model_idx[3] <- pair_bad$psa_params$model_idx[3] + 1L
expect_error_matching(check_price_pairing(pair_list, pair_bad, price_list, price_bio),
                      "model_idx identical.*first: draw 3")

pair_bad <- pair_bio
pair_bad$psa_params$u_np[5] <- 0.5
expect_error_matching(check_price_pairing(pair_list, pair_bad, price_list, price_bio),
                      "Non-price parameters.*u_np")

m_admin_negative <- m_admin
m_admin_negative[2, "crp"] <- -0.5
expect_error_matching(check_price_pairing(
  make_pair(m_cost_other, m_effect, m_admin_negative, price_list, v_model_idx, df_params),
  make_pair(m_cost_other, m_effect, m_admin_negative, price_bio, v_model_idx, df_params),
  price_list, price_bio), "non-negative")

m_admin_control <- m_admin
m_admin_control[4, "control"] <- 0.25
pair_ctl_list <- make_pair(m_cost_other, m_effect, m_admin_control, price_list, v_model_idx, df_params)
pair_ctl_bio <- make_pair(m_cost_other, m_effect, m_admin_control, price_bio, v_model_idx, df_params)
expect_error_matching(check_price_pairing(pair_ctl_list, pair_ctl_bio, price_list, price_bio,
                                          control_strategy = "control"),
                      "Control strategy receives no nivolumab")
expect_error_matching(check_price_pairing(pair_list, pair_bio, price_list, price_bio,
                                          max_administrations = 2),
                      "at most the scheduled number")

df_threshold_bad <- df_threshold
df_threshold_bad$discounted_administrations[2] <- df_threshold_bad$discounted_administrations[2] + 0.01
expect_error_matching(check_price_pairing(pair_list, pair_bio, price_list, price_bio,
                                          control_strategy = "control",
                                          df_threshold = df_threshold_bad),
                      "cea_threshold_price.csv")
expect_error_matching(check_price_pairing(pair_list, pair_bio, price_list, price_list),
                      "two different")
cat("PASS: QALY, model_idx, parameter, administration and table mismatches stop\n")

# ---------------------------------------------------------------------------
# 3. EVPI: dominance, a hand-computed 3-draw example, QALY scale
# ---------------------------------------------------------------------------
# Three draws, two strategies (A has no nivolumab). At price p:
#   cost_A = 100, cost_B = 50 + 10 p; QALYs A = (1, 1, 1), B = (0.5, 1.5, 1.2).
m_effect_3 <- cbind(A = c(1, 1, 1), B = c(0.5, 1.5, 1.2))
m_other_3 <- cbind(A = c(100, 100, 100), B = c(50, 50, 50))
m_admin_3 <- cbind(A = c(0, 0, 0), B = c(10, 10, 10))
df_params_3 <- data.frame(sim = 1:3, model_idx = 1:3)
pairing_3 <- check_price_pairing(
  make_pair(m_other_3, m_effect_3, m_admin_3, 20, 1:3, df_params_3),
  make_pair(m_other_3, m_effect_3, m_admin_3, 5, 1:3, df_params_3),
  20, 5, control_strategy = "A")
# WTP 1000, price 10: NMB_A = 900 each; NMB_B = (350, 1350, 1050).
# mean max = (900 + 1350 + 1050) / 3 = 1100; max mean = max(900, 916.67) = 916.67.
evpi_3 <- evpi_at_price(pairing_3, 10, 1000)
stopifnot(abs(evpi_3$evpi - (1100 - 2750 / 3)) < 1e-9,
          identical(evpi_3$optimal_strategy, "B"),
          isTRUE(all.equal(unname(evpi_3$prob_optimal), c(1 / 3, 2 / 3))))
# Price 200: B costs 2050 > every NMB gain, so A is optimal in every draw.
evpi_dominant <- evpi_at_price(pairing_3, 200, 1000)
stopifnot(evpi_dominant$evpi == 0, evpi_dominant$prob_optimal[["A"]] == 1)
# QALY scale: mean max = (1 + 1.5 + 1.2) / 3; max mean = max(1, 3.2 / 3).
stopifnot(abs(qaly_scale_evpi(pairing_3) - (3.7 / 3 - 3.2 / 3)) < 1e-12)

df_axis_3 <- run_price_axis_evpi(pairing_3, c(10, 200), c(1000, 2000),
                                 list(annual_incidence = 100, research_horizon = 2,
                                      discount_rate = 0),
                                 marked_prices = c(list_price = 200))
stopifnot(nrow(df_axis_3) == 4,
          abs(df_axis_3$evpi_per_patient[1] - evpi_3$evpi) < 1e-9,
          abs(df_axis_3$evpi_population_millions[1] - evpi_3$evpi * 200 / 1e6) < 1e-12,
          identical(df_axis_3$price_marker, c("", "list_price", "", "list_price")),
          all(abs(df_axis_3$evpi_qaly_scale - 0.5 / 3) < 1e-12),
          all(abs(df_axis_3$evpi_qaly_scale_population - 0.5 / 3 * 200) < 1e-9),
          abs(df_axis_3$evpi_qaly_equivalent[1] - evpi_3$evpi / 1000) < 1e-12)
cat("PASS: EVPI is zero under dominance and matches hand-computed values\n")

# ---------------------------------------------------------------------------
# 4. Price grid
# ---------------------------------------------------------------------------
v_grid <- price_axis_grid(13923, c(4641, 723.63, 4641), step = 250)
stopifnot(v_grid[1] == 0, utils::tail(v_grid, 1) == 13923, !is.unsorted(v_grid),
          !anyDuplicated(v_grid), all(c(4641, 723.63, 13750) %in% v_grid),
          length(v_grid) == length(seq(0, 13923, by = 250)) + 3)
expect_error_matching(price_axis_grid(13923, 20000), "within")
cat("PASS: price grid contains the marked prices\n")

# ---------------------------------------------------------------------------
# 5. Live check on the scenario cache
# ---------------------------------------------------------------------------
source(here::here("R", "cache_paths.R"))
scenario_file <- scenario_evppi_path("correct")
threshold_file <- here::here("outputs", "tables", "cea_threshold_price.csv")
s7_file <- here::here("outputs", "tables", "table_s7.csv")
if (!file.exists(scenario_file) || !file.exists(threshold_file) || !file.exists(s7_file)) {
  cat("SKIP: live price pairing (scenario cache or published tables absent)\n")
} else {
  scenario_cache <- readRDS(scenario_file)
  if (is.null(scenario_cache$fingerprint)) {
    cat("SKIP: live price pairing (", basename(scenario_file),
        " predates the #167 cache fingerprints; rerun analysis/12_scenario_EVPPIs.R)\n", sep = "")
  } else {
    df_threshold_live <- utils::read.csv(threshold_file)
    df_s7 <- utils::read.csv(s7_file)
    pairing_live <- scenario_price_pairing(scenario_cache, control_strategy = get_control_strategy(),
                                           df_threshold = df_threshold_live)
    stopifnot(all(pairing_live$checks$passed))
    df_axis_live <- run_price_axis_evpi(
      pairing_live, c(pairing_live$price_alt, pairing_live$price_ref), 51000,
      list(annual_incidence = 1, research_horizon = 1, discount_rate = 0))
    evpi_s7 <- df_s7$evpi[df_s7$scenario_id == "s1_decreased_nivo" & df_s7$group == "evpi"]
    evpi_list <- df_axis_live$evpi_per_patient[df_axis_live$c_drug_nivo == pairing_live$price_ref]
    evpi_bio <- df_axis_live$evpi_per_patient[df_axis_live$c_drug_nivo == pairing_live$price_alt]
    # Both prices against Table S7 (the base-price EVPI was 0 before #210).
    evpi_s7_base <- df_s7$evpi[df_s7$scenario_id == "base" & df_s7$group == "evpi"]
    stopifnot(length(evpi_s7) == 1, abs(evpi_bio - evpi_s7) < 1e-6,
              length(evpi_s7_base) == 1, abs(evpi_list - evpi_s7_base) < 1e-6)
    cat(sprintf("PASS: live pairing holds; EVPI at EUR %s and EUR 51,000 = %.6f (Table S7 %.6f)\n",
                format(pairing_live$price_alt, big.mark = ","), evpi_bio, evpi_s7))
  }
}

cat("PASS: test_voi_price_axis.R\n")
