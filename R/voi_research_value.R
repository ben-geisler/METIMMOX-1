# ===============================================================================
# RESEARCH VALUE: MAXIMUM JUSTIFIABLE COST, ENBS, RESEARCH THRESHOLD PRICE
# (issue #203, paper 3 plan step 6)
# ===============================================================================
# Inputs are the per-patient EVSI by design, n, price and evidence state from
# evsi_by_design() (R/evsi_functions.R) and the price-axis EVPI from the policy
# pairing of step 1b (policy_price_pairing(), evpi_at_price()). Nothing here
# re-estimates a value of information.
#
#   * population_discount_factor(start_year, horizon, rate): sum over years
#     y = start_year..horizon of (1 + rate)^-y; start_year = 1 reproduces the
#     annuity of calculate_population_evppi() (EVPI: perfect information now).
#   * A design reports `reporting_year` years after the research decision
#     (randomisation landmark + n / randomised rate + minimum follow-up;
#     R/voi_designs.R). The reporting year is ROUNDED UP to whole years and the
#     value starts the year after (design_value_start_year() = ceiling(R) + 1),
#     so patients treated before the report contribute nothing and a design
#     reporting after the horizon has zero value. Rounding up is conservative
#     (it drops the remainder of the reporting year) and keeps the annual
#     annuity convention of calculate_population_evppi().
#   * Population EVSI = per-patient EVSI x eligible patients per year x the
#     factor: the maximum justifiable research cost of the design at n.
#   * ENBS(n) = population EVSI(n) - trial cost(n). No trial budgets exist
#     (decision of 9 October 2026), so the trial cost is a scenario over a
#     cost band (trial_cost_bands()). Main text (PROPOSED in #213, revising
#     ruling (8) of 11 October 2026; to be confirmed by the user): the
#     Norwegian anchor, METIMMOX's own cost (NOK 50 million, nivolumab
#     donated; trial team, 15 October 2026), per patient, in 2023 euros at the
#     model's 10.90 NOK/EUR. Appendix: the oncology-specific band of Moore et
#     al. 2020 (the interquartile range for oncology pivotal trials) and the
#     full range of Speich et al. 2018, both converted to 2023 euros (US CPI-U
#     to 2023, then the ECB 2023 average rate). TRIAL_COST_BAND_MAIN selects
#     the main-text band. Each band is
#     indexed by u in [0, 1]: per-patient cost geometric between its ends,
#     fixed cost linear; u = 0 is the low end, 1 the high end. Nivolumab is donated (as in METIMMOX) or charged (administrations
#     per experimental patient x price). Screening tests are charged for every
#     registered patient (design_screening_cost()).
#   * ENBS-optimal n*: argmax over the evaluated n; NA when no n has ENBS > 0
#     (no trial). A probability that ENBS > 0 needs a distribution over the
#     trial cost, which does not exist; the share of the cost range (share of
#     the u grid) for which ENBS(n*(u)) > 0 is reported instead.
#   * Research threshold price: the highest price on the EVSI price grid at
#     which the best standalone design at its ENBS-optimal n has ENBS > 0,
#     i.e. its maximum justifiable cost exceeds the trial cost at the low and
#     at the high end of the cost range.
#   * The sequential value of a new trial after the METIMMOX-2 pooling (S4):
#     EVSI(S4 + Sx, n) - EVSI(S4), valued from the combined report year, with
#     the cost of Sx only (S4 is funded already).
#
# Requires R/evppi_functions.R (calculate_evpi_from_nmb(),
# calculate_population_evppi() for the check), R/voi_price_axis.R
# (evpi_at_price(), costs_at_price(), voi_output_dir()) and R/voi_designs.R.
# ===============================================================================

#' Source of the appendix trial-cost range (Speich et al. 2018)
TRIAL_COST_SOURCE <- paste0("Speich B et al. Systematic review on costs and resource use of randomized clinical ",
                            "trials. J Clin Epidemiol 2018;96:1-11. doi:10.1016/j.jclinepi.2017.12.018 ",
                            "(overall costs USD 43-103,254 per patient and USD 0.2-611.5 million per RCT, USD 2017)")
#' Overall RCT cost per randomised patient, USD 2017 (Speich et al. 2018,
#' abstract: range of the 16 articles reporting overall costs). TO VERIFY
#' against the full text before use in the manuscript.
TRIAL_COST_PER_PATIENT_USD2017 <- c(low = 43, high = 103254)
#' Fixed (set-up) cost per trial, USD 2017. Low end 0; high end the lowest
#' overall RCT cost reported by Speich et al. 2018 (USD 0.2 million), used as
#' a ceiling on the set-up cost of an academic cooperative-group trial. TO
#' VERIFY; ASSUMPTION.
TRIAL_COST_FIXED_USD2017 <- c(low = 0, high = 2e5)
#' EUR per USD 2017 (annual average reference rate 2017, about 1.13 USD per
#' EUR). TO VERIFY. Superseded for the conversion by usd_to_eur2023() (#203,
#' 11 October 2026: inflate in USD, convert at the 2023 rate); kept for the
#' provenance record of voi_research_parameters().
EUR_PER_USD2017 <- 1 / 1.1297
#' Oncology-specific estimated cost per patient of pivotal trials, USD (Moore
#' TJ, Heyward J, Anderson G, Alexander GC. Variation in the estimated costs
#' of pivotal clinical benefit trials supporting the US approval of new
#' therapeutic agents, 2015-2017: a cross-sectional study. BMJ Open
#' 2020;10:e038863. doi:10.1136/bmjopen-2020-038863, Table 3, oncology, 39
#' trials: median USD 100,271, IQR 80,880-155,714). Excludes the study drug
#' (sponsor costs for design, monitoring and medication are not included), so
#' nivolumab is added separately (donated or charged). "Current US dollars"
#' of the cost-estimating tool, taken as price year 2019 (the tool's manual
#' is dated 2019; the paper states no year). Industry pivotal trials, so the
#' band is high for an academic trial. PLACEHOLDER, TO VERIFY; an appendix
#' band since the Norwegian figure arrived (#213, proposed).
TRIAL_COST_ONCOLOGY_PER_PATIENT_USD2019 <- c(low = 80880, median = 100271, high = 155714)
#' Price year assumed for TRIAL_COST_ONCOLOGY_PER_PATIENT_USD2019. ASSUMPTION.
TRIAL_COST_ONCOLOGY_PRICE_YEAR <- 2019
#' US CPI-U, all items, U.S. city average, annual averages (1982-84 = 100;
#' BLS series CUUR0000SA0): 2017 245.120, 2019 255.657, 2023 304.702. TO
#' VERIFY against the BLS supplemental file.
US_CPI_U_ANNUAL <- c("2017" = 245.120, "2019" = 255.657, "2023" = 304.702)
#' USD per EUR, ECB euro reference rate, annual average 2023. TO VERIFY
#' against the ECB Data Portal (EXR.A.USD.EUR.SP00.A).
USD_PER_EUR_2023 <- 1.0813
#' Price year of every trial cost in the ENBS
TRIAL_COST_PRICE_YEAR <- 2023
#' METIMMOX's total trial cost, NOK, with nivolumab donated by BMS (trial
#' team, 15 October 2026, COST1). Taken as nominal NOK (the trial ran 2018 to
#' 2023; no cost year was given).
METIMMOX_TRIAL_COST_NOK <- 50e6
#' Patients the trial team divided METIMMOX's cost by (the analysed cohort;
#' about NOK 735,294 per patient) and patients enrolled (the low end of the
#' Norwegian band)
METIMMOX_TRIAL_COST_PATIENTS <- c(analysed = 68, enrolled = 80)
#' NOK per EUR of the model's costs (#36: the 2023 somatic-care unit price of
#' NOK 49,484 is EUR 4,541.60 in analysis/05_basecase_input_parameters.R)
NOK_PER_EUR_MODEL <- 10.90
#' Norwegian consumer price index, all items, annual averages (2015 = 100;
#' Statistics Norway, table 03013): 2018 108.4, 2023 130.4. TO VERIFY against
#' the SSB table.
NO_CPI_ANNUAL <- c("2018" = 108.4, "2023" = 130.4)
#' Cost band of the main text: "norway" (PROPOSED in #213, revising ruling (8)
#' of 11 October 2026, which had "oncology"; to be confirmed by the user)
TRIAL_COST_BAND_MAIN <- "norway"
#' Grid points of the trial-cost index u in [0, 1]
TRIAL_COST_GRID_POINTS <- 101L
#' Nivolumab options of the trial cost
TRIAL_NIVOLUMAB_OPTIONS <- c("donated", "charged")
#' Research cost of the METIMMOX-2 pooling (S4): the trial is funded and
#' running, so the incremental cost is taken as zero. ASSUMPTION.
S4_RESEARCH_COST <- 0

#' Present value of one patient-cohort per year from start_year to horizon
#'
#' @param start_year First year whose patients benefit (1 = from now; an
#'   integer >= 1; values above `horizon` give 0).
#' @param horizon Research horizon (years).
#' @param rate Annual discount rate.
#' @return sum_{y = start_year}^{horizon} (1 + rate)^-y (vectorised in
#'   start_year).
population_discount_factor <- function(start_year, horizon, rate) {
  if (!is.numeric(horizon) || length(horizon) != 1 || !is.finite(horizon) || horizon < 0) {
    stop("horizon must be one non-negative number")
  }
  if (!is.numeric(rate) || length(rate) != 1 || !is.finite(rate) || rate <= -1) stop("rate must exceed -1")
  if (any(!is.finite(start_year)) || any(start_year < 1) || any(start_year != round(start_year))) {
    stop("start_year must be whole years >= 1")
  }
  vapply(start_year, function(s) {
    if (s > horizon) return(0)
    sum((1 + rate)^-(seq(s, horizon)))
  }, numeric(1))
}

#' First year of benefit of a design that reports after `reporting_year`
#'
#' @param reporting_year Years from the decision to the report (continuous).
#' @return ceiling(reporting_year) + 1 (whole years; 1 for a report now).
design_value_start_year <- function(reporting_year) {
  if (any(!is.finite(reporting_year)) || any(reporting_year < 0)) stop("reporting_year must be non-negative")
  ceiling(reporting_year - 1e-9) + 1
}

#' Population scenarios of the research value
#'
#' @param eligible Named eligible patients per year (norway first, then the
#'   other populations).
#' @param accrual Named accrual scenario of each population (default: Norway's
#'   bounds use the Norwegian accrual, every other population its own).
#' @return Data frame with `population`, `eligible_per_year`, `accrual_scenario`.
population_scenarios <- function(eligible, accrual = NULL) {
  if (is.null(accrual)) {
    accrual <- stats::setNames(ifelse(grepl("^norway", names(eligible)), "norway", names(eligible)), names(eligible))
  }
  data.frame(population = names(eligible), eligible_per_year = unname(eligible),
             accrual_scenario = unname(accrual[names(eligible)]), stringsAsFactors = FALSE)
}

#' Convert US dollars of a price year to 2023 euros
#'
#' Inflates with the US CPI-U annual averages (US_CPI_U_ANNUAL) to 2023 and
#' converts at the ECB 2023 average rate (USD_PER_EUR_2023): inflation in the
#' currency of the estimate, then one conversion in the target year.
#'
#' @param x Amount(s) in USD of `price_year`.
#' @param price_year Price year (a name of US_CPI_U_ANNUAL).
#' @return EUR 2023 (names of `x` kept).
usd_to_eur2023 <- function(x, price_year) {
  key <- as.character(price_year)
  if (!key %in% names(US_CPI_U_ANNUAL)) stop("No CPI-U value for ", price_year)
  x * US_CPI_U_ANNUAL[["2023"]] / US_CPI_U_ANNUAL[[key]] / USD_PER_EUR_2023
}

#' Trial-cost bands of the ENBS (EUR 2023)
#'
#' `norway` (main text, proposed in #213): METIMMOX's cost per patient,
#' NOK 50 million with nivolumab donated, converted at NOK_PER_EUR_MODEL. The
#' point is the trial team's figure, NOK 50 million / 68 patients taken as
#' 2023 NOK (`per_patient_point`); the low end divides by the 80 enrolled
#' patients at face value, the high end inflates the per-analysed-patient cost
#' from 2018, the first inclusion year, to 2023 with the Norwegian CPI
#' (NO_CPI_ANNUAL), so the band spans both the denominator and the cost-year
#' uncertainty. No separate fixed cost (a whole-trial cost per patient).
#' `oncology` (appendix): per-patient cost over the oncology interquartile
#' range of Moore et al. 2020 (TRIAL_COST_ONCOLOGY_PER_PATIENT_USD2019), no
#' separate fixed cost (the estimates are whole-trial costs divided by the
#' patients). `speich` (appendix): the full range of Speich et al. 2018
#' (TRIAL_COST_PER_PATIENT_USD2017, TRIAL_COST_FIXED_USD2017). Both in 2023
#' euros (usd_to_eur2023()). Nivolumab is costed separately (donated or
#' charged) in trial_cost(). PLACEHOLDERS, TO VERIFY.
#'
#' @return Named list of bands (`norway`, `oncology`, `speich`), each
#'   list(`per_patient` = c(low, high), `fixed` = c(low, high), `price_year`,
#'   `currency`, `source`, `source_price_year`, `conversion`); `norway` also
#'   has `per_patient_point`.
trial_cost_bands <- function() {
  conversion <- sprintf("US CPI-U to 2023 (%s), then %.4f USD per EUR (ECB 2023 average)",
                        paste(names(US_CPI_U_ANNUAL), US_CPI_U_ANNUAL, sep = " ", collapse = ", "),
                        USD_PER_EUR_2023)
  v_nok <- METIMMOX_TRIAL_COST_NOK / METIMMOX_TRIAL_COST_PATIENTS
  list(
    norway = list(
      per_patient = c(low = v_nok[["enrolled"]],
                      high = v_nok[["analysed"]] * NO_CPI_ANNUAL[["2023"]] / NO_CPI_ANNUAL[["2018"]]) /
        NOK_PER_EUR_MODEL,
      per_patient_point = v_nok[["analysed"]] / NOK_PER_EUR_MODEL,
      fixed = c(low = 0, high = 0), price_year = TRIAL_COST_PRICE_YEAR, currency = "EUR",
      source = paste0("METIMMOX (NCT03388190): total cost NOK 50 million, nivolumab donated by Bristol-Myers ",
                      "Squibb (trial team, personal communication, 15 October 2026)"),
      source_price_year = "nominal NOK 2018-2023, taken as 2023",
      conversion = sprintf(paste0("NOK %.0f per patient (/ %d analysed) taken as 2023 NOK; low end / %d enrolled; ",
                                  "high end inflated 2018 to 2023 with the Norwegian CPI (%s); %.2f NOK per EUR ",
                                  "(the model's rate)"),
                           v_nok[["analysed"]], as.integer(METIMMOX_TRIAL_COST_PATIENTS[["analysed"]]),
                           as.integer(METIMMOX_TRIAL_COST_PATIENTS[["enrolled"]]),
                           paste(names(NO_CPI_ANNUAL), NO_CPI_ANNUAL, sep = " ", collapse = ", "),
                           NOK_PER_EUR_MODEL)),
    oncology = list(
      per_patient = usd_to_eur2023(TRIAL_COST_ONCOLOGY_PER_PATIENT_USD2019[c("low", "high")],
                                   TRIAL_COST_ONCOLOGY_PRICE_YEAR),
      fixed = c(low = 0, high = 0), price_year = TRIAL_COST_PRICE_YEAR, currency = "EUR",
      source = paste0("Moore TJ et al. BMJ Open 2020;10:e038863. doi:10.1136/bmjopen-2020-038863 ",
                      "(oncology pivotal trials, IQR of the estimated cost per patient, USD; drug excluded)"),
      source_price_year = TRIAL_COST_ONCOLOGY_PRICE_YEAR, conversion = conversion),
    speich = list(
      per_patient = usd_to_eur2023(TRIAL_COST_PER_PATIENT_USD2017, 2017),
      fixed = usd_to_eur2023(TRIAL_COST_FIXED_USD2017, 2017), price_year = TRIAL_COST_PRICE_YEAR,
      currency = "EUR", source = TRIAL_COST_SOURCE, source_price_year = 2017, conversion = conversion))
}

#' Trial-cost parameters over the cost range
#'
#' @param u Cost index in [0, 1].
#' @param band A trial_cost_bands() element (default the main-text band).
#' @param per_patient,fixed Ranges c(low, high) in EUR 2023 (default from
#'   `band`).
#' @return Data frame with `u`, `per_patient` and `fixed` (EUR 2023).
trial_cost_parameters <- function(u, band = trial_cost_bands()[[TRIAL_COST_BAND_MAIN]],
                                  per_patient = band$per_patient, fixed = band$fixed) {
  if (any(u < 0 | u > 1)) stop("u must lie in [0, 1]")
  data.frame(u = u,
             per_patient = per_patient[["low"]] * (per_patient[["high"]] / per_patient[["low"]])^u,
             fixed = fixed[["low"]] + u * (fixed[["high"]] - fixed[["low"]]))
}

#' Cost of a trial
#'
#' @param n_randomised Randomised patients (vector).
#' @param n_registered Registered (screened) patients (vector, same length).
#' @param screening_cost Test cost per registered patient.
#' @param p_experimental Allocation to the experimental arm.
#' @param administrations Nivolumab administrations per experimental patient.
#' @param price Nivolumab price.
#' @param nivolumab "donated" or "charged".
#' @param per_patient,fixed Per-patient and fixed trial cost (EUR).
#' @return EUR (vector).
trial_cost <- function(n_randomised, n_registered, screening_cost, p_experimental, administrations, price,
                       nivolumab = c("donated", "charged"), per_patient, fixed) {
  nivolumab <- match.arg(nivolumab)
  v_drug <- if (nivolumab == "charged") n_randomised * p_experimental * administrations * price else 0
  fixed + per_patient * n_randomised + screening_cost * n_registered + v_drug
}

#' Nivolumab administrations per experimental patient in the enrolled cells
#'
#' @param cells The `cells` component of a PSA object (psa_cells_v1).
#' @param cell_mix Named enrolled-cell shares (design_enrolment_fraction()).
#' @return Mean over draws of the mix-weighted discounted administrations
#'   under the experimental arm (progression-free occupancy weighted; the
#'   4% discounting over the 38-week schedule understates the count slightly).
administrations_per_patient <- function(cells, cell_mix) {
  m_admin <- matrix(cells$values[, names(cell_mix), "experimental", "nivo_admin"], ncol = length(cell_mix))
  sum(colMeans(m_admin) * cell_mix)
}

#' Drop policies dominated in every draw (exact for EVPI and EVSI)
#'
#' A policy whose NMB is at most that of one other policy in every draw has
#' a conditional mean that is never above that policy's for any data, so it
#' never attains the maximum in any value-of-information quantity. Column 1
#' (the reference) is always kept.
#'
#' @param m_nmb NMB matrix (draws x policies, named columns).
#' @return The matrix without the dominated columns; attribute `dropped`.
drop_dominated_policies <- function(m_nmb) {
  v_keep <- rep(TRUE, ncol(m_nmb))
  for (j in seq_len(ncol(m_nmb))[-1]) {
    for (k in seq_len(ncol(m_nmb))) {
      if (k == j || !v_keep[k]) next
      if (all(m_nmb[, k] >= m_nmb[, j])) {
        v_keep[j] <- FALSE
        break
      }
    }
  }
  m_out <- m_nmb[, v_keep, drop = FALSE]
  attr(m_out, "dropped") <- colnames(m_nmb)[!v_keep]
  m_out
}

#' Row indices of a sampling-importance resample
#'
#' @param log_weights Log importance weights (one per row).
#' @param n Size of the resample.
#' @param seed RNG seed.
#' @return List with `index` (rows, with replacement), `ess` and `weights`.
evidence_resample_index <- function(log_weights, n = length(log_weights), seed = 123L) {
  if (any(!is.finite(log_weights))) stop("log_weights must be finite")
  v_w <- exp(log_weights - max(log_weights))
  v_w <- v_w / sum(v_w)
  set.seed(seed)
  list(index = sample.int(length(v_w), n, replace = TRUE, prob = v_w), ess = 1 / sum(v_w^2), weights = v_w)
}

#' Rows of the external-evidence likelihood used for the E2 state
#'
#' One result per trial and endpoint: the populations matching
#' `population_pattern` (ITT), with a follow-up; maintenance trials only when
#' `maintenance_landmark_weeks` is given (it becomes `landmark_weeks`); per
#' trial and endpoint the row with the longest follow-up. Which rows to use is
#' open (#202); this rule is the default until it is decided.
#'
#' @param df_lik external_evidence_likelihood() output.
#' @param population_pattern Regular expression on `population`.
#' @param maintenance_landmark_weeks Landmark of maintenance designs (weeks)
#'   or NA to drop them.
#' @return Selected rows (with `landmark_weeks`).
select_external_evidence <- function(df_lik, population_pattern = "^ITT",
                                     maintenance_landmark_weeks = 16) {
  d <- df_lik[grepl(population_pattern, df_lik$population) & is.finite(df_lik$follow_up_weeks) &
                is.finite(df_lik$p_experimental), , drop = FALSE]
  if (is.na(maintenance_landmark_weeks)) d <- d[!d$maintenance, , drop = FALSE]
  d$landmark_weeks <- ifelse(d$maintenance, maintenance_landmark_weeks, 0)
  d <- d[order(d$trial, d$endpoint, -d$follow_up_weeks), , drop = FALSE]
  d <- d[!duplicated(paste(d$trial, d$endpoint)), , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' NMB matrices of the policy space at a set of prices
#'
#' @param pairing voi_price_pairing of the policies (policy_price_pairing()).
#' @param prices Nivolumab prices.
#' @param wtp Willingness-to-pay threshold.
#' @param rows Row index (e.g. an evidence resample), or NULL for all draws.
#' @return Named list (by price) of NMB matrices with dominated policies
#'   dropped (drop_dominated_policies()).
policy_nmb_at_prices <- function(pairing, prices, wtp, rows = NULL) {
  stats::setNames(lapply(prices, function(price) {
    m_nmb <- pairing$effect * wtp - costs_at_price(pairing, price)
    if (!is.null(rows)) m_nmb <- m_nmb[rows, , drop = FALSE]
    drop_dominated_policies(as.matrix(m_nmb))
  }), format(prices, scientific = FALSE, trim = TRUE))
}

#' EVSI of every design over prices (evsi_by_design.csv rows)
#'
#' @param l_nmb policy_nmb_at_prices() result.
#' @param theta_draws Joint coefficient draws row-aligned with the NMB.
#' @param designs Named list of evsi_design() objects (build_evsi_designs()).
#' @param meta Design metadata (build_evsi_designs()$meta): which n are evaluated.
#' @param evidence_state Label.
#' @param accrual_scenario Label.
#' @param methods,seed,B evsi_by_design() settings.
#' @return Data frame: design, n, price, evidence_state, accrual_scenario,
#'   evsi, mcse, method, evpi, error, n_role; n not evaluated (report after the
#'   horizon) carry NA with method "reports_after_horizon".
run_evsi_grid <- function(l_nmb, theta_draws, designs, meta, evidence_state, accrual_scenario,
                          methods = EVSI_METHODS, seed = 123L, B = 1000) {
  rows <- list()
  for (id in unique(meta$design)) {
    df_m <- meta[meta$design == id, , drop = FALSE]
    v_eval <- df_m$n[df_m$evaluated]
    for (k in seq_along(l_nmb)) {
      price <- as.numeric(names(l_nmb)[k])
      df <- if (length(v_eval) && !is.null(designs[[id]]) && ncol(l_nmb[[k]]) < 2) {
        # One undominated policy: it is optimal in every draw, EVPI is zero.
        data.frame(design = id, n = v_eval, evsi = 0, mcse = 0, method = "zero_evpi_bound", evpi = 0,
                   error = "", stringsAsFactors = FALSE)
      } else if (length(v_eval) && !is.null(designs[[id]])) {
        evsi_by_design(l_nmb[[k]], theta_draws, designs[[id]], v_eval, methods = methods, seed = seed, B = B)
      } else NULL
      v_skip <- setdiff(df_m$n, v_eval)
      if (length(v_skip)) {
        df <- rbind(df, data.frame(design = id, n = v_skip, evsi = NA_real_, mcse = NA_real_,
                                   method = "reports_after_horizon", evpi = calculate_evpi_from_nmb(l_nmb[[k]]),
                                   error = "reports after the longest research horizon; not evaluated",
                                   stringsAsFactors = FALSE))
      }
      df$price <- price
      rows[[length(rows) + 1]] <- df
    }
  }
  df_out <- do.call(rbind, rows)
  df_out$evidence_state <- evidence_state
  df_out$accrual_scenario <- accrual_scenario
  df_out$n_role <- meta$n_role[match(paste(df_out$design, df_out$n), paste(meta$design, meta$n))]
  rownames(df_out) <- NULL
  df_out[c("design", "n", "n_role", "price", "evidence_state", "accrual_scenario", "evsi", "mcse", "method",
           "evpi", "error")]
}

#' Incremental EVSI of a new design after the METIMMOX-2 pooling
#'
#' @param df_evsi run_evsi_grid() rows (one method).
#' @param s4 Id of the pooling design ("S4").
#' @return Rows with design "<S4>-><id>" (= EVSI(S4 + id, n) - EVSI(S4)),
#'   MCSE sqrt(mcse_c^2 + mcse_s4^2) (the two estimates are not independent;
#'   conservative), role "increment".
sequential_increments <- function(df_evsi, s4 = "S4") {
  v_comb <- unique(df_evsi$design[startsWith(df_evsi$design, paste0(s4, "+"))])
  df_s4 <- df_evsi[df_evsi$design == s4, , drop = FALSE]
  rows <- lapply(v_comb, function(id) {
    d <- df_evsi[df_evsi$design == id, , drop = FALSE]
    v_s4 <- match(paste(d$price, d$evidence_state, d$accrual_scenario, d$method),
                  paste(df_s4$price, df_s4$evidence_state, df_s4$accrual_scenario, df_s4$method))
    d$evsi <- d$evsi - df_s4$evsi[v_s4]
    d$mcse <- sqrt(d$mcse^2 + df_s4$mcse[v_s4]^2)
    d$design <- sub("+", "->", id, fixed = TRUE)
    d
  })
  do.call(rbind, rows)
}

#' Design role from its id
#'
#' @param design Design ids.
#' @return "standalone", "combined" ("S4+Sx") or "increment" ("S4->Sx").
design_role <- function(design) {
  ifelse(grepl("->", design, fixed = TRUE), "increment", ifelse(grepl("+", design, fixed = TRUE), "combined",
                                                                "standalone"))
}

#' Research value of the designs and EVPI (voi_research_value.csv)
#'
#' @param df_evsi EVSI rows (run_evsi_grid(), plus sequential_increments()).
#' @param meta Design metadata over accrual scenarios (build_evsi_designs()$meta
#'   with an `accrual_scenario` column).
#' @param df_evpi Data frame with `price`, `evidence_state`, `evpi` (EVPI per
#'   patient over the policy space).
#' @param populations population_scenarios().
#' @param horizons Research horizons (years).
#' @param rate Discount rate of research benefits.
#' @param method EVSI method reported (one of EVSI_METHODS).
#' @return Data frame per design (and "EVPI"), n, price, evidence state,
#'   population and horizon: per-patient value, reporting year, first year of
#'   benefit, discount factor, eligible population and the population value
#'   (`max_justifiable_cost`, EUR) with its MCSE.
research_value_table <- function(df_evsi, meta, df_evpi, populations, horizons, rate, method = "gam_regression") {
  d <- df_evsi[df_evsi$method %in% c(method, "zero_evpi_bound", "zero_information", "not_estimated",
                                     "reports_after_horizon"), , drop = FALSE]
  meta$base_id <- sub("->", "+", meta$design, fixed = TRUE)
  rows <- list()
  for (p in seq_len(nrow(populations))) {
    pop <- populations[p, ]
    d_p <- d[d$accrual_scenario == pop$accrual_scenario, , drop = FALSE]
    # Metadata are matched by design, n and accrual, and by evidence state when
    # they carry one (the power-based n and enrolment fraction follow the
    # state's mean draw).
    v_state_d <- if ("evidence_state" %in% names(meta)) d_p$evidence_state else ""
    v_state_m <- if ("evidence_state" %in% names(meta)) meta$evidence_state else ""
    v_meta <- match(paste(sub("->", "+", d_p$design, fixed = TRUE), d_p$n, d_p$accrual_scenario, v_state_d),
                    paste(meta$base_id, meta$n, meta$accrual_scenario, v_state_m))
    if (anyNA(v_meta)) stop("Design metadata missing for: ",
                            paste(unique(paste(d_p$design[is.na(v_meta)], d_p$n[is.na(v_meta)])), collapse = ", "))
    v_report <- meta$reporting_year[v_meta]
    for (h in horizons) {
      v_start <- design_value_start_year(v_report)
      v_factor <- population_discount_factor(v_start, h, rate)
      rows[[length(rows) + 1]] <- data.frame(
        design = d_p$design, role = design_role(d_p$design), n = d_p$n, n_role = d_p$n_role, price = d_p$price,
        evidence_state = d_p$evidence_state, population = pop$population,
        accrual_scenario = pop$accrual_scenario, horizon = h, evsi_per_patient = d_p$evsi, mcse = d_p$mcse,
        method = d_p$method, reporting_year = v_report, value_start_year = v_start, discount_factor = v_factor,
        eligible_per_year = pop$eligible_per_year,
        max_justifiable_cost = ifelse(v_factor == 0, 0, d_p$evsi * pop$eligible_per_year * v_factor),
        max_justifiable_cost_mcse = d_p$mcse * pop$eligible_per_year * v_factor,
        n_registered = meta$n_registered[v_meta], stringsAsFactors = FALSE)
      v_factor1 <- population_discount_factor(1, h, rate)
      rows[[length(rows) + 1]] <- data.frame(
        design = "EVPI", role = "evpi", n = NA_real_, n_role = NA_character_, price = df_evpi$price,
        evidence_state = df_evpi$evidence_state, population = pop$population,
        accrual_scenario = pop$accrual_scenario, horizon = h, evsi_per_patient = df_evpi$evpi, mcse = 0,
        method = "evpi", reporting_year = 0, value_start_year = 1, discount_factor = v_factor1,
        eligible_per_year = pop$eligible_per_year,
        max_justifiable_cost = df_evpi$evpi * pop$eligible_per_year * v_factor1, max_justifiable_cost_mcse = 0,
        n_registered = NA_real_, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

#' ENBS curves and optimal n over the trial-cost range
#'
#' @param df_value research_value_table() rows (designs only are used).
#' @param design_costs Data frame per design and accrual scenario: `design`,
#'   `accrual_scenario`, `screening_cost`, `p_experimental`,
#'   `administrations`, `fixed_research_cost` (NA for a design costed from
#'   the range; S4_RESEARCH_COST for S4). Combined ("S4+Sx") and increment
#'   ("S4->Sx") rows carry the cost parameters of Sx: only Sx is paid.
#' @param u_grid Cost-index grid.
#' @param nivolumab Nivolumab options.
#' @param cost_points Named u values written as curves (low, mid, high).
#' @param band trial_cost_bands() element of the trial cost.
#' @param cost_band Its name (column `cost_band` of both outputs).
#' @return List with `curves` (enbs_by_design.csv: one row per design, n,
#'   price, state, population, horizon, nivolumab option and cost point) and
#'   `optimal` (per group and u: n_star, enbs_star, max justifiable cost and
#'   trial cost at n_star, share of the cost range with ENBS > 0).
enbs_analysis <- function(df_value, design_costs, u_grid = seq(0, 1, length.out = TRIAL_COST_GRID_POINTS),
                          nivolumab = TRIAL_NIVOLUMAB_OPTIONS, cost_points = c(low = 0, mid = 0.5, high = 1),
                          band = trial_cost_bands()[[TRIAL_COST_BAND_MAIN]], cost_band = TRIAL_COST_BAND_MAIN) {
  d <- df_value[df_value$role != "evpi", , drop = FALSE]
  v_late <- d$method == "reports_after_horizon"
  d$max_justifiable_cost[v_late] <- 0
  d <- d[is.finite(d$max_justifiable_cost), , drop = FALSE]
  df_cost <- trial_cost_parameters(sort(unique(c(u_grid, cost_points))), band)
  v_group <- paste(d$design, d$price, d$evidence_state, d$population, d$horizon, sep = "|")
  curves <- list()
  optimal <- list()
  for (g in unique(v_group)) {
    d_g <- d[v_group == g, , drop = FALSE]
    d_g <- d_g[order(d_g$n), , drop = FALSE]
    dc <- design_costs[design_costs$design == d_g$design[1] & design_costs$accrual_scenario == d_g$accrual_scenario[1], ]
    if (nrow(dc) != 1) stop("No unique cost row for design ", d_g$design[1], " (", d_g$accrual_scenario[1], ")")
    for (niv in nivolumab) {
      m_enbs <- vapply(seq_len(nrow(df_cost)), function(j) {
        v_cost <- if (is.finite(dc$fixed_research_cost)) rep(dc$fixed_research_cost, nrow(d_g)) else
          trial_cost(d_g$n, d_g$n_registered, dc$screening_cost, dc$p_experimental, dc$administrations,
                     d_g$price, niv, df_cost$per_patient[j], df_cost$fixed[j])
        d_g$max_justifiable_cost - v_cost
      }, numeric(nrow(d_g)))
      m_enbs <- matrix(m_enbs, nrow = nrow(d_g))
      v_best <- apply(m_enbs, 2, which.max)
      v_star_enbs <- m_enbs[cbind(v_best, seq_len(ncol(m_enbs)))]
      v_in_grid <- df_cost$u %in% u_grid
      share <- mean(v_star_enbs[v_in_grid] > 0)
      v_n_power <- d_g$n[d_g$n_role %in% "power"]
      optimal[[length(optimal) + 1]] <- data.frame(
        design = d_g$design[1], role = d_g$role[1], price = d_g$price[1],
        evidence_state = d_g$evidence_state[1], population = d_g$population[1], horizon = d_g$horizon[1],
        nivolumab = niv, u = df_cost$u, per_patient_cost = df_cost$per_patient, fixed_cost = df_cost$fixed,
        n_star = ifelse(v_star_enbs > 0, d_g$n[v_best], NA_real_), enbs_star = pmax(v_star_enbs, 0),
        best_n_if_run = d_g$n[v_best], enbs_best_n = v_star_enbs,
        max_justifiable_cost_star = d_g$max_justifiable_cost[v_best],
        trial_cost_star = d_g$max_justifiable_cost[v_best] - v_star_enbs,
        share_cost_range_positive = share,
        n_power = if (length(v_n_power)) v_n_power[1] else NA_real_,
        enbs_power = if (length(v_n_power)) m_enbs[match(v_n_power[1], d_g$n), ] else NA_real_,
        stringsAsFactors = FALSE)
      for (cp in names(cost_points)) {
        j <- match(cost_points[[cp]], df_cost$u)
        curves[[length(curves) + 1]] <- data.frame(
          design = d_g$design, role = d_g$role, n = d_g$n, n_role = d_g$n_role, price = d_g$price,
          evidence_state = d_g$evidence_state, population = d_g$population, horizon = d_g$horizon,
          nivolumab = niv, cost_point = cp, per_patient_cost = df_cost$per_patient[j],
          fixed_cost = df_cost$fixed[j], max_justifiable_cost = d_g$max_justifiable_cost,
          trial_cost = d_g$max_justifiable_cost - m_enbs[, j], enbs = m_enbs[, j],
          is_n_star = seq_len(nrow(d_g)) == v_best[j] & v_star_enbs[j] > 0,
          share_cost_range_positive = share, stringsAsFactors = FALSE)
      }
    }
  }
  df_curves <- do.call(rbind, curves)
  df_optimal <- do.call(rbind, optimal)
  df_curves$cost_band <- cost_band
  df_optimal$cost_band <- cost_band
  rownames(df_curves) <- NULL
  rownames(df_optimal) <- NULL
  list(curves = df_curves, optimal = df_optimal)
}

#' Research threshold price (voi_threshold_price.csv)
#'
#' @param df_optimal enbs_analysis()$optimal.
#' @param cost_ends Named u values of the cost-range ends.
#' @param roles Design roles eligible as "best design" (standalone designs).
#' @return Data frame per cost band (when `df_optimal` has `cost_band`),
#'   evidence state, population, horizon, nivolumab option and cost end: `threshold_price` (highest grid price at which the
#'   best eligible design has ENBS > 0 at its optimal n; NA when none),
#'   `best_design`, `n_star`, `max_justifiable_cost`, `trial_cost`, `enbs`,
#'   `at_grid_top` (the threshold is the highest price evaluated) and
#'   `prices_evaluated`.
research_threshold_price <- function(df_optimal, cost_ends = c(low = 0, high = 1), roles = "standalone") {
  d <- df_optimal[df_optimal$role %in% roles & df_optimal$u %in% cost_ends, , drop = FALSE]
  if (!"cost_band" %in% names(d)) d$cost_band <- rep(NA_character_, nrow(d))
  v_group <- paste(d$cost_band, d$evidence_state, d$population, d$horizon, d$nivolumab, d$u, sep = "|")
  rows <- lapply(unique(v_group), function(g) {
    d_g <- d[v_group == g, , drop = FALSE]
    v_prices <- sort(unique(d_g$price))
    d_pos <- d_g[is.finite(d_g$n_star) & d_g$enbs_star > 0, , drop = FALSE]
    base <- data.frame(cost_band = d_g$cost_band[1], evidence_state = d_g$evidence_state[1],
                       population = d_g$population[1],
                       horizon = d_g$horizon[1], nivolumab = d_g$nivolumab[1],
                       cost_end = names(cost_ends)[match(d_g$u[1], cost_ends)], u = d_g$u[1],
                       per_patient_cost = d_g$per_patient_cost[1], fixed_cost = d_g$fixed_cost[1],
                       threshold_price = NA_real_, best_design = NA_character_, n_star = NA_real_,
                       max_justifiable_cost = NA_real_, trial_cost = NA_real_, enbs = NA_real_,
                       at_grid_top = NA, prices_evaluated = length(v_prices), stringsAsFactors = FALSE)
    if (!nrow(d_pos)) return(base)
    p_star <- max(d_pos$price)
    d_star <- d_pos[d_pos$price == p_star, , drop = FALSE]
    best <- d_star[which.max(d_star$enbs_star), ]
    base$threshold_price <- p_star
    base$best_design <- best$design
    base$n_star <- best$n_star
    base$max_justifiable_cost <- best$max_justifiable_cost_star
    base$trial_cost <- best$trial_cost_star
    base$enbs <- best$enbs_star
    base$at_grid_top <- p_star == max(v_prices)
    base
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

#' Write the research-value CSV outputs
#'
#' @param results List with `evsi` (evsi_by_design.csv), `designs`
#'   (voi_designs.csv), `research_value`, `enbs` (curves) and `threshold`;
#'   optional `power_sizes` (voi_power_sizes.csv) and `decision_effects`
#'   (voi_decision_effects.csv), written when present.
#' @param directory Output directory (created if needed; default data/output/voi).
#' @return Paths of the CSVs written, invisibly.
write_voi_research_csvs <- function(results, directory = voi_output_dir()) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  v_files <- c(evsi = "evsi_by_design.csv", designs = "voi_designs.csv",
               research_value = "voi_research_value.csv", enbs = "enbs_by_design.csv",
               threshold = "voi_threshold_price.csv")
  v_missing <- setdiff(names(v_files), names(results))
  if (length(v_missing)) stop("results lacks: ", paste(v_missing, collapse = ", "))
  v_optional <- c(power_sizes = "voi_power_sizes.csv", decision_effects = "voi_decision_effects.csv")
  v_present <- names(v_optional)[vapply(names(v_optional), function(k) !is.null(results[[k]]), logical(1))]
  v_files <- c(v_files, v_optional[v_present])
  v_paths <- file.path(directory, v_files)
  for (k in seq_along(v_files)) utils::write.csv(results[[names(v_files)[k]]], v_paths[k], row.names = FALSE)
  invisible(stats::setNames(v_paths, names(v_files)))
}

#' Trial-cost inputs of every evaluated design (input of enbs_analysis())
#'
#' @param specs candidate_designs().
#' @param params Parameter list with the test prices.
#' @param administrations Named nivolumab administrations per experimental
#'   patient, by randomised design id (administrations_per_patient()).
#' @param accrual Accrual scenario labels.
#' @param design_ids Design ids evaluated (candidates, "S4+Sx", "S4->Sx").
#' @param p_experimental Named experimental shares by design id
#'   (design_allocation()), overriding the specification's allocation;
#'   required for a stratified design, whose share follows from the cell mix.
#' @return Data frame per design and accrual scenario: `design`,
#'   `accrual_scenario`, `cost_design` (the design whose trial is paid),
#'   `screening_cost`, `p_experimental`, `administrations`,
#'   `fixed_research_cost` (S4_RESEARCH_COST for the single-arm pooling, NA
#'   otherwise).
design_cost_table <- function(specs, params, administrations, accrual, design_ids, p_experimental = NULL) {
  v_single <- names(specs)[vapply(specs, function(s) s$type == "single_arm", logical(1))]
  rows <- lapply(design_ids, function(id) {
    v_base <- id
    for (s4 in v_single) v_base <- sub(paste0("^", s4, "([+]|->)"), "", v_base)
    spec <- specs[[v_base]]
    if (is.null(spec)) stop("Unknown design ", id)
    single <- spec$type == "single_arm"
    p_exp <- if (v_base %in% names(p_experimental)) p_experimental[[v_base]] else spec$p_experimental
    if (!single && !is.finite(p_exp)) stop("No experimental share for design ", v_base)
    data.frame(design = id, accrual_scenario = accrual, cost_design = v_base,
               screening_cost = design_screening_cost(spec, params), p_experimental = p_exp,
               administrations = if (single) 0 else unname(administrations[[v_base]]),
               fixed_research_cost = if (single) S4_RESEARCH_COST else NA_real_, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

# -------------------------------------------------------------------------------
# Evidence state E2 and the step 5-6 driver (analysis/14_VOI.R)
# -------------------------------------------------------------------------------

#' Shared heterogeneity SD of the external trials (log hazard ratio scale),
#' a floor common to every trial of an endpoint. ASSUMPTION, elicitation
#' pending (I13).
E2_TAU_SHARED_DEFAULT <- 0.1
#' Trial-specific heterogeneity SD (other backbones, agents, schedules).
#' ASSUMPTION, elicitation pending (I13).
E2_TAU_TRIAL_DEFAULT <- 0.2
#' Scenario ranges of the two heterogeneity SDs (elicitation pending)
E2_TAU_RANGE <- list(shared = c(low = 0, high = 0.2), trial = c(low = 0.1, high = 0.4))
#' Randomisation landmark of the maintenance trials in the evidence table
#' (MODUL: up to 8 two-weekly induction cycles). ASSUMPTION.
E2_MAINTENANCE_LANDMARK_WEEKS <- 16
#' Smallest effective sample size of the E2 importance resample
E2_MIN_ESS <- 500

#' NMB rows and coefficient draws of an evidence state on the E0 PSA
#'
#' The E0 PSA rows resampled by the exact external likelihood of every row's
#' coefficient draw (the weights of build_evidence_state()$check, the
#' conjugate update's own check, mapped by `model_idx`). Since 12 October 2026
#' this is the check of E2, not its evaluation: E2's net benefit comes from a
#' PSA on its own draws (R/voi_evidence_psa.R; ruling (9) of 11 October 2026),
#' compared with this resample by evidence_state_comparison().
#'
#' @param log_weights Log weights over the component's draws.
#' @param model_idx PSA rows' cached-model indices.
#' @param seed RNG seed.
#' @param min_ess Smallest acceptable ESS (over the PSA rows).
#' @return evidence_resample_index() over the PSA rows; stops below min_ess.
evidence_state_rows <- function(log_weights, model_idx, seed, min_ess = E2_MIN_ESS) {
  if (any(model_idx < 1 | model_idx > length(log_weights))) stop("model_idx outside the component's draws")
  rs <- evidence_resample_index(log_weights[model_idx], length(model_idx), seed)
  if (rs$ess < min_ess) {
    stop("E2 importance resample has ESS ", round(rs$ess), " < ", min_ess,
         "; the external evidence is too strong for resampling the E0 PSA (run a PSA on the E2 component)")
  }
  rs
}

#' Evidence states script 14 can run: METIMMOX alone (E0, the primary state)
#' and METIMMOX plus the published trials (E2)
VOI14_EVIDENCE_STATES <- c("E0", "E2")

#' Validated evidence states of a script 14 run
#'
#' E0 is the primary state and always runs; E2 is optional. The
#' METIMMOX-only switch (user, 12 October 2026) is `evidence_states = "E0"`:
#' the preliminary analysis without the published-trials state, whose numbers
#' wait for a PSA on its own draws (ruling (9) of 11 October 2026).
#'
#' @param states Requested states (character; NULL or "" for both).
#' @return The states in VOI14_EVIDENCE_STATES order.
voi14_evidence_states <- function(states = VOI14_EVIDENCE_STATES) {
  states <- unique(trimws(as.character(states)))
  states <- states[nzchar(states)]
  if (!length(states)) return(VOI14_EVIDENCE_STATES)
  v_bad <- setdiff(states, VOI14_EVIDENCE_STATES)
  if (length(v_bad)) {
    stop("Unknown evidence state(s) ", paste(v_bad, collapse = ", "), "; choose from ",
         paste(VOI14_EVIDENCE_STATES, collapse = ", "), call. = FALSE)
  }
  if (!"E0" %in% states) stop("evidence_states must include E0, the primary state (METIMMOX alone)", call. = FALSE)
  VOI14_EVIDENCE_STATES[VOI14_EVIDENCE_STATES %in% states]
}

#' Evidence states requested through the environment (make.R --voi-states)
#'
#' make.R runs the analysis in a child process, so its flag `--voi-states E0`
#' reaches script 14 as the environment variable METIMMOX_VOI14_STATES
#' (comma-separated); options(metimmox.voi14_settings = list(evidence_states =
#' ...)) still takes precedence in an interactive session.
#'
#' @param value Value of the variable (default the environment).
#' @return voi14_evidence_states() of the value (both states when unset).
voi14_requested_states <- function(value = Sys.getenv("METIMMOX_VOI14_STATES")) {
  voi14_evidence_states(strsplit(value, ",", fixed = TRUE)[[1]])
}

#' Sample-size grid of the reduced first run of script 14 (ruling (10) of
#' 11 October 2026)
VOI14_N_GRID_REDUCED <- c(40, 80, 160, 320, 640, 1280)
#' Accrual scenarios of the full grid of script 14
VOI14_ACCRUAL_FULL <- c("norway", "nordic", "europe")

#' Grid settings of analysis/14_VOI.R
#'
#' Ruling (10) of 11 October 2026: script 14 runs on a reduced grid first,
#' until its runtime is known: the marked prices only (list, biosimilar,
#' adoption threshold), both evidence states (E0, E2), Norway only and the
#' coarse n grid VOI14_N_GRID_REDUCED. The full grid (the price axis, the
#' default n grid VOI_N_GRID, Norway, the Nordic countries and Europe) is
#' chosen with options(metimmox.voi14_settings = list(grid = "full")). The
#' decision-relevant effects of the power-based sizes are evaluated at the
#' marked prices in both grids. Either grid runs the evidence states
#' `evidence_states` (default both; "E0" is the METIMMOX-only switch).
#'
#' @param grid "reduced" (default) or "full".
#' @param marked_prices Named marked nivolumab prices (list, biosimilar,
#'   adoption threshold).
#' @param price_axis Price grid of the full run (price_axis_grid()).
#' @param evidence_states Evidence states (voi14_evidence_states()).
#' @return List with `grid`, `prices`, `n_grid`, `accrual`,
#'   `decision_prices` and `evidence_states`.
voi14_grid_settings <- function(grid = c("reduced", "full"), marked_prices, price_axis,
                                evidence_states = VOI14_EVIDENCE_STATES) {
  grid <- match.arg(grid)
  v_marked <- sort(unique(unname(marked_prices)))
  if (!length(v_marked) || any(!is.finite(v_marked))) stop("marked_prices must be finite")
  v_states <- voi14_evidence_states(evidence_states)
  if (grid == "reduced") {
    list(grid = grid, prices = v_marked, n_grid = VOI14_N_GRID_REDUCED, accrual = "norway",
         decision_prices = v_marked, evidence_states = v_states)
  } else {
    list(grid = grid, prices = sort(unique(c(price_axis, v_marked))), n_grid = VOI_N_GRID,
         accrual = VOI14_ACCRUAL_FULL, decision_prices = v_marked, evidence_states = v_states)
  }
}

#' Eligible populations of the accrual scenarios being run
#'
#' @param eligible Named eligible patients per year (norway, its bounds and
#'   the other regions).
#' @param accrual Accrual scenarios run.
#' @return `eligible` restricted to the populations whose accrual scenario
#'   (population_scenarios()) is run (Norway only keeps Norway and its bounds).
voi14_eligible <- function(eligible, accrual) {
  df_pop <- population_scenarios(eligible)
  eligible[df_pop$population[df_pop$accrual_scenario %in% accrual]]
}

#' Runtime estimate of steps 5-6 of script 14 on the configured grid
#'
#' Times one design_information() call (S1 at the smallest n of the grid,
#' Norwegian accrual), one EVSI evaluation (evsi_by_design() for that design
#' and n at the price with the largest EVPI, all configured methods) and, when
#' `params` can
#' run the model, one deterministic cell-path evaluation of the
#' decision-relevant effect; extrapolates to the grid: information calls
#' (states x accrual scenarios x randomised and stratified designs x n),
#' power searches (about `power_calls` information calls per target, method
#' and hypothesis), EVSI evaluations (states x accrual x (those designs +
#' their S4 combinations) x n x prices, plus S4 once per state and price) and the
#' decision-relevant roots (about `root_calls` evaluations per biomarker and
#' price). A planning figure, printed before the run; not a bound.
#'
#' @param pairing voi_price_pairing of the policies.
#' @param theta_draws Joint coefficient draws of the PSA rows.
#' @param specs candidate_designs().
#' @param component Joint sampling component (E0).
#' @param reference_population,weights Reference patients and weights.
#' @param params Parameter list (l_params_base), or NULL to skip the
#'   decision-relevant timing.
#' @param settings Script 14 settings (`prices`, `n_grid`, `accrual`,
#'   `decision_prices`, `methods`, `wtp`, `seed`, `B`, `reference_rate`,
#'   `include_power`).
#' @param n_states Evidence states.
#' @param power_calls,root_calls Information calls per power search and
#'   model evaluations per root search.
#' @param verbose Print the estimate.
#' @return Data frame of the components (`component`, `count`,
#'   `seconds_each`, `seconds`) with attribute `total_minutes`.
voi14_runtime_estimate <- function(pairing, theta_draws, specs, component, reference_population, weights = NULL,
                                   params = NULL, settings, n_states = 2, power_calls = 14, root_calls = 14,
                                   verbose = TRUE) {
  draw <- split_joint_draw(component_mean_draw(component), component)
  spec <- specs[[which(vapply(specs, function(s) s$type == "rct", logical(1)))[1]]]
  rate <- settings$reference_rate * design_enrolment_fraction(spec, reference_population, weights, component,
                                                              draw)$fraction
  trial <- as_trial_design(spec, rate)
  n1 <- min(settings$n_grid)
  t_info <- system.time(info <- design_information(
    trial, n1, draw, component, reference_population, weights,
    step = design_node_step(trial$min_follow_up_weeks + design_accrual_weeks(trial, n1))))[["elapsed"]]
  # The price with the largest EVPI: at zero EVPI the engine short-circuits.
  l_nmb <- lapply(settings$prices, function(p) as.matrix(pairing$effect * settings$wtp - costs_at_price(pairing, p)))
  m_nmb <- l_nmb[[which.max(vapply(l_nmb, calculate_evpi_from_nmb, numeric(1)))]]
  design <- evsi_design(spec$id, info$targets, stats::setNames(list(info$covariance), n1))
  t_evsi <- system.time(evsi_by_design(m_nmb, theta_draws, design, n1, methods = settings$methods,
                                       seed = settings$seed, B = settings$B))[["elapsed"]]
  t_eval <- 0
  if (!is.null(params) && !is.null(params$l_nivo)) {
    t_eval <- system.time(guided_policy_cell_outcomes(component_mean_draw(component), component, params,
                                                      reference_population, weights))[["elapsed"]]
  }
  v_rct <- vapply(specs, function(s) s$type == "rct", logical(1))
  # Randomised and stratified designs are evaluated on the n grid, alone and with S4.
  n_grid_designs <- sum(vapply(specs, function(s) s$type %in% c("rct", "stratified"), logical(1)))
  n_targets <- sum(vapply(specs[v_rct], function(s) {
    if (is.na(s$power_target)) 0L else if (s$power_target == "interaction") length(s$power_biomarkers) else 1L
  }, integer(1)))
  n_acc <- length(settings$accrual)
  n_dec <- if (isFALSE(settings$include_power)) 0 else length(settings$decision_prices)
  n_power <- if (isFALSE(settings$include_power)) 0 else
    n_targets * length(POWER_METHODS) * (n_states + n_dec) * n_acc * power_calls
  df <- data.frame(
    component = c("design information", "power-based sizes", "EVSI evaluations", "decision-relevant effects"),
    count = c(n_states * n_acc * n_grid_designs * length(settings$n_grid), n_power,
              n_states * (n_acc * 2 * n_grid_designs * length(settings$n_grid) + 1) * length(settings$prices),
              if (t_eval > 0) 2 * n_dec * root_calls else 0),
    seconds_each = c(t_info, t_info, t_evsi, t_eval), stringsAsFactors = FALSE)
  df$seconds <- df$count * df$seconds_each
  attr(df, "total_minutes") <- sum(df$seconds) / 60
  if (verbose) {
    cat(sprintf("Runtime estimate (%s grid: %d prices, %d n, %d accrual scenario(s), %d states):\n",
                if (is.null(settings$grid)) "configured" else settings$grid, length(settings$prices),
                length(settings$n_grid), n_acc, n_states))
    for (i in seq_len(nrow(df))) {
      cat(sprintf("  %-26s %6.0f x %6.2f s = %6.1f min\n", df$component[i], df$count[i], df$seconds_each[i],
                  df$seconds[i] / 60))
    }
    cat(sprintf("  total about %.0f minutes (one EVSI evaluation %.1f s, one information call %.2f s)\n",
                attr(df, "total_minutes"), t_evsi, t_info))
  }
  df
}

#' Paper 3 steps 5 and 6: EVSI by design, research value, ENBS, threshold price
#'
#' @param pairing voi_price_pairing of the 16 policies (policy_price_pairing()).
#' @param theta_draws Joint coefficient draws of the PSA rows (E0; aligned by
#'   model_idx).
#' @param states Named list of evidence states; the first is the primary
#'   state (E0). Each is list(`component` (for the design information and its
#'   mean draw), and either `rows` (row index of `pairing` and `theta_draws`,
#'   e.g. an importance resample; NULL for all) or the state's own PSA:
#'   `pairing` (its voi_price_pairing) with `theta` (its joint coefficient
#'   draws by PSA row; R/voi_evidence_psa.R)).
#' @param specs candidate_designs().
#' @param reference_population,weights Reference patients and weights.
#' @param params Parameter list with the test prices.
#' @param cells The PSA's `cells` component (administrations per patient).
#' @param eligible Named eligible populations (norway first; the names in
#'   `settings$accrual` double as accrual scenarios).
#' @param settings List: `wtp`, `prices`, `n_grid`, `horizons`,
#'   `primary_horizon`, `rate`, `methods`, `primary_method`, `seed`, `B`,
#'   `tau_m2`, `tau_soc` (default TAU_SOC_DEFAULT), `tau_m3` (default
#'   TAU_M3_DEFAULT), `include_power` (default
#'   TRUE), `accrual` (accrual scenario names), `reference_rate`
#'   (Norwegian registered patients per year); optional `decision_prices`
#'   (prices of the decision-relevant effects; default `prices`; computed
#'   for the primary state when `include_power` is TRUE and `params` can run
#'   the model) and `cost_bands` (trial_cost_bands() or a subset; default
#'   all bands, each with its own ENBS rows).
#' @param verbose Print progress.
#' @return List with `evsi`, `designs`, `research_value`, `enbs` (curves at
#'   the primary horizon), `enbs_optimal`, `threshold`, `evpi`, `meta`,
#'   `accrual`, `design_costs`, `power_sizes` (design_power_sizes() rows by
#'   state and accrual scenario, or NULL), `decision_effects`
#'   (design_decision_effects(), or NULL) and `cost_bands`.
run_voi_research_value <- function(pairing, theta_draws, states, specs, reference_population, weights = NULL,
                                   params, cells, eligible, settings, verbose = TRUE) {
  say <- function(...) if (verbose) cat(..., "\n", sep = "")
  df_accrual <- accrual_scenarios(eligible[c("norway", setdiff(settings$accrual, "norway"))], settings$reference_rate)
  populations <- population_scenarios(eligible)
  tau_soc <- if (is.null(settings$tau_soc)) TAU_SOC_DEFAULT else settings$tau_soc
  tau_m3 <- if (is.null(settings$tau_m3)) TAU_M3_DEFAULT else settings$tau_m3
  # The pooling and the single-centre stratified design do not depend on the
  # accrual scenario: their EVSI is computed once per state and copied.
  v_fixed <- accrual_independent_designs(specs)
  evsi_rows <- list(); meta_rows <- list(); evpi_rows <- list(); power_e0 <- list(); power_rows <- list()
  # Decision-relevant effects (ruling (4)): primary state, deterministic at the mean coefficients.
  df_decision <- NULL
  if (!isFALSE(settings$include_power) && !isFALSE(settings$decision_relevant) && !is.null(params$l_nivo)) {
    v_dec_prices <- if (is.null(settings$decision_prices)) settings$prices else settings$decision_prices
    say("  decision-relevant effects at ", length(v_dec_prices), " price(s)")
    df_decision <- design_decision_effects(specs, v_dec_prices, settings$wtp, states[[1]]$component, params,
                                           reference_population, weights)
  }
  for (state in names(states)) {
    st <- states[[state]]
    if (!is.null(st$pairing)) {
      # The state's own PSA (ruling (9) of 11 October 2026): its NMB and draws.
      if (!is.null(st$rows)) stop("State ", state, ": give either rows or its own pairing, not both")
      if (is.null(st$theta) || nrow(st$theta) != st$pairing$n_sim) {
        stop("State ", state, ": its own pairing needs theta with one row per PSA draw")
      }
      if (!identical(st$pairing$strategies, pairing$strategies)) stop("State ", state, ": different policies")
      l_nmb <- policy_nmb_at_prices(st$pairing, settings$prices, settings$wtp)
      theta <- st$theta
    } else {
      l_nmb <- policy_nmb_at_prices(pairing, settings$prices, settings$wtp, st$rows)
      theta <- if (is.null(st$rows)) theta_draws else theta_draws[st$rows, , drop = FALSE]
    }
    evpi_rows[[state]] <- data.frame(price = settings$prices, evidence_state = state,
                                     evpi = vapply(l_nmb, calculate_evpi_from_nmb, numeric(1)),
                                     n_policies_undominated = vapply(l_nmb, ncol, integer(1)))
    draw <- component_mean_draw(st$component)
    df_fixed <- NULL
    for (k in seq_len(nrow(df_accrual))) {
      sc <- df_accrual[k, ]
      say("  ", state, " | accrual ", sc$accrual_scenario, ": design information")
      built <- build_evsi_designs(specs, sc$registered_per_year, settings$n_grid, draw, st$component,
                                  reference_population, weights, max_reporting_year = max(settings$horizons),
                                  tau_m2 = settings$tau_m2, include_power = !isFALSE(settings$include_power),
                                  decision_effects = if (state == names(states)[1]) df_decision else NULL,
                                  tau_soc = tau_soc, tau_m3 = tau_m3)
      if (state == names(states)[1]) power_e0[[sc$accrual_scenario]] <- built$power
      if (!is.null(built$power_sizes)) {
        power_rows[[length(power_rows) + 1]] <- cbind(evidence_state = state, accrual_scenario = sc$accrual_scenario,
                                                      built$power_sizes, stringsAsFactors = FALSE)
      }
      meta <- built$meta
      meta$accrual_scenario <- sc$accrual_scenario
      meta$evidence_state <- state
      meta_rows[[length(meta_rows) + 1]] <- meta
      designs <- built$designs
      meta_eval <- meta
      # Accrual-independent designs: computed once per state.
      if (!is.null(df_fixed)) {
        designs <- designs[setdiff(names(designs), v_fixed)]
        meta_eval <- meta[!meta$design %in% v_fixed, , drop = FALSE]
      }
      say("  ", state, " | accrual ", sc$accrual_scenario, ": EVSI of ", length(designs), " designs at ",
          length(l_nmb), " prices")
      df <- run_evsi_grid(l_nmb, theta, designs, meta_eval, state, sc$accrual_scenario,
                          methods = settings$methods, seed = settings$seed, B = settings$B)
      if (is.null(df_fixed)) {
        df_fixed <- df[df$design %in% v_fixed, , drop = FALSE]
      } else {
        df_copy <- df_fixed
        df_copy$accrual_scenario <- sc$accrual_scenario
        df <- rbind(df, df_copy)
      }
      evsi_rows[[length(evsi_rows) + 1]] <- df
    }
  }
  df_evsi <- do.call(rbind, evsi_rows)
  df_meta <- do.call(rbind, meta_rows)
  df_evpi <- do.call(rbind, evpi_rows)
  rownames(df_evpi) <- NULL
  df_evsi_all <- rbind(df_evsi, sequential_increments(df_evsi))
  rownames(df_evsi_all) <- NULL
  df_value <- research_value_table(df_evsi_all, df_meta, df_evpi, populations, settings$horizons, settings$rate,
                                   method = settings$primary_method)
  # Trial costs: screening, allocation and the administrations per
  # experimental patient of each paid design's experimental cell mix (E0 mean draw).
  component0 <- states[[1]]$component
  draw0 <- component_mean_draw(component0)
  v_paid <- names(specs)[vapply(specs, function(s) s$type %in% c("rct", "stratified"), logical(1))]
  l_alloc <- lapply(specs[v_paid], function(spec) {
    design_allocation(spec, reference_population, weights, component0, split_joint_draw(draw0, component0))
  })
  v_admin <- vapply(l_alloc, function(a) administrations_per_patient(cells, a$experimental_mix), numeric(1))
  v_p_exp <- vapply(l_alloc, function(a) a$p_experimental, numeric(1))
  df_costs <- do.call(rbind, lapply(df_accrual$accrual_scenario, function(a)
    design_cost_table(specs, params, v_admin, a, unique(df_value$design[df_value$role != "evpi"]),
                      p_experimental = v_p_exp)))
  say("  research value, ENBS and threshold price")
  l_bands <- if (is.null(settings$cost_bands)) trial_cost_bands() else settings$cost_bands
  l_enbs <- lapply(names(l_bands), function(b) enbs_analysis(df_value, df_costs, band = l_bands[[b]], cost_band = b))
  enbs <- list(curves = do.call(rbind, lapply(l_enbs, `[[`, "curves")),
               optimal = do.call(rbind, lapply(l_enbs, `[[`, "optimal")))
  df_threshold <- research_threshold_price(enbs$optimal)
  df_power_sizes <- if (length(power_rows)) do.call(rbind, power_rows) else NULL
  if (!is.null(df_power_sizes)) rownames(df_power_sizes) <- NULL
  df_designs <- design_table(specs, df_accrual, params, reference_population, weights, component0, draw0,
                             power = power_e0, tau_m2 = settings$tau_m2, tau_soc = tau_soc, tau_m3 = tau_m3)
  df_designs$administrations_per_experimental <- unname(v_admin[sub("^S6_", "", df_designs$design)])
  list(evsi = df_evsi_all, designs = df_designs, research_value = df_value,
       enbs = enbs$curves[enbs$curves$horizon == settings$primary_horizon, , drop = FALSE],
       enbs_optimal = enbs$optimal, threshold = df_threshold, evpi = df_evpi, meta = df_meta,
       accrual = df_accrual, design_costs = df_costs, power_sizes = df_power_sizes,
       decision_effects = df_decision, cost_bands = l_bands)
}

#' Fingerprint of the step 5-6 results
#'
#' @param inputs Named list of everything the results depend on: `upstream`
#'   (voi_upstream_fingerprints() of the sampling cache, the PSA pair and the
#'   scenario cache), evidence file md5, settings, design specifications,
#'   named design and cost parameters, population inputs, reference
#'   population hash.
#' @return List with `fingerprint` and `inputs`; the implementation is the
#'   `voi` and `evppi` categories of calculation_identity() (#205).
voi_research_fingerprint <- function(inputs) {
  if (is.null(inputs$upstream$sampling) || is.null(inputs$upstream$psa)) {
    stop("voi_research_fingerprint() needs inputs$upstream with the sampling cache and the PSA pair ",
         "(voi_upstream_fingerprints())", call. = FALSE)
  }
  inputs <- c(list(schema = "voi_research_value_v2"), inputs,
              list(implementation = list(voi = calculation_identity("voi"),
                                         evppi = calculation_identity("evppi"))))
  list(fingerprint = cache_fingerprint(inputs), inputs = inputs)
}

#' Named design, cost and evidence parameters of steps 5 and 6
#'
#' @return Named list of every constant of R/voi_designs.R and this file that
#'   the results depend on (for the fingerprint and the provenance record).
voi_research_parameters <- function() {
  list(decision_date = VOI_DECISION_DATE, metimmox = list(METIMMOX_ACCRUAL_PATIENTS, METIMMOX_ACCRUAL_SITES,
                                                          METIMMOX_ACCRUAL_START, METIMMOX_ACCRUAL_END),
       metimmox2 = list(METIMMOX2_PATIENTS, METIMMOX2_ACCRUAL_START, METIMMOX2_ACCRUAL_END,
                        METIMMOX2_ANALYSIS_DATE, METIMMOX2_REPORT_LAG_YEARS, METIMMOX2_FULL_WEIGHT_WEEKS),
       tau_m2 = list(TAU_M2_DEFAULT, TAU_M2_RANGE), follow_up = DESIGN_MIN_FOLLOW_UP_YEARS,
       metimmox3 = list(METIMMOX3_N, METIMMOX3_ACCRUAL_START, METIMMOX3_ACCRUAL_END, METIMMOX3_MIN_FOLLOW_UP_YEARS,
                        METIMMOX3_STARTUP_WEEKS, METIMMOX3_DATA_CUTOFF, METIMMOX3_N_POSITIVE),
       tau_soc = list(TAU_SOC_DEFAULT, TAU_SOC_RANGE), tau_m3 = list(TAU_M3_DEFAULT, TAU_M3_RANGE),
       power = list(POWER_RATE_RATIO, POWER_LEVEL, POWER_ALPHA, POWER_ENDPOINT), max_nodes = DESIGN_MAX_NODES,
       cost = list(TRIAL_COST_PER_PATIENT_USD2017, TRIAL_COST_FIXED_USD2017, EUR_PER_USD2017,
                   TRIAL_COST_GRID_POINTS, S4_RESEARCH_COST, METIMMOX_TRIAL_COST_NOK, METIMMOX_TRIAL_COST_PATIENTS,
                   NOK_PER_EUR_MODEL, NO_CPI_ANNUAL, TRIAL_COST_BAND_MAIN),
       e2 = list(E2_TAU_SHARED_DEFAULT, E2_TAU_TRIAL_DEFAULT, E2_MAINTENANCE_LANDMARK_WEEKS, E2_MIN_ESS))
}
