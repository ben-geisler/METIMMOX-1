# Focused regression tests for the sampling-script B4 refactor. No trial data or
# sampling caches are required.

source("R/model_configs.R")
source("R/parameter_distributions.R")
source("R/psa_functions.R")

base_values <- list(
  c_drug_nivo = 100, c_drug_FLOX = 200,
  c_test_CT = 300, c_test_blood = 400,
  c_test_CRP = 500, c_test_NGS = 600,
  c_other_visit = 700, c_other_baseline = 800,
  c_other_follow = 900, c_other_last = 1000,
  u_np = 0.73, u_p = 0.59, u_decrement = 0.14,
  p_crp = 0.40, p_tmb_braf = 0.30,
  p_joint_00 = 0.42, p_joint_01 = 0.18, p_joint_10 = 0.28, p_joint_11 = 0.12,
  joint_counts = c("00" = 42, "01" = 18, "10" = 28, "11" = 12),
  # Later treatment sequences (issue #208): counts by arm and the shares.
  later_sequences = list(counts = data.frame(
    arm = c("control", "experimental"), n_patients = c(32, 36), n_seq3 = c(5, 13),
    n_seq4 = c(0, 4))),
  p_seq3_control = 5 / 32, p_seq3_experimental = 13 / 36,
  p_seq4_given_seq3_control = 0, p_seq4_given_seq3_experimental = 4 / 13
)
config <- configure_parameter_distributions(base_values)

# Unit prices are fixed in the PSA and varied in the DSA only (issue #154), so
# the PSA parameter set is the resource-use costs, the utilities and the
# prevalences. u_p is a derived column, not a draw.
v_expected_psa <- c(
  "c_other_visit", "c_other_baseline", "c_other_follow", "c_other_last",
  "u_np", "u_decrement", "u_p", "p_joint_00", "p_joint_01", "p_joint_10", "p_joint_11", "p_crp", "p_tmb_braf",
  LATER_SEQUENCE_SHARE_PARAMETERS
)
v_fixed_prices <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                  "c_test_CRP", "c_test_NGS")
stopifnot(
  identical(names(config$distributions), v_expected_psa),
  !any(v_fixed_prices %in% names(config$distributions)),
  identical(config$distributions$u_p$dist, "derived"),
  # Groups describe what is sampled: derived columns are excluded because they
  # are collinear with the parameters they come from.
  setequal(unlist(config$groups, use.names = FALSE),
           setdiff(v_expected_psa, c("u_p", "p_joint_11", "p_crp", "p_tmb_braf"))),
  identical(config$groups$utilities, c("u_np", "u_decrement")),
  # Only one cost group is still sampled, so "all_costs" would duplicate it.
  is.null(config$groups$all_costs),
  isTRUE(all.equal(config$distributions$c_other_visit$shape, 25)),
  isTRUE(all.equal(config$distributions$c_other_visit$rate, 25 / 700)),
  # Later-sequence shares (issue #208): beta from the counts, their own group,
  # sampled after every other parameter.
  identical(config$groups$later_sequences, LATER_SEQUENCE_SHARE_PARAMETERS),
  identical(tail(names(config$distributions), 4), LATER_SEQUENCE_SHARE_PARAMETERS),
  identical(unlist(config$distributions$p_seq3_control[c("shape1", "shape2")]),
            c(shape1 = 5, shape2 = 27)),
  identical(unlist(config$distributions$p_seq4_given_seq3_experimental[c("shape1", "shape2")]),
            c(shape1 = 4, shape2 = 9)),
  identical(unlist(config$distributions$p_seq4_given_seq3_control[c("shape1", "shape2")]),
            c(shape1 = 0, shape2 = 5))
)
# Adding the later-sequence rows leaves the draws of every other parameter
# unchanged; a zero count is a point mass; beta(0, 0) is a point mass at 0.
v_other <- setdiff(names(config$distributions), LATER_SEQUENCE_SHARE_PARAMETERS)
df_all <- generate_psa_samples(config$distributions, 2000, seed = 11L)
df_without <- generate_psa_samples(config$distributions[v_other], 2000, seed = 11L)
stopifnot(identical(df_all[c("sim", v_other)], df_without[c("sim", v_other)]),
          all(df_all$p_seq4_given_seq3_control == 0),
          abs(mean(df_all$p_seq3_experimental) - 13 / 36) < 0.01,
          identical(later_sequence_beta_shapes("p_seq4_given_seq3_control",
                                               data.frame(arm = "control", n_patients = 30,
                                                          n_seq3 = 0, n_seq4 = 0)),
                    c(shape1 = 0, shape2 = 1)))

# The DSA still varies the fixed prices, and never the PSA-only decrement.
df_dsa <- build_dsa_ranges(base_values, config$spec, mult = 0.2)
stopifnot(
  all(v_fixed_prices %in% df_dsa$pars),
  all(c("u_np", "u_p") %in% df_dsa$pars),
  !("u_decrement" %in% df_dsa$pars),
  isTRUE(all.equal(df_dsa$min[df_dsa$pars == "c_drug_nivo"], 80)),
  isTRUE(all.equal(df_dsa$max[df_dsa$pars == "c_drug_nivo"], 120)),
  # Utilities stay inside [0, 1].
  df_dsa$max[df_dsa$pars == "u_np"] <= 1,
  # Later-sequence shares are varied and stay inside [0, 1] (issue #208).
  all(LATER_SEQUENCE_SHARE_PARAMETERS %in% df_dsa$pars),
  all(df_dsa$min[df_dsa$pars %in% LATER_SEQUENCE_SHARE_PARAMETERS] >= 0),
  all(df_dsa$max[df_dsa$pars %in% LATER_SEQUENCE_SHARE_PARAMETERS] <= 1)
)

# No one-way DSA run has u_p > u_np (issue #190). Varying u_np holds the
# decrement fixed (u_p moves with it); varying u_p is capped at u_np. Checked
# for the CORRECT pair, the IPD pair (u_p near u_np, +20% of u_p above 1) and a
# large variation where the linked u_p would go negative.
v_utility_pairs <- list(c(0.73, 0.59), c(0.912, 0.897), c(0.30, 0.05))
for (v_pair in v_utility_pairs) {
  for (mult in c(0.2, 0.5)) {
    l_case <- modifyList(base_values, list(u_np = v_pair[1], u_p = v_pair[2],
                                           u_decrement = v_pair[1] - v_pair[2]))
    df_case <- build_dsa_ranges(l_case, config$spec, mult = mult)
    for (i in seq_len(nrow(df_case))) {
      for (side in c("min", "max")) {
        l_run <- dsa_endpoint_params(l_case, df_case, i, side)
        stopifnot(
          l_run$u_p <= l_run$u_np, l_run$u_p >= 0, l_run$u_np <= 1,
          isTRUE(all.equal(l_run$u_decrement, l_run$u_np - l_run$u_p))
        )
      }
    }
    # The u_np runs keep the base decrement unless u_p is floored at zero.
    i_np <- which(df_case$pars == "u_np")
    for (side in c("min", "max")) {
      l_run <- dsa_endpoint_params(l_case, df_case, i_np, side)
      stopifnot(isTRUE(all.equal(
        l_run$u_p, max(l_run$u_np - (v_pair[1] - v_pair[2]), 0)
      )))
    }
  }
}
# The CORRECT lower bound: u_np 0.584 with u_p 0.444, not the base 0.59.
i_np <- which(df_dsa$pars == "u_np")
stopifnot(
  isTRUE(all.equal(df_dsa$min[i_np], 0.584)),
  identical(df_dsa$linked_par[i_np], "u_p"),
  isTRUE(all.equal(df_dsa$linked_min[i_np], 0.444)),
  isTRUE(all.equal(df_dsa$linked_max[i_np], 0.736)),
  # Other parameters move nothing else.
  all(is.na(df_dsa$linked_par[-i_np]))
)

# The derived utility can never exceed the progression-free utility, including
# when the decrement draw is larger than u_np itself.
df_derived <- apply_derived_psa_parameters(data.frame(
  u_np = c(0.73, 0.50, 0.10),
  u_decrement = c(0.14, 0.60, 0.00)
))
stopifnot(
  all(df_derived$u_p <= df_derived$u_np),
  all(df_derived$u_p >= 0),
  isTRUE(all.equal(df_derived$u_p, c(0.59, 0.00, 0.10)))
)

# One-sided structural scenarios declare only the endpoint that differs.
scen <- dsa_structural_scenarios(base_dr = 0.04, base_horizon_years = 10)
stopifnot(
  setequal(names(scen), c("Discount_rate", "Time_horizon", "Post_progression_cost",
                          "Second_line_therapy", "Second_sequence")),
  identical(structural_scenario_sides(scen$Discount_rate), c("min", "max")),
  identical(structural_scenario_sides(scen$Post_progression_cost), "max"),
  identical(structural_scenario_sides(scen$Second_line_therapy), "min"),
  identical(structural_scenario_sides(scen$Second_sequence), "min")
)
# Issue #209: no second-line cost switches the pathway off; the flat rate per
# quarter replaces it rather than adding to it.
sl_base <- list(p_sl_treated = 0.8, c_other_pp = 0)
sl_off <- apply_structural_scenario("Second_line_therapy", 0, sl_base, 520)$params
sl_on <- apply_structural_scenario("Second_line_therapy", 1, sl_base, 520)$params
sl_flat <- apply_structural_scenario("Post_progression_cost", 5000, sl_base, 520)$params
stopifnot(sl_off$p_sl_treated == 0, sl_on$p_sl_treated == 0.8,
          sl_flat$p_sl_treated == 0, sl_flat$c_other_pp == 5000,
          identical(apply_structural_scenario("Post_progression_cost", 5000,
                                              list(c_other_pp = 0), 520)$params,
                    list(c_other_pp = 5000)))

# Second-line rows (issue #209): sampled shares and cycles with standard
# errors from params$psa_se, fixed drug costs and contacts, their own EVPPI
# group, sourced DSA bounds; a list without any second-line field keeps the
# former specification, a partial one stops.
sl_spec <- parameter_distribution_spec()
sl_rows <- sl_spec[sl_spec$group == "second_line", ]
v_sl_sampled <- sl_rows$parameter[sl_rows$psa]
v_sl_values <- c(p_pfs_progression = 0.98, p_sl_treated = 0.8, p_sl_folfiri = 0.5,
                 p_sl_ras_braf_wt = 0.362, p_sl_antiegfr_if_wt = 0.5,
                 p_sl_bevacizumab_if_no_antiegfr = 0.25, n_sl_cycles_none = 9.8,
                 n_sl_cycles_bevacizumab = 12, n_sl_cycles_panitumumab = 10.8,
                 n_sl_contacts_FOLFIRI = 2, n_sl_contacts_FLIRI = 2, n_sl_blood_tests = 1,
                 c_sl_drug_FOLFIRI = 822, c_sl_drug_FLIRI = 623, c_sl_drug_bevacizumab = 887,
                 c_sl_drug_panitumumab = 1729)
sl_values <- c(base_values, as.list(v_sl_values))
sl_values$psa_se <- c(p_pfs_progression = sqrt(0.98 * 0.02 / 51),
                      setNames(rep(0.05, length(v_sl_sampled) - 1), setdiff(v_sl_sampled, "p_pfs_progression")))
sl_values$dsa_bounds <- data.frame(pars = c("p_sl_treated", "n_sl_contacts_FOLFIRI"),
                                   min = c(0.54, 1), max = c(0.9, 2))
sl_config <- configure_parameter_distributions(sl_values)
stopifnot(
  setequal(sl_rows$parameter, names(v_sl_values)),
  identical(v_sl_sampled, c("p_pfs_progression", "p_sl_treated", "p_sl_folfiri", "p_sl_ras_braf_wt",
                            "p_sl_antiegfr_if_wt", "p_sl_bevacizumab_if_no_antiegfr",
                            "n_sl_cycles_none", "n_sl_cycles_bevacizumab", "n_sl_cycles_panitumumab")),
  # Appended after every earlier row: the earlier draws keep their RNG stream.
  identical(head(names(sl_config$distributions), length(v_expected_psa)), v_expected_psa),
  identical(sl_config$groups$second_line, v_sl_sampled),
  is.null(sl_config$groups$all_costs),
  all(vapply(sl_config$distributions[grep("^p_", v_sl_sampled, value = TRUE)],
             function(d) d$dist == "beta", logical(1))),
  all(vapply(sl_config$distributions[grep("^n_", v_sl_sampled, value = TRUE)],
             function(d) d$dist == "gamma", logical(1))),
  # Beta(49, 1) from the count-based standard error.
  isTRUE(all.equal(c(sl_config$distributions$p_pfs_progression$shape1,
                     sl_config$distributions$p_pfs_progression$shape2), c(49, 1))),
  !any(c("c_sl_drug_FOLFIRI", "n_sl_contacts_FOLFIRI") %in% names(sl_config$distributions)),
  identical(names(config$distributions), v_expected_psa),
  is.null(config$groups$second_line)
)
sl_ranges <- build_dsa_ranges(sl_values, mult = 0.2)
stopifnot(
  isTRUE(all.equal(unlist(sl_ranges[sl_ranges$pars == "p_sl_treated", c("min", "max")]), c(min = 0.54, max = 0.9))),
  isTRUE(all.equal(unlist(sl_ranges[sl_ranges$pars == "n_sl_contacts_FOLFIRI", c("min", "max")]), c(min = 1, max = 2))),
  # +/- 20% elsewhere, capped at 1 for shares.
  isTRUE(all.equal(sl_ranges$max[sl_ranges$pars == "p_pfs_progression"], 1)),
  isTRUE(all.equal(sl_ranges$min[sl_ranges$pars == "c_sl_drug_FLIRI"], 623 * 0.8)),
  !any(grepl("^(p_sl_|n_sl_|c_sl_|p_pfs_)", build_dsa_ranges(base_values, mult = 0.2)$pars))
)
sl_partial <- base_values
sl_partial$p_sl_treated <- 0.8
stopifnot(inherits(tryCatch(configure_parameter_distributions(sl_partial), error = identity), "error"))
sl_no_se <- sl_values
sl_no_se$psa_se <- NULL
stopifnot(grepl("No dispersion", tryCatch(configure_parameter_distributions(sl_no_se),
                                          error = conditionMessage)))
cat("PASS: second-line specification rows, EVPPI group, sourced DSA bounds and scenarios (#209)\n")

sampling_expressions <- parse("analysis/06_sampling.R")
#' The single top-level `name <- ...` expression of analysis/06_sampling.R, so a
#' production helper can be evaluated without running the script.
extract_assignment <- function(name) {
  v_matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(v_matches) == 1L)
  sampling_expressions[[which(v_matches)]]
}

test_env <- new.env(parent = baseenv())
eval(extract_assignment("resolve_best_distributions"), envir = test_env)
local_models <- list(best_fit = list(
  os_distribution = "gamma", pfs_distribution = "gengamma"
))
stopifnot(
  identical(
    suppressWarnings(test_env$resolve_best_distributions(local_models)),
    list(os = "gamma", pfs = "gengamma")
  ),
  identical(
    suppressWarnings(test_env$resolve_best_distributions(NULL)),
    list(os = "weibull", pfs = "weibull")
  )
)

test_env$weighted_survival_average <- function(curves, weights) as.vector(curves %*% (weights / sum(weights)))
test_env$extract_all_survival_probabilities <- function(prediction) {
  do.call(cbind, lapply(prediction$.pred, function(x) x$.pred_survival))
}
# Since issue #156 the prediction helper resolves the draw through the cache
# accessor; a legacy-shaped sample list is returned as-is here.
test_env$sampled_survival_models <- function(component, idx) component$samples[[idx]]
test_env$predict <- function(object, newdata, type, times) {
  v_rx_effect <- as.numeric(newdata$Rx == "experimental") * 10
  curves <- lapply(newdata$Age + v_rx_effect, function(value) {
    data.frame(.pred_survival = rep(value, length(times)))
  })
  data.frame(.pred = I(curves))
}
eval(
  extract_assignment("generate_psa_population_averaged_predictions"),
  envir = test_env
)

df_patient_data <- data.frame(
  Age = 1:4,
  Rx = factor(c("control", "experimental", "control", "experimental"),
              levels = c("control", "experimental")),
  crp = c(1, 1, 0, 0)
)
sampled_models <- list(samples = list(list(os = list(model = list()))))
v_control <- test_env$generate_psa_population_averaged_predictions(
  sampled_models, outcome = "os", data_original = df_patient_data,
  time_points = 0:2
)
subgroups <- test_env$generate_psa_population_averaged_predictions(
  sampled_models, biomarker_name = "crp", outcome = "os",
  data_original = df_patient_data, time_points = 0:2
)
# Control curve: every patient predicted with Rx = control (issue #151), so the
# mock's +10 experimental effect must NOT appear: mean(Age) = 2.5, not 7.5.
stopifnot(
  identical(v_control, rep(2.5, 3)),
  identical(subgroups$positive, rep(11.5, 3)),
  identical(subgroups$negative, rep(3.5, 3))
)

cat("Sampling B4 refactor regression tests passed.\n")
