# Focused regression tests for the sampling-script B4 refactor. No trial data or
# sampling caches are required.

source("scripts/R/functions/model_configs.R")
source("scripts/R/functions/parameter_distributions.R")

base_values <- list(
  c_drug_nivo = 100, c_drug_FLOX = 200,
  c_test_CT = 300, c_test_blood = 400,
  c_test_CRP = 500, c_test_NGS = 600,
  c_other_visit = 700, c_other_baseline = 800,
  c_other_follow = 900, c_other_last = 1000,
  u_np = 0.73, u_p = 0.59, u_decrement = 0.14,
  p_crp = 0.40, p_tmb_braf = 0.30
)
config <- configure_parameter_distributions(base_values)

# Unit prices are fixed in the PSA and varied in the DSA only (issue #154), so
# the PSA parameter set is the resource-use costs, the utilities and the
# prevalences. u_p is a derived column, not a draw.
expected_psa <- c(
  "c_other_visit", "c_other_baseline", "c_other_follow", "c_other_last",
  "u_np", "u_decrement", "u_p", "p_crp", "p_tmb_braf"
)
fixed_prices <- c("c_drug_nivo", "c_drug_FLOX", "c_test_CT", "c_test_blood",
                  "c_test_CRP", "c_test_NGS")
stopifnot(
  identical(names(config$distributions), expected_psa),
  !any(fixed_prices %in% names(config$distributions)),
  identical(config$distributions$u_p$dist, "derived"),
  # Groups describe what is sampled: derived columns are excluded because they
  # are collinear with the parameters they come from.
  setequal(unlist(config$groups, use.names = FALSE),
           setdiff(expected_psa, "u_p")),
  identical(config$groups$utilities, c("u_np", "u_decrement")),
  # Only one cost group is still sampled, so "all_costs" would duplicate it.
  is.null(config$groups$all_costs),
  isTRUE(all.equal(config$distributions$c_other_visit$shape, 25)),
  isTRUE(all.equal(config$distributions$c_other_visit$rate, 25 / 700))
)

# The DSA still varies the fixed prices, and never the PSA-only decrement.
dsa <- build_dsa_ranges(base_values, config$spec, mult = 0.2)
stopifnot(
  all(fixed_prices %in% dsa$pars),
  all(c("u_np", "u_p") %in% dsa$pars),
  !("u_decrement" %in% dsa$pars),
  isTRUE(all.equal(dsa$min[dsa$pars == "c_drug_nivo"], 80)),
  isTRUE(all.equal(dsa$max[dsa$pars == "c_drug_nivo"], 120)),
  # Utilities stay inside [0, 1].
  dsa$max[dsa$pars == "u_np"] <= 1
)

# The derived utility can never exceed the progression-free utility, including
# when the decrement draw is larger than u_np itself.
derived <- apply_derived_psa_parameters(data.frame(
  u_np = c(0.73, 0.50, 0.10),
  u_decrement = c(0.14, 0.60, 0.00)
))
stopifnot(
  all(derived$u_p <= derived$u_np),
  all(derived$u_p >= 0),
  isTRUE(all.equal(derived$u_p, c(0.59, 0.00, 0.10)))
)

# One-sided structural scenarios declare only the endpoint that differs.
scen <- dsa_structural_scenarios(base_dr = 0.04, base_horizon_years = 10)
stopifnot(
  setequal(names(scen), c("Discount_rate", "Time_horizon",
                          "Post_progression_cost", "Second_sequence")),
  identical(structural_scenario_sides(scen$Discount_rate), c("min", "max")),
  identical(structural_scenario_sides(scen$Post_progression_cost), "max"),
  identical(structural_scenario_sides(scen$Second_sequence), "min")
)

sampling_expressions <- parse("scripts/R/analysis/06_sampling.R")
extract_assignment <- function(name) {
  matches <- vapply(sampling_expressions, function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name(name))
  }, logical(1))
  stopifnot(sum(matches) == 1L)
  sampling_expressions[[which(matches)]]
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

test_env$extract_all_survival_probabilities <- function(prediction) {
  do.call(cbind, lapply(prediction$.pred, function(x) x$.pred_survival))
}
# Since issue #156 the prediction helper resolves the draw through the cache
# accessor; a legacy-shaped sample list is returned as-is here.
test_env$sampled_survival_models <- function(component, idx) component$samples[[idx]]
test_env$predict <- function(object, newdata, type, times) {
  rx_effect <- as.numeric(newdata$Rx == "experimental") * 10
  curves <- lapply(newdata$Age + rx_effect, function(value) {
    data.frame(.pred_survival = rep(value, length(times)))
  })
  data.frame(.pred = I(curves))
}
eval(
  extract_assignment("generate_psa_population_averaged_predictions"),
  envir = test_env
)

patient_data <- data.frame(
  Age = 1:4,
  Rx = factor(c("control", "experimental", "control", "experimental"),
              levels = c("control", "experimental")),
  crp = c(1, 1, 0, 0)
)
sampled_models <- list(samples = list(list(os = list(model = list()))))
control <- test_env$generate_psa_population_averaged_predictions(
  sampled_models, outcome = "os", data_original = patient_data,
  time_points = 0:2
)
subgroups <- test_env$generate_psa_population_averaged_predictions(
  sampled_models, biomarker_name = "crp", outcome = "os",
  data_original = patient_data, time_points = 0:2
)
# Control curve: every patient predicted with Rx = control (issue #151), so the
# mock's +10 experimental effect must NOT appear: mean(Age) = 2.5, not 7.5.
stopifnot(
  identical(control, rep(2.5, 3)),
  identical(subgroups$positive, rep(11.5, 3)),
  identical(subgroups$negative, rep(3.5, 3))
)

cat("Sampling B4 refactor regression tests passed.\n")
