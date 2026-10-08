# Regression tests for A13: configuration-driven strategy/biomarker coupling.

source("R/model_configs.R")
source("R/calculate_outcomes.R")
source("R/model_fun.R")
source("R/psa_functions.R")
source("R/evppi_functions.R")
source("R/report_format.R")

# Strategy metadata is joined by ID, not by position.
v_reordered <- c("tmb_braf", "control", "crp")
df_metadata <- get_strategy_metadata(v_reordered)
stopifnot(
  identical(df_metadata$id, v_reordered),
  identical(df_metadata$short_name, c("TMB/BRAF", "Control", "CRP")),
  identical(unname(strategy_labels[v_reordered]), df_metadata$report_label),
  identical(unname(strategy_display_name(v_reordered)), df_metadata$short_name),
  grepl("tumor mutation burden", df_metadata$name[1], fixed = TRUE),
  grepl("Standard of care", df_metadata$name[2], fixed = TRUE)
)

# Biomarker cost and prevalence parameter names are derived centrally.
stopifnot(
  identical(
    biomarker_cost_key(),
    c(crp = "c_test_CRP", tmb_braf = "c_test_NGS")
  ),
  identical(
    biomarker_prevalence_key(),
    c(crp = "p_crp", tmb_braf = "p_tmb_braf")
  )
)
cost_params <- list(
  c_test_CRP = 31,
  c_test_NGS = 4100,
  c_test_biomarker = list(crp = 0, tmb_braf = 0)
)
cost_params <- sync_biomarker_test_costs(cost_params)
stopifnot(
  cost_params$c_test_biomarker$crp == 31,
  cost_params$c_test_biomarker$tmb_braf == 4100
)

# Build a minimal valid model input with a control ordering violation.
#' Minimal 14-point model input with zero costs and unit utilities whose control
#' PFS exceeds OS at one time point; returns the parameter list.
make_model_params <- function() {
  v_biomarkers <- get_biomarkers()
  params <- list(
    dr_costs = 0,
    dr_effects = 0,
    u_np = 1,
    u_p = 1,
    c_drug_nivo = 0,
    c_drug_FLOX = 0,
    c_test_CT = 0,
    c_test_blood = 0,
    c_test_CRP = 0,
    c_test_NGS = 0,
    c_test_biomarker = setNames(as.list(rep(0, length(v_biomarkers))), v_biomarkers),
    c_other_visit = 0,
    c_other_baseline = 0,
    c_other_follow = 0,
    c_other_last = 0,
    l_nivo = rep(0, 14),
    l_FLOX_exp = rep(0, 14),
    l_FLOX_control = rep(0, 14),
    l_CT = rep(0, 14),
    l_blood = rep(0, 14),
    l_visit = rep(0, 14),
    p_os = list(control_OS = c(1, rep(0.8, 13))),
    p_pfs = list(control_PFS = c(1, 0.9, rep(0.7, 12)))
  )
  for (biomarker in v_biomarkers) {
    params[[biomarker_prevalence_key(biomarker)]] <- 0.5
    params$p_os[[paste0(biomarker, "_pos_OS")]] <- c(1, rep(0.8, 13))
    params$p_os[[paste0(biomarker, "_neg_OS")]] <- c(1, rep(0.8, 13))
    params$p_pfs[[paste0(biomarker, "_pos_PFS")]] <- c(1, rep(0.7, 13))
    params$p_pfs[[paste0(biomarker, "_neg_PFS")]] <- c(1, rep(0.7, 13))
  }
  params
}

ordering_warning <- NULL
invisible(withCallingHandlers(
  model_fun(
    make_model_params(),
    time_horizon = 13,
    cl = 1,
    determpsa = "curves",
    sim_idx = NULL
  ),
  survival_ordering_warning = function(w) {
    ordering_warning <<- w
    invokeRestart("muffleWarning")
  }
))
stopifnot(
  inherits(ordering_warning, "survival_ordering_warning"),
  identical(ordering_warning$strategy, get_control_strategy()),
  identical(ordering_warning$subgroup, "control"),
  identical(ordering_warning$n_violations, 1L),
  isTRUE(all.equal(ordering_warning$max_excess, 0.1))
)

# PSA aggregation routes on condition fields even when the message is opaque.
#' Stub model that raises a survival_ordering_warning with an opaque message and
#' returns zero costs and effects for every strategy.
model_fun <- function(params, ...) {
  warning(warningCondition(
    "opaque message",
    class = "survival_ordering_warning",
    strategy = "tmb_braf",
    subgroup = "positive",
    n_violations = 1L,
    max_excess = 0.1,
    sim_idx = 1L
  ))
  data.frame(Strategy = get_strategies(), Cost = 0, Effect = 0)
}
df_psa_params <- data.frame(sim = 1L)
v_output <- capture.output({
  psa_result <- run_psa_analysis(
    psa_params = df_psa_params,
    l_params_base = make_model_params(),
    param_distributions = list(),
    strategies = get_strategies(),
    time_horizon = 13,
    cl = 1,
    n_sim = 1
  )
})
stopifnot(
  identical(psa_result$pfs_os_violations$tmb_braf_positive, 1L),
  any(grepl("TMB/BRAF+: 1 iteration", v_output, fixed = TRUE))
)

# EVPPI interaction extraction and grouping follow every configured biomarker.
sampling_models <- setNames(lapply(get_biomarkers(), function(biomarker) {
  v_coefficients <- setNames(0.5, paste0(biomarker, "1:Rxexperimental"))
  list(samples = list(list(
    os = list(coefficients = v_coefficients),
    pfs = list(coefficients = v_coefficients)
  )))
}), get_biomarkers())
invisible(capture.output(
  df_interaction_coefficients <- extract_interaction_coefficients(sampling_models, 1L)
))
v_expected_interactions <- as.vector(outer(
  get_biomarkers(), c("os", "pfs"),
  function(biomarker, outcome) paste0("b_", biomarker, "_rx_", outcome)
))
stopifnot(
  setequal(names(df_interaction_coefficients), v_expected_interactions),
  setequal(get_interaction_evppi_params(), v_expected_interactions),
  all(paste0("interaction_", get_biomarkers()) %in%
        names(add_interaction_param_groups(list())))
)

psa_source <- paste(readLines("R/psa_functions.R"), collapse = "\n")
evppi_source <- paste(readLines("R/evppi_functions.R"), collapse = "\n")
stopifnot(
  !grepl("grepl(\"PFS > OS constraint enforced\"", psa_source, fixed = TRUE),
  !grepl("intersect(biomarkers, c(\"crp\", \"tmb_braf\"))",
         evppi_source, fixed = TRUE)
)

# Analysis scripts are sourced into the same global environment as the R/
# helpers, so a script-level assignment to a name that R/ defines silently
# replaces it (issue #187: 07_traces.R overwrote the id-keyed strategy_labels of
# R/report_format.R with an unnamed vector of trace labels). Static check of the
# assignment targets outside function bodies; needs no data.
#' Names assigned outside function bodies in an R file (top level and inside
#' if/for/braces blocks).
script_level_targets <- function(file) {
  v_out <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    head <- e[[1]]
    if ((identical(head, as.name("<-")) || identical(head, as.name("="))) &&
        is.symbol(e[[2]])) v_out <<- c(v_out, as.character(e[[2]]))
    if (any(vapply(c("<-", "=", "if", "for", "{"), function(op)
      identical(head, as.name(op)), logical(1)))) {
      for (k in as.list(e)[-1]) walk(k)
    }
  }
  for (e in parse(file, keep.source = FALSE)) walk(e)
  unique(v_out)
}
v_r_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
v_r_names <- unique(unlist(lapply(v_r_files, script_level_targets)))
v_shadowing <- unlist(lapply(list.files("analysis", pattern = "[.]R$", full.names = TRUE),
  function(f) {
    hit <- intersect(script_level_targets(f), v_r_names)
    if (length(hit)) paste0(basename(f), ": ", hit) else character(0)
  }))
if (length(v_shadowing) > 0) {
  stop("Analysis scripts overwrite globals defined in R/: ",
       paste(v_shadowing, collapse = "; "))
}
stopifnot("strategy_labels" %in% v_r_names,
          !"strategy_labels" %in% script_level_targets("analysis/07_traces.R"))

# The guard detects the defect it was written for.
v_bad_script <- tempfile(fileext = ".R")
writeLines(c("x <- 1", "if (TRUE) {", "  strategy_labels <- c('a', 'b')", "}"),
           v_bad_script)
stopifnot("strategy_labels" %in% intersect(script_level_targets(v_bad_script), v_r_names))
unlink(v_bad_script)

cat("Model-configuration coupling tests passed.\n")
