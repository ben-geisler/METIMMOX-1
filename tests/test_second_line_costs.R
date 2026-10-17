# ===============================================================================
# TEST: second-line therapy pathway reader and expected cost (issue #209, part 1)
# ===============================================================================
# 1. Vial costing: whole vials at the lowest cost without sharing, milligrams
#    at the lowest price per mg with full sharing, linear in between.
# 2. A synthetic table that meets the contract gives regimen shares that sum
#    to 1 and an expected cost equal to a hand calculation; contacts and blood
#    tests are costed at the unit costs passed in (`c_other_visit`,
#    `c_test_blood`), not at values in the table.
# 3. The low/high scenarios bracket the base case and equal the brute-force
#    minimum and maximum over every vertex of the synthetic parameter box;
#    one-way ranges are bracketed by them.
# 4. Every contract violation stops with a message naming it.
# 5. Live (when data/external/second_line_pathway.csv exists; it is not in
#    the public mirror): the sourced table meets the contract, the costs are
#    finite and ordered, with the model's unit costs read from the
#    l_params_base literal of analysis/05_basecase_input_parameters.R
#    (parsed, not run).
#
# Run: "C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/test_second_line_costs.R
# ===============================================================================

source(here::here("R", "second_line_costs.R"))

#' Stop unless `expr` raises an error whose message matches `pattern`.
expect_error_matching <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error), fixed = TRUE)) {
    stop("Expected an error matching '", pattern, "', got: ",
         if (inherits(error, "error")) conditionMessage(error) else "no error")
  }
  invisible(TRUE)
}

#' Stop unless two numbers agree within a relative tolerance.
expect_near <- function(x, y, tol = 1e-9) {
  if (!isTRUE(abs(x - y) <= tol * max(1, abs(y)))) stop("Expected ", y, ", got ", x)
  invisible(TRUE)
}

# --- 1. Vial costing -----------------------------------------------------------
# 450 mg from 100 mg (NOK 300) and 400 mg (NOK 1000) vials: 400 + 100 = 1300
# beats 5 x 100 = 1500; exact = 450 x 2.5 = 1125.
expect_near(second_line_vial_cost(450, c(100, 400), c(300, 1000)), 1300)
expect_near(second_line_vial_cost(450, c(100, 400), c(300, 1000), sharing = 1), 1125)
expect_near(second_line_vial_cost(450, c(100, 400), c(300, 1000), sharing = 0.5), 1212.5)
# 350 mg: 4 x 100 = 1200 against 1 x 400 = 1000.
expect_near(second_line_vial_cost(350, c(100, 400), c(300, 1000)), 1000)
expect_near(second_line_vial_cost(400, c(100, 400), c(300, 1000)), 1000)
expect_near(second_line_vial_cost(0, 100, 300), 0)
cat("PASS: whole-vial and shared-vial costs\n")

# --- 2. Synthetic table --------------------------------------------------------
#' One row of the synthetic table.
row <- function(component, parameter, value, lower = value, upper = value, unit,
                regimen = NA, drug = NA, vial_mg = NA, administrations_per_cycle = NA,
                cost_parameter = NA) {
  data.frame(component = component, parameter = parameter, regimen = regimen, drug = drug,
             value = value, lower = lower, upper = upper, unit = unit, vial_mg = vial_mg,
             administrations_per_cycle = administrations_per_cycle,
             cost_parameter = cost_parameter, source = "Synthetic", doi_or_url = "https://example.org",
             retrieved = "2026-10-17", note = "", stringsAsFactors = FALSE)
}
df_valid <- rbind(
  row("pathway", "share_treated", 0.8, 0.7, 0.9, unit = "proportion"),
  row("pathway", "share_folfiri", 0.5, 0.2, 0.8, unit = "proportion"),
  row("pathway", "share_ras_braf_wt", 0.4, 0.3, 0.5, unit = "proportion"),
  row("pathway", "share_antiegfr_if_wt", 0.5, unit = "proportion"),
  row("pathway", "share_bevacizumab_if_no_antiegfr", 0.5, 0.25, 0.75, unit = "proportion"),
  row("schedule", "cycle_length_weeks", 2, unit = "weeks"),
  row("patient", "bsa", 2, 1.8, 2.1, unit = "m2"),
  row("patient", "weight", 70, unit = "kg"),
  row("dose", "dose_FOLFIRI_irinotecan", 100, unit = "mg/m2", regimen = "FOLFIRI",
      drug = "irinotecan", administrations_per_cycle = 1),
  row("dose", "dose_FLIRI_irinotecan", 100, unit = "mg/m2", regimen = "FLIRI",
      drug = "irinotecan", administrations_per_cycle = 1),
  row("dose", "dose_FLIRI_folinate", 60, unit = "mg", regimen = "FLIRI",
      drug = "folinate", administrations_per_cycle = 2),
  row("dose", "dose_bevacizumab", 5, unit = "mg/kg", regimen = "bevacizumab",
      drug = "bevacizumab", administrations_per_cycle = 1),
  row("dose", "dose_panitumumab", 6, unit = "mg/kg", regimen = "panitumumab",
      drug = "panitumumab", administrations_per_cycle = 1),
  row("dose_intensity", "relative_dose_intensity", 1, 0.9, 1, unit = "proportion"),
  row("price", "price_irinotecan_100mg", 100, unit = "NOK per vial", drug = "irinotecan", vial_mg = 100),
  row("price", "price_folinate_100mg", 50, unit = "NOK per vial", drug = "folinate", vial_mg = 100),
  row("price", "price_bevacizumab_400mg", 4000, 4000, 6000, unit = "NOK per vial",
      drug = "bevacizumab", vial_mg = 400),
  row("price", "price_panitumumab_400mg", 10000, unit = "NOK per vial",
      drug = "panitumumab", vial_mg = 400),
  row("currency", "nok_per_eur", 10, unit = "NOK per EUR"),
  row("wastage", "vial_sharing", 0, 0, 1, unit = "proportion"),
  row("resource", "contacts_per_cycle_FOLFIRI", 2, 1, 2, unit = "contacts per cycle",
      cost_parameter = "c_other_visit"),
  row("resource", "contacts_per_cycle_FLIRI", 2, unit = "contacts per cycle",
      cost_parameter = "c_other_visit"),
  row("resource", "blood_tests_per_cycle", 1, unit = "tests per cycle", cost_parameter = "c_test_blood"),
  row("duration", "cycles_none", 6, 4, 8, unit = "cycles"),
  row("duration", "cycles_bevacizumab", 8, unit = "cycles"),
  row("duration", "cycles_panitumumab", 10, unit = "cycles"))
path <- tempfile(fileext = ".csv")
#' Write a variant of the synthetic table and return its path.
write_variant <- function(df) {
  utils::write.csv(df, path, row.names = FALSE)
  path
}
params <- list(c_other_visit = 200, c_test_blood = 20)
df_read <- read_second_line_pathway(write_variant(df_valid))
base <- second_line_cost_at(df_read, params)
df_reg <- base$regimens
expect_near(sum(df_reg$share), 1)
stopifnot(nrow(df_reg) == 6L)
# Hand calculation (EUR): irinotecan 200 mg = 2 vials = NOK 200 -> EUR 20;
# FLIRI folinate 60 mg twice = 2 vials = NOK 100 -> EUR 10; bevacizumab 350 mg
# = 1 vial = EUR 400; panitumumab 420 mg = 2 vials = EUR 2000. Contacts
# 2 x 200, test 20. Shares: anti-EGFR 0.2, bevacizumab 0.4, none 0.4, each
# split 0.5/0.5 between FOLFIRI and FLIRI.
v_cycle <- c(none = 420, bevacizumab = 820, panitumumab = 2420)
v_course_folfiri <- c(none = 6, bevacizumab = 8, panitumumab = 10) * (20 + v_cycle)
v_course_flir <- c(none = 6, bevacizumab = 8, panitumumab = 10) * (30 + v_cycle)
v_bio_share <- c(none = 0.4, bevacizumab = 0.4, panitumumab = 0.2)
per_treated <- sum(0.5 * v_bio_share * (v_course_folfiri + v_course_flir))
expect_near(base$per_treated, per_treated)
expect_near(base$per_progressing, 0.8 * per_treated)
expect_near(base$weeks_treated, 2 * sum(v_bio_share * c(6, 8, 10)))
expect_near(base$per_week_treated, per_treated / base$weeks_treated)
# Unit costs come from params, not from the table.
base_double <- second_line_cost_at(df_read, list(c_other_visit = 400, c_test_blood = 20))
expect_near(base_double$per_treated - base$per_treated,
            sum(v_bio_share * c(6, 8, 10)) * 2 * 200)
expect_error_matching(second_line_cost_at(df_read, list(c_other_visit = 200)),
                      "params$c_test_blood must be a finite nonnegative scalar")
cat("PASS: regimen shares and expected cost match the hand calculation\n")

# --- 3. Low/high scenarios against brute force ----------------------------------
result <- second_line_cost(df_read, params)
v_ranged <- second_line_ranged_parameters(df_read)
m_vertices <- as.matrix(expand.grid(rep(list(c("lower", "upper")), length(v_ranged)),
                                    stringsAsFactors = FALSE))
v_vertex_cost <- apply(m_vertices, 1, function(v_bound) second_line_cost_at(
  df_read, params, second_line_values(df_read, stats::setNames(v_bound, v_ranged)))$per_progressing)
expect_near(result$summary["low", "per_progressing"], min(v_vertex_cost))
expect_near(result$summary["high", "per_progressing"], max(v_vertex_cost))
stopifnot(result$summary["low", "per_progressing"] < result$summary["base", "per_progressing"],
          result$summary["base", "per_progressing"] < result$summary["high", "per_progressing"],
          identical(result$summary["base", "per_progressing"], base$per_progressing))
df_one_way <- second_line_one_way(df_read, params)
stopifnot(setequal(df_one_way$parameter, v_ranged),
          all(pmin(df_one_way$cost_lower, df_one_way$cost_upper) >=
                result$summary["low", "per_progressing"] - 1e-9),
          all(pmax(df_one_way$cost_lower, df_one_way$cost_upper) <=
                result$summary["high", "per_progressing"] + 1e-9),
          !is.unsorted(rev(df_one_way$swing)))
cat("PASS: low/high scenarios equal the brute-force extremes of the parameter box\n")

# --- 4. Contract violations ----------------------------------------------------
expect_error_matching(read_second_line_pathway(file.path(tempdir(), "absent.csv")), "File not found")
expect_error_matching(read_second_line_pathway(write_variant(df_valid[, names(df_valid) != "retrieved"])),
                      "lacks columns: retrieved")
expect_error_matching(read_second_line_pathway(write_variant(rbind(df_valid, df_valid[1, ]))),
                      "duplicated parameters: share_treated")
expect_error_matching(read_second_line_pathway(write_variant(df_valid[df_valid$parameter != "nok_per_eur", ])),
                      "missing currency parameters: nok_per_eur")
df_bad <- df_valid; df_bad$source[2] <- ""
expect_error_matching(read_second_line_pathway(write_variant(df_bad)), "rows without a source: share_folfiri")
df_bad <- df_valid; df_bad$doi_or_url[3] <- NA
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "rows without a DOI or URL: share_ras_braf_wt")
df_bad <- df_valid; df_bad$retrieved[4] <- "October 2026"
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "rows without an ISO retrieval date (YYYY-MM-DD): share_antiegfr_if_wt")
df_bad <- df_valid; df_bad$component[6] <- "timing"
expect_error_matching(read_second_line_pathway(write_variant(df_bad)), "unknown components: timing")
df_bad <- df_valid; df_bad$unit[7] <- "cm"
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "units not allowed for their component: bsa")
df_bad <- df_valid; df_bad$upper[1] <- NA
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "value, lower or upper not finite: share_treated")
df_bad <- df_valid; df_bad$lower[2] <- 0.6
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "bounds not ordered (lower <= value <= upper): share_folfiri")
df_bad <- df_valid; df_bad$upper[1] <- 1.2
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "proportions outside [0, 1]: share_treated")
df_bad <- df_valid; df_bad$lower[24] <- 0
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "quantities that must be positive: cycles_none")
df_bad <- df_valid; df_bad$regimen[9] <- "FOLFOX"
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "dose rows without a backbone or biologic regimen: dose_FOLFIRI_irinotecan")
df_bad <- df_valid; df_bad$administrations_per_cycle[11] <- 1.5
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "dose rows without a positive whole number of administrations per cycle: dose_FLIRI_folinate")
expect_error_matching(read_second_line_pathway(write_variant(df_valid[df_valid$parameter != "price_folinate_100mg", ])),
                      "dosed drugs without a price row: folinate")
df_bad <- df_valid; df_bad$vial_mg[15] <- 0
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "price rows without a positive vial size: price_irinotecan_100mg")
df_bad <- df_valid; df_bad$cost_parameter[21] <- "c_other_follow"
expect_error_matching(read_second_line_pathway(write_variant(df_bad)),
                      "resource rows without a model unit cost (cost_parameter): contacts_per_cycle_FOLFIRI")
unlink(path)
cat("PASS: every contract violation stops with a message naming it\n")

# --- 5. Live table -------------------------------------------------------------
#' Numeric literals of the l_params_base list in analysis/05 (parsed, not run).
model_unit_costs <- function(file = here::here("analysis", "05_basecase_input_parameters.R"),
                             names = SECOND_LINE_COST_PARAMETERS) {
  for (expr in parse(file, keep.source = FALSE)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("l_params_base"))) {
      v_args <- as.list(expr[[3]])[-1]
      return(lapply(stats::setNames(names, names), function(name) {
        x <- v_args[[name]]
        if (!is.numeric(x)) stop("l_params_base$", name, " is not a numeric literal in ", file)
        x
      }))
    }
  }
  stop("No l_params_base assignment in ", file)
}

live_path <- here::here(SECOND_LINE_PATHWAY_PATH)
if (!file.exists(live_path)) {
  cat("SKIP: live second-line pathway table absent (", SECOND_LINE_PATHWAY_PATH, ")\n", sep = "")
} else {
  df_live <- read_second_line_pathway(live_path)
  params_live <- model_unit_costs()
  live <- second_line_cost(df_live, params_live)
  v_cost <- live$summary$per_progressing
  stopifnot(all(is.finite(as.matrix(live$summary))), all(v_cost > 0),
            v_cost[2] <= v_cost[1], v_cost[1] <= v_cost[3])
  expect_near(sum(live$base$regimens$share), 1)
  cat(sprintf(paste0("PASS: live table: EUR %.0f per progressing patient (low %.0f, high %.0f); ",
                     "EUR %.0f per treated patient over %.1f weeks\n"),
              v_cost[1], v_cost[2], v_cost[3], live$base$per_treated, live$base$weeks_treated))
}

cat("PASS: test_second_line_costs.R\n")
