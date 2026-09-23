# Rebuild only the changed publication compilation from validated scenario
# results. PSA and EVPPI estimation functions are byte-for-byte equivalent as
# parsed R code. No expensive simulation/estimation input changed in #168.
source(here::here("R/report_setup.R"))
setup_report(sources = c("02", "03", "04", "05", "06"),
  funs = c("model_fun", "calculate_outcomes", "cea_helpers", "psa_functions",
           "evppi_functions", "scenario_analysis"), quiet_sources = TRUE)
expected <- current_scenario_fingerprint()
path <- scenario_evppi_path()
cache <- readRDS(path)
if (!identical(cache$fingerprint, expected$fingerprint)) {
  original <- system2("git", c("show", "b8e1b02:R/scenario_analysis.R"), stdout = TRUE)
  stopifnot(is.null(attr(original, "status")))
  before <- parse(text = original, keep.source = FALSE)
  after <- parse(here::here("R/scenario_analysis.R"), keep.source = FALSE)
  drop_compile <- function(exprs) {
    exprs[!vapply(exprs, function(e) is.call(e) &&
      identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("compile_evppi_results")), logical(1))]
  }
  stopifnot(identical(drop_compile(before), drop_compile(after)))
  definitions <- list()
  for (expr in before) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("function"))) {
      fn <- eval(expr[[3]], envir = baseenv())
      definitions[[as.character(expr[[2]])]] <- list(
        formals = paste(deparse(formals(fn)), collapse = "\n"),
        body = paste(deparse(body(fn)), collapse = "\n"))
    }
  }
  old_expected <- expected
  old_expected$inputs$implementation$code["R/scenario_analysis.R"] <- cache_fingerprint(
    list(source = paste(deparse(before), collapse = "\n"), definitions = definitions))
  old_expected$fingerprint <- cache_fingerprint(old_expected$inputs)
  cache <- load_scenario_cache(expected = old_expected, path = path)
  cat("Validated original scenario identity:", cache$fingerprint, "\n")
  cache$evppi_all_scenarios <- compile_evppi_results(cache$all_scenario_results)
  population <- function(x) calculate_population_evppi(x, annual_incidence_norway,
    research_horizon_years, discount_rate_research)
  cache$evppi_all_scenarios$evppi_population_millions <- population(cache$evppi_all_scenarios$evppi)
  cache$evppi_all_scenarios$evpi_population_millions <- population(cache$evppi_all_scenarios$evpi)
  cache$fingerprint <- expected$fingerprint
  cache$fingerprint_inputs <- expected$inputs
  cache$result_hash <- scenario_result_hash(cache)
  saveRDS(cache, path)
}
cache <- load_scenario_cache()
totals <- subset(cache$evppi_all_scenarios, parameter == "[EVPI]")
stopifnot(nrow(totals) == nrow(cache$scenarios))
cat("Validated rebuilt scenario identity:", cache$fingerprint, "\n")
print(totals[, c("scenario_id", "evpi")], row.names = FALSE)
