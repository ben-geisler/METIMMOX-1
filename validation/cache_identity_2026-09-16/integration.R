# Optional integration check for issues #167/#169. Requires local trial data.
# Runs 50 draws and six scenarios in a temporary cache directory (~6 minutes).
# Writes no production caches or rendered artifacts.
library(here)
library(dplyr)
library(ggplot2)
source('R/cache_paths.R')
live_paths <- list.files(cache_dir(), pattern = '\\.(rds|RData)$', full.names = TRUE)
live_before <- tools::md5sum(live_paths)
scratch <- tempfile('issues-167-169-')
dir.create(scratch)
original_options <- options(metimmox.cache_dir = scratch)
cat('Temporary cache directory:', scratch, '\n')
pdf(file.path(scratch, 'plots.pdf'))
for (s in c('02_setup_and_global_variables','03_biomarker_strategies',
            '04_parametric_survival_analysis','05_basecase_input_parameters')) {
  invisible(capture.output(source(paste0('analysis/', s, '.R'))))
}
n_samples <- n_sim <- 50L
for (s in c('06_sampling','10_PSA','11_EVPPIs')) {
  cat('RUN ', s, '\n')
  invisible(capture.output(source(paste0('analysis/', s, '.R'))))
}
cat('CHECK pipeline PSA and EVPPI loaders\n')
po <- load_current_psa_cache()
ef <- evppi_cache_fingerprint(po, WTP, current_evppi_config(), analysis_seed, current_population_inputs())
invisible(load_evppi_cache(po, ef, WTP))
cat('CHECK repeated sources are stable\n')
options(metimmox.sampling_allow_regenerate = FALSE)
invisible(capture.output(source('analysis/06_sampling.R')))
invisible(capture.output(source('analysis/10_PSA.R')))
stopifnot(psa_cached)
cat('RUN scenario producer\n')
invisible(capture.output(source('analysis/12_scenario_EVPPIs.R')))
sc <- load_scenario_cache()
stopifnot(length(sc$all_scenario_results) == 6L)
cat('CHECK report contracts against temporary caches\n')
source('tests/test_report_contracts.R')
cat('CHECK missing sampling cache under report setup\n')
source('R/report_setup.R')
n_samples <- 51L
blocked <- tryCatch({setup_report(sources = '06', packages = character(),
  set_knitr = FALSE, set_theme = FALSE, quiet_sources = TRUE); NULL}, error=identity)
stopifnot(inherits(blocked, 'error'), grepl('regeneration is disabled', conditionMessage(blocked)),
  !file.exists(sampling_cache_path(51)))

stopifnot(identical(live_before, tools::md5sum(live_paths)))
cat('PASS: 50-draw pipeline; PSA, EVPPI, scenarios and repeated loads; live caches unchanged\n')
dev.off()

options(original_options)
