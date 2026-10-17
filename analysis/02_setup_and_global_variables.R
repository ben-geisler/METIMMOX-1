# load faster binary format
rds_path <- here::here("data", "tidy", "METIMMOX.rds")
if (file.exists(rds_path)) {
  data <- readRDS(rds_path)
} else {
  stop("File not found: ", rds_path)
}
rm(rds_path)
for (nm in c("LastEvalwk", "CT1date")) {
  if (!nm %in% names(data)) stop("Missing ", nm, ": re-run analysis/01_data_prep.R.")
}

data$PFSwk <- data$`Days until progression`/7
data$OSwk <- data$`Days until death/last follow up`/7

# First on-treatment CT (issue #155): weeks from inclusion (the OS/PFS clock
# origin) to the scan at which target lesion reduction (TLR) is read. CT1date is
# the positional export column `Date...122`, the date column that follows the
# baseline target-lesion block (`TL LD...121`); 01_data_prep.R validates it
# against the export contract (R/data_contract.R, issue #174) together with the
# first on-treatment target-lesion sum `TL LD...130`. CT1wk is NA exactly when
# TLR is NA (no first on-treatment scan). It defines the per-patient landmark in
# tlr_landmark.R and the first-scan diagnostics of the week-9 landmark; the
# economic model ignores it.
data$CT1wk <- as.numeric(as.Date(data$CT1date) -
                           as.Date(data$`Date of inclusion`)) / 7

# Rename "Progression exit" to "Progression" for simplicity
names(data)[names(data) == "Progression exit"] <- "Progression"
names(data)[names(data) == "Sex 0female"] <- "sex"

# PFS endpoint (issue #181): two scheduled 8-week assessment intervals.
# Inf recovers the pre-181 death rule; LastEvalwk is the alternative anchor.
pfs_death_window_weeks <- 16
pfs_last_assessment <- "TTPwk"
source(here::here("R/pfs_endpoint.R"))
data <- derive_pfs_endpoint(data, death_window_weeks = pfs_death_window_weeks,
                            last_assessment = pfs_last_assessment)

# global variables
## time parameters
cl <- 1/52         # 1 week cycle (not accounting for leap years)
time_horizon = 520 # 10 years in weeks

## analysis parameters
WTP <- 51000       # CE threshold (in Euros)
DSA_mult <- 0.2    # +/- 20% variations as standard for DSA
n_samples <- 5000  # number of resampled models
n_sim <- n_samples # number of PSA simulations
analysis_seed <- 123L  # default passed explicitly to each stochastic block

## global discount rate
dr = 0.04

## population parameters for scaling EVPI/EVPPI to population level (issues
## #200, #201). Eligible patients per year (adults with unresectable metastatic
## MSS/pMMR CRC fit for first-line doublet chemotherapy) are read from the
## sourced table data/external/voi_incidence.csv (one eligible_per_year row per
## region with extreme low/high scenarios; derivation in the private
## repository's docs/manuscripts/3_value_of_information/inputs/
## incidence_derivation.md). Norway is primary; its bounds and the Nordic and
## European values are scenario inputs. The former unsourced placeholder of
## 1,500 per year matched the number of all MSS/pMMR metastatic patients per
## year in Norway before the fitness step (about 1,517). Since #212 the fitness
## share is the co-authors' 0.60 of that number (personal communication,
## A. H. Ree, 15 October 2026; range 0.50-0.70 an assumption pending the Cancer
## Registry figure), about 910 per year.
source(here::here("R/voi_incidence.R"))
voi_incidence <- read_voi_incidence()
annual_incidence_norway <- voi_eligible_per_year(voi_incidence, "Norway")
annual_incidence_norway_low <- voi_eligible_per_year(voi_incidence, "Norway", "lower")
annual_incidence_norway_high <- voi_eligible_per_year(voi_incidence, "Norway", "upper")
annual_incidence_scenarios <- c(
  norway_low = annual_incidence_norway_low, norway_high = annual_incidence_norway_high,
  nordic = voi_eligible_per_year(voi_incidence, "Nordic"),
  europe = voi_eligible_per_year(voi_incidence, "Europe"))
rm(voi_incidence)
research_horizon_years <- 10     # Time horizon for research value (years)
discount_rate_research <- 0.04   # Discount rate of research benefits, aligned with
                                 # the model's 4% (user decision, 9 October 2026)

## Utility source switch
## Controls which health state utility values are used as base case
## Options:
##   0 = "ipd"      - METIMMOX EQ-5D-5L IPD, derived in the pipeline by
##                     trial_eq5d_utility_uncertainty() (R/trial_utilities.R) with
##                     the EQ5D_VALUE_SET value set and patient-clustered
##                     bootstrap SEs for the PSA (issue #188). Danish set:
##                     u_np = 0.912, u_p = 0.897.
##   1 = "correct"  - EQ-5D from CORRECT trial (u_np = 0.73, u_p = 0.59)
##                     (Grothey et al. Lancet 2013, cited by Gourzoulidis et al. J Comp Effect Res 2018)
##
## Both can be set without editing this file, e.g.
## options(metimmox.utility_source = 0, metimmox.eq5d_value_set = "uk") before
## sourcing the scripts. IPD caches are labelled ipd_<value set>.
## IMPORTANT: When changed, PSA and EVPPI caches must be regenerated.
## Sampling cache (survival models) is NOT affected.
UTILITY_SOURCE <- getOption("metimmox.utility_source", 1)  # Default: CORRECT (Gourzoulidis et al. 2018)
EQ5D_VALUE_SET <- getOption("metimmox.eq5d_value_set", "dk")  # EQ-5D-5L value set for IPD utilities
if (!UTILITY_SOURCE %in% c(0, 1)) stop("UTILITY_SOURCE must be 0 (IPD) or 1 (CORRECT).")
if (!EQ5D_VALUE_SET %in% c("dk", "uk")) stop("EQ5D_VALUE_SET must be 'dk' or 'uk'.")
utility_source_label <- if (UTILITY_SOURCE == 0) paste0("ipd_", EQ5D_VALUE_SET) else "correct"

# ===============================================================================
# SOURCE MODEL CONFIGURATIONS (Single Source of Truth)
# ===============================================================================
# This sources the central economic model configuration file which defines all
# formulas and provides helper functions: get_strategies(), get_biomarkers(), etc.
source(here::here("R/model_configs.R"))
source(here::here("R/prediction_population.R"))

# Source cache-path helpers (utility label + cached-object file locations).
source(here::here("R/cache_paths.R"))

# Source the shared PSA/DSA/EVPPI parameter specification.
source(here::here("R/parameter_distributions.R"))

# Set global strategy and biomarker vectors for backward compatibility.
# Economic analyses include control, CRP, and TMB/BRAF only.
strategies <- get_strategies()
biomarkers <- get_biomarkers()
