# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

METIMMOX-1 is a cost-effectiveness analysis comparing biomarker-guided immunotherapy strategies for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer. The analysis uses a partitioned survival model implemented in R to evaluate three biomarker strategies (CRP, TLR, TMB/BRAF) against standard of care.

**Target Audience**: Health economists and researchers developing decision-analytic models in R.

## Essential Commands

### Running the Analysis

The numbered analysis scripts must be executed sequentially:

```r
# Core setup (run these first)
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")

# Survival analysis
source("scripts/R/analysis/06_parametric_survival analysis.R")
source("scripts/R/analysis/07_basecase_input_parameters.R")

# Survival resampling (generates cache - takes time on first run)
source("scripts/R/analysis/08_sampling.R")

# Model execution
source("scripts/R/analysis/09_traces.R")
source("scripts/R/analysis/10_basecase_analysis.R")

# Sensitivity analyses
source("scripts/R/analysis/11_DSA.R")  # Deterministic sensitivity analysis
source("scripts/R/analysis/12_PSA.R")  # Probabilistic sensitivity analysis
source("scripts/R/analysis/13_EVPPIs.R")  # Expected value of perfect partial information
```

### Package Installation

```r
# Install required packages using pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, readxl, dplyr, tableone, ggplot2, flexsurv,
               survival, survminer, gems, mstate, tidyverse, xtable,
               darthtools, dampack, mvtnorm, Matrix, here)

# Additional packages for Quarto reports
pacman::p_load(knitr, kableExtra, flextable, officer, scales, gridExtra, reshape2)
```

### Rendering Quarto Reports

The project includes comprehensive Quarto reports in `scripts/QMD/report/` that generate publication-ready PDF outputs:

```bash
# Render individual reports (from project root)
quarto render scripts/QMD/report/CEA.qmd
quarto render scripts/QMD/report/clinical_effectiveness.qmd
quarto render scripts/QMD/report/para_models.qmd
quarto render scripts/QMD/report/OWSA.qmd
quarto render scripts/QMD/report/EVPPIs.qmd
quarto render scripts/QMD/report/input_parameters.qmd
quarto render scripts/QMD/report/scenario_effect.qmd

# Render all reports at once
quarto render scripts/QMD/report/
```

**Prerequisites for rendering**:
- All analysis scripts (02-13) must be run first to generate required data objects
- Sampling cache must exist (from running [08_sampling.R](scripts/R/analysis/08_sampling.R))
- Results objects (e.g., `cea_results`, `owsa_results`, `psa_results`, `evppi_results`) must be in the R environment or saved as `.rds` files

### Complete Analysis & Reporting Workflow

To generate all analysis results and reports from scratch:

```r
# 1. Run all analysis scripts in order
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")
source("scripts/R/analysis/06_parametric_survival analysis.R")
source("scripts/R/analysis/07_basecase_input_parameters.R")
source("scripts/R/analysis/08_sampling.R")  # Takes time on first run
source("scripts/R/analysis/09_traces.R")
source("scripts/R/analysis/10_basecase_analysis.R")
source("scripts/R/analysis/11_DSA.R")
source("scripts/R/analysis/12_PSA.R")
source("scripts/R/analysis/13_EVPPIs.R")

# 2. Render reports (from terminal/command line)
# quarto render scripts/QMD/report/
```

Or render reports individually in the desired order:
```bash
quarto render scripts/QMD/report/para_models.qmd
quarto render scripts/QMD/report/input_parameters.qmd
quarto render scripts/QMD/report/clinical_effectiveness.qmd
quarto render scripts/QMD/report/CEA.qmd
quarto render scripts/QMD/report/OWSA.qmd
quarto render scripts/QMD/report/EVPPIs.qmd
quarto render scripts/QMD/report/scenario_effect.qmd
```

## Architecture & Key Concepts

### Model Structure

The analysis implements a **partitioned survival model (PSM)** with three health states:
- **Progression-free (PF)**: Patients alive without disease progression
- **Progressed (P)**: Patients alive with disease progression
- **Dead (D)**: Absorbing state

State occupancy is calculated from survival curves:
- PF state = PFS curve
- P state = max(OS - PFS, 0)
- D state = 1 - OS

### Critical Global Variables

Defined in [02_setup_and_global_variables.R](scripts/R/analysis/02_setup_and_global_variables.R):

```r
cl <- 1/52              # Cycle length: 1 week
time_horizon <- 520     # 10 years in weeks
WTP <- 51000           # Willingness-to-pay threshold (Euros)
DSA_mult <- 0.2        # ±20% variation for DSA
n_samples <- 5000      # Resampling/PSA sample size
dr <- 0.04             # Discount rate (4%)
USE_BOTH_MODELS <- 0   # 0 = full model only, 1 = both models
```

### Biomarker Strategies

Three biomarkers are evaluated (defined in [03_biomarker_strategies.R](scripts/R/analysis/03_biomarker_strategies.R)):

1. **CRP** (C-reactive protein): Binary variable, cut-off <5
2. **TLR** (Tumor lesion reduction): Binary variable, cut-off ≥10%
3. **TMB/BRAF**: Combined biomarker (TMB ≥9 mut/MB OR BRAF mutation)

Each strategy has:
- **Biomarker-positive subgroup**: Receives experimental treatment (alternating FLOX + nivolumab)
- **Biomarker-negative subgroup**: Receives standard treatment (FLOX only)
- **Population-level outcomes**: Weighted by biomarker prevalence

### Survival Resampling & Correlation

The model uses **correlated survival resampling** ([08_sampling.R](scripts/R/analysis/08_sampling.R:103-200)) to maintain the correlation between PFS and OS:

- Both PFS and OS models are fitted to the **same resampled patient cohort**
- Results are cached in `data/tidy/`
- Cache file naming: `sampling_models_n{n_samples}_{model_type}.rds`

**IMPORTANT**: The first run of `08_sampling.R` will take significant time (generates 5000 resampled models). Subsequent runs load from cache.

### Main Model Function

The core economic model is in [model_fun.R](scripts/R/functions/model_fun.R:9-381):

```r
model_fun(params, time_horizon = 520, cl = 1/52,
          determpsa = "det", return_traces = FALSE, sim_idx = NULL)
```

**Parameters**:
- `determpsa`: "det" for deterministic, "psa" for probabilistic
- `sim_idx`: Resampled model index (required for PSA mode)
- `return_traces`: If TRUE, returns state occupancy over time

**PSA Mode Behavior** (lines 17-177):
- When `determpsa = "psa"` and `sim_idx` is provided, survival curves are generated from resampled model #`sim_idx`
- PFS and OS always come from the **same** resampled model to maintain correlation
- If prediction from resampled model fails, falls back to base case curves

### Parameter Structure

The model expects a comprehensive parameter list `l_params_base` containing:

```r
list(
  # Time & discounting
  cl, time_horizon, dr_costs, dr_effects,

  # Utilities (beta distributions in PSA)
  u_np, u_p,

  # Drug costs (gamma distributions in PSA)
  c_drug_nivo, c_drug_FLOX,

  # Test costs
  c_test_CT, c_test_blood, c_test_NGS,

  # Other costs
  c_other_visit, c_other_baseline, c_other_follow, c_other_last,

  # Treatment schedules (vectors of length time_horizon+1)
  l_nivo, l_FLOX_exp, l_FLOX_control, l_CT, l_blood, l_visit,

  # Survival curves (lists of vectors)
  p_os = list(control_OS, crp_pos_OS, crp_neg_OS, ..., crp_weighted_OS, ...),
  p_pfs = list(control_PFS, crp_pos_PFS, crp_neg_PFS, ..., crp_weighted_PFS, ...),

  # Biomarker prevalence
  p_crp, p_tlr, p_tmb_braf
)
```

### Data Requirements

**Confidential Trial Data**: The clinical dataset is NOT included in the repository. The model expects:
- File location: `data/tidy/METIMMOX.rds`
- Required variables: `ID`, `PFSwk`, `Progression`, `OSwk`, `Death`, `Rx`, `Age`, `sex`, biomarker variables

**Biomarker derivation** (in [03_biomarker_strategies.R](scripts/R/analysis/03_biomarker_strategies.R:6-9)):
```r
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))
```

### Survival Model Formulas

**Control group** (age- and sex-adjusted):
```r
OS:  Surv(OSwk, Death) ~ Age + sex
PFS: Surv(PFSwk, Progression) ~ Age + sex
```

**Biomarker groups** (includes treatment interaction):
```r
OS:  Surv(OSwk, Death) ~ Age + sex + Rx + [biomarker]:Rx
PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + [biomarker]:Rx
```

Where `[biomarker]` is one of: `crp`, `tlr`, `tmb_braf`

### Key Functions

- **[calculate_outcomes.R](scripts/R/functions/calculate_outcomes.R)**: Calculates QALYs and costs from state occupancy traces
- **[bootstrap_survival_model.R](scripts/R/functions/bootstrap_survival_model.R)**: Alternative resampling approach (not used in main analysis)
- **[psa_functions.R](scripts/R/functions/psa_functions.R)**: PSA-related utilities
- **[evppi_functions.R](scripts/R/functions/evppi_functions.R)**: EVPPI calculation functions
- **[prediction_functions.R](scripts/R/functions/prediction_functions.R)**: Generate survival predictions from fitted models

### Treatment Schedules

Treatment administration is defined by binary vectors indicating weeks when treatments are given:

**Nivolumab** (experimental): Weeks 5, 7, 13, 15, 29, 31, 37, 39
**FLOX experimental**: Weeks 1, 3, 9, 11, 25, 27, 33, 35
**FLOX control**: All 16 time points from both nivolumab and FLOX experimental schedules

**Monitoring**:
- CT scans: Baseline + every 12 weeks
- Blood tests: Baseline + every 4 weeks
- Visits: Baseline + all treatment administration weeks

## Important Modeling Considerations

### Ensuring Non-Negative States

The model enforces `p_p = pmax(os - pfs, 0)` to prevent negative progressed state occupancy when PFS and OS curves cross (can happen in resampled models).

### Initial State Constraints

All cohorts start 100% progression-free:
```r
p_pf[1] <- 1.0
p_p[1] <- 0.0
p_d[1] <- 0.0
```

### Biomarker Prediction in PSA

When generating predictions from resampled models for biomarker strategies ([model_fun.R](scripts/R/functions/model_fun.R:73-100)), **all three biomarker variables must be included** in the newdata frame:
- Set target biomarker to 1 (positive) or 0 (negative)
- Set other biomarkers to reference values (modal values from data)
- Include age, sex, and treatment (Rx)

This is required because resampled models were fitted with all biomarkers as covariates.

### Discount Factor Application

Discount factors are applied using vector multiplication over the time horizon:
```r
v_dw_c <- 1 / (1 + dr_costs)^(seq(0, time_horizon) / 52)
v_dw_e <- 1 / (1 + dr_effects)^(seq(0, time_horizon) / 52)
```

## File Organization Principles

- **Numbered analysis scripts** (`scripts/R/analysis/`): Designed to run sequentially, building on previous steps
- **Functions directory** (`scripts/R/functions/`): Reusable components that are sourced by analysis scripts
- **Tests directory** (`scripts/R/tests/`): Validation and diagnostic scripts
- **Quarto reports** (`scripts/QMD/report/`): Publication-ready PDF reports with embedded R code
- **RMarkdown files** (`.Rmd`): Additional tables and visualizations for reporting

## Quarto Report Architecture

The Quarto reports in `scripts/QMD/report/` are self-contained documents that:
1. Set their working directory to the project root using `here::here()`
2. Load all required packages and helper functions
3. Source the necessary analysis scripts (02, 03, etc.) to recreate the analysis environment
4. Generate formatted tables, plots, and results
5. Output to PDF with consistent styling

### Report Dependencies & Execution Order

Each report has specific dependencies:

**[para_models.qmd](scripts/QMD/report/para_models.qmd)** - Parametric Survival Modeling
- **Sources**: 02, 03, 06
- **Shows**: Survival model fits, AIC/BIC comparisons, goodness-of-fit diagnostics
- **Note**: Forces `USE_BOTH_MODELS <- 1` to compare full and reduced models

**[clinical_effectiveness.qmd](scripts/QMD/report/clinical_effectiveness.qmd)** - Clinical Effectiveness Analysis
- **Sources**: 02, 03, 06, 07
- **Shows**: Baseline characteristics, survival curves, life-years gained

**[input_parameters.qmd](scripts/QMD/report/input_parameters.qmd)** - Input Parameters Summary
- **Sources**: 02, 03, 07
- **Shows**: All model input parameters (costs, utilities, prevalence rates, treatment schedules)

**[CEA.qmd](scripts/QMD/report/CEA.qmd)** - Cost-Effectiveness Analysis Report
- **Sources**: 02, 03, 05, 06, 07, 08, 09, 10
- **Requires**: Base case analysis results from script 10
- **Shows**: Incremental cost-effectiveness ratios (ICERs), cost-effectiveness plane, decision tables

**[OWSA.qmd](scripts/QMD/report/OWSA.qmd)** - One-Way Sensitivity Analysis (Deterministic)
- **Sources**: 02, 03, 06, 07, 08, 11
- **Requires**: DSA results from script 11
- **Shows**: Tornado diagrams, one-way sensitivity plots for all varied parameters

**[EVPPIs.qmd](scripts/QMD/report/EVPPIs.qmd)** - Value of Information Analysis
- **Sources**: 02, 03, 06, 07, 08, 12, 13
- **Requires**: PSA results (script 12) and EVPPI results (script 13)
- **Shows**: Expected value of perfect information (EVPI), expected value of perfect partial information (EVPPI) for parameter groups

**[scenario_effect.qmd](scripts/QMD/report/scenario_effect.qmd)** - Scenario Analysis
- **Sources**: 02, 03, 06, 07, 08, scenario analysis scripts
- **Shows**: Alternative scenario results (e.g., different time horizons, discount rates)

### Quarto Report Structure Pattern

All reports follow a consistent pattern:

```r
# 1. Setup chunk - global knitr options, working directory, theme
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, ...)
knitr::opts_knit$set(root.dir = here::here())

# 2. Load packages
pacman::p_load(here, knitr, kableExtra, ggplot2, ...)

# 3. Set ggplot theme for consistency
theme_set(theme_minimal() + theme(...))

# 4. Source analysis scripts
source(here("scripts/R/analysis/02_setup_and_global_variables.R"))
source(here("scripts/R/analysis/03_biomarker_strategies.R"))
# ... others as needed

# 5. Load or verify results objects
if (!exists("cea_results")) {
  # Either load from file or stop with error
}

# 6. Generate tables and figures using kableExtra/ggplot
```

### Key Quarto Report Features

**Self-Contained Execution**: Reports can be rendered independently if the prerequisite analysis scripts have been run, as they re-source all dependencies.

**Consistent Styling**: All reports use:
- A4 paper format
- 11pt font
- Numbered sections with TOC depth 3
- Figure captions below, table captions above
- Minimal ggplot2 theme with customized text sizes

**Error Handling**: Reports use `tryCatch()` blocks when loading data to provide informative error messages if prerequisites are missing.

**Output Location**: PDFs are generated in the same directory as the `.qmd` files (`scripts/QMD/report/`).

### Common Quarto Report Issues

1. **"Object not found" errors**: Run the required analysis scripts first (especially 02, 03, 06-13)
2. **Sampling cache missing**: Run [08_sampling.R](scripts/R/analysis/08_sampling.R) to generate sampling models
3. **Rendering hangs**: Some reports (especially CEA, EVPPI) may take minutes to render due to re-sourcing analysis scripts
4. **USE_BOTH_MODELS conflict**: The `para_models.qmd` report overrides this to 1 for comparison purposes; other reports respect the global setting

## Cache Management

The analysis uses two cache systems to speed up computation:

### 1. Sampling Cache (Survival Models)

**Location**: `data/tidy/sampling_models_n*.rds`
**Purpose**: Cached correlated PFS/OS survival model fits
**Generation**: Script [08_sampling.R](scripts/R/analysis/08_sampling.R) (~first run takes time)
**Size**: Hundreds of MB
**When to regenerate**: Delete cache file when:
- Survival model formulas change
- `n_samples` changes
- `USE_BOTH_MODELS` setting changes
- Clinical data is updated

**To regenerate**:
```r
# Delete cache file
cache_file <- here("data", "tidy",
                   paste0("sampling_models_n", n_samples, "_full.rds"))
file.remove(cache_file)
# Re-run 08_sampling.R
source("scripts/R/analysis/08_sampling.R")
```

### 2. PSA Cache (Analysis Results)

**Location**: `data/tidy/psa_obj.rds` and `psa_params.rds`
**Purpose**: Cached PSA simulation results (5000 runs)
**Generation**: Script [12_PSA.R](scripts/R/analysis/12_PSA.R) (~20-60 minutes first run)
**Size**: ~660 KB total
**When to regenerate**: Delete cache files when:
- Model structure changes ([model_fun.R](scripts/R/functions/model_fun.R))
- Base parameters change
- Sampling cache is regenerated
- Parameter distributions change

**To regenerate**:
```r
# Delete PSA cache files
file.remove(here("data", "tidy", "psa_obj.rds"))
file.remove(here("data", "tidy", "psa_params.rds"))
# Re-run 12_PSA.R
source("scripts/R/analysis/12_PSA.R")
```

### Cache Workflow

1. **First run**: [08_sampling.R](scripts/R/analysis/08_sampling.R) generates sampling cache
2. **PSA uses sampling cache**: [12_PSA.R](scripts/R/analysis/12_PSA.R) generates PSA cache
3. **Subsequent runs**: Both load from cache (fast)
4. **Manual invalidation**: Delete specific cache file(s) to regenerate

## Clinical Context

Based on the METIMMOX clinical trial (NCT03388190) comparing alternating chemotherapy + immunotherapy versus chemotherapy alone in MSS/pMMR metastatic colorectal cancer patients. The analysis focuses on identifying biomarker-based subgroups that may benefit from adding immunotherapy to standard chemotherapy.
