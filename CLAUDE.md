# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

METIMMOX-1 is a cost-effectiveness analysis comparing biomarker-guided immunotherapy strategies for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer. The analysis uses a partitioned survival model implemented in R to evaluate three biomarker strategies (CRP, TLR, TMB/BRAF) against standard of care.

**Target Audience**: Health economists and researchers developing decision-analytic models in R.

## Version Control

This repository uses Git for version control. **Important**: Claude Code should NOT commit changes or push to remote repositories. All Git operations (commits, pushes, branch management, pull requests) are the responsibility of the user.

Claude Code may:
- Read Git status and history for context
- Create or modify files as part of analysis workflows

Claude Code should NOT:
- Create commits
- Push changes to remote repositories
- Modify Git configuration
- Create or manage branches

## R Installation

R is installed at: `C:\Program Files\R\R-4.3.2\bin\x64\`

When running R scripts from the command line (e.g., via Bash tool), use:
```bash
"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" -e "your R code here"
```

**Note**: For interactive development, users typically work in RStudio rather than command line.

## GitHub CLI

GitHub CLI (`gh`) is installed and available for viewing issues, pull requests, and other repository information.

### Common GitHub CLI Commands

```bash
# View a specific issue (e.g., issue #42)
gh issue view 42

# List all open issues
gh issue list --state open

# List all issues (open and closed)
gh issue list --state all --limit 50

# Search for issues by keyword
gh issue list --search "bug" --state all

# View issue with comments
gh issue view 42 --comments

# View pull request
gh pr view 123

# Get repository information
gh repo view
```

**Important for Claude Code**: When a user references an issue number (e.g., "issue #42"), use `gh issue view <number>` to fetch the full issue description and context before proceeding with the fix.

## Essential Commands

### R Sourcing Patterns

R analysis scripts depend on functions defined in `scripts/R/functions/`. When running model code, you must source dependencies in the correct order:

#### Minimal Test Setup
```r
# For testing model_fun() and calculate_outcomes()
source("scripts/R/analysis/02_setup_and_global_variables.R")  # Global vars: time_horizon, cl, dr, etc.
source("scripts/R/analysis/03_biomarker_strategies.R")         # Biomarker definitions
source("scripts/R/analysis/06_parametric_survival_analysis.R") # Fit survival models
source("scripts/R/analysis/07_basecase_input_parameters.R")    # Parameter list: l_params_base

# Source required functions
source("scripts/R/functions/model_fun.R")
source("scripts/R/functions/calculate_outcomes.R")

# Now you can run the model
result <- model_fun(l_params_base, determpsa = "det", return_traces = FALSE)
```

#### Common Function Dependencies
- **model_fun.R** requires:
  - `calculate_outcomes.R` (cost/QALY calculations)
  - `prediction_functions.R` (for PSA mode with resampled models)

- **PSA/EVPPI scripts** (12, 13) require:
  - `model_fun.R`
  - `calculate_outcomes.R`
  - `psa_functions.R` and/or `evppi_functions.R`

#### From Command Line
```bash
# Test a fix by sourcing all dependencies
"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" -e "
  setwd('c:/Users/benjampg/git/METIMMOX-1');
  source('scripts/R/analysis/02_setup_and_global_variables.R');
  source('scripts/R/analysis/03_biomarker_strategies.R');
  source('scripts/R/analysis/06_parametric_survival_analysis.R');
  source('scripts/R/analysis/07_basecase_input_parameters.R');
  source('scripts/R/functions/model_fun.R');
  source('scripts/R/functions/calculate_outcomes.R');
  cat('Testing model_fun...\n');
  result <- model_fun(l_params_base, determpsa = 'det');
  cat('Success! Control cost:', result[['Cost']][1], '\n');
"
```

**Key principle**: Always check which analysis scripts source which function files (use `grep "source.*functions" scripts/R/analysis/*.R`) to understand dependencies.

### Running the Analysis

The numbered analysis scripts must be executed sequentially:

```r
# Core setup (run these first)
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")

# Survival analysis
source("scripts/R/analysis/06_parametric_survival_analysis.R")
source("scripts/R/analysis/07_basecase_input_parameters.R")

# Survival resampling (generates cache - takes time on first run)
source("scripts/R/analysis/08_sampling.R")

# Model execution
source("scripts/R/analysis/09_traces.R")
source("scripts/R/analysis/10_basecase_analysis.R")
source("scripts/R/analysis/10b_enriched_population_analysis.R")  # Optional: enriched population CEA

# Sensitivity analyses
source("scripts/R/analysis/11_DSA.R")  # Deterministic sensitivity analysis
source("scripts/R/analysis/12_PSA.R")  # Probabilistic sensitivity analysis
source("scripts/R/analysis/13_EVPPIs.R")  # Expected value of perfect partial information

# Extended analyses (optional)
source("scripts/R/analysis/14_scenario_EVPPIs.R")  # Scenario-based EVPPI analysis
source("scripts/R/analysis/14b_scenario_preview.R")  # Optional: quick scenario preview (500 iterations)
source("scripts/R/analysis/15_save_snapshot.R")    # Save results for impact assessment
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
quarto render scripts/QMD/report/biosimilar_scenario.qmd
quarto render scripts/QMD/report/enriched_population.qmd

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
source("scripts/R/analysis/06_parametric_survival_analysis.R")
source("scripts/R/analysis/07_basecase_input_parameters.R")
source("scripts/R/analysis/08_sampling.R")  # Takes time on first run
source("scripts/R/analysis/09_traces.R")
source("scripts/R/analysis/10_basecase_analysis.R")
source("scripts/R/analysis/10b_enriched_population_analysis.R")  # Optional: enriched population CEA
source("scripts/R/analysis/11_DSA.R")
source("scripts/R/analysis/12_PSA.R")
source("scripts/R/analysis/13_EVPPIs.R")
source("scripts/R/analysis/14_scenario_EVPPIs.R")  # Optional: scenario analysis
source("scripts/R/analysis/14b_scenario_preview.R")  # Optional: quick scenario preview (500 iterations)
source("scripts/R/analysis/15_save_snapshot.R")    # Optional: save for comparison

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
quarto render scripts/QMD/report/biosimilar_scenario.qmd
quarto render scripts/QMD/report/enriched_population.qmd
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

| Variable | Default | Description |
|----------|---------|-------------|
| `cl` | 1/52 | Cycle length (1 week) |
| `time_horizon` | 520 | 10 years in weeks |
| `WTP` | 51000 | Willingness-to-pay threshold (EUR) |
| `n_samples` | 5000 | Resampling/PSA sample size |
| `dr` | 0.04 | Discount rate (4%) |
| `USE_BOTH_MODELS` | 0 | 0=full model only, 1=compare both |
| `MODEL_STRUCTURE` | 0 | 0=joint (Model A), 1=focused (Model B), 2=separate (Model C) |
| `annual_incidence_norway` | 1500 | Annual eligible MSS/pMMR mCRC patients in Norway |
| `research_horizon_years` | 10 | Research value time horizon (years) for population EVPPI |
| `discount_rate_research` | 0.035 | Discount rate for research benefits (3.5%) |

**Model Structure Options:**
- **Model A (joint, 0)**: All biomarkers + all treatment interactions: `~ Age + sex + Rx + crp*Rx + tlr*Rx + tmb_braf*Rx`
- **Model B (focused, 1)**: Per-strategy formulas. CRP and TMB/BRAF strategies share: `~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx`. TLR strategy uses: `~ Age + sex + Rx + crp + tmb_braf + tlr:Rx`. All three biomarker strategies available.
- **Model C (separate, 2)**: One biomarker + its interaction only: `~ Age + sex + Rx + [biomarker]:Rx`

**When changed**: Regenerate sampling cache (USE_BOTH_MODELS or MODEL_STRUCTURE change) and PSA cache. Cache filenames encode these settings to prevent mixing results.

### Biomarker Strategies

Three biomarkers are evaluated (defined in [03_biomarker_strategies.R](scripts/R/analysis/03_biomarker_strategies.R)):

1. **CRP** (C-reactive protein): Binary variable, cut-off <5
2. **TLR** (Tumor lesion reduction): Binary variable, cut-off ≥10%
3. **TMB/BRAF**: Combined biomarker (TMB ≥9 mut/MB OR BRAF mutation)

**Note**: All three biomarker strategies are available in all model structures. Model B uses per-strategy formulas (CRP/TMB_BRAF share one formula, TLR has its own). Model C uses per-biomarker formulas (one biomarker + its interaction per strategy). See `model_configs.R` for the single source of truth.

Each strategy has:
- **Biomarker-positive subgroup**: Receives experimental treatment (alternating FLOX + nivolumab)
- **Biomarker-negative subgroup**: Receives standard treatment (FLOX only)
- **Analysis approach**: See "Survival Prediction Methodologies" section below for how population-level outcomes are calculated

### Survival Prediction Methodologies

The model uses different approaches for base case vs PSA, both properly accounting for patient heterogeneity:

**Base Case** (Issue #69): Population averaging across all patients using `generate_population_averaged_predictions()` in [prediction_functions.R](scripts/R/functions/prediction_functions.R:199-340). Predicts for all patients using their actual age, sex, and biomarker values, then averages.

**PSA** (Issue #70): Second-order Monte Carlo sampling one patient per iteration. See detailed methodology in [model_fun.R](scripts/R/functions/model_fun.R:17-148) comments. Over 5000 iterations, the distribution of sampled patients correctly propagates uncertainty from patient heterogeneity.

Both approaches are methodologically valid for their respective analytical purposes (deterministic point estimates vs. probabilistic uncertainty quantification).

### Survival Resampling & Correlation

The model uses **correlated survival resampling** ([08_sampling.R](scripts/R/analysis/08_sampling.R:103-200)) to maintain the correlation between PFS and OS:

- Both PFS and OS models are fitted to the **same resampled patient cohort**
- Results are cached in `data/tidy/`
- Cache file naming: `sampling_models_n{n_samples}_{full|both}_{joint|focused|separate}.rds`

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

**PSA Mode Behavior** (lines 17-148):
- Uses second-order Monte Carlo: samples one patient per iteration from resampled datasets (see "Survival Prediction Methodologies" section)
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

**Control group** (age- and sex-adjusted, same across all model structures):
```r
OS:  Surv(OSwk, Death) ~ Age + sex
PFS: Surv(PFSwk, Progression) ~ Age + sex
```

**Biomarker groups** (formulas vary by model structure):

- **Model A (joint)**: Single shared formula for all strategies: `~ Age + sex + Rx + crp*Rx + tlr*Rx + tmb_braf*Rx`
- **Model B (focused)**: Per-strategy formulas:
  - CRP / TMB_BRAF strategies: `~ Age + sex + Rx + crp + tmb_braf + crp:Rx + tmb_braf:Rx`
  - TLR strategy: `~ Age + sex + Rx + crp + tmb_braf + tlr:Rx`
- **Model C (separate)**: Per-biomarker: `~ Age + sex + Rx + [biomarker]:Rx`

See [model_configs.R](scripts/R/functions/model_configs.R) for the canonical formula definitions.

### Key Functions

**Model Configuration**:
- **[model_configs.R](scripts/R/functions/model_configs.R)**: Single source of truth for model structures, strategies, and formulas. Auto-sourced by `02_setup_and_global_variables.R`. Key functions: `get_model_configs()`, `get_strategies()`, `get_biomarkers()`, `get_strategy_formula()`, `get_control_formula()`, `uses_per_strategy_formulas()`, `get_model_type_label()`, `get_model_formulas()`.

**Core Model Functions**:
- **[model_fun.R](scripts/R/functions/model_fun.R)**: Main partitioned survival model with PSA support
- **[calculate_outcomes.R](scripts/R/functions/calculate_outcomes.R)**: Calculates QALYs and costs from state occupancy traces
- **[prediction_functions.R](scripts/R/functions/prediction_functions.R)**: Generate survival predictions from fitted models

**Sensitivity Analysis Functions**:
- **[psa_functions.R](scripts/R/functions/psa_functions.R)**: PSA-related utilities
- **[evppi_functions.R](scripts/R/functions/evppi_functions.R)**: EVPPI calculation functions, population-level EVPPI scaling (`calculate_population_evppi()`), and interaction coefficient EVPPI (`extract_interaction_coefficients()`, `get_interaction_evppi_params()`, `add_interaction_param_groups()`)
- **[scenario_analysis.R](scripts/R/functions/scenario_analysis.R)**: Scenario analysis framework

**Survival Modeling Functions**:
- **[para_model_fit.R](scripts/R/functions/para_model_fit.R)**: Parametric model fitting helper
- **[para_model_fit_table.R](scripts/R/functions/para_model_fit_table.R)**: Model fit summary tables
- **[survival_plots.R](scripts/R/functions/survival_plots.R)**: Survival curve visualization

**Visualization Functions**:
- **[create_tornado_plot.R](scripts/R/functions/create_tornado_plot.R)**: DSA tornado diagram generation

**Snapshot/Impact Assessment Functions**:
- **[snapshot_utils.R](scripts/R/functions/snapshot_utils.R)**: Snapshot management for bug fix impact assessment

**Archived Functions** (in `scripts/R/archive/`):
- `bootstrap_survival_model.R`: Alternative resampling approach (not used in main analysis)
- `ref_values_emm.R`: Reference value calculations (superseded)

### Multi-Model CEA Functions

The **[multi_model_cea.R](scripts/R/functions/multi_model_cea.R)** file provides functions for comparing cost-effectiveness across Model A/B/C structures:

**Core Functions:**
- `get_model_configs()` - Returns Model A (joint), B (focused), C (separate) configurations (canonical source is [model_configs.R](scripts/R/functions/model_configs.R))
- `run_all_basecase_analyses()` - Runs base case for all 3 model structures
- `run_basecase_for_structure()` - Runs base case for a single model structure (also used by `10b` and `15`)
- `load_all_psa_caches()` - Loads PSA caches for all models
- `format_multimodel_comparison_table()` - Creates side-by-side comparison tables
- `create_multimodel_ceac_plot()` - Multi-model CEAC visualization

**Workflow:** Used by CEA.qmd to generate cross-model comparisons. Requires PSA caches for each MODEL_STRUCTURE to be pre-generated via [12_PSA.R](scripts/R/analysis/12_PSA.R).

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

### Discount Factor Application

Discount factors are applied using vector multiplication over the time horizon:
```r
v_dw_c <- 1 / (1 + dr_costs)^(seq(0, time_horizon) / 52)
v_dw_e <- 1 / (1 + dr_effects)^(seq(0, time_horizon) / 52)
```

## File Organization Principles

- **Numbered analysis scripts** (`scripts/R/analysis/`): Designed to run sequentially, building on previous steps
  - Scripts 01-03: Core setup and data preparation
  - Scripts 04-05: Supplementary RMarkdown files (not part of main execution pipeline)
  - Scripts 06-13: Main analysis pipeline
  - Script 10b: Enriched population CEA (optional, runs across all model structures)
  - Scripts 14-15: Extended analyses (scenario EVPPIs, snapshot saving)
  - Script 14b: Quick scenario preview (500 iterations, Model A only)
- **Functions directory** (`scripts/R/functions/`): Reusable components that are sourced by analysis scripts
- **Archive directory** (`scripts/R/archive/`): Deprecated/unused code preserved for reference
- **Tests directory** (`scripts/R/tests/`): Validation and diagnostic scripts
- **Quarto reports** (`scripts/QMD/report/`): Publication-ready PDF reports with embedded R code
- **Technical docs** (`scripts/QMD/technical_docs/`): Bug fix impact reports and technical documentation

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
- **Models displayed**: Best-fit model (gamma) AND Weibull PH model for reference (issue #68)
- **Tables include**: Model parameter exponents for clinical interpretation
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

**[biosimilar_scenario.qmd](scripts/QMD/report/biosimilar_scenario.qmd)** - Biosimilar Nivolumab Pricing Scenario
- **Sources**: 02, 03, scenario cache (`scenario_evppi_results.rds` or `scenario_evppi_results_PREVIEW.rds`)
- **Shows**: ICER comparison for base case vs biosimilar pricing (EUR 13,923 vs EUR 4,641/dose) across Models A/B/C
- **Note**: Auto-detects full vs preview cache

**[enriched_population.qmd](scripts/QMD/report/enriched_population.qmd)** - Enriched Population Analysis
- **Sources**: 02, 03, 10b
- **Shows**: Enriched (biomarker-positive) ICERs vs base case across Models A/B/C
- **Note**: Not cached; re-runs on each render. Requires sampling cache.

**[biomarker_decomposition.qmd](scripts/QMD/report/biomarker_decomposition.qmd)** - Biomarker Effect Decomposition
- **Sources**: 02, 03, 06, 07, 08, 10, 12
- **Shows**: Decomposition of biomarker effects on cost-effectiveness outcomes

**[biomarker_distributions.qmd](scripts/QMD/report/biomarker_distributions.qmd)** - Biomarker Distributions
- **Sources**: 02, 03
- **Shows**: Biomarker prevalence and distribution analyses

**Technical Documentation** (in `scripts/QMD/technical_docs/`):
- **[age_effect_analysis.qmd](scripts/QMD/technical_docs/age_effect_analysis.qmd)**: Age effect on survival outcomes
- **[all_parametric_survival_models.qmd](scripts/QMD/technical_docs/all_parametric_survival_models.qmd)**: Full survival model diagnostics
- **[bug_fix_impact.qmd](scripts/QMD/technical_docs/bug_fix_impact.qmd)**: Bug fix impact documentation (v3.0, per-model comparisons with cross-model summary)

**Figure Vignettes** (in `scripts/QMD/vignettes/`):
- **figure1-4.qmd**: Publication-ready figures
- **suppl_figure_pfs_plots.qmd**: Supplementary PFS figures

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

### Report Formatting Standards

All Quarto reports must follow these formatting conventions for consistency.

**YAML Header Template** (complete example):
```yaml
---
title: "Clinical Effectiveness"
subtitle: "METIMMOX-1 Economic Evaluation"
author: "Ben Geisler"
date: "`r Sys.Date()`"
format:
  pdf:
    toc: true
    toc-depth: 3
    number-sections: true
    fig-width: 8
    fig-height: 6
    keep-tex: false
    geometry:
      - top=2.5cm
      - bottom=2.5cm
      - left=2.5cm
      - right=2.5cm
    fontsize: 11pt
    documentclass: article
    classoption: [a4paper]
execute:
  echo: false
  warning: false
  message: false
  cache: false
  fig-cap-location: bottom
  tbl-cap-location: top
---
```

**Title conventions**:
- **Title**: Simple, descriptive name matching the report content (e.g., "Clinical Effectiveness", "Cost-Effectiveness Analysis", "Input Parameters")
- **Subtitle**: Always "METIMMOX-1 Economic Evaluation"
- **No fancy LaTeX headers**: Do not use `include-in-header` with fancyhdr or custom footers

**Report signature** (at the very end of the document):
```markdown
---

**Report completed on:** `r Sys.Date()`\
**Repository:** ben-geisler/METIMMOX-1\
**Report version:** X.X
```

Key points for the signature:
- Horizontal rule (`---`) as separator before the signature
- Three lines: completion date, repository, version
- Use backslash (`\`) at end of each line for proper line breaks in PDF
- Left-justified (default markdown alignment)
- No trailing content after the signature

**Footer versioning**: Reports include version numbers (e.g., 2.0, 3.0) corresponding to major methodology updates. Increment version when making significant changes to analysis methodology.

**Special characters**: Avoid Unicode Greek letters (e.g., α, β, λ) in text that will render to PDF. Instead, spell out the word (e.g., "alpha = 0" instead of "α = 0") or use LaTeX math mode (`$\alpha$`) if mathematical formatting is needed.

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
5. **Results changed after methodological updates**: If cost-effectiveness results differ from earlier versions, check if methodological fixes were applied. Issues #69 and #70 (Nov 2024) changed survival prediction methodology from reference patient to population averaging/individual sampling. This **should** change results - it's a methodological improvement. Regenerate both sampling cache and PSA cache after these fixes. See issue #64 for impact documentation approach.

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
# Delete cache file (example for joint model structure)
cache_file <- here("data", "tidy",
                   paste0("sampling_models_n", n_samples, "_full_joint.rds"))
file.remove(cache_file)
# Re-run 08_sampling.R
source("scripts/R/analysis/08_sampling.R")
```

### 2. PSA Cache (Analysis Results)

**Location**: `data/tidy/psa_obj_{joint|focused|separate}.rds` and `psa_params_{joint|focused|separate}.rds`
**Purpose**: Cached PSA simulation results (5000 runs)
**Generation**: Script [12_PSA.R](scripts/R/analysis/12_PSA.R) (~20-60 minutes first run)
**Size**: ~660 KB total
**When to regenerate**: Delete cache files when:
- Model structure changes ([model_fun.R](scripts/R/functions/model_fun.R) or [calculate_outcomes.R](scripts/R/functions/calculate_outcomes.R))
- **Prediction methodology changes** (e.g., issues #69, #70 fixes to survival curve generation)
- Base parameters change (costs, utilities, time horizon, discount rates)
- **Sampling cache is regenerated** (PSA depends on specific resampled models - always regenerate PSA after regenerating sampling cache)
- Parameter distributions change (distributional assumptions, means, SDs, correlations)

**To regenerate**:
```r
# Delete PSA cache files (example for joint model structure)
file.remove(here("data", "tidy", "psa_obj_joint.rds"))
file.remove(here("data", "tidy", "psa_params_joint.rds"))
# Re-run 12_PSA.R
source("scripts/R/analysis/12_PSA.R")
```

### 3. EVPPI Cache

**Location**: `data/tidy/evppi_results.RData`
**Purpose**: Cached EVPPI results for parameter groups
**Generation**: Script [13_EVPPIs.R](scripts/R/analysis/13_EVPPIs.R)
**When to regenerate**: Delete cache file when PSA cache is regenerated or EVPPI parameter groupings change

### 4. Scenario EVPPI Cache

**Location**: `data/tidy/scenario_evppi_results.RData`
**Purpose**: Cached scenario-based EVPPI analysis results
**Generation**: Script [14_scenario_EVPPIs.R](scripts/R/analysis/14_scenario_EVPPIs.R)
**When to regenerate**: Delete cache file when PSA cache is regenerated or scenario definitions change

### 5. Scenario Preview Cache

**Location**: `data/tidy/scenario_evppi_results_PREVIEW.rds`
**Purpose**: Quick-test scenario results (500 iterations, Model A only)
**Generation**: Script [14b_scenario_preview.R](scripts/R/analysis/14b_scenario_preview.R)
**Note**: Used by `biosimilar_scenario.qmd` as fallback when full scenario cache is absent

### Cache Workflow

1. **First run**: [08_sampling.R](scripts/R/analysis/08_sampling.R) generates sampling cache
2. **PSA uses sampling cache**: [12_PSA.R](scripts/R/analysis/12_PSA.R) generates PSA cache
3. **EVPPI uses PSA cache**: [13_EVPPIs.R](scripts/R/analysis/13_EVPPIs.R) generates EVPPI cache
4. **Subsequent runs**: All load from cache (fast)
5. **Manual invalidation**: Delete specific cache file(s) to regenerate

## Snapshot System (Bug Fix Impact Assessment)

The repository includes a snapshot comparison system for assessing the impact of bug fixes on model results.

### Components

- **[15_save_snapshot.R](scripts/R/analysis/15_save_snapshot.R)**: Saves multi-model snapshots (all 3 model structures A/B/C) using `run_basecase_for_structure()` from `multi_model_cea.R`
- **[snapshot_utils.R](scripts/R/functions/snapshot_utils.R)**: Utility functions for snapshot management, including `is_multimodel_snapshot()` and `get_snapshot_model_names()`
- **[compare_snapshots.R](scripts/R/tests/compare_snapshots.R)**: Compares before/after snapshots to quantify changes
- **[bug_fix_impact.qmd](scripts/QMD/technical_docs/bug_fix_impact.qmd)**: Report documenting bug fix impacts

### Snapshot Storage

**Location**: `data/output/snapshots/`
**Format**: `snapshot_NN_[baseline/fixed]_HASH.rds` and `psa_NN_[baseline/fixed]_HASH.rds`
**Content**: Base case results, PSA results, metadata (git commit, timestamp, issue number)

### Workflow

1. Save "baseline" snapshot before applying a fix
2. Apply fix
3. Save "fixed" snapshot after applying a fix
4. Run `compare_snapshots.R` to quantify the impact
5. Render `bug_fix_impact.qmd` to document changes

## Test Files

The test suite in `scripts/R/tests/` includes:

- **[test_sex_variable_fix.R](scripts/R/tests/test_sex_variable_fix.R)**: Validates Issue #72 fix (sex variable type consistency in resampled models)
- **[compare_snapshots.R](scripts/R/tests/compare_snapshots.R)**: Compares before/after snapshots for bug fix impact assessment
- **[para_models.Rmd](scripts/R/tests/para_models.Rmd)**: Parametric model fit validation and diagnostics
- **[snapshot.R](scripts/R/tests/snapshot.R)**: Helper script for running snapshot saves

**Diagnostic Scripts** (for troubleshooting PSA/sampling issues):
- **diagnose_prediction_failures.R**: Analyzes why PSA iterations fail
- **test_psa_error_rate.R**: Quantifies PSA iteration error rates
- **test_sampling_convergence_rate.R**: Tests sampling convergence

## Clinical Context

Based on the METIMMOX clinical trial (NCT03388190) comparing alternating chemotherapy + immunotherapy versus chemotherapy alone in MSS/pMMR metastatic colorectal cancer patients. The analysis focuses on identifying biomarker-based subgroups that may benefit from adding immunotherapy to standard chemotherapy.
