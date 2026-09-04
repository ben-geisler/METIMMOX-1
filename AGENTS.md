## Project Overview

METIMMOX-1 is a cost-effectiveness analysis comparing biomarker-guided immunotherapy strategies for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer. The analysis uses a partitioned survival model implemented in R to evaluate two pre-immunotherapy biomarker strategies (week-4 CRP and baseline TMB/BRAF) against standard of care. TLR is retained only for clinical effectiveness, DAG, and biomarker distribution analyses because it is a post-randomization mediator measured on treatment rather than a treatment-selection biomarker.

**Target Audience**: Health economists and researchers developing decision-analytic models in R.

## Version Control

This repository uses Git for version control. The agent may create local commits as part of analysis workflows. Pushing to remote repositories remains the responsibility of the user unless the user explicitly requests it.

The agent may:
- Read Git status and history for context
- Create or modify files as part of analysis workflows
- Create local commits

The agent should NOT:
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

**Important for the agent**: When a user references an issue number (e.g., "issue #42"), use `gh issue view <number>` to fetch the full issue description and context before proceeding with the fix.

## Essential Commands

### R Sourcing Patterns

R analysis scripts depend on functions defined in `scripts/R/functions/`. When running model code, you must source dependencies in the correct order:

#### Minimal Test Setup
```r
# For testing model_fun() and calculate_outcomes()
source("scripts/R/analysis/02_setup_and_global_variables.R")  # Global vars: time_horizon, cl, dr, etc.
source("scripts/R/analysis/03_biomarker_strategies.R")         # Biomarker definitions
source("scripts/R/analysis/04_parametric_survival_analysis.R") # Fit survival models
source("scripts/R/analysis/05_basecase_input_parameters.R")    # Parameter list: l_params_base

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
  source('scripts/R/analysis/04_parametric_survival_analysis.R');
  source('scripts/R/analysis/05_basecase_input_parameters.R');
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
source("scripts/R/analysis/04_parametric_survival_analysis.R")
source("scripts/R/analysis/05_basecase_input_parameters.R")

# Survival resampling (generates cache - takes time on first run)
source("scripts/R/analysis/06_sampling.R")

# Model execution
source("scripts/R/analysis/07_traces.R")
source("scripts/R/analysis/08_basecase_analysis.R")
source("scripts/R/analysis/08b_enriched_population_analysis.R")  # Optional: enriched population CEA for economic biomarkers

# Sensitivity analyses
source("scripts/R/analysis/09_DSA.R")  # Deterministic sensitivity analysis
source("scripts/R/analysis/10_PSA.R")  # Probabilistic sensitivity analysis
source("scripts/R/analysis/11_EVPPIs.R")  # Expected value of perfect partial information

# Extended analyses (optional)
source("scripts/R/analysis/12_scenario_EVPPIs.R")  # Scenario-based EVPPI analysis
# 13_save_snapshot.R: optional, interactive, standalone -- run via Rscript with <issue#> <baseline|fixed>, not sourced
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

Output format depends on report type, and this changes the render command:

- **Clinical/DAG/descriptive** reports (`clinical_effectiveness`, `dag`, `dag_associations`, `biomarker_distributions`, `survival_model_specification`) declare `format:` with both `pdf:` and `gfm:` and render cleanly to **both** a `.pdf` and a readable `.md` (e.g. `scripts/QMD/report/clinical_effectiveness.md`).
- **Economic** reports (`CEA`, `OWSA`, `EVPPIs`, `scenario_effect`, `biosimilar_scenario`, `enriched_population`, `biomarker_decomposition`, `input_parameters`, `para_models`) plus the survival technical docs use kableExtra HTML tables, so their **GFM pass fails** (`Functions that produce HTML output found in document targeting commonmark output`) and aborts the whole render, leaving a STALE `.pdf`. Render these with `--to pdf` and verify the text with `pdftotext` (no `.md` is produced).

```bash
# Clinical/DAG report -> PDF + MD (read the .md directly)
quarto render scripts/QMD/report/clinical_effectiveness.qmd

# Economic report -> PDF only, then verify content as text
quarto render scripts/QMD/report/CEA.qmd --to pdf
pdftotext scripts/QMD/report/CEA.pdf - | grep -ciE '\btlr\b'   # economic reports: expect 0
```

- Do NOT run `quarto render scripts/QMD/report/` (whole directory) — it fails on every economic report's GFM pass. Render economic reports one at a time with `--to pdf`.
- `pdftotext` (poppler) is at `/mingw64/bin`; `pdftoppm` (needed for the Read tool's visual PDF rendering) is NOT installed — verify PDFs with `pdftotext`, not by reading them directly.

**Prerequisites for rendering**:
- Run the analysis scripts required by the report (up to 11 for the full core analysis; script 12 additionally generates scenario outputs)
- Sampling cache must exist (from running [06_sampling.R](scripts/R/analysis/06_sampling.R))
- Results objects (e.g., `cea_results`, `owsa_results`, `psa_results`, `evppi_results`) must be in the R environment or saved as `.rds` files

### Complete Analysis & Reporting Workflow

To generate all analysis results and reports from scratch:

```r
# 1. Run all analysis scripts in order
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")
source("scripts/R/analysis/04_parametric_survival_analysis.R")
source("scripts/R/analysis/05_basecase_input_parameters.R")
source("scripts/R/analysis/06_sampling.R")  # Takes time on first run
source("scripts/R/analysis/07_traces.R")
source("scripts/R/analysis/08_basecase_analysis.R")
source("scripts/R/analysis/08b_enriched_population_analysis.R")  # Optional: enriched population CEA for economic biomarkers
source("scripts/R/analysis/09_DSA.R")
source("scripts/R/analysis/10_PSA.R")
source("scripts/R/analysis/11_EVPPIs.R")
source("scripts/R/analysis/12_scenario_EVPPIs.R")  # Optional: scenario analysis
# 13_save_snapshot.R: optional, interactive, standalone -- NOT sourced here.
# Run separately: Rscript scripts/R/analysis/13_save_snapshot.R <issue#> <baseline|fixed>

# 2. Render reports (from terminal/command line)
# quarto render scripts/QMD/report/
```

**Running the pipeline from a clean / non-interactive session** (e.g. driving it with `Rscript`):
- Attach packages first — only `01_data_prep.R` calls `p_load`, so `pacman::p_load(...)` the full set above before sourcing, or `07_traces.R` fails with `could not find function "ggplot"`.
- `07_traces.R` assumes `model_fun()` is already loaded (scripts 08-13 source it themselves); source `model_fun.R`, `calculate_outcomes.R`, `prediction_functions.R` before it.
- `06_sampling.R` loads `sampling_models_n{n}_full.rds` if present (fast); only `10_PSA.R` (~20 min) and `12_scenario_EVPPIs.R` (~40 min) are slow at n=5000.
- Scenario reports and publication vignettes require the full cache generated by `12_scenario_EVPPIs.R`.

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
| `UTILITY_SOURCE` | 1 | 0=IPD-derived (u_np=0.9077, u_p=0.9005), 1=CORRECT trial (u_np=0.73, u_p=0.59) |
| `annual_incidence_norway` | 1500 | Annual eligible MSS/pMMR mCRC patients in Norway |
| `research_horizon_years` | 10 | Research value time horizon (years) for population EVPPI |
| `discount_rate_research` | 0.035 | Discount rate for research benefits (3.5%) |

**Utility-source provenance**: Reported economic results were generated with `UTILITY_SOURCE = 1` (CORRECT trial utilities) and correspond to the `*_correct.rds` / `*_correct.RData` cache files; the manuscript's IPD values (approximately 0.91/0.90) are retained as a different, historical utility base case.

**Single economic survival model**:
```r
OS:  Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
Control arm: the same joint models predicted with Rx = control for every patient
```

There is no model-structure switch or multi-structure comparison layer. The economic model is the joint CRP + TMB/BRAF formula defined in [model_configs.R](scripts/R/functions/model_configs.R). There is no separate control-arm model: in both the base case and the PSA the standard-of-care curves are the joint model predicted with treatment set to control (issue #151). The sampling cache stores one joint bootstrap under each biomarker key; a legacy `control` component in an older cache is ignored.

**Parametric distribution selection**: OS and PFS distributions are selected jointly from the nine candidate families in [04_parametric_survival_analysis.R](scripts/R/analysis/04_parametric_survival_analysis.R). The selected pair is the minimum-combined-AIC pair that preserves OS >= PFS for control and every economic biomarker subgroup at every modeled weekly time point. With the current data and 10-year horizon, the ordering-constrained selection is **gamma for OS and gamma for PFS**.

**When changed**: Regenerate sampling cache when survival formulas, the economic strategy/biomarker set, `n_samples`, or clinical data change. Regenerate PSA and EVPPI caches after regenerating sampling cache or changing economic parameters, distributions, prediction methodology, or `UTILITY_SOURCE`.

### Biomarker Strategies

Two pre-immunotherapy biomarkers are evaluated economically (defined in [03_biomarker_strategies.R](scripts/R/analysis/03_biomarker_strategies.R) and selected by [model_configs.R](scripts/R/functions/model_configs.R)):

1. **CRP** (C-reactive protein): Binary variable, cut-off <5 mg/L, **week-4 value** (`CRP1cat`: cycle 3 day 1, trial visit 3), measured after the two FLOX cycles that both arms receive and before the first nivolumab dose
2. **TMB/BRAF**: Combined biomarker (TMB >=9 mut/MB OR BRAF mutation), from baseline NGS

**CRP timing (issue #150)**: CRP is not a baseline (pre-randomization) measurement, but it is available before the immunotherapy decision, because every patient receives two FLOX cycles before the first nivolumab dose. Describe it as "week-4 CRP, before the first nivolumab dose", never as "baseline" or "pre-treatment" CRP. The true baseline value (`CRP0cat`, cycle 1 day 1) is derived in `01_data_prep.R` but unused. Consequences: (i) week-4 CRP prevalence differs by arm (17/36 experimental vs 7/35 control, Fisher p = 0.023) whereas baseline CRP is balanced (p = 1); the DAG report reads this as chance arm imbalance on a week-4 measurement in a small trial, not as a randomisation failure, and the DAG has no `T -> CRP` edge because nivolumab has not been given at week 4 and both arms receive identical FLOX until then; (ii) the 3 patients without a week-4 CRP are early deaths (weeks 2.4, 15.7, 20.9) and are dropped from the complete-case data; (iii) the model clock starts at randomization (t = 0) while the CRP decision is made at week 4; this has no economic consequence because both arms receive identical FLOX, monitoring, and visits up to that point, and the week-4 blood test is already in the monitoring schedule. No model code or cache depends on this framing.

**TLR** (tumor lesion reduction) is intentionally excluded from all economic analyses because it is a post-randomization mediator measured on treatment, not a treatment-selection biomarker. `data$tlr` and `p_tlr` may still be created for clinical effectiveness, DAG, and biomarker distribution reports that analyze TLR directly from the trial data.

Each economic biomarker strategy has:
- **Biomarker-positive subgroup**: Receives experimental treatment (alternating FLOX + nivolumab)
- **Biomarker-negative subgroup**: Receives standard treatment (FLOX only)
- **Analysis approach**: See "Survival Prediction Methodologies" section below for how population-level outcomes are calculated

### Survival Prediction Methodologies

The model uses different approaches for base case vs PSA, both properly accounting for patient heterogeneity:

**Base Case** (Issue #69): Population averaging across all patients using `generate_population_averaged_predictions()` in [prediction_functions.R](scripts/R/functions/prediction_functions.R:199-340). Predicts for all patients using their actual age, sex, and biomarker values, then averages.

**PSA** (Issue #70): Second-order Monte Carlo sampling one patient per iteration. See detailed methodology in [model_fun.R](scripts/R/functions/model_fun.R:17-148) comments. Over 5000 iterations, the distribution of sampled patients correctly propagates uncertainty from patient heterogeneity.

Both approaches are methodologically valid for their respective analytical purposes (deterministic point estimates vs. probabilistic uncertainty quantification).

### Survival Resampling & Correlation

The model uses **correlated survival resampling** ([06_sampling.R](scripts/R/analysis/06_sampling.R:103-200)) to maintain the correlation between PFS and OS:

- Both PFS and OS models are fitted to the **same resampled patient cohort**
- One joint-model bootstrap serves every strategy: control, biomarker-positive, and biomarker-negative curves in a PSA draw all come from the same resampled joint fit (control = joint fit with `Rx = control`), so control and biomarker effects stay correlated (issue #151)
- Results are cached in `data/tidy/`
- Cache file naming: `sampling_models_n{n_samples}_full.rds`

**IMPORTANT**: The first run of `06_sampling.R` will take significant time (generates 5000 resampled models). Subsequent runs load from cache.

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
  p_crp, p_tmb_braf
)
```

### Data Requirements

**Confidential Trial Data**: The clinical dataset is NOT included in the repository. The model expects:
- File location: `data/tidy/METIMMOX.rds`
- Required variables: `ID`, `PFSwk`, `Progression`, `OSwk`, `Death`, `Rx`, `Age`, `sex`, biomarker variables

**PFS endpoint derivation** (issue #149, in [pfs_endpoint.R](scripts/R/functions/pfs_endpoint.R), called by [02_setup_and_global_variables.R](scripts/R/analysis/02_setup_and_global_variables.R)): the trial export records exit-for-progression only (`Progression exit`) and time to that exit (`Days until progression`), so patients who died without a recorded progression would be censored for PFS. `derive_pfs_endpoint()` recodes `Progression` as progression OR death and `PFSwk` as time to progression if progressed, time to death if died progression-free, and follow-up time otherwise, matching the trial definition (Ree et al. 2024) and the partitioned-survival-model requirement that the progression-free state means alive and progression-free. The raw values are kept as `ProgressionExit` and `TTPwk`. Any report that reads `data/tidy/METIMMOX.rds` directly (currently `biomarker_distributions.qmd`) must call `derive_pfs_endpoint()` after renaming `Progression exit`.

**Biomarker derivation** (in [03_biomarker_strategies.R](scripts/R/analysis/03_biomarker_strategies.R:6-9)):
```r
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))
```

`data$tlr` is kept for clinical effectiveness, DAG, and biomarker distribution reports only. It is not included in economic strategies, economic survival formulas, PSA parameters, EVPPI groups, DSA, scenario analyses, enriched-population CEA, or snapshots.

### Adverse-Event Scope

Adverse-event/toxicity costs and disutilities are not modeled separately. This is a deliberate scope choice based on the intended tolerability of the alternating short-course FLOX-nivolumab regimen and the lack of sufficiently robust treatment-specific trial data on adverse-event incidence, resource use, and utility decrements for economic parameterization.

### Survival Model Formulas

**Economic model** (shared by the control, CRP-guided, and TMB/BRAF-guided strategies):
```r
OS:  Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
```

**Control arm**: not a separate model. Standard-of-care curves are the joint OS and PFS models predicted with `Rx` set to the control level for every patient in `data_complete`, then averaged. This holds in the base case (`generate_population_averaged_predictions()`) and in every PSA draw (`generate_psa_population_averaged_predictions()` with `biomarker_name = NULL`), so PSA means are centred on the base case. An earlier age/sex-only control bootstrap fitted on the 36 control-arm patients was removed in issue #151.

See [model_configs.R](scripts/R/functions/model_configs.R) for the canonical formula definitions.

### Key Functions

**Model Configuration**:
- **[model_configs.R](scripts/R/functions/model_configs.R)**: Single source of truth for the economic strategies, biomarkers, and formulas. Auto-sourced by `02_setup_and_global_variables.R`. Key functions: `get_model_configs()`, `get_current_model_config()`, `get_strategies()`, `get_biomarkers()`, `get_strategy_formula()`, `get_model_formulas()`. (`get_control_formula()` was removed in issue #151; there is no separate control-arm formula.)
- **[pfs_endpoint.R](scripts/R/functions/pfs_endpoint.R)**: Single source of truth for the composite PFS endpoint (progression or death). Auto-sourced by `02_setup_and_global_variables.R`. Key function: `derive_pfs_endpoint()`.
- **[cache_paths.R](scripts/R/functions/cache_paths.R)**: Single source of truth for the utility-source label and cached-object file paths (sampling, PSA, EVPPI, scenario). Auto-sourced by `02_setup_and_global_variables.R`. Key functions: `resolve_util_label()`, `sampling_cache_path()`, `psa_obj_path()`, `psa_params_path()`, `evppi_path()`, `scenario_evppi_path()`.
- **[report_setup.R](scripts/R/functions/report_setup.R)**: One-call Quarto report setup (knitr options, package loading, shared ggplot theme, and sourcing of analysis scripts/function files), used to remove duplicated setup boilerplate across the economic reports. Key function: `setup_report(sources, funs, packages, set_theme)`.

**Core Model Functions**:
- **[model_fun.R](scripts/R/functions/model_fun.R)**: Main partitioned survival model with PSA support
- **[calculate_outcomes.R](scripts/R/functions/calculate_outcomes.R)**: Calculates QALYs and costs from state occupancy traces
- **[prediction_functions.R](scripts/R/functions/prediction_functions.R)**: Generate survival predictions from fitted models
- **[cea_helpers.R](scripts/R/functions/cea_helpers.R)**: Single-model CEA execution and summary helpers wrapping dampack (`run_basecase()`, `load_psa_cache()`, `create_ceac_plot()`, `create_psa_summary_table()`). Renamed from the legacy `multi_model_cea.R`.
- **[eq5d5l_utility.R](scripts/R/functions/eq5d5l_utility.R)**: Vectorized Danish and UK EQ-5D-5L value-set functions retained from the archived QALY notebook

**Shared Report and Clinical-Analysis Helpers**:
- **[assoc_tests.R](scripts/R/functions/assoc_tests.R)**: Shared categorical/continuous association tests and formatted results used by DAG reports and vignettes
- **[cox_extract.R](scripts/R/functions/cox_extract.R)**: Shared Cox-model fitting and tidy coefficient extraction helpers
- **[dag_helpers.R](scripts/R/functions/dag_helpers.R)**: Shared DAG construction, styling, validation, and rendering helpers
- **[report_format.R](scripts/R/functions/report_format.R)**: Shared strategy/biomarker labels and economic-result number formatting

**Sensitivity Analysis Functions**:
- **[psa_functions.R](scripts/R/functions/psa_functions.R)**: PSA-related utilities
- **[evppi_functions.R](scripts/R/functions/evppi_functions.R)**: EVPPI calculation functions, population-level EVPPI scaling (`calculate_population_evppi()`), and interaction coefficient EVPPI (`extract_interaction_coefficients()`, `get_interaction_evppi_params()`, `add_interaction_param_groups()`)
- **[scenario_analysis.R](scripts/R/functions/scenario_analysis.R)**: Scenario analysis framework

**Survival Modeling Functions**:
- **[para_model_fit.R](scripts/R/functions/para_model_fit.R)**: Parametric model fitting helper
- **[survival_plots.R](scripts/R/functions/survival_plots.R)**: Survival curve visualization

**Visualization Functions**:
- **[create_tornado_plot.R](scripts/R/functions/create_tornado_plot.R)**: DSA tornado diagram generation

**Snapshot/Impact Assessment Functions**:
- **[snapshot_utils.R](scripts/R/functions/snapshot_utils.R)**: Snapshot management for bug fix impact assessment

**Archived Functions** (in `scripts/R/archive/`):
- `bootstrap_survival_model.R`: Alternative resampling approach (not used in main analysis)
- `ref_values_emm.R`: Reference value calculations (superseded)

### Treatment Schedules

Treatment administration is defined by binary vectors aligned to `time_points <- seq(0, time_horizon, by = 1)`. Because R vectors are one-indexed while the model grid starts at week 0, vector position `i` represents modeled week `i - 1`; schedule subscripts in the code are positions, not week numbers.

**Nivolumab** (experimental): Modeled weeks 4, 6, 12, 14, 28, 30, 36, 38 (R positions 5, 7, 13, 15, 29, 31, 37, 39)
**FLOX experimental**: Modeled weeks 0, 2, 8, 10, 24, 26, 32, 34 (R positions 1, 3, 9, 11, 25, 27, 33, 35)
**FLOX control**: Modeled weeks 0, 2, 4, 6, 8, 10, 12, 14, 24, 26, 28, 30, 32, 34, 36, 38 (the union of the two experimental-arm position sets)

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
  - Scripts 04-11: Main analysis pipeline
  - Script 08b: Enriched population CEA for economic biomarker-positive populations (optional)
  - Script 12: Extended scenario EVPPI analysis (optional)
  - Script 13: Standalone snapshot utility (optional; run via `Rscript`, not sourced)
- **Functions directory** (`scripts/R/functions/`): Reusable components that are sourced by analysis scripts
- **Archive directory** (`scripts/R/archive/`): Deprecated/unused code preserved for reference, including the former baseline-characteristics and QALY exploratory notebooks
- **Tests directory** (`scripts/R/tests/`): Validation and diagnostic scripts
- **Quarto reports** (`scripts/QMD/report/`): Publication-ready PDF reports with embedded R code
- **Technical docs** (`scripts/QMD/technical_docs/`): Bug fix impact reports and technical documentation

### Retired Files

The B1 cleanup deleted obsolete or superseded files rather than leaving dead entry points in the live tree:

- `scripts/R/analysis/14b_scenario_preview.R` (the full scenario cache from `12_scenario_EVPPIs.R` is the only supported scenario output)
- `scripts/R/functions/para_model_fit_table.R`
- `scripts/R/tests/diagnose_prediction_failures.R`
- `scripts/R/tests/test_psa_error_rate.R`
- `scripts/R/tests/test_sampling_convergence_rate.R`
- `scripts/R/tests/test_sex_variable_fix.R`
- `scripts/R/tests/tlr.R`

`scripts/R/tests/fit_independent_biomarker_models.R` was moved to `scripts/R/archive/fit_independent_biomarker_models.R`; it was archived, not deleted.

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
- **Sources**: 02, 03, 04, 05
- **Shows**: Survival model fits, AIC/BIC comparisons, goodness-of-fit diagnostics
- **Models displayed**: Ordering-constrained, minimum-combined-AIC pair (currently gamma for OS and gamma for PFS) AND Weibull PH models for reference (issue #68)
- **Tables include**: Model parameter exponents for clinical interpretation
- **Note**: Economic survival fits use the single joint CRP + TMB/BRAF model

**[clinical_effectiveness.qmd](scripts/QMD/report/clinical_effectiveness.qmd)** - Clinical Effectiveness Analysis
- **Sources**: 02, 03 (with report-specific clinical models fitted in the document)
- **Shows**: Baseline characteristics, clinical survival curves, life-years gained; may include TLR clinical analyses independent of the economic model

**[input_parameters.qmd](scripts/QMD/report/input_parameters.qmd)** - Input Parameters Summary
- **Sources**: 02, 03, 04, 05
- **Shows**: All model input parameters (costs, utilities, prevalence rates, treatment schedules)

**[CEA.qmd](scripts/QMD/report/CEA.qmd)** - Cost-Effectiveness Analysis Report
- **Sources**: 02, 03, 04, 05; runs the base case through shared helpers
- **Requires**: PSA cache from script 10 for probabilistic outputs
- **Shows**: Incremental cost-effectiveness ratios (ICERs), cost-effectiveness plane, decision tables

**[OWSA.qmd](scripts/QMD/report/OWSA.qmd)** - One-Way Sensitivity Analysis (Deterministic)
- **Sources**: 02, 03, 04, 05, then 09
- **Requires**: DSA results from script 09
- **Shows**: Tornado diagrams, one-way sensitivity plots for all varied parameters

**[EVPPIs.qmd](scripts/QMD/report/EVPPIs.qmd)** - Value of Information Analysis
- **Sources**: 02, 03
- **Requires**: PSA results (script 10) and EVPPI results (script 11)
- **Shows**: Expected value of perfect information (EVPI), expected value of perfect partial information (EVPPI) for parameter groups

**[scenario_effect.qmd](scripts/QMD/report/scenario_effect.qmd)** - Scenario Analysis
- **Sources**: 02, 03; loads scenario results generated by script 12
- **Shows**: Alternative scenario results (e.g., different time horizons, discount rates)

**[biosimilar_scenario.qmd](scripts/QMD/report/biosimilar_scenario.qmd)** - Biosimilar Nivolumab Pricing Scenario
- **Sources**: 02, 03, 04, 05
- **Shows**: ICER comparison for base case vs biosimilar pricing (EUR 13,923 vs EUR 4,641/dose) for the single joint economic model

**[enriched_population.qmd](scripts/QMD/report/enriched_population.qmd)** - Enriched Population Analysis
- **Sources**: 02, 03, 08b
- **Shows**: Enriched (biomarker-positive) ICERs vs base case for CRP and TMB/BRAF
- **Note**: Not cached; re-runs on each render. Requires sampling cache.

**[biomarker_decomposition.qmd](scripts/QMD/report/biomarker_decomposition.qmd)** - Biomarker Effect Decomposition
- **Sources**: 02, 03, 04, 05
- **Shows**: Decomposition of biomarker effects on cost-effectiveness outcomes

**[biomarker_distributions.qmd](scripts/QMD/report/biomarker_distributions.qmd)** - Biomarker Distributions
- **Sources**: 02, 03
- **Shows**: Biomarker prevalence and distribution analyses, including clinical-only TLR summaries

**[survival_model_specification.qmd](scripts/QMD/report/survival_model_specification.qmd)** - Parametric Survival Model Specification
- **Sources**: 02, 03, 04 (uses the `models`, `os_candidates`/`pfs_candidates`, and `models$ordered_selection` objects created by script 04; no refitting)
- **Shows**: (1) regression coefficients of the single joint OS and PFS models for all nine candidate distributions (selected gamma/gamma first, then the previously used Weibull and log-normal), with exp(coefficient) interpreted per family (time ratio / hazard ratio / gamma rate ratio); (2) the current gamma/gamma pair versus the previous unconstrained Weibull-OS/log-normal-PFS pair on the population-averaged subgroup curves, with OS >= PFS violation counts and a PFS-minus-OS gap plot; (3) a gallery of all candidate distributions per subgroup (OS | PFS facets), unlabelled first and then labelled, over the Kaplan-Meier curves
- **Render**: dual-format (`pdf` + `gfm`). Render `--to pdf` first and `--to gfm` second: the PDF pass deletes the `_files/` figure directory, so a combined render leaves the `.md` with dangling image links

**Technical Documentation** (in `scripts/QMD/technical_docs/`):
- **[age_effect_analysis.qmd](scripts/QMD/technical_docs/age_effect_analysis.qmd)**: Age effect on survival outcomes
- **[all_parametric_survival_models.qmd](scripts/QMD/technical_docs/all_parametric_survival_models.qmd)**: Full survival model diagnostics
- **[bug_fix_impact.qmd](scripts/QMD/technical_docs/bug_fix_impact.qmd)**: Bug fix impact documentation for the single economic model

**Figure Vignettes** (in `scripts/QMD/vignettes/`):
- **figure1-4.qmd**: Publication-ready figures
- **suppl_figure_pfs_plots.qmd**: Supplementary PFS figures
- **figure_pfs_os_curves.qmd**: OS and PFS in the same panel, one PNG per strategy (`figs/figure_pfs_os_curves_{soc,crp,tmb_braf}.png`). Overlays Kaplan-Meier step curves (complete-case data) on the base-case parametric curves from the `predictions` object; biomarker panels show positive (experimental arm) and negative (control arm) subgroups. Sources 02/03/04 only; unnumbered for now.

**Clinical effectiveness vignettes** (in `scripts/QMD/vignettes/`): Files for the clinical effectiveness paper use the prefix `clin_effect_` followed by the paper figure/table label. Each vignette is self-contained (HTML, `embed-resources: true`) and saves its output to `figs/` (figures, via `ggsave`) or `tables/` (CSV, via `write.csv`). Naming examples:
- `clin_effect_figure1.qmd` → `figs/clin_effect_figure1.png` (six-panel Kaplan-Meier grid: OS left / PFS right columns; CRP, TMB/BRAF, and TLR rows, with the TLR row on the week-9 landmark cohort)
- `clin_effect_figure_s1.qmd` → `figs/clin_effect_figure_s1.png` (full DAG)
- `clin_effect_table_s2.qmd` → `tables/clin_effect_table_s2.csv` (DAG association consolidated summary)
- `clin_effect_figure_simplified_dag.qmd` → `figs/clin_effect_figure_simplified_dag.png` (simplified DAG; unpublished, off the numbered figure sequence)
- `clin_effect_figure_sensitivity_dag.qmd` → `figs/clin_effect_figure_sensitivity_dag.png` (sensitivity DAG; unpublished, off the numbered figure sequence)

**Poster vignette** (in `scripts/QMD/vignettes/`): `poster_biomarker_correlation.qmd` supports the SMDM poster "Modeling Multiple Biomarkers in Precision Oncology: How Ignoring Biomarker Correlations Can Reverse Clinical Conclusions". Self-contained (HTML, `embed-resources: true`). **Standalone and intentionally off-pipeline**: it re-introduces TLR as a third biomarker-guided strategy purely to illustrate the methodological point, even though the main economic model deliberately excludes TLR as a post-randomization mediator. It sources 02/03/06/07 + `model_fun`/`calculate_outcomes`/`prediction_functions`, then overrides `get_strategies()`/`get_biomarkers()` locally to add `tlr`; it does not modify the main pipeline. Outputs:
- `figs/poster_forest.png` (treatment x biomarker interaction, single-biomarker vs joint correlation-adjusted model, parametric AFT time ratios, PFS + OS)
- `figs/poster_ceac.png` (two-panel cost-effectiveness acceptability curves: single-biomarker vs joint modelling, at biosimilar nivolumab pricing; built via a self-contained parametric-bootstrap PSA reusing the cached economic-parameter draws)
- `tables/poster_cea.csv` (three-strategy CEA), `tables/poster_forest_data.csv` (forest interaction estimates), and `tables/poster_ceac_data.csv` (CEAC probabilities)

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

1. **"Object not found" errors**: Run the required analysis scripts first (especially 02-11, in sequence)
2. **Sampling cache missing**: Run [06_sampling.R](scripts/R/analysis/06_sampling.R) to generate sampling models
3. **Rendering hangs**: Some reports (especially CEA, EVPPI) may take minutes to render due to re-sourcing analysis scripts
4. **Stale caches after model changes**: Economic survival model or strategy-set changes require regeneration of sampling, PSA, EVPPI, scenario, and snapshot outputs.
5. **Results changed after methodological updates**: If cost-effectiveness results differ from earlier versions, check if methodological fixes were applied. Issues #69 and #70 (Nov 2024) changed survival prediction methodology from reference patient to population averaging/individual sampling. This **should** change results - it's a methodological improvement. Regenerate both sampling cache and PSA cache after these fixes. See issue #64 for impact documentation approach.

## Cache Management

The analysis uses several cache systems to speed up computation:

### 1. Sampling Cache (Survival Models)

**Location**: `data/tidy/sampling_models_n*.rds`
**Purpose**: Cached correlated PFS/OS survival model fits
**Generation**: Script [06_sampling.R](scripts/R/analysis/06_sampling.R) (~first run takes time)
**Size**: Hundreds of MB
**When to regenerate**: Delete cache file when:
- Survival model formulas change
- Economic strategy/biomarker set changes
- `n_samples` changes
- Clinical data is updated

**To regenerate**:
```r
# Delete cache file (example for full single-model cache)
cache_file <- here("data", "tidy",
                   paste0("sampling_models_n", n_samples, "_full.rds"))
file.remove(cache_file)
# Re-run 06_sampling.R
source("scripts/R/analysis/06_sampling.R")
```

### 2. PSA Cache (Analysis Results)

**Location**: `data/tidy/psa_obj_{ipd|correct}.rds` and `psa_params_{ipd|correct}.rds`
**Purpose**: Cached PSA simulation results (5000 runs)
**Generation**: Script [10_PSA.R](scripts/R/analysis/10_PSA.R) (~20-60 minutes first run)
**Size**: ~660 KB total
**When to regenerate**: Delete cache files when:
- Economic model logic changes ([model_fun.R](scripts/R/functions/model_fun.R), [calculate_outcomes.R](scripts/R/functions/calculate_outcomes.R), or [model_configs.R](scripts/R/functions/model_configs.R))
- **Prediction methodology changes** (e.g., issues #69, #70 fixes to survival curve generation)
- Base parameters change (costs, time horizon, discount rates)
- **UTILITY_SOURCE changes** (cache filenames encode utility source to prevent mixing)
- **Sampling cache is regenerated** (PSA depends on specific resampled models - always regenerate PSA after regenerating sampling cache)
- Parameter distributions change (distributional assumptions, means, SDs, correlations)

**To regenerate**:
```r
# Delete PSA cache files (example for IPD utilities)
file.remove(here("data", "tidy", "psa_obj_ipd.rds"))
file.remove(here("data", "tidy", "psa_params_ipd.rds"))
# Re-run 10_PSA.R
source("scripts/R/analysis/10_PSA.R")
```

### 3. EVPPI Cache

**Location**: `data/tidy/evppi_results_{ipd|correct}.RData`
**Purpose**: Cached EVPPI results for parameter groups
**Generation**: Script [11_EVPPIs.R](scripts/R/analysis/11_EVPPIs.R)
**When to regenerate**: Delete cache file when PSA cache is regenerated or EVPPI parameter groupings change

### 4. Scenario EVPPI Cache

**Location**: `data/tidy/scenario_evppi_results_{ipd|correct}.rds`
**Purpose**: Cached scenario-based EVPPI analysis results
**Generation**: Script [12_scenario_EVPPIs.R](scripts/R/analysis/12_scenario_EVPPIs.R)
**When to regenerate**: Delete cache file when PSA cache is regenerated or scenario definitions change

### Cache Workflow

1. **First run**: [06_sampling.R](scripts/R/analysis/06_sampling.R) generates sampling cache
2. **PSA uses sampling cache**: [10_PSA.R](scripts/R/analysis/10_PSA.R) generates PSA cache
3. **EVPPI uses PSA cache**: [11_EVPPIs.R](scripts/R/analysis/11_EVPPIs.R) generates EVPPI cache
4. **Subsequent runs**: All load from cache (fast)
5. **Manual invalidation**: Delete specific cache file(s) to regenerate

## Snapshot System (Bug Fix Impact Assessment)

The repository includes a snapshot comparison system for assessing the impact of bug fixes on model results.

### Components

- **[13_save_snapshot.R](scripts/R/analysis/13_save_snapshot.R)**: Saves single-model snapshots for the joint economic survival model. **Interactive and standalone, NOT part of the cache pipeline** — it reads the issue number and `baseline|fixed` from stdin, so run it as `Rscript scripts/R/analysis/13_save_snapshot.R <issue#> <baseline|fixed>` (sourcing it non-interactively just errors on the empty prompt). Optional; not required to regenerate report caches.
- **[snapshot_utils.R](scripts/R/functions/snapshot_utils.R)**: Utility functions for snapshot management
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

- **[test_survival_ordering.R](scripts/R/tests/test_survival_ordering.R)**: Verifies deterministic OS >= PFS ordering and PSA-only handling of resampled crossings
- **[test_biomarker_test_cost_mapping.R](scripts/R/tests/test_biomarker_test_cost_mapping.R)**: Verifies data-driven diagnostic-test cost assignment
- **[test_canonical_prevalence.R](scripts/R/tests/test_canonical_prevalence.R)**: Verifies weighted curves use canonical full-cohort biomarker prevalence
- **[test_psa_fallback_reporting.R](scripts/R/tests/test_psa_fallback_reporting.R)**: Verifies PSA fallback-rate reporting and its failure threshold
- **[test_sim_idx_validation.R](scripts/R/tests/test_sim_idx_validation.R)**: Verifies PSA resampling-index validation
- **[test_psa_basecase_alignment.R](scripts/R/tests/test_psa_basecase_alignment.R)**: Verifies structurally that the PSA control curve is the joint model with `Rx = control`, then checks that PSA strategy means (costs and QALYs) sit within 5 Monte Carlo standard errors of the base case and prints the incremental comparison (issue #151); skips the numerical check when the PSA cache is absent. The numerical criterion is a known failure at n_sim = 5000 (see the test protocol table)
- **[compare_snapshots.R](scripts/R/tests/compare_snapshots.R)**: Compares before/after snapshots for bug fix impact assessment
- **[para_models.Rmd](scripts/R/tests/para_models.Rmd)**: Parametric model fit validation and diagnostics
- **[snapshot.R](scripts/R/tests/snapshot.R)**: Helper script for running snapshot saves

### Test protocol and results

Focused executable regression tests live in [`scripts/R/tests/`](scripts/R/tests/). The durable record of the broader 12 July 2026 black-box and extreme-value run lives in [`report_Opus/findings/findings_blackbox.json`](report_Opus/findings/findings_blackbox.json); that run used an external scratch harness, so the JSON record, rather than the harness itself, is retained in this repository. Unless a row specifies otherwise, numeric comparisons use a tolerance of `1e-10`; any non-finite value, unexpected warning/error, or unmet criterion is a failure.

| Invariant / boundary test | Predefined pass criterion | Location | Recorded result |
|---|---|---|---|
| OS/PFS ordering | OS >= PFS at all 521 modeled time points for control and all four economic biomarker subgroups; an injected deterministic crossing must error, while a PSA crossing must be reported and capped | [`test_survival_ordering.R`](scripts/R/tests/test_survival_ordering.R) | Pass after ordering-constrained distribution selection |
| Cohort conservation | For every strategy and cycle, PF + P + D = 1 and each occupancy is within [0, 1] | Black-box record (`BB-TR0/1`) in [`findings_blackbox.json`](report_Opus/findings/findings_blackbox.json) | Pass for control, CRP, and TMB/BRAF |
| Utilities = 1 | With both state utilities set to 1, discounted QALYs equal discounted life-years | Black-box record (`BB-U1`) | Pass |
| Utilities = 0 | With both state utilities set to 0, total QALYs equal 0 for every strategy | Black-box record (`BB-U0`) | Pass |
| Costs = 0 | With every drug, test, visit, follow-up, and end-of-life unit cost set to 0, total cost equals 0 for every strategy | Black-box record (`BB-C0`) | Pass |
| No mortality | With OS and PFS fixed at 1, death occupancy remains 0 and undiscounted life-years equal the stated 10-year horizon | Black-box record (`BB-M0`) | Fail (minor): 521 weekly grid points produce 10.019 rather than 10.000 life-years; equal across strategies |
| Near-certain mortality | With survival forced near 0 after baseline, more than 99% of the cohort is dead by cycle 3 | Black-box record (`BB-M1`) | Pass |
| Determinism | Two deterministic runs with identical inputs produce bitwise-identical costs and QALYs | Black-box record (`BB-RE`) | Pass |
| PSA centred on base case | For every strategy, the PSA mean cost and mean QALYs are within 5 Monte Carlo standard errors of the deterministic base-case values; the control PSA curve must be the joint model with `Rx = control` (structural check, no cache needed) | [`test_psa_basecase_alignment.R`](scripts/R/tests/test_psa_basecase_alignment.R) | Structural check: Pass. Numerical check: **Fail (known)** after the issue #151 fix at n_sim = 5000: every strategy's mean QALYs sit +0.025 to +0.033 above the base case (6.3 to 9.6 SE, all the same sign), the bootstrap-mean bias of the extrapolated gamma curves that is shared by all arms because they come from one joint fit; costs are within 4.7 SE. The increments versus control, which the shared bias cancels out of, are within 1.2 SE (QALYs) and 3.4 SE (cost). Before the fix the control QALY mean was 21 SE *below* the base case and the CRP incremental QALYs were 0.122 versus 0.020. |

Run focused tests from the repository root with `"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" scripts/R/tests/<test-file>.R`. A test passes only if it exits with status 0 and all documented assertions succeed. Update the recorded result whenever model logic or the corresponding acceptance criterion changes; do not overwrite a known failure with a looser criterion.

## Clinical Context

Based on the METIMMOX clinical trial (NCT03388190) comparing alternating chemotherapy + immunotherapy versus chemotherapy alone in MSS/pMMR metastatic colorectal cancer patients. The analysis focuses on identifying biomarker-based subgroups that may benefit from adding immunotherapy to standard chemotherapy.
