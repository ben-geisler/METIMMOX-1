## Project Overview

METIMMOX-1 is a cost-effectiveness analysis comparing biomarker-guided immunotherapy strategies for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer. The analysis uses a partitioned survival model implemented in R to evaluate two pre-immunotherapy biomarker strategies (week-4 CRP and baseline TMB/BRAF) against standard of care. TLR is retained only for clinical effectiveness, DAG, and biomarker distribution analyses because it is a post-randomization mediator measured on treatment rather than a treatment-selection biomarker.

**Target Audience**: Health economists and researchers developing decision-analytic models in R.

## Repository Layout

Research-compendium layout (September 2026; the earlier `scripts/R/...` and `scripts/QMD/...` tree is gone):

```
R/                  model and helper functions (not an R package; sourced with here::here("R", ...))
analysis/           numbered pipeline scripts 01-13 (run in order; 13 is a standalone Rscript)
tests/              executable regression tests (run from the repo root with Rscript tests/<file>.R)
archive/            superseded code kept for reference
reports/            Quarto reports: each .qmd renders to .pdf + .md (+ <name>_files/ images)
reports/technical/  technical documentation (bug-fix impact, survival model specification, ...)
outputs/vignettes/  figure/table/poster generators (HTML); they write to outputs/figs and outputs/tables
outputs/figs/, outputs/tables/   publication figures and tables (tracked)
data/               confidential trial data and caches (ignored) except data/output/snapshots/
docs/manuscript/, docs/references/   manuscript sources and bibliography (currently placeholders)
validation/         external validation reports (validateHE_Opus_2026-07-12: report.qmd + findings/*.json)
publish/            publish_reports.R renders every report and publishes the PDFs
```

Git policy for renders: the `.md` files and their `_files/` assets are tracked so the reports are readable on GitHub and by LLMs; `.pdf`, `.html` and LaTeX intermediates are ignored and published by `Rscript publish/publish_reports.R --push` (gh-pages branch) or `--release` (GitHub release). No CI renders anything because rendering needs the confidential data. The validateHE skill writes its report to a sibling directory of the repo; copy a finished run into `validation/validateHE_<model>_<date>/` to keep it.

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

R analysis scripts depend on functions defined in `R/`. When running model code, you must source dependencies in the correct order:

#### Minimal Test Setup
```r
# For testing model_fun() and calculate_outcomes()
source("analysis/02_setup_and_global_variables.R")  # Global vars: time_horizon, cl, dr, etc.
source("analysis/03_biomarker_strategies.R")         # Biomarker definitions
source("analysis/04_parametric_survival_analysis.R") # Fit survival models
source("analysis/05_basecase_input_parameters.R")    # Parameter list: l_params_base

# Source required functions
source("R/model_fun.R")
source("R/calculate_outcomes.R")

# Now you can run the model
result <- model_fun(l_params_base, determpsa = "det", return_traces = FALSE)
```

#### Common Function Dependencies
- **model_fun.R** requires:
  - `calculate_outcomes.R` (cost/QALY calculations)
  - `prediction_functions.R` (population-averaged survival predictions; the draw-specific PSA predictions come from `generate_psa_population_averaged_predictions()` in script 06)

- **PSA/EVPPI scripts** (12, 13) require:
  - `model_fun.R`
  - `calculate_outcomes.R`
  - `psa_functions.R` and/or `evppi_functions.R`

#### From Command Line
```bash
# Test a fix by sourcing all dependencies
"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" -e "
  setwd(here::here());  # run from the repository root
  source('analysis/02_setup_and_global_variables.R');
  source('analysis/03_biomarker_strategies.R');
  source('analysis/04_parametric_survival_analysis.R');
  source('analysis/05_basecase_input_parameters.R');
  source('R/model_fun.R');
  source('R/calculate_outcomes.R');
  cat('Testing model_fun...\n');
  result <- model_fun(l_params_base, determpsa = 'det');
  cat('Success! Control cost:', result[['Cost']][1], '\n');
"
```

**Key principle**: Always check which analysis scripts source which function files (use `grep "source.*functions" analysis/*.R`) to understand dependencies.

### Running the Analysis

The numbered analysis scripts must be executed sequentially:

```r
# Core setup (run these first)
source("analysis/01_data_prep.R")
source("analysis/02_setup_and_global_variables.R")
source("analysis/03_biomarker_strategies.R")

# Survival analysis
source("analysis/04_parametric_survival_analysis.R")
source("analysis/05_basecase_input_parameters.R")

# Survival coefficient draws (paired covariance bootstrap; several minutes)
source("analysis/06_sampling.R")

# Model execution
source("analysis/07_traces.R")
source("analysis/08_basecase_analysis.R")
source("analysis/08b_enriched_population_analysis.R")  # Optional: enriched population CEA for economic biomarkers

# Sensitivity analyses
source("analysis/09_DSA.R")  # Deterministic sensitivity analysis
source("analysis/10_PSA.R")  # Probabilistic sensitivity analysis
source("analysis/11_EVPPIs.R")  # Expected value of perfect partial information

# Extended analyses (optional)
source("analysis/12_scenario_EVPPIs.R")  # Scenario-based EVPPI analysis
# 13_save_snapshot.R: optional, interactive, standalone -- run via Rscript with <issue#> <baseline|fixed>, not sourced
```

### Package Installation

```r
# Install required packages using pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, readxl, dplyr, tableone, ggplot2, flexsurv,
               survival, survminer, gems, mstate, tidyverse, xtable,
               darthtools, dampack, mvtnorm, Matrix, here, voi)

# Additional packages for Quarto reports
pacman::p_load(knitr, kableExtra, flextable, officer, scales, gridExtra, reshape2)
```

### Rendering Quarto Reports

Every report and technical doc declares `format:` with both `pdf:` and `gfm:` and renders to a `.pdf` (people) and a `.md` (LLMs, GitHub). Render PDF first and GFM second: the PDF pass deletes the `<name>_files/` image directory, so the reverse order leaves the `.md` with dangling links.

```bash
quarto render reports/CEA.qmd --to pdf && quarto render reports/CEA.qmd --to gfm
Rscript publish/publish_reports.R            # all reports, both formats, then _site/ for publishing
Rscript publish/publish_reports.R --only CEA,OWSA
```

- Tables: `knitr::kable()` piped into the `tbl_*` wrappers of [report_tables.R](R/report_tables.R) (`tbl_style`, `tbl_column_spec`, `tbl_row_spec`, `tbl_header_above`, `tbl_pack_rows`, `tbl_footnote`, `tbl_landscape`; `format = report_table_format()` where a format is needed). They apply kableExtra under LaTeX only; `setup_report()` also sets `kableExtra.auto_format = FALSE` and pins `knitr.table.format` to `pipe` outside LaTeX, because attaching kableExtra otherwise switches every `knitr::kable()` to HTML in the GFM render. Calling kableExtra directly makes the GFM pass fail with `Functions that produce HTML output found in document targeting commonmark output`. Under GFM, group headers and spanners are dropped and footnotes become a paragraph below the table.
- Verify content from the `.md`; `pdftotext` (poppler, `/mingw64/bin`) is the fallback for PDF-only checks (`pdftoppm` is not installed, so PDFs cannot be read visually).

**Prerequisites for rendering**:
- Run the analysis scripts required by the report (up to 11 for the full core analysis; script 12 additionally generates scenario outputs)
- Sampling cache must exist (from running [06_sampling.R](analysis/06_sampling.R))
- Results objects (e.g., `cea_results`, `owsa_results`, `psa_results`, `evppi_results`) must be in the R environment or saved as `.rds` files

### Complete Analysis & Reporting Workflow

To generate all analysis results and reports from scratch:

```r
# 1. Run all analysis scripts in order
source("analysis/01_data_prep.R")
source("analysis/02_setup_and_global_variables.R")
source("analysis/03_biomarker_strategies.R")
source("analysis/04_parametric_survival_analysis.R")
source("analysis/05_basecase_input_parameters.R")
source("analysis/06_sampling.R")  # Takes time on first run
source("analysis/07_traces.R")
source("analysis/08_basecase_analysis.R")
source("analysis/08b_enriched_population_analysis.R")  # Optional: enriched population CEA for economic biomarkers
source("analysis/09_DSA.R")
source("analysis/10_PSA.R")
source("analysis/11_EVPPIs.R")
source("analysis/12_scenario_EVPPIs.R")  # Optional: scenario analysis
# 13_save_snapshot.R: optional, interactive, standalone -- NOT sourced here.
# Run separately: Rscript analysis/13_save_snapshot.R <issue#> <baseline|fixed>

# 2. Render reports (from terminal/command line)
# quarto render reports/
```

**Running the pipeline from a clean / non-interactive session** (e.g. driving it with `Rscript`):
- Attach packages first — only `01_data_prep.R` calls `p_load`, so `pacman::p_load(...)` the full set above before sourcing, or `07_traces.R` fails with `could not find function "ggplot"`.
- `07_traces.R` assumes `model_fun()` is already loaded (scripts 08-13 source it themselves); source `model_fun.R`, `calculate_outcomes.R`, `prediction_functions.R` before it.
- `06_sampling.R` loads `sampling_models_n{n}_full.rds` if present (fast); regenerating it takes several minutes for 1,000 paired bootstrap fits. `10_PSA.R` (~20 min) and `12_scenario_EVPPIs.R` (~40 min) are the longest steps at n=5000.
- Scenario reports and publication vignettes require the full cache generated by `12_scenario_EVPPIs.R`.

Or render reports individually in the desired order:
```bash
quarto render reports/para_models.qmd
quarto render reports/input_parameters.qmd
quarto render reports/clinical_effectiveness.qmd
quarto render reports/CEA.qmd
quarto render reports/OWSA.qmd
quarto render reports/EVPPIs.qmd
quarto render reports/scenario_effect.qmd
quarto render reports/biosimilar_scenario.qmd
quarto render reports/enriched_population.qmd
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

Defined in [02_setup_and_global_variables.R](analysis/02_setup_and_global_variables.R):

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

There is no model-structure switch or multi-structure comparison layer. The economic model is the joint CRP + TMB/BRAF formula defined in [model_configs.R](R/model_configs.R). There is no separate control-arm model: in both the base case and the PSA the standard-of-care curves are the joint model predicted with treatment set to control (issue #151). The sampling cache stores one joint coefficient-draw component; a legacy `control` component in an older cache is ignored.

**Parametric distribution selection**: OS and PFS distributions are selected jointly from the nine candidate families in [04_parametric_survival_analysis.R](analysis/04_parametric_survival_analysis.R). The selected pair is the minimum-combined-AIC pair that preserves OS >= PFS for control and every economic biomarker subgroup at every modeled weekly time point. With the current data and 10-year horizon, the ordering-constrained selection is **gamma for OS and gamma for PFS**.

**When changed**: Regenerate sampling cache when survival formulas, the economic strategy/biomarker set, `n_samples`, or clinical data change. Regenerate PSA and EVPPI caches after regenerating sampling cache or changing economic parameters, distributions, prediction methodology, or `UTILITY_SOURCE`. Regenerate the EVPPI and scenario-EVPPI caches (scripts 11 and 12) after changing the EVPPI estimator or parameter groups.

### Value of Information (EVPPI) Estimator

EVPPI is estimated by nonparametric regression (Strong, Oakley & Brennan 2014) through `voi::evppi()` in [evppi_functions.R](R/evppi_functions.R) (issue #152): the incremental NMB of each strategy versus control is regressed on the parameter(s) with a GAM; groups of up to four parameters use voi's tensor-product cubic regression spline, larger groups use additive cubic regression splines because NMB is linear and additive in unit costs. Since issue #154 the sampled parameter set no longer contains the unit prices, so the cost side of the EVPPI table is the single `other_costs` group (resource-use costs) and there is no `drug_costs`, `test_costs` or `all_costs` row; the `utilities` group is `{u_np, u_decrement}`, with the derived `u_p` excluded from groups but still reported as a single-parameter row. Every estimate carries a Monte Carlo standard error (`evppi_se`, the SD of the EVPPI over 1,000 draws of the regression coefficients, seeded with `analysis_seed`). Failed estimates are `NA`, never zero, with the reason in the `error` column. `run_evppi_analysis()` asserts that every `[GROUP]` row is at least its largest member within two combined SEs (floor 1% of EVPI) and attaches the check as attribute `group_consistency`; `add_interaction_param_groups()` adds per-biomarker interaction groups plus the joint `interaction_all` group that Figure 5 and Table S7 report as "Biomarker-treatment interaction". Figure 5 and Table S7 use `[GROUP]` rows only; EVPPI is not additive and rows must never be summed. After issue #159 (5,000 joint OS/PFS normal draws), EVPI is **EUR 0** per patient at WTP EUR 51,000: standard of care has the largest NMB in every draw. The existing low-EVPI short-circuit stores no EVPPI rows; all 301 report/cache contracts pass, with one inapplicable row-check block skipped. Before #159, the regenerated issue #166 PSA cache (independent OS/PFS normal draws and shared joint-population weights) gave EVPI **EUR 0.22** per patient at WTP EUR 51,000. For that #166 cache, all EVPPI point estimates were zero; the joint interaction group had SE 0.25. That cache contained the full parameter/group table, and all 296 report/cache contract checks passed without skips. Before issue #166, EVPI was EUR 0 and the EVPPI cache was empty. Under the pre-#156 bootstrap cache total EVPI was EUR 22.13 with every economic-parameter EVPPI 0 (SE 0) and only the interaction coefficients non-zero (`[GROUP] interaction_all` 1.18, SE 1.94); that value came from the bootstrap's extreme resamples, which the normal approximation does not produce. The previous kNN estimator (1,000 nearest of 5,000 draws at 500 grid points, capped at EVPI) reported neighbourhood re-weighting bias as value of information and was removed.

### Biomarker Strategies

Two pre-immunotherapy biomarkers are evaluated economically (defined in [03_biomarker_strategies.R](analysis/03_biomarker_strategies.R) and selected by [model_configs.R](R/model_configs.R)):

1. **CRP** (C-reactive protein): Binary variable, cut-off <5 mg/L, **week-4 value** (`CRP1cat`: cycle 3 day 1, trial visit 3), measured after the two FLOX cycles that both arms receive and before the first nivolumab dose
2. **TMB/BRAF**: Combined biomarker (TMB >=9 mut/MB OR BRAF mutation), from baseline NGS

**CRP timing (issue #150)**: CRP is not a baseline (pre-randomization) measurement, but it is available before the immunotherapy decision, because every patient receives two FLOX cycles before the first nivolumab dose. Describe it as "week-4 CRP, before the first nivolumab dose", never as "baseline" or "pre-treatment" CRP. The true baseline value (`CRP0cat`, cycle 1 day 1) is derived in `01_data_prep.R` but unused. Consequences: (i) week-4 CRP prevalence differs by arm (17/36 experimental vs 7/35 control, Fisher p = 0.023) whereas baseline CRP is balanced (p = 1); the DAG report reads this as chance arm imbalance on a week-4 measurement in a small trial, not as a randomisation failure, and the DAG has no `T -> CRP` edge because nivolumab has not been given at week 4 and both arms receive identical FLOX until then; (ii) the 3 patients without a week-4 CRP are early deaths (weeks 2.4, 15.7, 20.9) and are dropped from the complete-case data; (iii) the model clock starts at randomization (t = 0) while the CRP decision is made at week 4; this has no economic consequence because both arms receive identical FLOX, monitoring, and visits up to that point, and the week-4 blood test is already in the monitoring schedule. No model code or cache depends on this framing.

**TLR** (tumor lesion reduction) is intentionally excluded from all economic analyses because it is a post-randomization mediator measured on treatment, not a treatment-selection biomarker. `data$tlr` and `p_tlr` may still be created for clinical effectiveness, DAG, and biomarker distribution reports that analyze TLR directly from the trial data.

**TLR landmark and DAG validation (issue #155)**: TLR is read at the first on-treatment CT, whose date is the positional column `Date...122` of the trial export; `02_setup_and_global_variables.R` derives `CT1wk` (weeks from inclusion to that scan, NA exactly when TLR is NA) and `03_biomarker_strategies.R` keeps `CT1wk`, `ProgressionExit` and `TTPwk` in `data` for the clinical and DAG reports (no economic script reads them). Scans fell at weeks 6.7 to 12.1 (median 8.6), and 8 patients progressed at that scan (all TLR-negative, because progression and TLR are read from the same image). [tlr_landmark.R](R/tlr_landmark.R) is the single definition of the landmark cohorts (`build_tlr_landmark_cohorts()`, `build_all_tlr_landmark_cohorts()`, endpoint time strictly after the landmark, clock reset to the landmark) shared by `clinical_effectiveness.qmd`, `clin_effect_figure1.qmd` and the DAG association tests. The **primary landmark stays fixed at week 9** (user decision); it retains 3 of the 8 first-scan progressors (landmark times 0.3 to 1.0 weeks) and 19 patients whose TLR was read after week 9, and every report states both counts. The **per-patient scan-date landmark** and a **fixed week-12 landmark** are sensitivity analyses in `clinical_effectiveness.qmd`; both exclude all 8 first-scan progressors (PFS cohort 60 instead of 63 in the DAG report's 68-patient TLR subset, 57 instead of 60 in the clinical report's 65-patient subset) and leave the direction of the TLR x Rx PFS estimate unchanged. Never describe the week-9 cut as guaranteeing that every retained patient reached the scan. [dag_association_tests.R](R/dag_association_tests.R) derives the edge list and the 19 implied conditional independencies from the canonical `dag` object (`dag_edge_coverage()`, `dag_ci_coverage()`, `run_dag_edge_tests()`, `run_dag_ci_tests()`), stops when a DAG statement has no test or a test has no DAG statement, tests both `TxCRP` edges, and replaces two tautological tests: `PFS -> OS` uses a time-dependent progression indicator (`survival::tmerge` on `ProgressionExit`/`TTPwk`, Firth Cox with counting-process time; HR 3.75 for death after versus before progression), and `TLR -> PFS` is estimated on the week-9 landmark PFS cohort. Edges from `U` and into the interaction nodes are listed as untested; the two statements relating `TxCRP` and `TxTMB` given `T` are reported as holding by construction. `dag_associations.qmd` no longer defines its own graph.

Each economic biomarker strategy has:
- **Biomarker-positive subgroup**: Receives experimental treatment (alternating FLOX + nivolumab)
- **Biomarker-negative subgroup**: Receives standard treatment (FLOX only)
- **Analysis approach**: See "Survival Prediction Methodologies" section below for how population-level outcomes are calculated

### Common Economic Target Population (issue #166)

`economic_prediction_population()` in `R/prediction_population.R` defines one complete-case trial population (68 patients), including both economic biomarkers and both endpoints. Script 03 derives economic prevalences from this population (CRP 23/68, TMB/BRAF 30/68); clinical-only TLR still uses its own available-case cohort. `l_params_base$prediction_population` is passed explicitly to all PSA predictions and to the enriched analysis. `population_curves` stores the individual curve basis for deterministic reweighting; `set_population_predictions()` must be used whenever replacing a curve set, so its population and weights remain aligned.

Prevalence uncertainty represents uncertainty in the trial population's **joint** biomarker distribution, not transport. The PSA samples Dirichlet(n00, n01, n10, n11) joint cell masses (Bayesian-bootstrap cell weights), keeping the within-cell age/sex distribution fixed. `p_joint_00`, `p_joint_01`, and `p_joint_10` are the free EVPPI coordinates; `p_joint_11`, `p_crp`, and `p_tmb_braf` are derived. Population weights and survival coefficients are drawn independently. A one-way change of either marginal rakes the observed joint distribution to the requested marginal while retaining the other marginal and the joint odds ratio. **Every population change reweights control and all guided strategies together.** The identical-treatment regression covers ordered synthetic curves with heterogeneous prognosis, both marginal changes, joint draws, deterministic and PSA paths, and an endpoint-missing patient.

### Survival Prediction Methodologies

The model uses different approaches for base case vs PSA, both properly accounting for patient heterogeneity:

**Base Case** (Issue #69): Population averaging across all patients using `generate_population_averaged_predictions()` in [prediction_functions.R](R/prediction_functions.R:199-340). Predicts for all patients using their actual age, sex, and biomarker values, then averages.

**PSA** (issues #70, #156): second-order Monte Carlo over the survival-model coefficients, with the same population averaging as the base case inside every draw. `06_sampling.R` (`sample_survival_coefficients()`) draws `n_samples` combined OS/PFS coefficient vectors from a joint multivariate normal on flexsurvreg's optimisation scale (issue #159). Each endpoint retains its fitted mean and marginal covariance; paired patient bootstrap fits estimate cross-endpoint dependence. PSA draw `sim_idx` takes row `sim_idx` of both draw matrices, rebuilds the two flexsurvreg objects (`sampled_survival_models()`), and `model_fun()` predicts every patient of `data_complete` (control: all patients with `Rx = control`; biomarker-positive: the positives with `Rx = experimental`; biomarker-negative: the negatives with `Rx = control`) and averages the predicted curves (`generate_psa_population_averaged_predictions()` in script 06). Parameter uncertainty therefore varies between draws while patient heterogeneity is integrated out within each draw, so EVPI measures the value of reducing parameter uncertainty. No patient is sampled per iteration; the earlier "one patient per iteration" description and the nonparametric bootstrap it referred to no longer apply (see "Survival Parameter Sampling" below).

Both approaches are methodologically valid for their respective analytical purposes (deterministic point estimates vs. probabilistic uncertainty quantification).

### Survival Parameter Sampling (issues #156, #159)

Survival-model parameter uncertainty enters the PSA through **joint multivariate-normal coefficient draws** in [06_sampling.R](analysis/06_sampling.R), with the covariance helpers in [joint_survival_sampling.R](R/joint_survival_sampling.R). The two gamma models retain their original maximum-likelihood estimates and observed-information marginal covariance matrices on flexsurvreg's optimisation scale. Issue #159 replaces the independent OS/PFS draws of #156 using its first proposed fix:

- **Paired covariance bootstrap.** Both endpoints are refitted on each of 1,000 identical, unstratified patient resamples. A pair is excluded if either fit fails convergence, covariance or coefficient-set checks; no original estimate is substituted and no coefficient-tail or crossing rejection is applied. More than 10% failed pairs or inadequate covariance rank stops generation. Failure counts and reasons are stored.
- **Marginal-preserving cross-covariance.** The empirical paired bootstrap covariance is block-whitened and recoloured to the original fits' covariance matrices using symmetric matrix roots. This preserves positive semidefiniteness and exactly retains both marginal covariance blocks. Only the calibrated cross-block changes. Bootstrap coefficients are never used as PSA draws. The [technical specification](reports/technical/joint_survival_sampling.qmd) gives the formula and limitations.
- **One joint component.** The cache is `list(joint = <component>, biomarkers, fingerprint, fingerprint_inputs, creation_time)` plus canonical input/runtime diagnostics. The component has `method = "mvn_joint_v2"`, `draws$os`, `draws$pfs`, the original fits, coefficient summaries, the seed, and `joint_covariance` (joint covariance, raw bootstrap covariance, retained bootstrap coefficient matrices, attempted/successful/failed counts and failure reasons). The normal approximation produces no failed coefficient draws (`n_failed = 0`); this is distinct from `joint_covariance$n_failed`. The extreme interaction flag remains `log(10)` and is diagnostic only. Cache generation now takes several minutes because it estimates covariance by refitting.
- **Accessors unchanged.** `get_joint_sampling_models()` resolves the component; `sampled_survival_models(component, idx)` rebuilds the two flexsurvreg objects with `build_sampled_flexsurv_model()`. Legacy bootstrap fixtures still resolve through `$samples`. Every strategy within one PSA draw shares the same endpoint vectors. Population weights and economic draws remain independent of survival coefficients.
- **Validity is a fingerprint.** `calculation_identity()` includes `R/joint_survival_sampling.R` in both sampling and PSA dependencies. Sampling method, helper implementation, formulas, data, distributions, sample size and seed invalidate the chain. The active file remains `data/tidy/sampling_models_n{n_samples}_full.rds`.
- **Different outcome estimands.** Joint covariance does not force outcome means to equal the deterministic result: averaging nonlinear survival curves differs from evaluating them at fitted coefficients. With fixed marginal distributions, dependence alone cannot change expected raw-curve QALYs. It can change crossings, the clamp correction and the net-benefit distribution. No recentering or relaxed five-SE threshold is used. `test_psa_zero_uncertainty.R` verifies equality of every strategy's cost and QALYs when covariance is zero and economic parameters/population weights are fixed; `test_joint_survival_sampling.R` checks paired refits, marginal covariance preservation, dependence and reproducibility on synthetic data.

### Main Model Function

The core economic model is in [model_fun.R](R/model_fun.R:9-381):

```r
model_fun(params, time_horizon = 520, cl = 1/52,
          determpsa = "det", return_traces = FALSE, sim_idx = NULL)
```

**Parameters**:
- `determpsa`: "det" for deterministic, "psa" for probabilistic
- `sim_idx`: Coefficient-draw index, the row of the OS and PFS draw matrices in the sampling cache (required for PSA mode)
- `return_traces`: If TRUE, returns state occupancy over time

**PSA Mode Behavior** (lines 17-148):
- Uses second-order Monte Carlo: population-averaged predictions from the coefficient draw `sim_idx` of the joint model (see "Survival Prediction Methodologies" section)
- The OS and PFS curves of a draw come from the same joint component (row `sim_idx` of the OS and PFS draw matrices, drawn jointly since issue #159)
- If prediction from the sampled model fails, the draw is flagged `fallback_used` and the PSA loop re-pairs the economic-parameter draw with up to 10 other cached models; the model actually used is recorded per row as `model_idx` (issue #156)

### Parameter Structure

The model expects a comprehensive parameter list `l_params_base` containing:

```r
list(
  # Time & discounting
  cl, time_horizon, dr_costs, dr_effects,

  # Utilities: u_np beta in PSA; u_p derived as u_np - u_decrement (gamma)
  u_np, u_p, u_decrement,

  # Drug costs (FIXED in PSA, varied in DSA only)
  c_drug_nivo, c_drug_FLOX,

  # Test costs (FIXED in PSA, varied in DSA only)
  c_test_CT, c_test_blood, c_test_NGS,

  # Other costs (gamma distributions in PSA); c_other_pp = 0 in the base case
  c_other_visit, c_other_baseline, c_other_follow, c_other_pp, c_other_last,

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

**PFS endpoint derivation** (issue #149, in [pfs_endpoint.R](R/pfs_endpoint.R), called by [02_setup_and_global_variables.R](analysis/02_setup_and_global_variables.R)): the trial export records exit-for-progression only (`Progression exit`) and time to that exit (`Days until progression`), so patients who died without a recorded progression would be censored for PFS. `derive_pfs_endpoint()` recodes `Progression` as progression OR death and `PFSwk` as time to progression if progressed, time to death if died progression-free, and follow-up time otherwise, matching the trial definition (Ree et al. 2024) and the partitioned-survival-model requirement that the progression-free state means alive and progression-free. The raw values are kept as `ProgressionExit` and `TTPwk`. Any report that reads `data/tidy/METIMMOX.rds` directly (currently `biomarker_distributions.qmd`) must call `derive_pfs_endpoint()` after renaming `Progression exit`.

**Biomarker derivation** (in [03_biomarker_strategies.R](analysis/03_biomarker_strategies.R:6-9)):
```r
data$crp <- as.numeric(data$CRP1cat == 1)
data$tlr <- as.numeric(data$TLRcat == 1)
data$tmb_braf <- as.numeric((data$TMBcat == 1) | (data$Mutation == "BRAF"))
```

`data$tlr` is kept for clinical effectiveness, DAG, and biomarker distribution reports only. It is not included in economic strategies, economic survival formulas, PSA parameters, EVPPI groups, DSA, scenario analyses, enriched-population CEA, or snapshots.

### Adverse-Event Scope

Adverse-event/toxicity costs and disutilities are not modeled separately. This is a deliberate scope choice based on the intended tolerability of the alternating short-course FLOX-nivolumab regimen and the lack of sufficiently robust treatment-specific trial data on adverse-event incidence, resource use, and utility decrements for economic parameterization.

### Cost and Scope Decisions (issue #154)

**Unit prices are fixed in the PSA.** `c_drug_nivo`, `c_drug_FLOX`, `c_test_CT`, `c_test_blood` and the biomarker test costs (`c_test_CRP`, `c_test_NGS`) are published tariffs or assumed list prices, not quantities research could resolve. They are `psa = FALSE, dsa = TRUE` in `parameter_distribution_spec()`: excluded from PSA sampling and from every EVPPI group, still varied +/-20% in the one-way DSA. Under the pre-#152 kNN estimator, sampling them put `c_drug_FLOX` (44% of EVPI) and `c_drug_nivo` (34%) at the top of the value-of-information ranking, which is what issue #154 reported; the regression estimator already gave both an EVPPI of 0, but their variance still inflated total EVPI, which fell from 45.20 to 22.13 when they were fixed. Their deterministic influence is undiminished: `c_drug_nivo` remains the second-largest one-way driver for both guided strategies. Visit, baseline, follow-up and end-of-life costs stay probabilistic because they bundle genuine resource-use uncertainty. Consequences: there is no `drug_costs` or `test_costs` EVPPI group, and `all_costs` is emitted only when more than one cost group is still sampled (currently it is not, so `create_parameter_groups()` returns `other_costs` alone) — Figure 5 reports `[GROUP] other_costs` as "Resource-use costs".

**The nivolumab price is an assumption, not a tariff.** EUR 13,923 per administration is an assumed Norwegian hospital acquisition cost, roughly 40% below list, giving about EUR 111,000 over eight administrations. Norwegian hospital prices are set by confidential LIS tender and cannot be cited; international list prices are higher (a US list price is about USD 8,100 per 240 mg vial). It is the largest single cost driver and is bounded by the one-way DSA and the biosimilar scenario (EUR 4,641 per administration).

**Utilities are sampled jointly, never independently.** The PSA draws `u_np` (beta) and a non-negative decrement `u_decrement` (gamma, mean `u_np - u_p`), then derives `u_p = pmax(u_np - u_decrement, 0)` in `apply_derived_psa_parameters()`. Independent beta draws put 15.9% of draws (48.5% under `UTILITY_SOURCE = 0`) in the impossible region `u_p > u_np`. `u_p` is a `derived` spec row: it is a model input so it takes part in the DSA, but it is excluded from EVPPI groups because it is exactly collinear with its inputs; the `utilities` group is `{u_np, u_decrement}`. `l_params_base$u_decrement` is derived in script 05 and asserted positive. `CEA.qmd` prints the realised reversal fraction (zero by construction).

**No post-progression treatment cost.** `c_other_pp` (base case 0) is charged on the quarterly follow-up cycles against `p_p` in `calculate_outcomes()`. Second-line systemic therapy, post-progression imaging and post-progression visits are not modeled, so the progressed state accrues only `c_other_follow` and `c_other_last`. The strategies differ in time spent progressed, so the omission is differential. The parameter has a degenerate one-way range at zero and is instead varied by the `Post_progression_cost` structural scenario (EUR 5,000 per quarter, an illustrative upper bound rather than a costed Norwegian pathway).

**The second treatment sequence is given to all progression-free patients.** Schedule positions 25-39 (modeled weeks 24-38) repeat the first eight-cycle sequence in both arms, for everyone still progression-free; METIMMOX re-treated on progression during the break. A partitioned survival model has no on-treatment substate, so re-treatment cannot be triggered on progression without restructuring the model. The `Second_sequence` structural scenario removes the second sequence (drug and visit costs; monitoring and survival unchanged), giving a cost-side bound.

**Structural scenarios can be one-sided.** `dsa_structural_scenarios()` entries carry a `sides` field read through `structural_scenario_sides()`; `Post_progression_cost` runs `"max"` only and `Second_sequence` runs `"min"` only, because their other endpoint is the base case. `summarise_param_ranges()` treats an absent endpoint as the base case (NMB difference 0) rather than dropping the parameter, so one-sided scenarios still get a tornado bar spanning base case to scenario. `09_DSA.R` groups DSA rows through `parameter_group_lookup()` (the full spec) rather than the PSA groups, so fixed prices still appear as drug/test costs in the tornado plots.

### Survival Model Formulas

**Economic model** (shared by the control, CRP-guided, and TMB/BRAF-guided strategies):
```r
OS:  Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx
```

**Control arm**: not a separate model. Standard-of-care curves are the joint OS and PFS models predicted with `Rx` set to the control level for every patient in `data_complete`, then averaged. This holds in the base case (`generate_population_averaged_predictions()`) and in every PSA draw (`generate_psa_population_averaged_predictions()` with `biomarker_name = NULL`), so the deterministic and PSA paths agree when uncertainty is removed; nonlinear averaging can still shift uncertain PSA outcome means. An earlier age/sex-only control bootstrap fitted on the 36 control-arm patients was removed in issue #151.

See [model_configs.R](R/model_configs.R) for the canonical formula definitions.

### Key Functions

**Model Configuration**:
- **[model_configs.R](R/model_configs.R)**: Single source of truth for the economic strategies, biomarkers, and formulas. Auto-sourced by `02_setup_and_global_variables.R`. Key functions: `get_model_configs()`, `get_current_model_config()`, `get_strategies()`, `get_biomarkers()`, `get_strategy_formula()`, `get_model_formulas()`. (`get_control_formula()` was removed in issue #151; there is no separate control-arm formula.)
- **[pfs_endpoint.R](R/pfs_endpoint.R)**: Single source of truth for the composite PFS endpoint (progression or death). Auto-sourced by `02_setup_and_global_variables.R`. Key function: `derive_pfs_endpoint()`.
- **[cache_paths.R](R/cache_paths.R)**: Single source of truth for the utility-source label, cached-object file paths (sampling, PSA, EVPPI, scenario) and cache fingerprints (issue #156). Auto-sourced by `02_setup_and_global_variables.R`. Key functions: `resolve_util_label()`, `sampling_cache_path()`, `psa_obj_path()`, `psa_params_path()`, `evppi_path()`, `scenario_evppi_path()`, `cache_fingerprint()` (rlang hash of a canonicalised object; formulas deparsed, factors as character, list and column order ignored), `sampling_cache_fingerprint()` (formulas, fitting data, distributions, n_samples, seed, method), `psa_cache_fingerprint()` (sampling fingerprint, the whole `l_params_base`, the PSA distributions, strategies, n_sim, seed, horizon, cycle length), `cache_file_provenance()` (md5 and mtime of a file).
- **[report_setup.R](R/report_setup.R)**: One-call Quarto report setup (knitr options, package loading, shared ggplot theme, and sourcing of analysis scripts/function files), used to remove duplicated setup boilerplate across the economic reports. Key function: `setup_report(sources, funs, packages, set_theme)`.
- **[report_tables.R](R/report_tables.R)**: Format-aware table wrappers (`tbl_style()`, `tbl_column_spec()`, `tbl_row_spec()`, `tbl_header_above()`, `tbl_pack_rows()`, `tbl_footnote()`, `tbl_landscape()`, `report_table_format()`) that apply kableExtra under LaTeX only, so every report renders to PDF and GFM. Sourced by `report_setup.R`.

**Core Model Functions**:
- **[model_fun.R](R/model_fun.R)**: Main partitioned survival model with PSA support
- **[calculate_outcomes.R](R/calculate_outcomes.R)**: Calculates QALYs and costs from state occupancy traces
- **[prediction_functions.R](R/prediction_functions.R)**: Generate survival predictions from fitted models; sampling-cache accessors `get_joint_sampling_models()`, `sampled_survival_models()` and `build_sampled_flexsurv_model()` (issue #156)
- **[cea_helpers.R](R/cea_helpers.R)**: Single-model CEA execution and summary helpers wrapping dampack (`run_basecase()`, `load_psa_cache()`, `create_ceac_plot()`, `create_psa_summary_table()`, `calculate_pairwise_icers()`, `frontier_status()`, `frontier_label()`). Renamed from the legacy `multi_model_cea.R`. **Reports stop on missing or stale caches (issue #157)**: `validate_psa_cache()` mirrors the acceptance checks of `10_PSA.R` (strategy set, requested `n_sim`, failed-draw policy, `model_idx`, sampling fingerprint, full `psa_cache_fingerprint()` with the differing input named); `load_psa_cache()` and `load_psa_params_cache()` take `required = TRUE` plus the expected fingerprint, strategies, `n_sim`, seed and sampling fingerprint, and `CEA.qmd` and `EVPPIs.qmd` pass them (both source script 06 to compute the fingerprint) instead of rendering "--" placeholders; `EVPPIs.qmd` also stops when the EVPPI cache was built from a different PSA fingerprint. Callers without expectations (snapshot script, tests, `table_5.qmd`) behave as before. **Dominance conventions (issue #153)**: `CEA.qmd` reports dampack *frontier* ICERs (status ND/D/ED across all three strategies). The biosimilar, enriched-population and Table S6 outputs report *pairwise* ICERs of each guided strategy versus standard of care (`Status` values "Reference", "Pairwise ICER vs SoC", "Dominated by SoC", "Cost-saving vs SoC") and always print the dampack `Frontier_Status` alongside, with captions stating "pairwise versus standard of care". A strategy can be frontier-dominated (more costly and less effective than the other guided strategy) yet have a finite pairwise ICER; percentage changes between pairwise ICERs are only reported for frontier strategies with finite, positive ICERs in both analyses (`Pct_Change` in `08b_enriched_population_analysis.R`).
- **[parameter_distributions.R](R/parameter_distributions.R)**: Single source of truth for the PSA/DSA/EVPPI parameter set (`parameter_distribution_spec()`, `configure_parameter_distributions()`, `parameter_group_lookup()`, `apply_derived_psa_parameters()`), the one-way DSA bounds (`build_dsa_ranges()`: base +/- `DSA_mult`, capped at 0 for costs and prevalences and at 1 for utilities and prevalences) and the deterministic structural scenarios (`dsa_structural_scenarios()`, `structural_scenario_sides()`: discount rate 0% and 8% applied to costs and QALYs, time horizon 5 and 20 years, post-progression cost EUR 5,000 per quarter, second sequence omitted). The spec carries `psa`, `dsa` and `derived` flags, so PSA membership and DSA membership are no longer the same set: unit prices are `psa = FALSE, dsa = TRUE`, and `u_p` is `derived` (issue #154). `build_dsa_ranges()` takes the spec, not the distribution list. `09_DSA.R`, `figure2.qmd` and `table_1.qmd` all read these helpers, so the published ranges in Table 1 are exactly the ranges run (issue #153). The scenario code is also here (issue #157): `apply_structural_scenario(scenario_name, value, base_params, base_horizon, models, strategies_df, data_complete)` returns `list(params, time_horizon)` for `model_fun()`, using `build_horizon_params()` (re-predicts the base-case curves on the new weekly grid from `models$best_fit` and rebuilds the schedule vectors with the script-05 rules) and `drop_second_sequence()`; `structural_scenario_label()` names an endpoint. `09_DSA.R` (`run_structural_scenario()`) and `table_s8.qmd` both call it, so the DSA and Table S8 run identical scenarios; scenario rows enter `dsa_results` with `group = "structural"` and appear in the tornado plots and in the OWSA structural-scenario table.
- **[eq5d5l_utility.R](R/eq5d5l_utility.R)**: Vectorized Danish and UK EQ-5D-5L value-set functions retained from the archived QALY notebook

**Shared Report and Clinical-Analysis Helpers**:
- **[assoc_tests.R](R/assoc_tests.R)**: Shared categorical/continuous association tests and formatted results used by DAG reports and vignettes
- **[cox_extract.R](R/cox_extract.R)**: Shared Cox-model fitting and tidy coefficient extraction helpers
- **[dag_helpers.R](R/dag_helpers.R)**: Shared DAG construction, styling, validation, and rendering helpers
- **[dag_association_tests.R](R/dag_association_tests.R)**: DAG-derived edge and conditional-independence test lists with coverage assertions, the time-dependent `PFS -> OS` test and the landmark `TLR -> PFS` test, used by `dag_associations.qmd` and `clin_effect_table_s2.qmd` (issue #155)
- **[tlr_landmark.R](R/tlr_landmark.R)**: Single definition of the TLR landmark cohorts (fixed week 9 primary; per-patient scan date and fixed week 12 sensitivity) with first-scan-progressor diagnostics (issue #155)
- **[report_format.R](R/report_format.R)**: Shared strategy/biomarker labels and economic-result number formatting. `format_icer()` prints "Cost-saving" for a "Cost-saving vs SoC"/"Dominant" status and "--" for a negative ratio with no informative status (a negative ICER is never printed as a cost per QALY); `restricted_mean_survival(y, cycle_length)` is the trapezoidal restricted mean of a weekly survival curve in years, restricted to the curve's horizon (520 weeks) and used by `biomarker_decomposition.qmd` and `age_effect_analysis.qmd` (issue #157)

**Sensitivity Analysis Functions**:
- **[psa_functions.R](R/psa_functions.R)**: PSA-related utilities, including `generate_psa_samples()` (which skips `derived` distributions and then calls `apply_derived_psa_parameters()`) and `utility_reversal_fraction()`, the `u_p > u_np` diagnostic printed by `CEA.qmd`. **Replaced draws keep their parameter row but not their model index (issue #156)**: `run_psa_analysis()` returns `model_idx`, the cached survival model behind each retained row (equal to the draw number except for replaced draws, which record the model that finally succeeded); `build_psa_obj()` stores it as `psa_obj$model_idx` and `psa_params$model_idx` and binds `additional_params` (the interaction coefficients in the scenario PSA) by `model_idx` after the run, never by draw number. `11_EVPPIs.R` indexes the interaction coefficients by `psa_params$model_idx`. The pre-#156 cache had five replaced rows (708, 716, 1598, 3400, 3552) whose `b_*` columns came from the failed draw's model rather than the model that produced the outcomes.
- **[evppi_functions.R](R/evppi_functions.R)**: Regression-based EVPPI via `voi::evppi()` (`calculate_evppi_regression()`, `run_evppi_analysis()`, `evppi_gam_formula()`), group-consistency checks (`check_evppi_group_consistency()`, `assert_evppi_group_consistency()`), population-level EVPPI scaling (`calculate_population_evppi()`), and interaction coefficient EVPPI (`extract_interaction_coefficients()`, `get_interaction_evppi_params()`, `add_interaction_param_groups()`)
- **[scenario_analysis.R](R/scenario_analysis.R)**: Scenario analysis framework. `define_scenarios()` (WTP 51,000/100,000/150,000 by list-price and decreased nivolumab price) is the PSA/EVPPI scenario set of script 12; `define_utility_scenarios()` (issue #157) holds the deterministic literature-utility scenario of Table S8 (u_np 0.80, u_p 0.65) with a `source` column that is still a **citation placeholder**: the reference must be supplied before publication, and the placeholder text is printed in the Table S8 footnote so it cannot be overlooked

**Survival Modeling Functions**:
- **[para_model_fit.R](R/para_model_fit.R)**: Parametric model fitting helper
- **[survival_plots.R](R/survival_plots.R)**: Survival curve visualization

**Visualization Functions**:
- **[create_tornado_plot.R](R/create_tornado_plot.R)**: DSA tornado diagram generation. `summarise_param_ranges()` and `create_tornado_plot()` take `measure = "NMB_diff"` (change in the strategy's own NMB) or `"INMB_diff"` (change in incremental NMB versus control, the decision-sensitivity measure that Figure 2 plots; issue #156). `09_DSA.R` adds `INMB` and `INMB_diff` to `dsa_results`, keeps a second frame `dsa_results_inmb` whose survival-model endpoints are chosen on incremental NMB, and produces `incremental_impact_summary` alongside `impact_summary`. Under the NMB measure `u_np` and `c_other_last` top every tornado because they move all strategies equally; under the incremental measure the survival family, the second treatment sequence, the nivolumab price and the CRP prevalence lead for the CRP-guided strategy.

**Snapshot/Impact Assessment Functions**:
- **[snapshot_utils.R](R/snapshot_utils.R)**: Snapshot management for bug fix impact assessment, including cache provenance (`collect_cache_provenance()`, `get_git_commit_time()`, `describe_psa_provenance()`; issue #156)

**PSA Ordering-Violation Diagnostics**:
- **[pfs_os_violation_diagnostics.R](R/pfs_os_violation_diagnostics.R)**: Per-draw PFS > OS diagnostics behind `reports/technical/pfs_os_violations.qmd`. `run_pfs_os_violation_diagnostics()` re-predicts the five subgroup curves (control; CRP and TMB/BRAF positive and negative) for each retained PSA row, using its explicit prediction population, joint-population weights, sampled utilities and model_idx, and records, per row and curve, the number of violating weekly points, first/last violating week, maximum excess, the PFS-over-OS excess area (integral of max(PFS - OS, 0), trapezoidal, patient-years; the report calls the shaded region in the figures the trough), the progressed-state area raw and clamped, and the discounted-QALY overstatement of the raw curves (`sum((PFS - OS)^+ (u_np - u_p) cl w_t)`, which equals the raw-minus-clamped QALYs of `model_fun()` exactly). Curves come from `direct_gamma_pop_avg()`, a closed-form gamma AFT evaluation that is about twenty times faster than `flexsurv::predict()`; `validate_direct_curves()` checks it against `generate_psa_population_averaged_predictions()` (max difference 1e-16 on the current cache) and stops above 1e-8. The result is cached in `data/tidy/pfs_os_violations_n{n_samples}.rds` (`pfs_os_violation_cache_path()`) with a fingerprint of the sampling-cache fingerprint, explicit prediction population, PSA parameter table, base parameters, utilities, discount rate, horizon and cycle length; a full pass over 5,000 draws takes about 10 minutes. Gamma only: the closed form stops on any other distribution.

**Archived Functions** (in `archive/`):
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

Positions 25-39 (modeled weeks 24-38) are the **second treatment sequence**, and the model gives it to every patient still progression-free at that point in both arms. See "Cost and Scope Decisions (issue #154)" above; the `Second_sequence` structural scenario zeroes those positions and rebuilds `l_visit`.

### Report and Vignette Housekeeping (issue #157)

Issue #157 changed no model logic and no headline number; the base case, PSA and EVPPI caches are unchanged. What it fixed, so that it is not undone:

- **Enriched-population screening cost.** `08b_enriched_population_analysis.R` charges the test-and-treat arm the expected screening cost per identified positive, `c_test / prevalence` (complete-case target prevalence from `strategies_df` since issue #166), in place of one test, and the standard-of-care arm no test; both are explicit columns (`Enr_Test_Cost`, `Enr_Screening_Cost`) in `enriched_population.qmd`, Table 6 and Table S6. Enriched CRP pairwise ICER moved from 1,587,067 to 1,578,523; TMB/BRAF stays dominated.
- **Restricted mean survival** is the trapezoidal rule over the 521-point weekly grid, restricted to the 10-year horizon and stated as such (`restricted_mean_survival()`); the earlier left Riemann sum counted t = 0 as a full cycle.
- **Hard-coded numbers.** HRs, CIs, p-values, ridge values and counts in `clinical_effectiveness.qmd`, `dag_associations.qmd`, `biomarker_decomposition.qmd` and `table_s2.qmd` are inline R from the objects the document computes. The only typed numbers left are the baseline-CRP counts (8/38 vs 7/35, p = 1) in the clinical report, because `CRP0cat` is not carried into `data`.
- **Table S5** is the full result set with the pairwise comparison against standard of care, the frontier ICER and the frontier status; before #157 it was a byte-identical copy of Table 4.
- **Tornado (Figure 2)** has a lower-/upper-bound legend and omits parameters whose bounds leave the increment unchanged (`create_tornado_plot(drop_empty = TRUE)`, the omitted names are listed under the figure). **CEAC grids** contain the base-case WTP as an evaluated point (Figure 4, poster). **Population EVPPI** in Figure 5 and Table S7 comes from `calculate_population_evppi()` (years 1 to 10); there is no local annuity fallback.
- **Table S1** uses the 68-patient economic cohort and reports censored outcomes descriptively (see the vignette list); **Table S8** is built from the pipeline's scenario definitions, and its literature-utility scenario still needs its citation (see `scenario_analysis.R`).
- **Cycle length and horizon** print as "1 week (1/52 year = 0.0192 years)" and "520 weeks (10 years)" in `input_parameters.qmd`.
- `test_sim_idx_validation.R` and the AGENTS.md command example use `here::here()`; no file in the live tree carries a machine-specific path.

## Important Modeling Considerations

### Ensuring Non-Negative States

The model enforces `p_p = pmax(os - pfs, 0)` to prevent negative progressed state occupancy when PFS and OS curves cross (can happen for a sampled coefficient draw, because joint coefficient draws do not guarantee OS >= PFS).

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

- **Numbered analysis scripts** (`analysis/`): Designed to run sequentially, building on previous steps
  - Scripts 01-03: Core setup and data preparation
  - Scripts 04-11: Main analysis pipeline
  - Script 08b: Enriched population CEA for economic biomarker-positive populations (optional)
  - Script 12: Extended scenario EVPPI analysis (optional)
  - Script 13: Standalone snapshot utility (optional; run via `Rscript`, not sourced)
- **Functions directory** (`R/`): Reusable components that are sourced by analysis scripts
- **Archive directory** (`archive/`): Deprecated/unused code preserved for reference, including the former baseline-characteristics and QALY exploratory notebooks
- **Tests directory** (`tests/`): Validation and diagnostic scripts
- **Quarto reports** (`reports/`): Publication-ready reports with embedded R code, rendered to PDF and Markdown
- **Technical docs** (`reports/technical/`): Bug fix impact reports and technical documentation
- **Vignettes and outputs** (`outputs/vignettes/`, `outputs/figs/`, `outputs/tables/`): figure and table generators and their products
- **Validation and publishing** (`validation/`, `publish/`): external validation reports; the render-and-publish script (see "Repository Layout")

### Retired Files

The B1 cleanup deleted obsolete or superseded files rather than leaving dead entry points in the live tree:

- `analysis/14b_scenario_preview.R` (the full scenario cache from `12_scenario_EVPPIs.R` is the only supported scenario output)
- `R/para_model_fit_table.R`
- `tests/diagnose_prediction_failures.R`
- `tests/test_psa_error_rate.R`
- `tests/test_sampling_convergence_rate.R`
- `tests/test_sex_variable_fix.R`
- `tests/tlr.R`

Issue #157 deleted `outputs/figs/forrest-plot.png`, an orphaned June 2026 image that no vignette wrote (user decision); every file in `outputs/figs/` and `outputs/tables/` now has a generating vignette listed under "Figure Vignettes".

`tests/fit_independent_biomarker_models.R` was moved to `archive/fit_independent_biomarker_models.R`; it was archived, not deleted.

## Quarto Report Architecture

The Quarto reports in `reports/` are self-contained documents that:
1. Set their working directory to the project root using `here::here()`
2. Load all required packages and helper functions
3. Source the necessary analysis scripts (02, 03, etc.) to recreate the analysis environment
4. Generate formatted tables, plots, and results
5. Output to PDF and GitHub-flavoured Markdown with consistent styling

### Report Dependencies & Execution Order

Each report has specific dependencies:

**[para_models.qmd](reports/para_models.qmd)** - Parametric Survival Modeling
- **Sources**: 02, 03, 04, 05
- **Shows**: Survival model fits, AIC/BIC comparisons, goodness-of-fit diagnostics
- **Models displayed**: Ordering-constrained, minimum-combined-AIC pair (currently gamma for OS and gamma for PFS) AND Weibull PH models for reference (issue #68)
- **Tables include**: Model parameter exponents for clinical interpretation
- **Note**: Economic survival fits use the single joint CRP + TMB/BRAF model
- **Selection reporting (issue #153)**: the selected distributions are read from `models$best_fit$os_distribution` / `pfs_distribution` and marked in a "Selected" column of the marginal-AIC tables (the lowest marginal-AIC row need not be the selected one); the joint pair audit `models$ordered_selection$pairs` (combined AIC, OS >= PFS flag, violating time points) is displayed; coefficient tables carry a `Scale` column (baseline parameters on the natural scale, covariates on the log location scale). The same rule applies to `table_s3.qmd` (also writes `outputs/tables/table_s3_pairs.csv`), `table_s4.qmd`, `figure_s1.qmd` and `figure_s2.qmd`, which no longer hard-code `dist = "gamma"` or "Gamma (selected)".

**[clinical_effectiveness.qmd](reports/clinical_effectiveness.qmd)** - Clinical Effectiveness Analysis
- **Sources**: 02, 03 (with report-specific clinical models fitted in the document)
- **Shows**: Baseline characteristics, clinical survival curves, life-years gained; may include TLR clinical analyses independent of the economic model
- **Cohorts (issue #153)**: the primary Firth/standard Cox models, PH diagnostics, ridge sensitivity and the Overall/CRP/TMB-BRAF columns of Table 1 use `data_complete`, the same 68-patient complete-case cohort as the economic survival models (complete on Age, sex, Rx, CRP, TMB/BRAF and both endpoints). TLR completeness is required only for the TLR-complete subset `data_tlr` (65 patients; the 3 missing are early control-arm deaths before the first on-treatment CT), which feeds the TLR columns/row of Table 1, the Cramer's V pairs involving TLR, and the exploratory TLR responder and landmark analyses.
- **Sex coding**: `sex` is the renamed trial column "Sex 0female", so level "0" is female and level "1" is male (`SEX_FEMALE_LEVEL`/`SEX_MALE_LEVEL` in the data-preparation chunk); never take "the last factor level" as female.
- **TLR landmark (issue #155)**: cohorts come from `build_all_tlr_landmark_cohorts(data_tlr)` in `tlr_landmark.R`; `LANDMARK_WK` is read from the week-9 spec, not typed. The primary week-9 landmark section reports how many first-scan progressors and after-landmark scans it retains, and a "Landmark Sensitivity" subsection tabulates the Firth TLR x Rx HR under the week-9, scan-date and week-12 landmarks (`lm_sensitivity`); the Discussion and Summary quote the scan-date and week-12 PFS HRs from that table.
- **Ridge sensitivity**: `glmnet::cv.glmnet(family = "cox", alpha = 0)` on an explicit design matrix with pre-built `crp_x_rx` and `tmb_braf_x_rx` product terms, penalty factor 1 for the CRP and TMB/BRAF main effects and 0 for Age, sex, Rx and the interactions, lambda = `lambda.min` from seeded 10-fold CV; ridge results are point estimates (no SE). The earlier `survival::ridge()` formula with `crp_num:Rx` coded the interaction per arm level (a within-arm CRP slope, not the CRP x Rx contrast) and was removed. Narrative wording about Firth-to-ridge shifts is computed (`shift_word()`, 1.5-fold HR threshold), not hard-coded.

**[input_parameters.qmd](reports/input_parameters.qmd)** - Input Parameters Summary
- **Sources**: 02, 03, 04, 05
- **Shows**: All model input parameters (costs, utilities, prevalence rates, treatment schedules), which costs are sampled versus fixed in the PSA, the nivolumab price provenance, and the post-progression-cost and second-sequence scope limitations (issue #154)

**[CEA.qmd](reports/CEA.qmd)** - Cost-Effectiveness Analysis Report
- **Sources**: 02, 03, 04, 05; runs the base case through shared helpers
- **Requires**: PSA cache from script 10 for probabilistic outputs
- **Shows**: Incremental cost-effectiveness ratios (ICERs), cost-effectiveness plane, decision tables, the utility-ordering diagnostic (fraction of PSA draws with u_p > u_np, zero by construction), and a Scope Limitations section (issue #154)

**[OWSA.qmd](reports/OWSA.qmd)** - One-Way Sensitivity Analysis (Deterministic)
- **Sources**: 02, 03, 04, 05, then 09
- **Requires**: DSA results from script 09
- **Shows**: Tornado diagrams on the strategy's own NMB and, for the guided strategies, on incremental NMB versus standard of care (issue #156), the corresponding impact rankings, the survival-distribution structural sensitivity, and (issues #153, #154) the discount-rate (0%, 8%), time-horizon (5, 20 years), post-progression-cost (EUR 5,000 per quarter) and no-second-sequence structural scenarios with NMB per strategy and the optimal strategy; failed scenarios are listed from `dsa_scenario_status`
- **Survival-family sensitivity (issue #156)**: `09_DSA.R` evaluates every candidate family with the same family for OS and PFS; a family whose OS/PFS pair violates OS >= PFS (from `models$ordered_selection$pairs`) is not run through the model, because the deterministic model stops on an ordering violation, but is listed in `dsa_distribution_status` with its violation count. The report states "7 of 9 candidate families evaluated" in the table caption and tabulates the skipped ones (gengamma, 814 violating curve-time points; genf, 825), instead of silently showing seven rows. The base-case anchor is the ordering-constrained selected pair (`models$best_fit$os_distribution`/`pfs_distribution`, labelled "gamma" or "os/pfs" when they differ), taken from the base-case rows, not the OS family alone

**[EVPPIs.qmd](reports/EVPPIs.qmd)** - Value of Information Analysis
- **Sources**: 02, 03
- **Requires**: PSA results (script 10) and EVPPI results (script 11)
- **Shows**: Expected value of perfect information (EVPI); single-parameter and joint group EVPPI with Monte Carlo standard errors (regression estimator, `voi::evppi()`); the group-versus-largest-member consistency table; an estimator description

**[scenario_effect.qmd](reports/scenario_effect.qmd)** - Scenario Analysis
- **Sources**: 02, 03; loads scenario results generated by script 12
- **Shows**: PSA and EVPPI results under the scenarios of `define_scenarios()` in [scenario_analysis.R](R/scenario_analysis.R), which vary only the willingness-to-pay threshold (EUR 51,000, 100,000, 150,000) and the nivolumab price (base case versus biosimilar EUR 4,641 per administration). Discount-rate and time-horizon variations are not scenarios of this report: they are the deterministic structural scenarios of `dsa_structural_scenarios()` reported in `OWSA.qmd`

**[biosimilar_scenario.qmd](reports/biosimilar_scenario.qmd)** - Biosimilar Nivolumab Pricing Scenario
- **Sources**: 02, 03, 04, 05
- **Shows**: ICER comparison for base case vs biosimilar pricing (EUR 13,923 vs EUR 4,641/dose) for the single joint economic model

**[enriched_population.qmd](reports/enriched_population.qmd)** - Enriched Population Analysis
- **Sources**: 02, 03, 08b
- **Shows**: Enriched (biomarker-positive) ICERs vs base case for CRP and TMB/BRAF
- **Note**: Not cached; re-runs on each render. Requires sampling cache.

**[biomarker_decomposition.qmd](reports/biomarker_decomposition.qmd)** - Biomarker Effect Decomposition
- **Sources**: 02, 03, 04, 05
- **Shows**: Decomposition of biomarker effects on cost-effectiveness outcomes

**[biomarker_distributions.qmd](reports/biomarker_distributions.qmd)** - Biomarker Distributions
- **Sources**: 02, 03
- **Shows**: Biomarker prevalence and distribution analyses, including clinical-only TLR summaries
- **Stratum labels (issue #153)**: arm-by-biomarker strata are built by `make_arm_strata()`, which uses `interaction(..., lex.order = TRUE)` and relabels by level name; never relabel `interaction()` output positionally (its default order is Control/-, Exp/-, Control/+, Exp/+, which swapped two columns in every stratified table)

**[survival_model_specification.qmd](reports/technical/survival_model_specification.qmd)** - Parametric Survival Model Specification
- **Sources**: 02, 03, 04 (uses the `models`, `os_candidates`/`pfs_candidates`, and `models$ordered_selection` objects created by script 04; no refitting)
- **Shows**: (1) regression coefficients of the single joint OS and PFS models for all nine candidate distributions (selected gamma/gamma first, then the previously used Weibull and log-normal), with exp(coefficient) interpreted per family (time ratio / hazard ratio / gamma rate ratio); (2) the current gamma/gamma pair versus the previous unconstrained Weibull-OS/log-normal-PFS pair on the population-averaged subgroup curves, with OS >= PFS violation counts and a PFS-minus-OS gap plot; (3) a gallery of all candidate distributions per subgroup (OS | PFS facets), unlabelled first and then labelled, over the Kaplan-Meier curves
- **Render**: PDF first, then GFM, like every report (see "Rendering Quarto Reports")

**Technical Documentation** (in `reports/technical/`):
- **[age_effect_analysis.qmd](reports/technical/age_effect_analysis.qmd)**: Age effect on survival outcomes
- **[all_parametric_survival_models.qmd](reports/technical/all_parametric_survival_models.qmd)**: Full survival model diagnostics
- **[bug_fix_impact.qmd](reports/technical/bug_fix_impact.qmd)**: Bug fix impact documentation for the single economic model
- **[pfs_os_violations.qmd](reports/technical/pfs_os_violations.qmd)**: PFS > OS ordering violations in the PSA. Sources 02-06 and requires the PSA cache (for the draw-to-row `model_idx` mapping). Reports the number and share of the 5,000 coefficient draws, and of retained PSA rows, with PFS above OS on at least one of the five subgroup curves (overall, per curve, and by number of curves affected); where in the horizon the curves cross and by how much (share of draws violating at each week, extent table, worst-case draw per curve with the excess region shaded); the PFS-over-OS excess area before the clamp (zero after it) and the progressed-state and progression-free areas before and after; and the discounted-QALY consequence of the clamp per curve and per strategy, including its effect on the increments versus standard of care compared with the deterministic base-case increment. Reads the diagnostics cache of `pfs_os_violation_diagnostics.R` (about 10 minutes on first render) and stops unless the sampling, PSA and diagnostics caches carry the same sampling fingerprint and the diagnostics file is younger than the sampling cache; a provenance table lists the three files with dates and fingerprints. Tables use short curve labels (SoC, CRP+, CRP-, TMB/BRAF+, TMB/BRAF-) and are split rather than scaled down; bump the report version on every change. Historical results before issue #166 (see the regenerated report for current values): 2,422 of 5,000 draws (48.4%) violate on at least one curve (control 14.9%, CRP+ 28.2%, CRP- 10.5%, TMB/BRAF+ 33.1%, TMB/BRAF- 9.9%); the share violating rises through the horizon and is highest at week 520; the median violating draw has a maximum excess of 0.004 and an excess area of 0.011 patient-years (maximum 1.28); the mean QALY overstatement of the raw curves is 0.0006 (control) to 0.0014 (TMB/BRAF) per patient, and the clamp shifts the mean CRP increment versus standard of care by only 0.0003 QALYs, so it does not explain the PSA-versus-base-case increment gap recorded in the test protocol table.

**Figure Vignettes** (in `outputs/vignettes/`):
- Each vignette is self-contained (HTML, `embed-resources: true`) and writes its outputs to `outputs/figs/` or `outputs/tables/`; every file listed here is regenerated by the named vignette and by nothing else:
  - `figure1.qmd` → `outputs/figs/figure1.png` (CRP-guided strategy survival curves: `models$best_fit` predicted for the plotted subsets, the control-arm patients and the observed CRP-guided cohort, with 95% percentile ribbons from 500 sampling-cache draws predicted for the same subsets; the shared legend is extracted with a local `extract_legend()` because `cowplot::get_legend()` returns an empty grob under ggplot2 3.5; issue #157)
  - `figure2.qmd` → `outputs/figs/figure2.png` (one-way tornado of the incremental NMB, CRP-guided versus standard of care, with a lower-/upper-bound legend; parameters whose bounds leave the increment unchanged are omitted and listed under the figure; issue #157)
  - `figure3.qmd` → `outputs/figs/figure3.png` and `outputs/figs/figure_s3.png` (Figure S3: the same PSA cost-effectiveness plane drawn with smaller, fainter points)
  - `figure4.qmd` → `outputs/figs/figure4.png` (CEACs, base case and biosimilar pricing, faceted) and `outputs/figs/figure_s4.png` (Figure S4: base-case CEAC alone with a WTP 100,000 reference line); the WTP grid contains the base-case threshold as an evaluated point and strategy labels are looked up by name (issue #157)
  - `figure5.qmd` → `outputs/figs/figure5.png` (population EVPPI by parameter group and scenario, `[GROUP]` rows only)
  - `figure_s1.qmd` → `outputs/figs/figure_s1.png` (Kaplan-Meier versus parametric fits by biomarker subgroup); `figure_s2.qmd` → `outputs/figs/figure_s2.png` (extrapolations under alternative distributions, predicted for the control-arm patients that the Kaplan-Meier reference curve uses; issue #157); `figure_s5.qmd` → `outputs/figs/figure_s5.png` (base-case cost-effectiveness frontier)
  - `table_1.qmd` → `outputs/tables/table_1.csv`; `table_4.qmd` → `outputs/tables/table_4.csv` (dampack frontier results) and `outputs/tables/table_s5.csv` (Table S5: full results with the pairwise comparison against standard of care, the frontier ICER and the frontier status; before issue #157 it was a byte-identical copy of `table_4.csv`); `table_5.qmd` → `outputs/tables/table_5.csv`; `table_6.qmd` → `outputs/tables/table_6.csv` and `table_s6.qmd` → `outputs/tables/table_s6.csv` (enriched-population results from script 08b, including the screening cost per identified positive; issue #157); `table_s1.qmd`, `table_s2.qmd`, `table_s4.qmd`, `table_s7.qmd`, `table_s8.qmd` → the matching `outputs/tables/table_s*.csv`; `table_s3.qmd` → `outputs/tables/table_s3.csv` and `outputs/tables/table_s3_pairs.csv`
  - `table_s1.qmd` → `outputs/tables/table_s1.csv` (patient characteristics on the 68-patient economic complete-case cohort `data_complete`, TLR shown with a "Not assessed" level for the 3 early deaths; follow-up as the reverse Kaplan-Meier median (IQR), deaths and progression-or-death counts, and Kaplan-Meier medians with 95% CI reported without tests; `n` is computed; issue #157)
  - `table_s8.qmd` → `outputs/tables/table_s8.csv` (deterministic scenario table built only from pipeline definitions: `define_scenarios()` for WTP and nivolumab price, `dsa_structural_scenarios()` applied with `apply_structural_scenario()` for discount rate, horizon, post-progression cost and second sequence, and `define_utility_scenarios()` for the literature utilities; issue #157)
- **figure_pfs_os_curves.qmd**: OS and PFS in the same panel, one PNG per strategy (`outputs/figs/figure_pfs_os_curves_{soc,crp,tmb_braf}.png`). Overlays Kaplan-Meier step curves of a patient subset (control arm; biomarker-positive patients from the experimental arm; biomarker-negative patients from the control arm) on the selected joint models (`models$best_fit`) predicted for exactly those patients and averaged (issue #157). Sources 02/03/04 only; unnumbered for now.

**KM overlay convention (issue #157)**: a Kaplan-Meier curve of a subset is only ever overlaid on a parametric curve predicted for that same subset (each patient's own covariates and actual arm, then averaged), so the overlay is a fit check. The economic model's curves (`predictions`, `l_params_base`) are standardised to all 68 complete-case patients with `Rx` set by strategy; they are never drawn over a subset KM curve, because the covariate mix differs (the control arm has 6/32 CRP-positive patients against 23/68 in the full cohort) and the gap would read as lack of fit. Captions and notes in `figure1.qmd`, `figure_s2.qmd` and `figure_pfs_os_curves.qmd` say which curves are shown.

**Clinical effectiveness vignettes** (in `outputs/vignettes/`): Files for the clinical effectiveness paper use the prefix `clin_effect_` followed by the paper figure/table label. Each vignette is self-contained (HTML, `embed-resources: true`) and saves its output to `outputs/figs/` (figures, via `ggsave`) or `outputs/tables/` (CSV, via `write.csv`). Naming examples:
- `clin_effect_figure1.qmd` → `outputs/figs/clin_effect_figure1.png` and `outputs/figs/clin_effect_figure1.eps` (the same figure in PNG and in EPS for journal submission; six-panel Kaplan-Meier grid: OS left / PFS right columns; CRP, TMB/BRAF, and TLR rows, with the TLR row on the week-9 landmark cohorts from `tlr_landmark.R`; every panel carries its own time-axis title because the time origin differs between rows, and the caption, built as `fig1_caption` and passed through `!expr`, states both landmark cohort sizes, the retained first-scan progressors and the after-landmark scans; issue #155)
- `clin_effect_figure_s1.qmd` → `outputs/figs/clin_effect_figure_s1.png` (full DAG)
- `clin_effect_table_s2.qmd` → `outputs/tables/clin_effect_table_s2.csv` (DAG association consolidated summary from `dag_association_tests.R`: 18 tested edges including both `TxCRP` edges, 19 implied independencies with two marked as holding by construction, landmark `TLR -> PFS` and time-dependent `PFS -> OS`; issue #155)
- `clin_effect_figure_simplified_dag.qmd` → `outputs/figs/clin_effect_figure_simplified_dag.png` (simplified DAG; unpublished, off the numbered figure sequence)
- `clin_effect_figure_sensitivity_dag.qmd` → `outputs/figs/clin_effect_figure_sensitivity_dag.png` (sensitivity DAG; unpublished, off the numbered figure sequence)

**Poster vignette** (in `outputs/vignettes/`): `poster_biomarker_correlation.qmd` supports the SMDM poster "Modeling Multiple Biomarkers in Precision Oncology: How Ignoring Biomarker Correlations Can Reverse Clinical Conclusions". Self-contained (HTML, `embed-resources: true`). **Standalone and intentionally off-pipeline**: it re-introduces TLR as a third biomarker-guided strategy purely to illustrate the methodological point, even though the main economic model deliberately excludes TLR as a post-randomization mediator. It sources 02 to 06 (06 for the sampling cache and `sample_survival_coefficients()`) + `model_fun`/`calculate_outcomes`/`prediction_functions`, then overrides `get_strategies()`/`get_biomarkers()` locally to add `tlr`; it does not modify the main pipeline. **Stated simplifications (issue #157)**: TLR is modelled as if it were known at treatment selection at no separate test cost (it is read from the first on-treatment CT, which is in the monitoring schedule), its prevalence is the experimental-arm response rate (19/36, not the pooled 41/65, because it is measured on treatment and differs by arm), and all prevalences are those of the 65-patient TLR-complete cohort; the table notes say so, and the conclusion sentence is computed from the dampack `Status` column, not typed. The CEAC uses the pipeline's survival sampler: joint-model draws come from the sampling cache aligned to the cached economic-parameter draws through `model_idx`, single-biomarker models are drawn with `sample_survival_coefficients()` (same distributions, 1,000 draws), the standard-of-care curve is the joint model with `Rx = control` in both panels (issue #151), patient curves are predicted with `predict()` and averaged with the same draw-specific joint-population weights for control and both guided strategies (issue #166), and each draw runs `model_fun(determpsa = "psa")`, which caps PFS > OS crossings instead of stopping; the WTP grid contains the base-case threshold and the caption reports the probabilities at that threshold from the curves. Rendering takes about 11 minutes. Outputs:
- `outputs/figs/poster_forest.png` (treatment x biomarker interaction, single-biomarker vs joint correlation-adjusted model, parametric AFT time ratios, PFS + OS)
- `outputs/figs/poster_ceac.png` (two-panel cost-effectiveness acceptability curves: single-biomarker vs joint modelling, at biosimilar nivolumab pricing)
- `outputs/tables/poster_cea.csv` (three-strategy CEA), `outputs/tables/poster_forest_data.csv` (forest interaction estimates), and `outputs/tables/poster_ceac_data.csv` (CEAC probabilities)

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
source(here("analysis/02_setup_and_global_variables.R"))
source(here("analysis/03_biomarker_strategies.R"))
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

**Output Location**: PDFs and `.md` files are generated next to the `.qmd` files (`reports/`); PDFs are gitignored and published with `publish/publish_reports.R`.

### Common Quarto Report Issues

1. **"Object not found" errors**: Run the required analysis scripts first (especially 02-11, in sequence)
2. **Sampling cache missing**: Run [06_sampling.R](analysis/06_sampling.R) to generate sampling models
3. **Rendering hangs**: Some reports (especially CEA, EVPPI) may take minutes to render due to re-sourcing analysis scripts
4. **Stale caches after model changes**: Economic survival model or strategy-set changes require regeneration of sampling, PSA, EVPPI, scenario, and snapshot outputs.
5. **Results changed after methodological updates**: If cost-effectiveness results differ from earlier versions, check if methodological fixes were applied. Issues #69 and #70 (Nov 2024) changed survival prediction methodology from reference patient to population averaging/individual sampling. This **should** change results - it's a methodological improvement. Regenerate both sampling cache and PSA cache after these fixes. See issue #64 for impact documentation approach.

## Cache Management

The analysis uses several cache systems to speed up computation:

### 1. Sampling Cache (Survival Models)

**Location**: `data/tidy/sampling_models_n*.rds`
**Purpose**: Joint multivariate-normal OS/PFS coefficient draws with paired-bootstrap cross-endpoint dependence (one `joint` component; issue #159)
**Generation**: Script [06_sampling.R](analysis/06_sampling.R) (several minutes for the paired covariance bootstrap)
**Size**: Compact coefficient matrices and two original fits, including covariance bootstrap diagnostics; no stored fit per bootstrap replicate
**Validity (issues #156, #169)**: script 06 compares a versioned fingerprint of the joint formulas, selected distributions, sample size, seed, sampling implementation and runtime dependency versions. Its data identity covers only formula variables (including both endpoints) and `ID`, retaining row order, factor levels and contrasts. Canonical values are hashed using serialization version 2 and SHA-256, which avoids dependence on ALTREP materialisation. The ignored local cache stores `canonical_data`, per-column hashes in `fingerprint_inputs$data_columns`, and `runtime` (R, packages, library paths and locale) for diagnosis. Logs print column names and hashes, never patient values. The review-session data hash (`6cc5a549...`) was reproduced during this fix, but the alternate pre-review hash (`13f4cb24...`) and its cause were not. `rlang::hash()` does distinguish compact from materialised vectors with identical values; this was independently demonstrated, without establishing it as the historical cause.

**Read-only consumers**: `setup_report()` temporarily sets `options(metimmox.sampling_allow_regenerate = FALSE)` and restores the caller's option on exit. Script 06 also refuses regeneration while `knitr.in.progress` is true. Missing, legacy or stale sampling caches stop report setup with a pipeline instruction. Ordinary pipeline sourcing retains regeneration by default. Tests that need caches must set `options(metimmox.cache_dir = <temporary directory>)`; the index-validation test uses synthetic functions/data and never sources the pipeline.

### 2. PSA Cache (Analysis Results)

**Location**: `data/tidy/psa_obj_{ipd|correct}.rds` and `psa_params_{ipd|correct}.rds`.
**Generation**: script [10_PSA.R](analysis/10_PSA.R).

**Validity (issue #167)**: the PSA input fingerprint covers sampling provenance, base parameters, schedules/curves, distributions, strategies, time settings, seed, the common prediction population, parsed calculation source and loaded function definitions, relevant package versions, R version, contrasts and RNG kind. Uncommitted source changes and interactive function-body edits invalidate it automatically. Code dependencies are enumerated in `calculation_identity()` in [cache_provenance.R](R/cache_provenance.R); extend that list when introducing a calculation dependency. Comments and line-ending changes do not invalidate results.

The existing two-file layout is retained. `bind_psa_pair()` attaches identical `pair_metadata` to the outcome object and parameter-table attributes: a content-addressed generation ID, input fingerprint, ordered parameter hash and outcome hash. The outcome object also records `sim`. `validate_psa_pair()` checks hashes, unique positive draw IDs, row counts and exact `sim`/`model_idx` alignment. Mixed generations after an interrupted save or partial restore are rejected. Both `load_psa_cache()` and `load_psa_params_cache()` validate the pair; report consumers use `load_current_psa_cache()` (or the equivalent explicit current expectations). Never stamp legacy cache metadata to make it pass: rebuild with script 10.

### 3. EVPPI Cache

**Location**: `data/tidy/evppi_results_{ipd|correct}.RData`.
**Generation**: script [11_EVPPIs.R](analysis/11_EVPPIs.R), which first reloads a validated PSA pair before deriving interaction columns.

`evppi_fingerprint` and `evppi_fingerprint_inputs` identify the exact PSA generation/result, WTP, full parameter/group configuration, estimator code/settings/dependencies, seed and population-scaling inputs. `load_evppi_cache()` validates this identity and recalculates EVPI cheaply from the current PSA and WTP. The retained `evppi_psa_fingerprint` supports older provenance displays but is insufficient by itself. Changing WTP requires EVPPI regeneration, not PSA regeneration.

### 4. Scenario EVPPI Cache

**Location**: `data/tidy/scenario_evppi_results_{ipd|correct}.rds`.
**Generation**: script [12_scenario_EVPPIs.R](analysis/12_scenario_EVPPIs.R), after validating the current base PSA.

The scenario-specific fingerprint covers the full expected PSA identity, scenario definitions (including WTPs and prices), parameter/groups, estimator/scenario implementation, seed and population inputs. Each scenario's embedded PSA outcome/parameter pair is bound and validated. `scenario_effect.qmd`, Figures 4/5 and Table S7 use `load_scenario_cache()` with current expectations. Figure 3, Table 5 and the poster's input PSA also use validated loaders; Figure 1 receives the sampling cache validated by script 06 in read-only report setup.

### 5. PFS/OS Violation Diagnostics Cache

**Location**: `data/tidy/pfs_os_violations_n{n_samples}.rds`.
**Generation**: `run_pfs_os_violation_diagnostics()` in [pfs_os_violation_diagnostics.R](R/pfs_os_violation_diagnostics.R), called by the technical report (about 10 minutes for 5,000 draws).
**Validity**: sampling fingerprint, utilities, effect discount rate, horizon, cycle length and number of draws. The report's sampling input must already be current; rendering cannot rebuild it.

**Migration for #167/#169, completed with #166 (18 September 2026)**: regenerated sampling, the 5,000-draw PSA, EVPPI, and both 5,000-draw scenario configurations from the repaired calculation sources. Cleared the active root caches and installed these freshly generated, fingerprint-validated results before taking the fixed snapshot. All 296 live report/cache contract checks pass without skips. The legacy PSA-versus-base-case numerical criterion still fails, as recorded below; it is distinct from cache validity.

### Cache Workflow

1. **First run**: [06_sampling.R](analysis/06_sampling.R) generates sampling cache
2. **PSA uses sampling cache**: [10_PSA.R](analysis/10_PSA.R) generates PSA cache
3. **EVPPI uses PSA cache**: [11_EVPPIs.R](analysis/11_EVPPIs.R) generates EVPPI cache
4. **Subsequent runs**: scripts 06 and 10 reuse valid caches; scripts 11 and 12 regenerate their respective analyses when explicitly sourced.
5. **Manual invalidation**: Delete specific cache file(s) to regenerate

## Snapshot System (Bug Fix Impact Assessment)

The repository includes a snapshot comparison system for assessing the impact of bug fixes on model results.

### Components

- **[13_save_snapshot.R](analysis/13_save_snapshot.R)**: Saves single-model snapshots for the joint economic survival model. **Interactive and standalone, NOT part of the cache pipeline** — it reads the issue number and `baseline|fixed` from stdin, so run it as `Rscript analysis/13_save_snapshot.R <issue#> <baseline|fixed>` (sourcing it non-interactively just errors on the empty prompt). Optional; not required to regenerate report caches.
- **[snapshot_utils.R](R/snapshot_utils.R)**: Utility functions for snapshot management
- **[compare_snapshots.R](tests/compare_snapshots.R)**: Compares before/after snapshots to quantify changes
- **[bug_fix_impact.qmd](reports/technical/bug_fix_impact.qmd)**: Report documenting bug fix impacts

**Issue #159 snapshot provenance**: baseline and fixed snapshots precede the requested commit and therefore both carry HEAD `6ca19ab`; their PSA fingerprints and result hashes differ. The baseline was saved before any source change. The active sampling, PSA pair, EVPPI, scenario and ordering-diagnostic caches were then cleared and regenerated, and reports/publication outputs re-rendered before the fixed snapshot. The cumulative impact report selects snapshots by their saved timestamps, because #159 was implemented after #166.

**Issue #166 snapshot provenance**: both baseline and fixed snapshots were taken before the requested fix commit, so both filenames carry HEAD 96d0710; their PSA fingerprints and result hashes differ. The baseline PSA was freshly generated under the original model. Its writer omitted pair metadata, so that newly generated pair was bound and EVPPI regenerated before applying the model fix; metadata records `snapshot_writer_pair_repaired = TRUE`. The writer now calls `bind_psa_pair()` before saving new PSA results. Historical caches were not relabelled as current.

### Snapshot Storage

**Location**: `data/output/snapshots/`
**Format**: `snapshot_NN_[baseline/fixed]_HASH.rds` and `psa_NN_[baseline/fixed]_HASH.rds`
**Content**: Base case results, PSA results, metadata (git commit, timestamp, issue number) and, from issue #152 onwards, the EVPPI cache (`evppi$evpi`, `evppi$evppi_results`) when it exists at snapshot time; `bug_fix_impact.qmd` adds an EVPPI before/after table whenever both snapshots of a pair carry it

**PSA provenance (issue #156)**: the `psa_NN_*.rds` file is a copy of the PSA cache, so before issue #156 a "fixed" snapshot could carry the pre-fix PSA unchanged (the psa files of #146 baseline, #146 fixed and #147 baseline are byte-identical, as are #147 fixed, #149 baseline and the live cache at the time) and the impact report showed "no PSA change" by construction. `metadata$caches` now records md5 and mtime of the sampling, PSA object, PSA parameter, EVPPI and scenario caches, the PSA fingerprint, `git_commit_time`, and the flags `psa_stale_reason`, `psa_regenerated`, `evppi_regenerated`, `scenario_evppi_stale`. A **fixed** snapshot regenerates the PSA (about 20 minutes) and then the EVPPI cache (script 11) when the PSA cache's fingerprint differs from the current inputs or its mtime predates the HEAD commit (user decision, issue #156); the scenario-EVPPI cache is not regenerated (40 minutes) but its staleness is recorded and printed. `compare_snapshots.R` and `bug_fix_impact.qmd` print `describe_psa_provenance()`, which says outright when the two PSA files are the same cache.

### Workflow

1. Save "baseline" snapshot before applying a fix
2. Apply fix
3. Save "fixed" snapshot after applying a fix
4. Run `compare_snapshots.R` to quantify the impact
5. Render `bug_fix_impact.qmd` to document changes

## Test Files

The test suite in `tests/` includes:

- **[test_survival_ordering.R](tests/test_survival_ordering.R)**: Verifies deterministic OS >= PFS ordering and PSA-only handling of resampled crossings
- **[test_biomarker_test_cost_mapping.R](tests/test_biomarker_test_cost_mapping.R)**: Verifies data-driven diagnostic-test cost assignment
- **[test_canonical_prevalence.R](tests/test_canonical_prevalence.R)**: Verifies a common target population, zero identical-treatment increments under marginal and joint prevalence changes, and no PSA endpoint-missingness leak
- **[test_prediction_population_live.R](tests/test_prediction_population_live.R)**: Confirms zero identical-treatment increments at every prevalence DSA endpoint on trial data, and that weighted diagnostic curves reproduce PSA QALYs using sampled utilities and repeated model indices. Uses temporary caches only.
- **[test_psa_fallback_reporting.R](tests/test_psa_fallback_reporting.R)**: Verifies PSA fallback-rate reporting and its failure threshold
- **[test_cache_provenance.R](tests/test_cache_provenance.R)**: Synthetic cache contracts for #167/#169: fresh-session value hashes, read-only sampling guards, interactive code edits, mixed/reordered PSA files, draw indices, independent EVPPI/scenario inputs, cheap EVPI reconciliation and forced index-test failure status. Uses temporary caches only.
- **[test_sim_idx_validation.R](tests/test_sim_idx_validation.R)**: Verifies PSA resampling-index validation with a synthetic fixture; sources no analysis scripts and exits nonzero on assertion failure
- **[test_psa_basecase_alignment.R](tests/test_psa_basecase_alignment.R)**: Verifies structurally that the PSA control curve is the joint model with `Rx = control`, then checks that PSA strategy means (costs and QALYs) sit within 5 Monte Carlo standard errors of the base case and prints the incremental comparison (issue #151); skips the numerical check when the PSA cache is absent. The numerical criterion is a known failure at n_sim = 5000 (see the test protocol table)
- **[test_psa_zero_uncertainty.R](tests/test_psa_zero_uncertainty.R)**: With zero joint coefficient covariance and fixed economic parameters/population weights, all deterministic and PSA costs agree within EUR 1e-8 and QALYs within 1e-12 (issue #159). Uses trial data, no production cache.
- **[test_joint_survival_sampling.R](tests/test_joint_survival_sampling.R)**: Synthetic tests of paired endpoint refits, calibrated marginal covariance preservation, joint dependence, RNG reproducibility and invalid covariance guards (issue #159).
- **[test_evppi_estimator.R](tests/test_evppi_estimator.R)**: Verifies the regression EVPPI estimator on a synthetic problem with a closed-form answer, a pure-noise parameter, group-versus-member consistency (including the additive formula for more than four parameters), seed-reproducible standard errors, and NA-not-zero failure reporting (issue #152)
- **[test_sampling_rework.R](tests/test_sampling_rework.R)**: Parameter-specification contracts — PSA membership excludes the fixed unit prices, `u_p` is derived rather than drawn, EVPPI groups exclude derived parameters and drop `all_costs` when only one cost group is sampled, the DSA still covers the fixed prices but never `u_decrement`, the derived utility never exceeds `u_np`, and one-sided structural scenarios declare only their differing endpoint (issue #154); plus the script-06 distribution-resolution helper
- **[test_report_contracts.R](tests/test_report_contracts.R)**: Cache and report wiring contracts; for the EVPPI cache it requires the interaction rows with finite standard errors, the joint `interaction_all` group, and no group below its largest member; since issue #156 also that the PSA cache carries `model_idx` and a fingerprint aligned with `psa_params`, that the EVPPI cache was built from that PSA fingerprint, and that the sampling cache is a single fingerprinted `joint` multivariate-normal component that the PSA and scenario caches reference
- **[test_mvn_sampling.R](tests/test_mvn_sampling.R)**: Multivariate-normal sampling and draw/model alignment (issue #156) — a flexsurvreg object rebuilt from a drawn coefficient vector reproduces the closed-form gamma survival to 1e-10 on simulated data, the accessor returns the legacy sample shape and rejects out-of-range indices, `get_joint_sampling_models()` resolves both cache layouts, `build_psa_obj()` records `model_idx = c(1, 3, 3, 4)` when draw 2 is re-run with model 3 and binds the model-derived column by that index, and the sampling/PSA fingerprints change with every input and not with column or list order
- **[test_sampling_failure_fallback.R](tests/test_sampling_failure_fallback.R)**: The sampler produces finite, independent OS and PFS draw matrices with `n_failed = 0` for converged fits and stops, naming the outcome, on a fit without a covariance matrix (issue #156; previously tested the bootstrap's original-fit substitution)
- **[test_rng_reproducibility.R](tests/test_rng_reproducibility.R)**: PSA draws and coefficient draws are identical under the same seed regardless of prior RNG use and differ under another seed; the cache seed check resolves the joint component
- **[test_dag_landmark_contracts.R](tests/test_dag_landmark_contracts.R)**: DAG-validation and landmark contracts (issue #155) — every DAG-implied edge and conditional independence is tested or classified latent/definitional, both `TxCRP` edges and all six `TxCRP` independencies are covered, a DAG edit that adds an untested edge stops `run_dag_edge_tests()`; landmark cohorts contain only endpoint times strictly after the landmark, the scan-date and week-12 landmarks retain no first-scan progressor, the week-9 landmark reports the ones it retains; `PFS -> OS` is estimated on counting-process rows with a time-dependent progression flag (HR > 1) and `TLR -> PFS` on the landmark cohort. Structural checks run without the trial data; data checks are skipped when `data/tidy/METIMMOX.rds` is absent
- **[compare_snapshots.R](tests/compare_snapshots.R)**: Compares before/after snapshots for bug fix impact assessment
- **[para_models.Rmd](tests/para_models.Rmd)**: Parametric model fit validation and diagnostics
- **[snapshot.R](tests/snapshot.R)**: Helper script for running snapshot saves

### Test protocol and results

Focused executable regression tests live in [`tests/`](tests/). The durable record of the broader 12 July 2026 black-box and extreme-value run lives in [`validation/validateHE_Opus_2026-07-12/findings/findings_blackbox.json`](validation/validateHE_Opus_2026-07-12/findings/findings_blackbox.json); that run used an external scratch harness, so the JSON record, rather than the harness itself, is retained in this repository. Unless a row specifies otherwise, numeric comparisons use a tolerance of `1e-10`; any non-finite value, unexpected warning/error, or unmet criterion is a failure.

| Invariant / boundary test | Predefined pass criterion | Location | Recorded result |
|---|---|---|---|
| OS/PFS ordering | OS >= PFS at all 521 modeled time points for control and all four economic biomarker subgroups; an injected deterministic crossing must error, while a PSA crossing must be reported and capped | [`test_survival_ordering.R`](tests/test_survival_ordering.R) | Pass after ordering-constrained distribution selection |
| Cohort conservation | For every strategy and cycle, PF + P + D = 1 and each occupancy is within [0, 1] | Black-box record (`BB-TR0/1`) in [`findings_blackbox.json`](validation/validateHE_Opus_2026-07-12/findings/findings_blackbox.json) | Pass for control, CRP, and TMB/BRAF |
| Utilities = 1 | With both state utilities set to 1, discounted QALYs equal discounted life-years | Black-box record (`BB-U1`) | Pass |
| Utilities = 0 | With both state utilities set to 0, total QALYs equal 0 for every strategy | Black-box record (`BB-U0`) | Pass |
| Costs = 0 | With every drug, test, visit, follow-up, and end-of-life unit cost set to 0, total cost equals 0 for every strategy | Black-box record (`BB-C0`) | Pass |
| No mortality | With OS and PFS fixed at 1, death occupancy remains 0 and undiscounted life-years equal the stated 10-year horizon | Black-box record (`BB-M0`) | Fail (minor): 521 weekly grid points produce 10.019 rather than 10.000 life-years; equal across strategies |
| Near-certain mortality | With survival forced near 0 after baseline, more than 99% of the cohort is dead by cycle 3 | Black-box record (`BB-M1`) | Pass |
| Determinism | Two deterministic runs with identical inputs produce bitwise-identical costs and QALYs | Black-box record (`BB-RE`) | Pass |
| PSA centred on base case | For every strategy, the PSA mean cost and mean QALYs are within 5 Monte Carlo standard errors of the deterministic base-case values; the control PSA curve must be the joint model with `Rx = control` (structural check, no cache needed) | [`test_psa_basecase_alignment.R`](tests/test_psa_basecase_alignment.R) | Structural check: Pass. Numerical check: **Fail (known)** after #159 at n_sim = 5000: four of six level comparisons exceed the unchanged five-SE criterion. QALY differences are +0.03682, +0.02255 and +0.03035 (8.75, 5.92 and 7.78 SE) for control, CRP and TMB/BRAF; cost differences are -0.73, -4.75 and -6.11 SE. CRP incremental QALYs are 0.00628 versus 0.02054 deterministic (-5.34 SE); TMB/BRAF -0.01800 versus -0.01153 (-2.61 SE). No recentering or threshold relaxation. The zero-uncertainty regression passes: nonlinear averaging explains why uncertain outcome means need not equal fitted-parameter results. Crossings decreased from 2,422/5,000 (48.44%) to 1,198/5,000 (23.96%); the clamp narrows the CRP gap by about 0.0002 QALYs and does not explain it. Before #159, after #166: CRP incremental QALYs 0.00175 (-7.45 SE), TMB/BRAF -0.02256 (-4.90 SE); QALY level differences 6.61 to 10.49 SE. |
| Zero-uncertainty PSA | Zero joint coefficient covariance, fixed economic parameters and population weights reproduce all deterministic costs within EUR 1e-8 and QALYs within 1e-12, without fallback | [`test_psa_zero_uncertainty.R`](tests/test_psa_zero_uncertainty.R) | Pass (#159) |
| Joint OS/PFS covariance | Paired endpoint refits; exact fitted marginal covariance blocks; realised joint covariance within 5 MCSE; seed reproducibility | [`test_joint_survival_sampling.R`](tests/test_joint_survival_sampling.R), [`test_report_contracts.R`](tests/test_report_contracts.R) | Pass (#159): 996/1,000 bootstrap pairs retained; corresponding endpoint correlations 0.36 to 0.68; maximum realised covariance deviation 2.46 MCSE |
| EVPPI estimator | On a synthetic two-strategy problem with incremental NMB = theta + noise, theta ~ N(300, 1000^2): EVPPI(theta) within 4 SE + 2% of the closed-form value E[max(theta,0)] - max(E[theta],0) evaluated on the realised draws, within 10% of the analytic value, and below EVPI; a pure-noise parameter within max(4 SE, 2% EVPI) of 0; the groups {theta, noise} and the five-parameter additive group within 4 combined SE + 2% EVPI of the single-parameter value; no group-consistency violation, while a constructed shortfall of 50% of EVPI must raise an error; identical point estimates and SEs under the same seed, different SEs under another seed; missing or constant parameters return NA with an error message | [`test_evppi_estimator.R`](tests/test_evppi_estimator.R) | Pass (issue #152): analytic 266.8, realised-draw value 273.3, estimate 270.0 (SE 2.7), EVPI 285.1; noise 0.000 (SE 0.000); group 270.3 (SE 2.6); five-parameter additive 270.3 (SE 2.5) |

Run focused tests from the repository root with `"C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe" tests/<test-file>.R`. A test passes only if it exits with status 0 and all documented assertions succeed. Update the recorded result whenever model logic or the corresponding acceptance criterion changes; do not overwrite a known failure with a looser criterion.

## Clinical Context

Based on the METIMMOX clinical trial (NCT03388190) comparing alternating chemotherapy + immunotherapy versus chemotherapy alone in MSS/pMMR metastatic colorectal cancer patients. The analysis focuses on identifying biomarker-based subgroups that may benefit from adding immunotherapy to standard chemotherapy.
