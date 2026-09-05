# METIMMOX-1

Cost-Effectiveness Analysis of Biomarker-Guided Immunotherapy in Metastatic MSS/pMMR Colorectal Cancer

## Overview

This repository contains the R code for a cost-effectiveness analysis comparing two pre-immunotherapy biomarker strategies that guide the addition of immunotherapy (PD1/PDL1 inhibitor) to standard of care treatment for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer patients receiving first-line treatment.

**METIMMOX** stands for: **Colorectal Cancer METastasis - Shaping Anti-tumor IMMunity by OXaliplatin**

### Biomarker Strategies Evaluated

1. **CRP Strategy**: C-reactive protein <5 mg/L, measured at week 4 (cycle 3 day 1) after two FLOX cycles common to both arms and before the first nivolumab dose
2. **TMB/BRAF Strategy**: Tumor mutation burden ≥9 mut/MB or presence of a BRAF mutation (baseline next-generation sequencing)

Both biomarkers are available before the decision to add immunotherapy is made. CRP is not a baseline (pre-randomization) measurement: the trial gives two cycles of FLOX to every patient before the first nivolumab dose, and the CRP used here is the value at that decision point (issue #150). The baseline (cycle 1 day 1) CRP is retained in the data as `CRP0` but is not used.

**Clinical-only TLR analysis**: Tumor lesion reduction (TLR) is retained in DAG, clinical effectiveness, and biomarker distribution reports, but is excluded from the economic model because it is a post-randomization mediator measured on treatment rather than a treatment-selection biomarker.

Each economic strategy is compared against standard of care alone (platinum-based Nordic FLOX regimen without immunotherapy).

The deployed analysis uses one joint economic survival model (`Age + sex + Rx + crp*Rx + tmb_braf*Rx`) for both biomarker-guided strategies.

## Target Audience

This repository is designed for **health economists** and researchers developing decision-analytic models in R. The code provides a framework that can be adopted and adapted for similar cost-effectiveness analyses.

## Repository Structure

```
METIMMOX-1/
├── scripts/
│   ├── R/
│   │   ├── analysis/      # Numbered analysis scripts (01-13)
│   │   ├── functions/     # Reusable model functions
│   │   ├── tests/         # Validation and diagnostic scripts
│   │   └── archive/       # Deprecated code (for reference)
│   └── QMD/
│       ├── report/        # Publication-ready Quarto reports (PDF)
│       ├── vignettes/     # Figure generation (Figures 1-4, supplemental)
│       └── technical_docs/# Technical documentation (bug impacts, methodological analyses)
├── data/                  # Data files (not included - confidential)
│   ├── tidy/             # Processed data and caches
│   └── output/           # Analysis outputs and snapshots
└── METIMMOX-1.Rproj      # RStudio project file
```

### Key Files

- **Analysis scripts** (`scripts/R/analysis/`): Numbered R scripts (01-13) containing the core decision-analytic model workflow and optional extended analyses
- **Functions** (`scripts/R/functions/`): Reusable functions including `model_fun.R` (main model), `calculate_outcomes.R`, `cea_helpers.R` (single-model CEA execution and summary helpers), and sensitivity analysis utilities
- **Quarto reports** (`scripts/QMD/report/`): Publication-ready PDF reports covering cost-effectiveness, clinical effectiveness, sensitivity analyses, biomarker decomposition, and more
- **Vignettes** (`scripts/QMD/vignettes/`): Publication figure generation scripts (Figures 1-4 and supplemental plots)
- **Technical docs** (`scripts/QMD/technical_docs/`): Bug fix impact assessments and methodological analyses
- **Tests** (`scripts/R/tests/`): Validation, diagnostic, and convergence testing scripts

The shared report and clinical-analysis helper modules introduced in B2 are:

- `assoc_tests.R`: association tests and result formatting used by DAG reports and vignettes
- `cox_extract.R`: Firth-corrected Cox fitting and tidy coefficient extraction
- `dag_helpers.R`: shared DAG specifications, preparation, and plotting
- `report_format.R`: shared strategy/biomarker labels and economic-result number formatting

## Installation

### Prerequisites

This project is built in R. The following packages are used throughout the analysis:

```r
# Core packages for survival analysis and modeling
- survival
- flexsurv
- survminer
- gems
- mstate

# Economic evaluation packages
- dampack
- darthtools

# Data manipulation and visualization
- dplyr
- tidyverse
- ggplot2
- readxl
- scales
- gridExtra
- reshape2

# Report generation
- knitr
- kableExtra
- flextable
- officer

# Other utilities
- pacman (for package management)
- mvtnorm
- Matrix
- here
```

### Setup

```r
# Clone the repository
git clone https://github.com/ben-geisler/METIMMOX-1.git

# Set working directory
setwd("METIMMOX-1")

# Install required packages using pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools, readxl, dplyr, tableone, ggplot2, flexsurv,
               survival, survminer, gems, mstate, tidyverse, xtable,
               darthtools, dampack, mvtnorm, Matrix, here, voi,
               knitr, kableExtra, flextable, officer, scales, gridExtra, reshape2)
```

## Usage

Run the numbered analysis scripts in the `scripts/R/analysis/` folder in sequential order:

```r
# Core setup
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")

# Survival analysis
source("scripts/R/analysis/04_parametric_survival_analysis.R")
source("scripts/R/analysis/05_basecase_input_parameters.R")

# Survival resampling (first run generates cache - takes time)
source("scripts/R/analysis/06_sampling.R")

# Model execution
source("scripts/R/analysis/07_traces.R")
source("scripts/R/analysis/08_basecase_analysis.R")

# Sensitivity analyses
source("scripts/R/analysis/09_DSA.R")   # Deterministic sensitivity analysis
source("scripts/R/analysis/10_PSA.R")   # Probabilistic sensitivity analysis
source("scripts/R/analysis/11_EVPPIs.R") # Expected value of perfect partial information
```

**Notes**:
- The obsolete `04_baseline_characteristics.Rmd` and `05_QALYs.Rmd` exploratory notebooks are retained in `scripts/R/archive/`; their reusable EQ-5D-5L scoring functions are in `scripts/R/functions/eq5d5l_utility.R`
- The clinical trial dataset is confidential and not included in this repository
- First runs of scripts 06 (sampling) and 10 (PSA) generate caches and may take significant time; subsequent runs are much faster

### Treatment-Schedule Semantics

Treatment schedule vectors are aligned to the zero-origin weekly grid `time_points <- seq(0, time_horizon, by = 1)`. R vector position `i` therefore represents modeled week `i - 1`; the numeric subscripts used to construct a schedule are vector positions, not week labels.

- Nivolumab in the experimental arm: modeled weeks 4, 6, 12, 14, 28, 30, 36, and 38 (R positions 5, 7, 13, 15, 29, 31, 37, and 39)
- FLOX in the experimental arm: modeled weeks 0, 2, 8, 10, 24, 26, 32, and 34 (R positions 1, 3, 9, 11, 25, 27, 33, and 35)
- FLOX in the control arm: modeled weeks 0, 2, 4, 6, 8, 10, 12, 14, 24, 26, 28, 30, 32, 34, 36, and 38

### Optional and Extended Analyses

Additional analysis scripts provide extended functionality:

```r
# Extended analyses (optional)
source("scripts/R/analysis/12_scenario_EVPPIs.R")  # Scenario-based EVPPI analysis

# Standalone snapshot utility (run from a terminal, not with source())
# Rscript scripts/R/analysis/13_save_snapshot.R <issue_number> <baseline|fixed>
```

**When to use**:
- `12_scenario_EVPPIs.R`: For scenario-specific value of information analysis
- `13_save_snapshot.R`: For documenting model state before/after bug fixes or methodological changes

### Retired Files

The B1 cleanup deleted these obsolete or superseded files:

- `scripts/R/analysis/14b_scenario_preview.R`
- `scripts/R/functions/para_model_fit_table.R`
- `scripts/R/tests/diagnose_prediction_failures.R`
- `scripts/R/tests/test_psa_error_rate.R`
- `scripts/R/tests/test_sampling_convergence_rate.R`
- `scripts/R/tests/test_sex_variable_fix.R`
- `scripts/R/tests/tlr.R`

The full scenario cache produced by `12_scenario_EVPPIs.R` is the only supported scenario output. The former `scripts/R/tests/fit_independent_biomarker_models.R` was moved to `scripts/R/archive/fit_independent_biomarker_models.R` rather than deleted.

### Generating Reports

The project includes comprehensive Quarto reports in `scripts/QMD/report/` that generate publication-ready PDFs:

```bash
# Render individual reports (from project root)
quarto render scripts/QMD/report/para_models.qmd
quarto render scripts/QMD/report/input_parameters.qmd
quarto render scripts/QMD/report/clinical_effectiveness.qmd
quarto render scripts/QMD/report/CEA.qmd
quarto render scripts/QMD/report/OWSA.qmd
quarto render scripts/QMD/report/EVPPIs.qmd
quarto render scripts/QMD/report/scenario_effect.qmd
quarto render scripts/QMD/report/biomarker_decomposition.qmd
quarto render scripts/QMD/report/biomarker_distributions.qmd
```

Clinical/DAG/descriptive reports may be rendered to all declared formats. Economic reports use HTML-producing tables and must be rendered individually with `--to pdf` (for example, `quarto render scripts/QMD/report/CEA.qmd --to pdf`). Do not render the whole report directory in one command.

**Report Descriptions**:
1. **para_models.qmd** - Parametric survival model fits and diagnostics
2. **input_parameters.qmd** - Model input parameters summary
3. **clinical_effectiveness.qmd** - Survival outcomes and life-years gained
4. **CEA.qmd** - Cost-effectiveness analysis with ICERs
5. **OWSA.qmd** - One-way deterministic sensitivity analysis (tornado diagrams)
6. **EVPPIs.qmd** - Value of information analysis
7. **scenario_effect.qmd** - Scenario analysis results
8. **biomarker_decomposition.qmd** - Biomarker effect decomposition analysis
9. **biomarker_distributions.qmd** - Biomarker distribution and prevalence sensitivity

**Prerequisites**: Run the scripts required by each report; the full core analysis runs scripts 02-11 in sequence, while scenario outputs additionally require script 12.

## Key Features

### Core Modeling
- **Partitioned survival model** for cost-effectiveness analysis
- **Biomarker-guided treatment strategies** comparing two pre-immunotherapy biomarkers (week-4 CRP and baseline TMB/BRAF) against standard of care
- **Microsatellite-stable (MSS) colorectal cancer** focus
- **Parametric survival modeling** using multiple distributions (Weibull, exponential, gamma, etc.)

### Sensitivity and Uncertainty Analysis
- **Deterministic sensitivity analysis (DSA)** with tornado diagrams
- **Probabilistic sensitivity analysis (PSA)** with second-order Monte Carlo simulation
- **Expected Value of Perfect Partial Information (EVPPI)** analysis
- **Cost-effectiveness acceptability curves** and analysis
- **Scenario analysis framework** for alternative assumptions

### Advanced Analytical Features
- **Biomarker effect decomposition** - quantifying direct biomarker effects vs. treatment interactions
- **Biomarker distribution and prevalence sensitivity** analysis
- **Age effect analysis** - investigating age as prognostic factor
- **Correlated survival resampling** - maintaining PFS/OS correlation in PSA

### Quality Assurance and Reporting
- **Snapshot system** for bug fix impact assessment and model validation
- **Quarto-based report generation** - publication-ready PDFs with embedded R code
- **Comprehensive test suite** - validation, diagnostic, and convergence tests
- **Based on real-world clinical trial data** (METIMMOX trial, NCT03388190)

## Model Structure

The analysis implements a **partitioned survival model (PSM)** with three health states:
- **Progression-free (PF)**: Patients alive without disease progression
- **Progressed (P)**: Patients alive with disease progression
- **Dead (D)**: Absorbing state

State occupancy is derived from parametric survival curves:
- PF state = PFS curve
- P state = OS - PFS (bounded at 0)
- D state = 1 - OS

### Adverse-event scope

Adverse-event/toxicity costs and disutilities are not modeled separately. This is a deliberate scope choice based on the intended tolerability of the alternating short-course FLOX-nivolumab regimen and the lack of sufficiently robust treatment-specific trial data on adverse-event incidence, resource use, and utility decrements for economic parameterization.

### Cost and scope decisions

- **No post-progression treatment cost.** After progression the model charges only the quarterly follow-up contact and the one-time end-of-life cost. Second-line systemic therapy, post-progression imaging and post-progression visits are outside the modeled scope. Because the strategies differ in time spent in the progressed state, the omission is differential rather than a common offset. The parameter `c_other_pp` makes it explicit (base case 0) and is varied in a structural scenario at EUR 5,000 per quarter.
- **Second treatment sequence given to every progression-free patient.** Both arms receive a second eight-cycle sequence at modeled weeks 24-38, applied to everyone still progression-free at that point; METIMMOX re-treated on progression during the treatment break. A partitioned survival model has no on-treatment substate, so re-treatment cannot be triggered on progression without restructuring the model. A structural scenario removes the second sequence, with survival held at the trial estimate, giving a cost-side bound.
- **Unit prices are fixed in the probabilistic analysis.** Drug and test unit costs are published tariffs or assumed list prices, so they are varied deterministically only and carry no value of information. Visit, baseline, follow-up and end-of-life costs remain probabilistic because they bundle genuine resource-use uncertainty.
- **The nivolumab price is an assumption.** EUR 13,923 per administration is an assumed Norwegian hospital acquisition cost (roughly 40% below list), not a citable tariff: Norwegian hospital prices are set by confidential LIS tender. It is the largest single cost driver, bounded by the one-way analysis and the biosimilar scenario.
- **Health-state utilities are sampled jointly.** The PSA draws the progression-free utility and a non-negative decrement, deriving the progressed utility as the difference, so no draw places the progressed utility above the progression-free utility.

## Publication Figures

The repository includes Quarto vignettes (`scripts/QMD/vignettes/`) for generating publication-ready figures:

- **figure1.qmd** - Figure 1 (model structure/patient flow)
- **figure2.qmd** - Figure 2 (survival curves)
- **figure3.qmd** - Figure 3 (cost-effectiveness results)
- **figure4.qmd** - Figure 4 (sensitivity analysis)
- **suppl_figure_pfs_plots.qmd** - Supplemental PFS plots

Render individual figures or all at once:
```bash
quarto render scripts/QMD/vignettes/figure1.qmd
# Or render all figures
quarto render scripts/QMD/vignettes/
```

## Technical Documentation

The `scripts/QMD/technical_docs/` directory contains methodological documentation:

- **bug_fix_impact.qmd** - Template for documenting bug fix impacts on model results
- **age_effect_analysis.qmd** - Analysis of age as prognostic factor in survival models
- **all_parametric_survival_models.qmd** - Comparison of all candidate parametric model fits

These reports support model validation and document methodological decisions. For comprehensive guidance on using this codebase, see [CLAUDE.md](CLAUDE.md).

## Related Publications

This work is based on and related to the METIMMOX clinical trial (ClinicalTrials.gov Identifier: NCT03388190). Please cite the following publications:

1. **Primary Trial Results**:  
   Ree AH, Šaltytė Benth J, Hamre HM, et al. First-line oxaliplatin-based chemotherapy and nivolumab for metastatic microsatellite-stable colorectal cancer-the randomised METIMMOX trial. *Br J Cancer*. 2024;130(12):1921-1928. [doi:10.1038/s41416-024-02696-6](https://doi.org/10.1038/s41416-024-02696-6)

2. **Radiologic Response**:  
   Meltzer S, Negård A, Bakke KM, et al. Early radiologic signal of responsiveness to immune checkpoint blockade in microsatellite-stable/mismatch repair-proficient metastatic colorectal cancer. *Br J Cancer*. 2022;127(12):2227-2233. [doi:10.1038/s41416-022-02004-0](https://doi.org/10.1038/s41416-022-02004-0)

3. **CRP Predictive Value**:  
   Meltzer S, Berg JP, Hamre HM, et al. 632P Predictive value of C-reactive protein (CRP) in microsatellite-stable (MSS) metastatic colorectal cancer (mCRC) patients given first-line alternating short-course oxaliplatin-based chemotherapy (FLOX) and nivolumab. *Annals of Oncology*. 2023;34:S449. [doi:10.1016/j.annonc.2023.09.1822](https://doi.org/10.1016/j.annonc.2023.09.1822)

4. **TMB, BRAF, and CRP Biomarkers**:  
   Ree AH, Bousquet PA, Nilsen HL, et al. 543P Tumor mutational burden (TMB), BRAF status, and C-reactive protein (CRP) predict response to first-line alternating oxaliplatin-based chemotherapy and nivolumab in metastatic microsatellite-stable (MSS) colorectal cancer (CRC). *Annals of Oncology*. 2024;35:S453. [doi:10.1016/j.annonc.2024.08.612](https://doi.org/10.1016/j.annonc.2024.08.612)

## Citation

If you use this code or adapt it for your own research, please cite:

```bibtex
@software{metimmox_cea,
  author = {Geisler, Ben and [Co-authors]},
  title = {METIMMOX-1: Cost-Effectiveness Analysis of Biomarker-Guided Immunotherapy},
  year = {2025},
  url = {https://github.com/ben-geisler/METIMMOX-1}
}
```

*Note: A research paper describing this analysis is currently in preparation. This citation will be updated once published.*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

This repository contains research code for a specific clinical trial analysis. If you have suggestions or find issues, please open an issue in the repository.

## Contact

- **Repository Owner**: [ben-geisler](https://github.com/ben-geisler)
- **Clinical Trial**: METIMMOX (ClinicalTrials.gov: NCT03388190)

## Acknowledgments

This work is based on the METIMMOX clinical trial data and related research. We gratefully acknowledge:

- All investigators and clinical staff involved in the METIMMOX trial
- The patients who participated in the trial
- **Frederick Thielen** for contributions to the survival modeling code

---

**Disclaimer**: This repository contains research code for academic purposes. The clinical trial data is confidential and not included in this repository.
