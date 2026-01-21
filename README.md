# METIMMOX-1

Cost-Effectiveness Analysis of Biomarker-Guided Immunotherapy in Metastatic MSS/pMMR Colorectal Cancer

## Overview

This repository contains the R code for a cost-effectiveness analysis comparing three biomarker strategies that guide the addition of immunotherapy (PD1/PDL1 inhibitor) to standard of care treatment for metastatic microsatellite-stable (MSS)/mismatch repair-proficient (pMMR) colorectal cancer patients receiving first-line treatment.

**METIMMOX** stands for: **Colorectal Cancer METastasis - Shaping Anti-tumor IMMunity by OXaliplatin**

### Biomarker Strategies Evaluated

1. **CRP Strategy**: C-reactive protein levels
2. **TLR Strategy**: Tumor lesion reduction
3. **TMB/BRAF Strategy**: Tumor mutation burden ≥9 mut/MB or presence of a BRAF mutation

Each strategy is compared against standard of care alone (platinum-based Nordic FLOX regimen without immunotherapy).

## Target Audience

This repository is designed for **health economists** and researchers developing decision-analytic models in R. The code provides a framework that can be adopted and adapted for similar cost-effectiveness analyses.

## Repository Structure

```
METIMMOX-1/
├── scripts/
│   ├── R/
│   │   ├── analysis/      # Numbered analysis scripts (01-15, main workflow)
│   │   ├── functions/     # Reusable model functions
│   │   ├── tests/         # Validation and testing scripts
│   │   └── archive/       # Deprecated code (for reference)
│   └── QMD/
│       ├── report/        # Publication-ready Quarto reports (PDF)
│       ├── vignettes/     # Figure generation scripts
│       └── technical_docs/# Technical documentation
├── data/                  # Data files (not included - confidential)
│   ├── tidy/             # Processed data and caches
│   └── output/           # Analysis outputs and snapshots
└── METIMMOX-1.Rproj      # RStudio project file
```

### Key Files

- **Analysis scripts** (`scripts/R/analysis/`): Numbered R scripts (01-15) containing the core decision-analytic model workflow
- **Functions** (`scripts/R/functions/`): Reusable functions including `model_fun.R` (main model), `calculate_outcomes.R`, and sensitivity analysis utilities
- **Quarto reports** (`scripts/QMD/report/`): Publication-ready PDF reports (CEA, clinical effectiveness, sensitivity analyses)
- **Tests** (`scripts/R/tests/`): Validation and testing scripts

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
               darthtools, dampack, mvtnorm, Matrix, here)
```

## Usage

Run the numbered analysis scripts in the `scripts/R/analysis/` folder in sequential order:

```r
# Core setup
source("scripts/R/analysis/01_data_prep.R")
source("scripts/R/analysis/02_setup_and_global_variables.R")
source("scripts/R/analysis/03_biomarker_strategies.R")

# Survival analysis
source("scripts/R/analysis/06_parametric_survival_analysis.R")
source("scripts/R/analysis/07_basecase_input_parameters.R")

# Survival resampling (first run generates cache - takes time)
source("scripts/R/analysis/08_sampling.R")

# Model execution
source("scripts/R/analysis/09_traces.R")
source("scripts/R/analysis/10_basecase_analysis.R")

# Sensitivity analyses
source("scripts/R/analysis/11_DSA.R")   # Deterministic sensitivity analysis
source("scripts/R/analysis/12_PSA.R")   # Probabilistic sensitivity analysis
source("scripts/R/analysis/13_EVPPIs.R") # Expected value of perfect partial information
```

**Notes**:
- Scripts 04-05 are supplementary RMarkdown files, not part of the main pipeline
- The clinical trial dataset is confidential and not included in this repository

## Key Features

- Decision-analytic model for cost-effectiveness analysis
- Biomarker-guided treatment strategy comparison
- Microsatellite-stable (MSS) colorectal cancer focus
- Parametric survival modeling using multiple distributions
- Deterministic sensitivity analysis (DSA)
- Probabilistic sensitivity analysis (PSA)
- Expected Value of Perfect Partial Information (EVPPI) analysis
- Cost-effectiveness acceptability analysis
- Quarto-based report generation (publication-ready PDFs)
- Based on real-world clinical trial data (METIMMOX trial, NCT03388190)

## Model Structure

The analysis implements a **partitioned survival model (PSM)** with three health states:
- **Progression-free (PF)**: Patients alive without disease progression
- **Progressed (P)**: Patients alive with disease progression
- **Dead (D)**: Absorbing state

State occupancy is derived from parametric survival curves:
- PF state = PFS curve
- P state = OS - PFS (bounded at 0)
- D state = 1 - OS

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
