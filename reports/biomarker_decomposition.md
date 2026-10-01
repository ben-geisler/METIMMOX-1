# Biomarker Effect Decomposition
Ben Geisler
2026-10-01

- [Overview](#overview)
- [Biomarker Prevalence](#biomarker-prevalence)
- [Restricted Mean Survival](#restricted-mean-survival)
- [Survival Decomposition](#survival-decomposition)
- [Interpretation](#interpretation)
  - [Validation warnings](#validation-warnings)

# Overview

This report decomposes the survival contribution of the two economic
biomarker strategies in the single joint model. Each strategy is shown
as biomarker-positive, biomarker-negative, prevalence-weighted, and
standard-of-care survival curves.

# Biomarker Prevalence

| Biomarker | Prevalence |
|:----------|-----------:|
| CRP       |      33.8% |
| TMB/BRAF  |      44.1% |

Economic biomarker prevalences

# Restricted Mean Survival

Restricted mean survival is the area under each survival curve over the
modelled horizon of 520 weeks (10 years), integrated with the
trapezoidal rule on the weekly model grid; it is therefore bounded by 10
years and is not an estimate of unrestricted mean survival.

| Biomarker | Endpoint                  | Curve    | Years |
|:----------|:--------------------------|:---------|------:|
| CRP       | Overall survival          | Positive |  2.81 |
| CRP       | Overall survival          | Negative |  1.87 |
| CRP       | Overall survival          | Weighted |  2.19 |
| CRP       | Overall survival          | Control  |  2.19 |
| CRP       | Progression-free survival | Positive |  1.81 |
| CRP       | Progression-free survival | Negative |  0.93 |
| CRP       | Progression-free survival | Weighted |  1.23 |
| CRP       | Progression-free survival | Control  |  0.94 |
| TMB/BRAF  | Overall survival          | Positive |  2.22 |
| TMB/BRAF  | Overall survival          | Negative |  2.09 |
| TMB/BRAF  | Overall survival          | Weighted |  2.15 |
| TMB/BRAF  | Overall survival          | Control  |  2.19 |
| TMB/BRAF  | Progression-free survival | Positive |  1.45 |
| TMB/BRAF  | Progression-free survival | Negative |  0.82 |
| TMB/BRAF  | Progression-free survival | Weighted |  1.10 |
| TMB/BRAF  | Progression-free survival | Control  |  0.94 |

Restricted mean survival (years) over the 520-week (10-year) horizon,
trapezoidal rule, by biomarker subgroup and strategy curve

# Survival Decomposition

<img
src="biomarker_decomposition_files/figure-commonmark/decomposition-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Survival decomposition by economic biomarker" />

# Interpretation

The prevalence-weighted strategy curve combines biomarker-positive
patients assigned to experimental treatment and biomarker-negative
patients assigned to standard care. The distance between the weighted
strategy curve and the control curve is the survival component of the
economic strategy effect.

------------------------------------------------------------------------

**Report completed on:** 2026-10-01  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.0

## Validation warnings

No warnings recorded during rendering.
