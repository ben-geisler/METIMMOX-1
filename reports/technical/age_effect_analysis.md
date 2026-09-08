# Age Effect Analysis
Ben Geisler
2026-09-08

- [Overview](#overview)
- [Age Distribution](#age-distribution)
- [Restricted Mean Survival
  Comparison](#restricted-mean-survival-comparison)
- [Interpretation](#interpretation)

# Overview

This technical note evaluates how population averaging over observed
ages affects survival predictions for the single joint economic model.
The comparison holds all patients at the mean observed age while
retaining their observed sex and biomarker values.

# Age Distribution

| Metric      | Value |
|:------------|------:|
| Patients    |    68 |
| Mean age    |  64.4 |
| Median age  |  65.0 |
| Minimum age |  38.0 |
| Maximum age |  80.0 |

Age distribution in the economic survival complete-case data

# Restricted Mean Survival Comparison

| Strategy | Endpoint | Base Case Years | Mean-Age Years | Difference |
|:---|:---|---:|---:|---:|
| Standard of Care | Overall survival | 2.187 | 2.185 | -0.001 |
| Standard of Care | Progression-free survival | 1.122 | 1.121 | -0.001 |
| CRP-guided | Overall survival | 2.187 | 2.186 | -0.001 |
| CRP-guided | Progression-free survival | 1.285 | 1.284 | -0.002 |
| TMB/BRAF-guided | Overall survival | 2.149 | 2.148 | -0.001 |
| TMB/BRAF-guided | Progression-free survival | 1.208 | 1.207 | -0.002 |

Restricted mean survival (years) over the 520-week (10-year) horizon,
trapezoidal rule, under observed-age and mean-age predictions

<img src="age_effect_analysis_files/figure-commonmark/age-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Observed-age and mean-age survival predictions" />

# Interpretation

The base case averages predictions over the observed age distribution.
The mean-age scenario is a diagnostic check that removes age
heterogeneity while keeping the same fitted survival model and strategy
definitions.

------------------------------------------------------------------------

**Report completed on:** 2026-09-08  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.0
