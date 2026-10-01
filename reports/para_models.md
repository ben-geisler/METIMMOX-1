# Parametric Survival Models
Ben Geisler
2026-10-01

- [Overview](#overview)
- [Distribution Selection](#distribution-selection)
- [Coefficients](#coefficients)
- [Strategy Survival Curves](#strategy-survival-curves)
- [Summary](#summary)
  - [Validation warnings](#validation-warnings)

# Overview

This report documents the parametric survival models used in the
economic evaluation. A single joint model is fitted for both endpoints:

| Endpoint | Formula |
|:---|:---|
| Overall survival | Surv(OSwk, Death) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |
| Progression-free survival | Surv(PFSwk, Progression) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |

Economic survival formulas

# Distribution Selection

The OS and PFS distributions are selected **jointly** from the nine
candidate families: the selected pair is the one with the lowest
combined AIC among all pairs whose population-averaged curves satisfy OS
\>= PFS for control and every economic biomarker subgroup at every
weekly time point of the 10-year horizon. Under this rule the selected
pair is **OS gamma / PFS gamma** (combined AIC 1167.47), which is 4.43
AIC above the unconstrained minimum-AIC pair. The lowest marginal-AIC
distribution of an endpoint taken on its own is therefore not
necessarily the selected one; the “Selected” column marks the
distributions in use.

| Distribution |    AIC |    BIC | delta_AIC | delta_BIC | Selected |
|:-------------|-------:|-------:|----------:|----------:|:--------:|
| weibull      | 666.98 | 686.96 |      0.00 |      0.00 |          |
| weibullph    | 666.98 | 686.96 |      0.00 |      0.00 |          |
| gamma        | 667.31 | 687.28 |      0.32 |      0.32 |   Yes    |
| llogis       | 668.86 | 688.83 |      1.87 |      1.87 |          |
| gengamma     | 668.97 | 691.17 |      1.99 |      4.21 |          |
| lognormal    | 670.97 | 690.95 |      3.99 |      3.99 |          |
| genf         | 670.98 | 695.39 |      3.99 |      8.43 |          |
| gompertz     | 673.66 | 693.63 |      6.67 |      6.67 |          |
| exponential  | 686.14 | 703.89 |     19.15 |     16.93 |          |

Overall survival distribution comparison (ordered by marginal AIC)

| Distribution |    AIC |    BIC | delta_AIC | delta_BIC | Selected |
|:-------------|-------:|-------:|----------:|----------:|:--------:|
| gengamma     | 496.05 | 518.24 |      0.00 |      0.23 |          |
| lognormal    | 498.04 | 518.02 |      1.99 |      0.00 |          |
| genf         | 498.07 | 522.48 |      2.02 |      4.46 |          |
| gamma        | 500.16 | 520.13 |      4.11 |      2.11 |   Yes    |
| weibull      | 501.29 | 521.26 |      5.24 |      3.24 |          |
| weibullph    | 501.29 | 521.26 |      5.24 |      3.24 |          |
| llogis       | 501.56 | 521.53 |      5.51 |      3.52 |          |
| gompertz     | 505.31 | 525.29 |      9.26 |      7.27 |          |
| exponential  | 509.48 | 527.24 |     13.43 |      9.22 |          |

Progression-free survival distribution comparison (ordered by marginal
AIC)

2 candidate fits raised a warning while fitting, listed below. None of
them is a selected distribution, so the economic model and the PSA
covariance are unaffected; the warnings concern the reliability of the
standard errors of rejected candidates.

| Endpoint | Distribution | Warning |
|:---|:---|:---|
| PFS | gengamma | Hessian not positive definite: smallest eigenvalue is -1.0e+00 (threshold: -1.0e-05). |
| PFS | genf | Hessian not positive definite: smallest eigenvalue is -1.0e+00 (threshold: -1.0e-05). |

Warnings raised while fitting candidate distributions

| Rank | OS distribution | PFS distribution | Combined AIC | OS \>= PFS | Violating points | Selected |
|---:|:---|:---|---:|:--:|---:|:--:|
| 1 | weibull | gengamma | 1163.03 | No | 1614 |  |
| 2 | weibullph | gengamma | 1163.03 | No | 1614 |  |
| 3 | gamma | gengamma | 1163.36 | No | 1590 |  |
| 4 | llogis | gengamma | 1164.90 | No | 1480 |  |
| 5 | gengamma | gengamma | 1165.02 | No | 1614 |  |
| 6 | weibull | lognormal | 1165.03 | No | 761 |  |
| 7 | weibullph | lognormal | 1165.03 | No | 761 |  |
| 8 | weibull | genf | 1165.05 | No | 1590 |  |
| 9 | weibullph | genf | 1165.05 | No | 1590 |  |
| 10 | gamma | lognormal | 1165.35 | No | 507 |  |
| 11 | gamma | genf | 1165.37 | No | 1564 |  |
| 12 | llogis | lognormal | 1166.90 | No | 8 |  |
| 20 | gamma | gamma | 1167.47 | Yes | 0 | Yes |

Joint (OS, PFS) pair audit: top 12 of 81 pairs by combined AIC plus the
selected pair

*Note:* Violating points: number of (subgroup, week) points at which the
population-averaged PFS curve exceeds the OS curve over the model
horizon. The selected pair is the lowest combined-AIC pair with no
violation.

# Coefficients

Coefficient scales differ by row type. Rows named after a distribution
parameter (for example `shape` and `rate` for the gamma family, `shape`
and `scale` for Weibull PH) are the baseline parameters on their natural
scale. All covariate rows (Age, sex, Rx, biomarkers and their
interactions) are additive effects on the log of the distribution’s
location parameter (`log(rate)` for gamma, `log(scale)` for Weibull PH,
so that exp(coefficient) is a rate ratio for gamma and a hazard ratio
for Weibull PH). The `Scale` column states this for every row.

| Outcome | Distribution | Parameter | Scale | Estimate | SE | 95% CI |
|:---|:---|:---|:---|---:|---:|---:|
| OS | gamma | shape | baseline parameter (natural scale) | 2.4303 | 0.4236 | \[1.7270, 3.4198\] |
| OS | gamma | rate | baseline parameter (natural scale) | 0.0185 | 0.0123 | \[0.0051, 0.0678\] |
| OS | gamma | Age | covariate effect on log(rate) | 0.0036 | 0.0092 | \[-0.0145, 0.0217\] |
| OS | gamma | sex1 | covariate effect on log(rate) | 0.1796 | 0.1675 | \[-0.1486, 0.5079\] |
| OS | gamma | RxExperimental arm | covariate effect on log(rate) | 0.2215 | 0.2365 | \[-0.2420, 0.6850\] |
| OS | gamma | crp1 | covariate effect on log(rate) | -0.3896 | 0.3306 | \[-1.0375, 0.2583\] |
| OS | gamma | tmb_braf1 | covariate effect on log(rate) | -0.0663 | 0.2474 | \[-0.5513, 0.4187\] |
| OS | gamma | RxExperimental arm:crp1 | covariate effect on log(rate) | -0.1711 | 0.4137 | \[-0.9820, 0.6397\] |
| OS | gamma | RxExperimental arm:tmb_braf1 | covariate effect on log(rate) | -0.0982 | 0.3451 | \[-0.7746, 0.5782\] |
| PFS | gamma | shape | baseline parameter (natural scale) | 1.8980 | 0.3348 | \[1.3432, 2.6818\] |
| PFS | gamma | rate | baseline parameter (natural scale) | 0.0644 | 0.0511 | \[0.0136, 0.3045\] |
| PFS | gamma | Age | covariate effect on log(rate) | -0.0069 | 0.0112 | \[-0.0288, 0.0150\] |
| PFS | gamma | sex1 | covariate effect on log(rate) | 0.1142 | 0.2114 | \[-0.3002, 0.5287\] |
| PFS | gamma | RxExperimental arm | covariate effect on log(rate) | 0.4970 | 0.2963 | \[-0.0837, 1.0777\] |
| PFS | gamma | crp1 | covariate effect on log(rate) | 0.0190 | 0.4703 | \[-0.9028, 0.9408\] |
| PFS | gamma | tmb_braf1 | covariate effect on log(rate) | -0.2555 | 0.3378 | \[-0.9175, 0.4066\] |
| PFS | gamma | RxExperimental arm:crp1 | covariate effect on log(rate) | -0.9574 | 0.5642 | \[-2.0633, 0.1485\] |
| PFS | gamma | RxExperimental arm:tmb_braf1 | covariate effect on log(rate) | -0.3255 | 0.4553 | \[-1.2178, 0.5669\] |
| OS | weibullph | shape | baseline parameter (natural scale) | 1.7474 | 0.1883 | \[1.4147, 2.1583\] |
| OS | weibullph | scale | baseline parameter (natural scale) | 0.0001 | 0.0002 | \[0.0000, 0.0023\] |
| OS | weibullph | Age | covariate effect on log(scale) | 0.0075 | 0.0154 | \[-0.0227, 0.0378\] |
| OS | weibullph | sex1 | covariate effect on log(scale) | 0.3440 | 0.2744 | \[-0.1938, 0.8818\] |
| OS | weibullph | RxExperimental arm | covariate effect on log(scale) | 0.4496 | 0.3718 | \[-0.2791, 1.1783\] |
| OS | weibullph | crp1 | covariate effect on log(scale) | -0.8269 | 0.5524 | \[-1.9097, 0.2559\] |
| OS | weibullph | tmb_braf1 | covariate effect on log(scale) | -0.1029 | 0.4000 | \[-0.8868, 0.6810\] |
| OS | weibullph | RxExperimental arm:crp1 | covariate effect on log(scale) | -0.0989 | 0.6809 | \[-1.4333, 1.2356\] |
| OS | weibullph | RxExperimental arm:tmb_braf1 | covariate effect on log(scale) | -0.1275 | 0.5628 | \[-1.2305, 0.9756\] |
| PFS | weibullph | shape | baseline parameter (natural scale) | 1.4630 | 0.1602 | \[1.1804, 1.8132\] |
| PFS | weibullph | scale | baseline parameter (natural scale) | 0.0066 | 0.0082 | \[0.0006, 0.0758\] |
| PFS | weibullph | Age | covariate effect on log(scale) | -0.0106 | 0.0161 | \[-0.0422, 0.0210\] |
| PFS | weibullph | sex1 | covariate effect on log(scale) | 0.1245 | 0.3025 | \[-0.4683, 0.7173\] |
| PFS | weibullph | RxExperimental arm | covariate effect on log(scale) | 0.6989 | 0.4160 | \[-0.1164, 1.5141\] |
| PFS | weibullph | crp1 | covariate effect on log(scale) | 0.0415 | 0.6663 | \[-1.2645, 1.3475\] |
| PFS | weibullph | tmb_braf1 | covariate effect on log(scale) | -0.4120 | 0.4852 | \[-1.3629, 0.5390\] |
| PFS | weibullph | RxExperimental arm:crp1 | covariate effect on log(scale) | -1.2657 | 0.8155 | \[-2.8641, 0.3328\] |
| PFS | weibullph | RxExperimental arm:tmb_braf1 | covariate effect on log(scale) | -0.5361 | 0.6688 | \[-1.8468, 0.7747\] |

Selected (OS gamma, PFS gamma) and Weibull PH coefficient estimates

*Note:* Baseline parameters are on their natural scale; covariate rows
are additive effects on the log location parameter (gamma: log(rate),
exp(coef) = rate ratio; Weibull PH: log(scale), exp(coef) = hazard
ratio).

# Strategy Survival Curves

<img src="para_models_files/figure-commonmark/survival-curves-1.png"
style="width:100.0%" data-fig-align="center"
alt="Predicted survival curves for economic biomarker strategies" />

# Summary

The economic model uses one joint formula for both CRP-guided and
TMB/BRAF-guided strategies. The OS and PFS distributions are selected
jointly as the minimum combined-AIC pair that keeps OS \>= PFS for every
subgroup over the model horizon (currently OS gamma, PFS gamma), not by
marginal AIC per endpoint; Weibull PH is retained as a familiar
reference fit.

The comparisons in this report are dependent validation: the fits are
checked against the data they were fitted to. Model-versus-Kaplan-Meier
survival at 1, 2, 3 and 5 years, and the comparison of the extrapolated
standard-of-care curves with published first-line mCRC survival, are in
the technical report `survival_external_validation.qmd` (issue \#162).

------------------------------------------------------------------------

**Report completed on:** 2026-10-01  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.3

## Validation warnings

No warnings recorded during rendering.
