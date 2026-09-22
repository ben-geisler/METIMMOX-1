# Parametric Survival Model Specification
Ben Geisler
2026-09-22

- [Model specification](#model-specification)
- [Regression coefficients](#regression-coefficients)
  - [Fit statistics](#fit-statistics)
  - [Interpretation of the exponentiated
    coefficients](#interpretation-of-the-exponentiated-coefficients)
  - [Coefficient tables by
    distribution](#coefficient-tables-by-distribution)
    - [Gamma (selected for OS and PFS)](#gamma-selected-for-os-and-pfs)
    - [Weibull (AFT) (previous OS
      distribution)](#weibull-aft-previous-os-distribution)
    - [Weibull (PH)](#weibull-ph)
    - [Log-normal (previous PFS
      distribution)](#log-normal-previous-pfs-distribution)
    - [Exponential](#exponential)
    - [Log-logistic](#log-logistic)
    - [Gompertz](#gompertz)
    - [Generalised gamma](#generalised-gamma)
    - [Generalised F](#generalised-f)
- [Current versus previous distribution
  pair](#current-versus-previous-distribution-pair)
- [All candidate distributions](#all-candidate-distributions)
  - [Unlabelled curves](#unlabelled-curves)
  - [Labelled curves](#labelled-curves)
- [Notes](#notes)

# Model specification

The economic model uses **one joint parametric regression model per
endpoint**: one for overall survival (OS) and one for progression-free
survival (PFS). Both share the same linear predictor:

    OS : Surv(OSwk, Death) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx
    PFS: Surv(PFSwk, Progression) ~ Age + sex + Rx + crp * Rx + tmb_braf * Rx

The three strategies are **not three models**. They are three ways of
setting the treatment indicator `Rx` in the same fitted model before
predicting each patient’s survival curve and averaging over the patients
in the relevant subgroup (population averaging; see
`generate_population_averaged_predictions()`):

| Strategy | Rx set per patient | Non-zero treatment terms |
|:---|:---|:---|
| Standard of care strategy | Control for every patient | None |
| CRP-guided strategy | Experimental if CRP+, else control | Rx, Rx x CRP+ (and Rx x TMB/BRAF+ if TMB/BRAF+) |
| TMB/BRAF-guided strategy | Experimental if TMB/BRAF+, else control | Rx, Rx x TMB/BRAF+ (and Rx x CRP+ if CRP+) |

How the single joint model is used by each strategy

Because every patient enters the prediction with their observed age,
sex, and both biomarker values, the strategy only determines the value
of `Rx` (and therefore which interaction terms are non-zero) for each
patient; the age, sex, and biomarker main effects are always active.
Prevalence-weighted strategy curves are then formed from the
biomarker-positive and biomarker-negative subgroup curves.

Nine parametric families were fitted to each endpoint (Gamma, Weibull
(AFT), Weibull (PH), Log-normal, Exponential, Log-logistic, Gompertz,
Generalised gamma, Generalised F) on the complete-case data set (n =
68). The base case uses the pair with the smallest combined AIC among
all pairs that satisfy OS \>= PFS at every weekly time point for the
control group and all four biomarker subgroups: **Gamma for OS and Gamma
for PFS**.

# Regression coefficients

## Fit statistics

| Distribution      | OS AIC | PFS AIC | OS dAIC | PFS dAIC |
|:------------------|-------:|--------:|--------:|---------:|
| Gamma             | 667.31 |  636.86 |    0.32 |     1.90 |
| Weibull (AFT)     | 666.98 |  638.25 |    0.00 |     3.29 |
| Weibull (PH)      | 666.98 |  638.25 |    0.00 |     3.29 |
| Log-normal        | 670.97 |  634.96 |    3.99 |     0.00 |
| Exponential       | 686.14 |  645.75 |   19.15 |    10.80 |
| Log-logistic      | 668.86 |  637.39 |    1.87 |     2.44 |
| Gompertz          | 673.66 |  640.68 |    6.67 |     5.73 |
| Generalised gamma | 668.97 |  636.96 |    1.99 |     2.00 |
| Generalised F     | 670.98 |  638.96 |    3.99 |     4.00 |

AIC by distribution and endpoint (dAIC = difference from the
best-fitting distribution for that endpoint). Weibull (AFT) and Weibull
(PH) are reparametrisations of the same model and therefore share the
same AIC and fitted curves.

## Interpretation of the exponentiated coefficients

`flexsurv` places covariates on the location parameter of each
distribution. What the exponentiated coefficient (`exp(Est.)` in the
tables below) means therefore depends on the family (AFT = accelerated
failure time; PH = proportional hazards; TR = time ratio; RR = rate
ratio):

| Distribution      | Family | Location parameter | exp(coefficient) is a    |
|:------------------|:-------|:-------------------|:-------------------------|
| Gamma             | AFT    | rate               | Rate ratio (TR = 1 / RR) |
| Weibull (AFT)     | AFT    | scale              | Time ratio               |
| Weibull (PH)      | PH     | scale              | Hazard ratio             |
| Log-normal        | AFT    | meanlog            | Time ratio               |
| Exponential       | PH     | rate               | Hazard ratio             |
| Log-logistic      | AFT    | scale              | Time ratio               |
| Gompertz          | PH     | rate               | Hazard ratio             |
| Generalised gamma | AFT    | mu                 | Time ratio               |
| Generalised F     | AFT    | mu                 | Time ratio               |

Location parameter and interpretation of exp(coefficient) by
distribution

- **Time ratio** (AFT families): values above 1 lengthen survival times
  (Weibull AFT, log-logistic, log-normal, generalised gamma, generalised
  F).
- **Hazard ratio** (PH families): values above 1 increase the hazard
  (exponential, Weibull PH, Gompertz).
- **Rate ratio** (gamma): covariates act on the gamma rate parameter, so
  values above 1 *shorten* survival times; the equivalent time ratio is
  the reciprocal.

Baseline parameters (`shape`, `rate`, `scale`, `meanlog`, `sdlog`, `mu`,
`sigma`, `Q`, `P`) are shown on their natural scale and are not
exponentiated (blank `exp(Est.)` column). Covariate coefficients are
shown on the estimation scale of the location parameter. `Rx` is the
experimental arm (alternating FLOX + nivolumab) versus the control arm
(FLOX); `CRP+` and `TMB/BRAF+` are the biomarker-positive indicators;
`Rx x CRP+` and `Rx x TMB/BRAF+` are the treatment-by-biomarker
interactions that carry the biomarker-guided treatment benefit. Values
with an absolute magnitude of 10,000 or more are shown in scientific
notation; this only occurs for the generalised F fits, whose `P`
parameter is estimated at its lower boundary (P = 0, where the
generalised F reduces to the generalised gamma), so its confidence
interval is not meaningful.

## Coefficient tables by distribution

### Gamma (selected for OS and PFS)

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | shape          |    2.430 | 0.424 |  \[1.727, 3.420\] |           |
| OS       | rate           |    0.019 | 0.012 |  \[0.005, 0.068\] |           |
| OS       | Age            |    0.004 | 0.009 | \[-0.014, 0.022\] |     1.004 |
| OS       | Sex (male)     |    0.180 | 0.167 | \[-0.149, 0.508\] |     1.197 |
| OS       | Rx             |    0.222 | 0.236 | \[-0.242, 0.685\] |     1.248 |
| OS       | CRP+           |   -0.390 | 0.331 | \[-1.038, 0.258\] |     0.677 |
| OS       | TMB/BRAF+      |   -0.066 | 0.247 | \[-0.551, 0.419\] |     0.936 |
| OS       | Rx x CRP+      |   -0.171 | 0.414 | \[-0.982, 0.640\] |     0.843 |
| OS       | Rx x TMB/BRAF+ |   -0.098 | 0.345 | \[-0.775, 0.578\] |     0.906 |
| PFS      | shape          |    1.790 | 0.292 |  \[1.300, 2.466\] |           |
| PFS      | rate           |    0.047 | 0.034 |  \[0.011, 0.197\] |           |
| PFS      | Age            |   -0.004 | 0.011 | \[-0.025, 0.017\] |     0.996 |
| PFS      | Sex (male)     |    0.055 | 0.199 | \[-0.334, 0.445\] |     1.057 |
| PFS      | Rx             |    0.473 | 0.276 | \[-0.068, 1.014\] |     1.605 |
| PFS      | CRP+           |   -0.296 | 0.413 | \[-1.105, 0.512\] |     0.743 |
| PFS      | TMB/BRAF+      |   -0.189 | 0.296 | \[-0.770, 0.392\] |     0.828 |
| PFS      | Rx x CRP+      |   -0.608 | 0.508 | \[-1.603, 0.387\] |     0.545 |
| PFS      | Rx x TMB/BRAF+ |   -0.298 | 0.420 | \[-1.121, 0.525\] |     0.743 |

Gamma: OS and PFS regression coefficients; exp(Est.) is a rate ratio (tr
= 1 / rr)

### Weibull (AFT) (previous OS distribution)

| Endpoint | Parameter      | Estimate |     SE |              95% CI | exp(Est.) |
|:---------|:---------------|---------:|-------:|--------------------:|----------:|
| OS       | shape          |    1.747 |  0.188 |    \[1.415, 2.158\] |           |
| OS       | scale          |  158.034 | 94.104 | \[49.191, 507.707\] |           |
| OS       | Age            |   -0.004 |  0.009 |   \[-0.022, 0.013\] |     0.996 |
| OS       | Sex (male)     |   -0.197 |  0.156 |   \[-0.502, 0.108\] |     0.821 |
| OS       | Rx             |   -0.257 |  0.211 |   \[-0.672, 0.157\] |     0.773 |
| OS       | CRP+           |    0.473 |  0.313 |   \[-0.141, 1.087\] |     1.605 |
| OS       | TMB/BRAF+      |    0.059 |  0.229 |   \[-0.390, 0.508\] |     1.061 |
| OS       | Rx x CRP+      |    0.057 |  0.390 |   \[-0.707, 0.821\] |     1.058 |
| OS       | Rx x TMB/BRAF+ |    0.073 |  0.322 |   \[-0.559, 0.705\] |     1.076 |
| PFS      | shape          |    1.391 |  0.138 |    \[1.145, 1.691\] |           |
| PFS      | scale          |   40.319 | 28.758 |  \[9.963, 163.166\] |           |
| PFS      | Age            |    0.004 |  0.011 |   \[-0.017, 0.025\] |     1.004 |
| PFS      | Sex (male)     |   -0.021 |  0.198 |   \[-0.408, 0.366\] |     0.979 |
| PFS      | Rx             |   -0.431 |  0.264 |   \[-0.949, 0.087\] |     0.650 |
| PFS      | CRP+           |    0.286 |  0.402 |   \[-0.503, 1.075\] |     1.331 |
| PFS      | TMB/BRAF+      |    0.228 |  0.287 |   \[-0.335, 0.791\] |     1.256 |
| PFS      | Rx x CRP+      |    0.527 |  0.499 |   \[-0.451, 1.505\] |     1.694 |
| PFS      | Rx x TMB/BRAF+ |    0.294 |  0.419 |   \[-0.527, 1.114\] |     1.342 |

Weibull (AFT): OS and PFS regression coefficients; exp(Est.) is a time
ratio

### Weibull (PH)

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | shape          |    1.747 | 0.188 |  \[1.415, 2.158\] |           |
| OS       | scale          |    0.000 | 0.000 |  \[0.000, 0.002\] |           |
| OS       | Age            |    0.008 | 0.015 | \[-0.023, 0.038\] |     1.008 |
| OS       | Sex (male)     |    0.344 | 0.274 | \[-0.194, 0.882\] |     1.411 |
| OS       | Rx             |    0.450 | 0.372 | \[-0.279, 1.178\] |     1.568 |
| OS       | CRP+           |   -0.827 | 0.552 | \[-1.910, 0.256\] |     0.437 |
| OS       | TMB/BRAF+      |   -0.103 | 0.400 | \[-0.887, 0.681\] |     0.902 |
| OS       | Rx x CRP+      |   -0.099 | 0.681 | \[-1.433, 1.236\] |     0.906 |
| OS       | Rx x TMB/BRAF+ |   -0.127 | 0.563 | \[-1.231, 0.976\] |     0.880 |
| PFS      | shape          |    1.391 | 0.138 |  \[1.145, 1.691\] |           |
| PFS      | scale          |    0.006 | 0.007 |  \[0.001, 0.053\] |           |
| PFS      | Age            |   -0.006 | 0.015 | \[-0.035, 0.024\] |     0.994 |
| PFS      | Sex (male)     |    0.029 | 0.275 | \[-0.509, 0.568\] |     1.030 |
| PFS      | Rx             |    0.600 | 0.370 | \[-0.124, 1.325\] |     1.823 |
| PFS      | CRP+           |   -0.398 | 0.560 | \[-1.496, 0.700\] |     0.672 |
| PFS      | TMB/BRAF+      |   -0.317 | 0.401 | \[-1.104, 0.470\] |     0.728 |
| PFS      | Rx x CRP+      |   -0.733 | 0.694 | \[-2.093, 0.627\] |     0.480 |
| PFS      | Rx x TMB/BRAF+ |   -0.409 | 0.583 | \[-1.552, 0.735\] |     0.664 |

Weibull (PH): OS and PFS regression coefficients; exp(Est.) is a hazard
ratio

### Log-normal (previous PFS distribution)

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | meanlog        |    4.878 | 0.695 |  \[3.517, 6.239\] |           |
| OS       | sdlog          |    0.741 | 0.070 |  \[0.616, 0.892\] |           |
| OS       | Age            |   -0.008 | 0.010 | \[-0.027, 0.012\] |     0.992 |
| OS       | Sex (male)     |   -0.187 | 0.185 | \[-0.549, 0.174\] |     0.829 |
| OS       | Rx             |   -0.142 | 0.283 | \[-0.696, 0.412\] |     0.867 |
| OS       | CRP+           |    0.074 | 0.365 | \[-0.642, 0.789\] |     1.076 |
| OS       | TMB/BRAF+      |    0.122 | 0.280 | \[-0.427, 0.671\] |     1.130 |
| OS       | Rx x CRP+      |    0.576 | 0.460 | \[-0.325, 1.477\] |     1.779 |
| OS       | Rx x TMB/BRAF+ |    0.086 | 0.384 | \[-0.668, 0.839\] |     1.090 |
| PFS      | meanlog        |    3.845 | 0.770 |  \[2.336, 5.355\] |           |
| PFS      | sdlog          |    0.826 | 0.074 |  \[0.693, 0.985\] |           |
| PFS      | Age            |   -0.002 | 0.011 | \[-0.024, 0.019\] |     0.998 |
| PFS      | Sex (male)     |   -0.183 | 0.205 | \[-0.585, 0.220\] |     0.833 |
| PFS      | Rx             |   -0.587 | 0.318 | \[-1.210, 0.035\] |     0.556 |
| PFS      | CRP+           |    0.188 | 0.434 | \[-0.662, 1.037\] |     1.206 |
| PFS      | TMB/BRAF+      |    0.125 | 0.320 | \[-0.503, 0.753\] |     1.133 |
| PFS      | Rx x CRP+      |    1.023 | 0.536 | \[-0.027, 2.073\] |     2.782 |
| PFS      | Rx x TMB/BRAF+ |    0.221 | 0.433 | \[-0.628, 1.070\] |     1.247 |

Log-normal: OS and PFS regression coefficients; exp(Est.) is a time
ratio

### Exponential

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | rate           |    0.006 | 0.006 |  \[0.001, 0.042\] |           |
| OS       | Age            |    0.007 | 0.015 | \[-0.022, 0.036\] |     1.007 |
| OS       | Sex (male)     |    0.214 | 0.268 | \[-0.310, 0.739\] |     1.239 |
| OS       | Rx             |    0.287 | 0.370 | \[-0.438, 1.012\] |     1.332 |
| OS       | CRP+           |   -0.534 | 0.550 | \[-1.612, 0.544\] |     0.586 |
| OS       | TMB/BRAF+      |   -0.135 | 0.397 | \[-0.915, 0.644\] |     0.873 |
| OS       | Rx x CRP+      |   -0.132 | 0.679 | \[-1.463, 1.200\] |     0.876 |
| OS       | Rx x TMB/BRAF+ |   -0.121 | 0.554 | \[-1.207, 0.965\] |     0.886 |
| PFS      | rate           |    0.026 | 0.025 |  \[0.004, 0.171\] |           |
| PFS      | Age            |   -0.004 | 0.014 | \[-0.032, 0.025\] |     0.996 |
| PFS      | Sex (male)     |    0.072 | 0.267 | \[-0.452, 0.596\] |     1.074 |
| PFS      | Rx             |    0.471 | 0.367 | \[-0.248, 1.191\] |     1.602 |
| PFS      | CRP+           |   -0.324 | 0.548 | \[-1.399, 0.750\] |     0.723 |
| PFS      | TMB/BRAF+      |   -0.236 | 0.396 | \[-1.012, 0.540\] |     0.790 |
| PFS      | Rx x CRP+      |   -0.618 | 0.679 | \[-1.949, 0.713\] |     0.539 |
| PFS      | Rx x TMB/BRAF+ |   -0.298 | 0.566 | \[-1.408, 0.811\] |     0.742 |

Exponential: OS and PFS regression coefficients; exp(Est.) is a hazard
ratio

### Log-logistic

| Endpoint | Parameter      | Estimate |     SE |              95% CI | exp(Est.) |
|:---------|:---------------|---------:|-------:|--------------------:|----------:|
| OS       | shape          |    2.431 |  0.266 |    \[1.962, 3.012\] |           |
| OS       | scale          |   92.445 | 63.536 | \[24.036, 355.551\] |           |
| OS       | Age            |   -0.001 |  0.010 |   \[-0.021, 0.018\] |     0.999 |
| OS       | Sex (male)     |   -0.154 |  0.178 |   \[-0.502, 0.195\] |     0.857 |
| OS       | Rx             |   -0.162 |  0.274 |   \[-0.699, 0.374\] |     0.850 |
| OS       | CRP+           |    0.223 |  0.401 |   \[-0.563, 1.009\] |     1.250 |
| OS       | TMB/BRAF+      |   -0.031 |  0.291 |   \[-0.602, 0.540\] |     0.970 |
| OS       | Rx x CRP+      |    0.332 |  0.477 |   \[-0.602, 1.266\] |     1.393 |
| OS       | Rx x TMB/BRAF+ |    0.206 |  0.373 |   \[-0.525, 0.937\] |     1.228 |
| PFS      | shape          |    2.083 |  0.218 |    \[1.696, 2.558\] |           |
| PFS      | scale          |   49.980 | 40.511 | \[10.206, 244.757\] |           |
| PFS      | Age            |   -0.002 |  0.012 |   \[-0.025, 0.021\] |     0.998 |
| PFS      | Sex (male)     |   -0.206 |  0.211 |   \[-0.620, 0.208\] |     0.814 |
| PFS      | Rx             |   -0.730 |  0.333 |  \[-1.382, -0.078\] |     0.482 |
| PFS      | CRP+           |    0.321 |  0.440 |   \[-0.540, 1.183\] |     1.379 |
| PFS      | TMB/BRAF+      |    0.004 |  0.331 |   \[-0.645, 0.653\] |     1.004 |
| PFS      | Rx x CRP+      |    0.912 |  0.556 |   \[-0.179, 2.002\] |     2.488 |
| PFS      | Rx x TMB/BRAF+ |    0.404 |  0.448 |   \[-0.475, 1.282\] |     1.497 |

Log-logistic: OS and PFS regression coefficients; exp(Est.) is a time
ratio

### Gompertz

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | shape          |    0.010 | 0.002 |  \[0.005, 0.015\] |           |
| OS       | rate           |    0.002 | 0.002 |  \[0.000, 0.019\] |           |
| OS       | Age            |    0.010 | 0.016 | \[-0.021, 0.040\] |     1.010 |
| OS       | Sex (male)     |    0.344 | 0.275 | \[-0.196, 0.883\] |     1.410 |
| OS       | Rx             |    0.444 | 0.374 | \[-0.289, 1.176\] |     1.558 |
| OS       | CRP+           |   -0.971 | 0.570 | \[-2.087, 0.145\] |     0.379 |
| OS       | TMB/BRAF+      |   -0.136 | 0.401 | \[-0.921, 0.649\] |     0.872 |
| OS       | Rx x CRP+      |    0.054 | 0.686 | \[-1.292, 1.399\] |     1.055 |
| OS       | Rx x TMB/BRAF+ |   -0.119 | 0.563 | \[-1.223, 0.985\] |     0.888 |
| PFS      | shape          |    0.009 | 0.003 |  \[0.003, 0.016\] |           |
| PFS      | rate           |    0.018 | 0.018 |  \[0.003, 0.130\] |           |
| PFS      | Age            |   -0.004 | 0.015 | \[-0.033, 0.026\] |     0.996 |
| PFS      | Sex (male)     |    0.025 | 0.273 | \[-0.511, 0.561\] |     1.025 |
| PFS      | Rx             |    0.579 | 0.369 | \[-0.144, 1.302\] |     1.784 |
| PFS      | CRP+           |   -0.348 | 0.557 | \[-1.441, 0.745\] |     0.706 |
| PFS      | TMB/BRAF+      |   -0.371 | 0.405 | \[-1.164, 0.422\] |     0.690 |
| PFS      | Rx x CRP+      |   -0.778 | 0.700 | \[-2.150, 0.594\] |     0.459 |
| PFS      | Rx x TMB/BRAF+ |   -0.421 | 0.588 | \[-1.573, 0.731\] |     0.656 |

Gompertz: OS and PFS regression coefficients; exp(Est.) is a hazard
ratio

### Generalised gamma

| Endpoint | Parameter      | Estimate |    SE |            95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|------------------:|----------:|
| OS       | mu             |    4.955 | 0.627 |  \[3.725, 6.185\] |           |
| OS       | sigma          |    0.590 | 0.110 |  \[0.408, 0.851\] |           |
| OS       | Q              |    0.905 | 0.504 | \[-0.082, 1.892\] |           |
| OS       | Age            |   -0.003 | 0.009 | \[-0.021, 0.014\] |     0.997 |
| OS       | Sex (male)     |   -0.193 | 0.161 | \[-0.509, 0.124\] |     0.825 |
| OS       | Rx             |   -0.249 | 0.222 | \[-0.683, 0.186\] |     0.780 |
| OS       | CRP+           |    0.457 | 0.334 | \[-0.199, 1.112\] |     1.579 |
| OS       | TMB/BRAF+      |    0.059 | 0.234 | \[-0.399, 0.516\] |     1.060 |
| OS       | Rx x CRP+      |    0.077 | 0.423 | \[-0.753, 0.907\] |     1.080 |
| OS       | Rx x TMB/BRAF+ |    0.083 | 0.330 | \[-0.563, 0.730\] |     1.087 |
| PFS      | mu             |    3.841 | 0.818 |  \[2.238, 5.444\] |           |
| PFS      | sigma          |    0.826 | 0.076 |  \[0.690, 0.988\] |           |
| PFS      | Q              |    0.009 | 0.581 | \[-1.130, 1.148\] |           |
| PFS      | Age            |   -0.002 | 0.013 | \[-0.027, 0.023\] |     0.998 |
| PFS      | Sex (male)     |   -0.181 | 0.235 | \[-0.642, 0.280\] |     0.834 |
| PFS      | Rx             |   -0.586 | 0.326 | \[-1.225, 0.052\] |     0.556 |
| PFS      | CRP+           |    0.190 | 0.459 | \[-0.709, 1.089\] |     1.209 |
| PFS      | TMB/BRAF+      |    0.125 | 0.320 | \[-0.503, 0.753\] |     1.133 |
| PFS      | Rx x CRP+      |    1.017 | 0.658 | \[-0.273, 2.308\] |     2.765 |
| PFS      | Rx x TMB/BRAF+ |    0.223 | 0.450 | \[-0.659, 1.104\] |     1.250 |

Generalised gamma: OS and PFS regression coefficients; exp(Est.) is a
time ratio

### Generalised F

| Endpoint | Parameter      | Estimate |    SE |              95% CI | exp(Est.) |
|:---------|:---------------|---------:|------:|--------------------:|----------:|
| OS       | mu             |    4.959 | 0.628 |    \[3.728, 6.190\] |           |
| OS       | sigma          |    0.587 | 0.111 |    \[0.405, 0.852\] |           |
| OS       | Q              |    0.914 | 0.512 |   \[-0.089, 1.917\] |           |
| OS       | P              |    0.003 | 0.073 | \[0.000, 1.21e+22\] |           |
| OS       | Age            |   -0.003 | 0.009 |   \[-0.020, 0.014\] |     0.997 |
| OS       | Sex (male)     |   -0.194 | 0.161 |   \[-0.510, 0.122\] |     0.824 |
| OS       | Rx             |   -0.249 | 0.221 |   \[-0.682, 0.185\] |     0.780 |
| OS       | CRP+           |    0.457 | 0.334 |   \[-0.198, 1.111\] |     1.579 |
| OS       | TMB/BRAF+      |    0.059 | 0.233 |   \[-0.398, 0.515\] |     1.060 |
| OS       | Rx x CRP+      |    0.075 | 0.423 |   \[-0.753, 0.903\] |     1.078 |
| OS       | Rx x TMB/BRAF+ |    0.082 | 0.329 |   \[-0.563, 0.727\] |     1.086 |
| PFS      | mu             |    3.842 | 0.818 |    \[2.238, 5.445\] |           |
| PFS      | sigma          |    0.825 | 0.076 |    \[0.690, 0.988\] |           |
| PFS      | Q              |    0.008 | 0.582 |   \[-1.132, 1.148\] |           |
| PFS      | P              |    0.000 | 0.014 | \[0.000, 3.77e+30\] |           |
| PFS      | Age            |   -0.002 | 0.013 |   \[-0.027, 0.023\] |     0.998 |
| PFS      | Sex (male)     |   -0.181 | 0.235 |   \[-0.642, 0.280\] |     0.834 |
| PFS      | Rx             |   -0.586 | 0.326 |   \[-1.225, 0.052\] |     0.556 |
| PFS      | CRP+           |    0.190 | 0.459 |   \[-0.709, 1.089\] |     1.209 |
| PFS      | TMB/BRAF+      |    0.124 | 0.320 |   \[-0.503, 0.752\] |     1.132 |
| PFS      | Rx x CRP+      |    1.017 | 0.659 |   \[-0.274, 2.308\] |     2.766 |
| PFS      | Rx x TMB/BRAF+ |    0.224 | 0.450 |   \[-0.658, 1.105\] |     1.250 |

Generalised F: OS and PFS regression coefficients; exp(Est.) is a time
ratio

# Current versus previous distribution pair

The previous base case used the unconstrained minimum-combined-AIC pair,
**Weibull (AFT) for OS and Log-normal for PFS**. That pair fits the
observed data marginally better but produces population-averaged PFS
curves that exceed the OS curves in parts of the horizon, which is
impossible in a partitioned survival model (the progressed-state
occupancy OS - PFS would be negative). The current base case, **Gamma /
Gamma**, is the best-fitting pair that respects OS \>= PFS everywhere.

| Pair     | OS            | PFS        |     AIC | Ordered | Violations | Max gap |
|:---------|:--------------|:-----------|--------:|--------:|-----------:|--------:|
| Current  | Gamma         | Gamma      | 1304.16 |     Yes |   0 / 2605 |   0.000 |
| Previous | Weibull (AFT) | Log-normal | 1301.94 |      No | 865 / 2605 |   0.025 |

Current and previous distribution pairs. AIC is the combined OS + PFS
AIC (endpoint AICs are in the fit-statistics table); Ordered = OS \>=
PFS at every time point; Violations are counted over 521 weekly time
points x 5 subgroups; Max gap is the largest PFS - OS difference. The
ordering constraint costs 2.22 AIC points.

| Subgroup                     |   n | Violating weeks | Max gap |
|:-----------------------------|----:|----------------:|--------:|
| Standard of care (FLOX)      |  32 |              30 |   0.001 |
| CRP+ (FLOX + nivolumab)      |  17 |             272 |   0.025 |
| CRP- (FLOX)                  |  26 |             245 |   0.003 |
| TMB/BRAF+ (FLOX + nivolumab) |  16 |             267 |   0.014 |
| TMB/BRAF- (FLOX)             |  18 |              51 |   0.001 |

Where the previous Weibull (AFT) / Log-normal pair violates OS \>= PFS,
by subgroup (of 521 weekly time points each)

<img
src="survival_model_specification_files/figure-commonmark/pair-curves-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Population-averaged OS and PFS curves under the current and previous distribution pairs for each modeled subgroup, overlaid on the Kaplan-Meier estimates." />

<img
src="survival_model_specification_files/figure-commonmark/gap-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Difference between the PFS and OS curves under each distribution pair. The previous pair crosses above zero in every subgroup; the current pair does not." />

# All candidate distributions

Each figure shows one subgroup with OS (left) and PFS (right); every
candidate distribution is drawn together with the observed Kaplan-Meier
curve. The unlabelled versions are presented first so that visual fit
can be judged before the identity of each curve is revealed. Because
Weibull (AFT) and Weibull (PH) are the same model, they coincide exactly
and only eight distinct curves are visible in each panel.

## Unlabelled curves

<img
src="survival_model_specification_files/figure-commonmark/blind-soc-1.png"
style="width:100.0%" data-fig-align="center"
alt="Unlabelled candidate distributions, standard of care." />

<img
src="survival_model_specification_files/figure-commonmark/blind-crp-pos-1.png"
style="width:100.0%" data-fig-align="center"
alt="Unlabelled candidate distributions, CRP-positive patients on FLOX + nivolumab." />

<img
src="survival_model_specification_files/figure-commonmark/blind-crp-neg-1.png"
style="width:100.0%" data-fig-align="center"
alt="Unlabelled candidate distributions, CRP-negative patients on FLOX." />

<img
src="survival_model_specification_files/figure-commonmark/blind-tmb-pos-1.png"
style="width:100.0%" data-fig-align="center"
alt="Unlabelled candidate distributions, TMB/BRAF-positive patients on FLOX + nivolumab." />

<img
src="survival_model_specification_files/figure-commonmark/blind-tmb-neg-1.png"
style="width:100.0%" data-fig-align="center"
alt="Unlabelled candidate distributions, TMB/BRAF-negative patients on FLOX." />

## Labelled curves

<img
src="survival_model_specification_files/figure-commonmark/lab-soc-1.png"
style="width:100.0%" data-fig-align="center"
alt="Labelled candidate distributions, standard of care." />

<img
src="survival_model_specification_files/figure-commonmark/lab-crp-pos-1.png"
style="width:100.0%" data-fig-align="center"
alt="Labelled candidate distributions, CRP-positive patients on FLOX + nivolumab." />

<img
src="survival_model_specification_files/figure-commonmark/lab-crp-neg-1.png"
style="width:100.0%" data-fig-align="center"
alt="Labelled candidate distributions, CRP-negative patients on FLOX." />

<img
src="survival_model_specification_files/figure-commonmark/lab-tmb-pos-1.png"
style="width:100.0%" data-fig-align="center"
alt="Labelled candidate distributions, TMB/BRAF-positive patients on FLOX + nivolumab." />

<img
src="survival_model_specification_files/figure-commonmark/lab-tmb-neg-1.png"
style="width:100.0%" data-fig-align="center"
alt="Labelled candidate distributions, TMB/BRAF-negative patients on FLOX." />

# Notes

- All parametric curves are population-averaged predictions over the
  complete-case patients in each subgroup, using each patient’s observed
  age, sex, and biomarker values, with `Rx` set as the strategy
  dictates. They are the same curves that
  `04_parametric_survival_analysis.R` uses for the ordering-constrained
  selection; for the selected pair they are identical to the base-case
  curves in `l_params_base`.
- Kaplan-Meier curves are estimated on the same complete-case patients
  (control arm for standard of care and the biomarker-negative
  subgroups; experimental arm for the biomarker-positive subgroups).
- The previous pair referred to in Section 3 is the unconstrained
  minimum-combined-AIC pair, which was the base case before the ordering
  constraint was introduced.
- Coefficient estimates, standard errors, and confidence intervals are
  taken directly from `flexsurvreg` (`fit$res`).

------------------------------------------------------------------------

**Report completed on:** 2026-09-22  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.1
