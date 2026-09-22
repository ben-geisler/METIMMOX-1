# Clinical Effectiveness
Ben Geisler
2026-09-22

- [Overview](#overview)
- [Methodological Notes](#methodological-notes)
  - [Software Packages](#software-packages)
  - [Firth-Corrected Cox Regression](#firth-corrected-cox-regression)
  - [Ridge Regression Sensitivity
    Analysis](#ridge-regression-sensitivity-analysis)
- [Data Preparation](#data-preparation)
  - [Sample Characteristics](#sample-characteristics)
  - [Biomarker Correlations](#biomarker-correlations)
- [Proportional Hazards Assumption](#proportional-hazards-assumption)
  - [Schoenfeld Residual Tests](#schoenfeld-residual-tests)
  - [Log-Log Survival Plots](#log-log-survival-plots)
- [Primary Analysis: Firth-Corrected Cox
  Models](#primary-analysis-firth-corrected-cox-models)
  - [Standard Cox vs Firth](#standard-cox-vs-firth)
  - [Unified Model Results](#unified-model-results)
    - [Overall Survival](#overall-survival)
    - [Progression-Free Survival](#progression-free-survival)
  - [Forest Plot](#forest-plot)
- [Sensitivity Analysis: Ridge
  Regression](#sensitivity-analysis-ridge-regression)
- [Exploratory: TLR as a Predictive Biomarker (Responder
  Analysis)](#exploratory-tlr-as-a-predictive-biomarker-responder-analysis)
  - [Descriptive Context](#descriptive-context)
    - [TLR Prevalence by Treatment
      Arm](#tlr-prevalence-by-treatment-arm)
    - [Kaplan-Meier Curves Stratified by
      TLR](#kaplan-meier-curves-stratified-by-tlr)
  - [Standard Cox vs Firth](#standard-cox-vs-firth-1)
  - [Firth-Corrected Cox Model](#firth-corrected-cox-model)
    - [Overall Survival](#overall-survival-1)
    - [Progression-Free Survival](#progression-free-survival-1)
  - [Ridge Sensitivity Analysis (TLR Terms
    Penalized)](#ridge-sensitivity-analysis-tlr-terms-penalized)
- [Landmark Analysis: TLR at Week 9](#landmark-analysis-tlr-at-week-9)
  - [Landmark Cohort and Attrition](#landmark-cohort-and-attrition)
  - [Descriptive Context (Landmark
    Cohort)](#descriptive-context-landmark-cohort)
    - [Characteristics of the Week-9 Landmark Cohorts by TLR
      Status](#characteristics-of-the-week-9-landmark-cohorts-by-tlr-status)
    - [TLR Prevalence by Treatment Arm at Week
      9](#tlr-prevalence-by-treatment-arm-at-week-9)
    - [Landmark Kaplan-Meier Curves Stratified by
      TLR](#landmark-kaplan-meier-curves-stratified-by-tlr)
  - [Standard Cox vs Firth (Landmark)](#standard-cox-vs-firth-landmark)
  - [Firth-Corrected Cox Model
    (Landmark)](#firth-corrected-cox-model-landmark)
    - [Overall Survival](#overall-survival-2)
    - [Progression-Free Survival](#progression-free-survival-2)
  - [Ridge Sensitivity Analysis (Landmark, TLR Terms
    Penalized)](#ridge-sensitivity-analysis-landmark-tlr-terms-penalized)
    - [Landmark Sensitivity: Scan-Date and Week-12
      Landmarks](#landmark-sensitivity-scan-date-and-week-12-landmarks)
- [Discussion](#discussion)
  - [Primary interaction results](#primary-interaction-results)
  - [Proportional hazards](#proportional-hazards)
  - [Ridge regression sensitivity](#ridge-regression-sensitivity)
  - [TLR responder analysis
    (exploratory)](#tlr-responder-analysis-exploratory)
  - [TLR landmark analysis
    (exploratory)](#tlr-landmark-analysis-exploratory)
- [Summary](#summary)

# Overview

This report presents clinical effectiveness results from the METIMMOX-1
trial, evaluating biomarker-guided treatment strategies for selecting
immunotherapy candidates in metastatic MSS/pMMR colorectal cancer.

The analysis uses **Firth-corrected Cox proportional hazards
regression** with profile likelihood confidence intervals to assess
treatment effect heterogeneity across two biomarkers that are available
before the decision to add immunotherapy, identified by the causal DAG
analysis:

- **CRP** (C-reactive protein): Low CRP (\<5 mg/L) at week 4 (cycle 3
  day 1), measured after the two FLOX cycles that both arms receive and
  before the first nivolumab dose — a prognostic and potentially
  predictive marker
- **TMB/BRAF**: High tumor mutational burden (≥9 mut/MB) or BRAF
  mutation — a baseline genomic marker

**CRP timing.** CRP is not a baseline (pre-randomization) measurement:
the value used throughout is the week-4 CRP at the start of nivolumab,
matching the post-hoc CRP analysis of the trial paper. Because nivolumab
has not yet been given at week 4 and both arms receive identical FLOX
until then, the DAG carries no `T -> CRP` edge. Two consequences follow.
First, week-4 CRP-positivity differs by arm (experimental 17/36, control
7/35; Fisher p = 0.023) whereas baseline CRP was balanced (8/38 vs 7/35;
p = 1; `CRP0cat` in `01_data_prep.R`, which is not part of the analysis
data); this is treated as chance imbalance on a week-4 measurement in a
small trial, not as a failure of randomization, and is what the `CRP:Rx`
interaction is estimated against. Second, the 3 patients with no week-4
CRP are early deaths (weeks 2.4, 15.7, 20.9) and are excluded from the
complete-case data, so the CRP subgroups are conditional on surviving to
the week-4 measurement. The clock remains randomization; a landmark
analysis at week 4 would drop no further patients from the CRP
subgroups.

**Model specification is motivated by the causal DAG** (see `dag.qmd`
and `dag_associations.qmd`). The DAG encodes treatment (T) as randomized
with no parents, so no confounding adjustment is needed for causal
identification; Age and Sex are included as precision covariates. TLR
(Tumor Lesion Reduction) is encoded as a post-randomization intermediate
variable (T → TLR → PFS): it is not a pre-treatment baseline
characteristic and conditioning on it as a covariate would block part of
the treatment effect pathway. TLR is therefore excluded from the primary
survival models and examined separately as a descriptive endpoint.

The **unified DAG-informed model** is:

$$\text{Surv}(\text{time}, \text{event}) \sim \text{Age} + \text{sex} + \text{Rx} + \text{CRP} + \text{TMB/BRAF} + \text{CRP} \times \text{Rx} + \text{TMB/BRAF} \times \text{Rx}$$

Ridge regression (L2-penalized Cox) is applied as a sensitivity analysis
to assess the stability of the interaction estimates under shrinkage.
TLR is reported in a separate exploratory section.

# Methodological Notes

## Software Packages

This analysis uses the following R packages:

- `coxphf` version 1.13.4: Firth-corrected Cox regression with profile
  likelihood CIs
- `survival` version 3.7.0: Cox regression and Schoenfeld residual
  diagnostics
- `glmnet` version 4.1.8: Ridge regression sensitivity analysis
- `vcd` version 1.4.13: Cramér’s V for biomarker correlation assessment

## Firth-Corrected Cox Regression

We use Firth’s penalized maximum likelihood method for Cox regression,
which addresses:

1.  **Monotone likelihood**: When a covariate perfectly separates events
    from non-events
2.  **Small sample bias**: Reduces bias in coefficient estimates with
    sparse data
3.  **Unstable confidence intervals**: Profile likelihood CIs are more
    reliable than Wald-based CIs in small samples

Statistical inference is based on **penalized likelihood ratio tests
(PLRT)**, which are more appropriate than Wald tests for testing
treatment-biomarker interactions in this setting.

## Ridge Regression Sensitivity Analysis

We use Ridge regression as a sensitivity analysis to assess the
stability of our interaction estimates. Ridge regression (also called
L2-penalized regression or Tikhonov regularization) belongs to the
family of regularized or penalized regression methods. In survival
analysis, this is implemented as a penalized Cox model.

**How Ridge regression works:**

- Ridge adds a penalty term that discourages large coefficients. The
  bigger a coefficient tries to get, the more the model resists, pulling
  estimates toward zero (HR toward 1.0).
- Unlike LASSO (alpha = 1), which sets some coefficients exactly to zero
  and performs variable selection, Ridge (alpha = 0) shrinks all
  coefficients but retains all predictors.
- The degree of shrinkage is controlled by the tuning parameter lambda,
  which is selected via **cross-validation (CV)**: the data are
  repeatedly split into training and validation sets, and the lambda
  that produces the best predictive performance is chosen. Here this is
  `glmnet::cv.glmnet(family = "cox", alpha = 0)` with 10-fold CV of the
  partial-likelihood deviance and `lambda.min`; fold assignment is
  seeded so the result is reproducible.
- **Implementation detail.** The model is fitted on an explicit design
  matrix (Age, sex, Rx, CRP, TMB/BRAF, CRP x Rx, TMB/BRAF x Rx) in which
  the product terms are constructed before fitting. This matters because
  `survival::ridge()` with a factor treatment variable and a formula
  interaction such as `crp:Rx` codes the interaction separately for each
  arm (the CRP slope within the control arm and within the experimental
  arm) rather than as the single treatment-by-biomarker contrast; the
  pre-built product term recovers the same contrast that the Firth model
  estimates. In the primary model only the CRP and TMB/BRAF main effects
  are penalized (penalty factor 1); Age, sex, Rx and the two interaction
  terms have penalty factor 0. glmnet does not provide standard errors,
  so ridge results are point estimates.

**Why not LASSO or Elastic Net?**

We determined that LASSO and elastic net are **not appropriate** for
this analysis:

1.  **Correlated biomarkers**: LASSO arbitrarily selects one of
    correlated predictors while setting others to zero, which does not
    reflect the clinical question of comparing ALL biomarker strategies.

2.  **Hierarchy principle violation**: LASSO may select interaction
    terms (e.g., `crp:Rx`) while setting the main effect (`crp`) to
    zero, producing models that are difficult to interpret clinically.

3.  **Clinical objective**: Our goal is to estimate treatment effect
    heterogeneity across ALL biomarkers simultaneously, not to select a
    “winning” biomarker.

Elastic net (alpha = 0.1–0.3) represents a compromise but does not fully
address these concerns.

**Interpreting Ridge results:**

Ridge regression is not replacing the Firth estimates—it provides a
reality check. When Ridge substantially shrinks an interaction
coefficient toward zero while Firth estimates a large effect, this
indicates the effect is too unstable to trust for practical use. The
CV-selected lambda reflects how much regularization is needed for
reasonable out-of-sample prediction given the sample size and
collinearity structure.

# Data Preparation

## Sample Characteristics

| Characteristic | Overall | Positive | Negative | Positive | Negative | Positive | Negative |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| Patients, n | 68 | 23 | 45 | 41 | 24 | 30 | 38 |
| Age, mean (SD) | 64.4 (10.0) | 64.5 (10.1) | 64.3 (10.0) | 65.0 (9.9) | 62.2 (10.2) | 64.9 (10.4) | 63.9 (9.7) |
| Male | 36 (52.9%) | 11 (47.8%) | 25 (55.6%) | 22 (53.7%) | 13 (54.2%) | 15 (50.0%) | 21 (55.3%) |
| Female | 32 (47.1%) | 12 (52.2%) | 20 (44.4%) | 19 (46.3%) | 11 (45.8%) | 15 (50.0%) | 17 (44.7%) |
| Control (FLOX alone) | 32 (47.1%) | 6 (26.1%) | 26 (57.8%) | 22 (53.7%) | 7 (29.2%) | 14 (46.7%) | 18 (47.4%) |
| Experimental (FLOX/nivolumab) | 36 (52.9%) | 17 (73.9%) | 19 (42.2%) | 19 (46.3%) | 17 (70.8%) | 16 (53.3%) | 20 (52.6%) |
| Deaths (OS events) | 59 (86.8%) | 17 (73.9%) | 42 (93.3%) | 32 (78.0%) | 24 (100.0%) | 23 (76.7%) | 36 (94.7%) |
| PFS events (progression or death) | 63 (92.6%) | 19 (82.6%) | 44 (97.8%) | 36 (87.8%) | 24 (100.0%) | 25 (83.3%) | 38 (100.0%) |
| CRP-positive (low CRP) | 23 (33.8%) | 23 (100.0%) | 0 (0.0%) | 18 (43.9%) | 4 (16.7%) | 11 (36.7%) | 12 (31.6%) |
| TLR-positive (early response) | 41 (63.1%) | 18 (81.8%) | 23 (53.5%) | 41 (100.0%) | 0 (0.0%) | 18 (62.1%) | 23 (63.9%) |
| TMB/BRAF-positive | 30 (44.1%) | 11 (47.8%) | 19 (42.2%) | 18 (43.9%) | 11 (45.8%) | 30 (100.0%) | 0 (0.0%) |

Sample Characteristics by Biomarker Status

*Note:* Percentages within each column use the column subgroup as
denominator. CRP/TLR/TMB-BRAF rows in the corresponding subgroup column
equal 100% by definition. Overall, CRP and TMB/BRAF columns are the
primary complete-case cohort (n = 68); the TLR columns and the
TLR-positive row use the 65 patients with a first on-treatment CT scan
(3 early deaths have no TLR value). Sex: level 0 of the trial variable
is female.

## Biomarker Correlations

| Biomarker Pair  | Cram\<U+00E9\>r’s V | Interpretation |
|:----------------|--------------------:|:---------------|
| CRP vs TLR      |               0.278 | Weak           |
| CRP vs TMB/BRAF |               0.053 | Negligible     |
| TLR vs TMB/BRAF |               0.019 | Negligible     |

Biomarker Correlations (Cram\<U+00E9\>r’s V)

**Note:** Cramér’s V interpretation: \<0.1 = Negligible, 0.1–0.3 = Weak,
\>0.3 = Moderate/Strong. The observed correlations are weak to
negligible. While modest, any correlation among candidate biomarkers
supports our decision not to use LASSO, which would arbitrarily select
among correlated predictors.



# Proportional Hazards Assumption

## Schoenfeld Residual Tests

The proportional hazards assumption was tested using scaled Schoenfeld
residuals. A significant p-value indicates potential violation of the PH
assumption for that covariate.

| Covariate | \<U+03C7\>\<U+00B2\> (OS) | p-value | \<U+03C7\>\<U+00B2\> (PFS) | p-value |
|:---|---:|:---|---:|:---|
| Age | 0.68 | 0.411 | 0.07 | 0.796 |
| sex | 0.12 | 0.729 | 0.16 | 0.689 |
| Rx | 4.83 | 0.028 | 0.01 | 0.926 |
| crp_num | 0.38 | 0.536 | 2.29 | 0.130 |
| tmb_braf_num | 0.94 | 0.332 | 2.98 | 0.084 |
| Rx:crp_num | 3.08 | 0.079 | 1.20 | 0.272 |
| Rx:tmb_braf_num | 0.97 | 0.324 | 1.62 | 0.203 |
| GLOBAL | 8.98 | 0.254 | 11.48 | 0.119 |

Schoenfeld Residual Test for Proportional Hazards Assumption

**Note:** GLOBAL is the omnibus test for all covariates. Individual
covariate tests help identify specific PH violations.

## Log-Log Survival Plots

<img
src="clinical_effectiveness_files/figure-commonmark/loglog-plots-1.png"
style="width:100.0%" data-fig-align="center" />

**Interpretation:** Approximately parallel lines on the log-log plot
support the proportional hazards assumption. Crossing or converging
lines suggest time-varying effects.



# Primary Analysis: Firth-Corrected Cox Models

## Standard Cox vs Firth

The standard Cox model is the unpenalized maximum-likelihood reference.
Divergence between standard Cox and Firth estimates flags small-sample
or near-separation instability that Firth corrects.

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| CRP x Rx          | 0.71 (0.18-2.72)         | 0.65 (0.19-2.60)  |
| TMB/BRAF x Rx     | 0.97 (0.32-2.92)         | 0.95 (0.32-2.82)  |
| Rx (Experimental) | 1.51 (0.73-3.15)         | 1.53 (0.74-3.18)  |
| CRP               | 0.57 (0.19-1.68)         | 0.63 (0.19-1.62)  |
| TMB/BRAF          | 0.89 (0.40-1.94)         | 0.91 (0.41-1.96)  |
| Age               | 1.01 (0.98-1.04)         | 1.01 (0.98-1.04)  |
| Sex               | 1.44 (0.84-2.49)         | 1.44 (0.84-2.48)  |

Unified Model - Overall Survival: Standard Cox vs Firth

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| CRP x Rx          | 0.44 (0.11-1.74)         | 0.42 (0.11-1.68)  |
| TMB/BRAF x Rx     | 0.66 (0.21-2.12)         | 0.66 (0.21-2.05)  |
| Rx (Experimental) | 1.79 (0.87-3.72)         | 1.79 (0.88-3.73)  |
| CRP               | 0.75 (0.25-2.26)         | 0.81 (0.26-2.21)  |
| TMB/BRAF          | 0.77 (0.35-1.71)         | 0.79 (0.36-1.72)  |
| Age               | 1.00 (0.97-1.03)         | 1.00 (0.97-1.03)  |
| Sex               | 1.07 (0.62-1.84)         | 1.06 (0.62-1.84)  |

Unified Model - Progression-Free Survival: Standard Cox vs Firth

## Unified Model Results

The unified DAG-informed model includes CRP (week 4, before the first
nivolumab dose) and TMB/BRAF (baseline) as pre-immunotherapy biomarkers
with their treatment interactions, adjusted for Age and Sex.

**Formula:**
`Surv(time, event) ~ Age + sex + Rx + CRP + TMB/BRAF + CRP:Rx + TMB/BRAF:Rx`

### Overall Survival

| Term                            | HR (95% CI)      | PLRT p-value |
|:--------------------------------|:-----------------|:-------------|
| RxExperimental arm:crp_num      | 0.65 (0.19-2.60) | 0.520        |
| RxExperimental arm:tmb_braf_num | 0.95 (0.32-2.82) | 0.919        |
| Age                             | 1.01 (0.98-1.04) | 0.691        |
| sex1                            | 1.44 (0.84-2.48) | 0.187        |
| RxExperimental arm              | 1.53 (0.74-3.18) | 0.250        |
| crp_num                         | 0.63 (0.19-1.62) | 0.355        |
| tmb_braf_num                    | 0.91 (0.41-1.96) | 0.815        |

Unified Model: Overall Survival - Full Coefficient Table

### Progression-Free Survival

| Term                            | HR (95% CI)      | PLRT p-value |
|:--------------------------------|:-----------------|:-------------|
| RxExperimental arm:crp_num      | 0.42 (0.11-1.68) | 0.213        |
| RxExperimental arm:tmb_braf_num | 0.66 (0.21-2.05) | 0.469        |
| Age                             | 1.00 (0.97-1.03) | 0.826        |
| sex1                            | 1.06 (0.62-1.84) | 0.821        |
| RxExperimental arm              | 1.79 (0.88-3.73) | 0.108        |
| crp_num                         | 0.81 (0.26-2.21) | 0.687        |
| tmb_braf_num                    | 0.79 (0.36-1.72) | 0.550        |

Unified Model: Progression-Free Survival - Full Coefficient Table



## Forest Plot

<img
src="clinical_effectiveness_files/figure-commonmark/forest-plot-1.png"
style="width:100.0%" data-fig-align="center" />

**Interpretation:** HR \< 1 indicates that biomarker-positive patients
derive greater benefit from experimental treatment (reduced hazard of
death/progression). The dashed line at HR = 1 represents no differential
treatment effect.



# Sensitivity Analysis: Ridge Regression

Ridge regression (L2-penalized Cox, alpha = 0) shrinks coefficients
toward zero but does not eliminate any, providing a sensitivity check
for the stability of the interaction estimates. The same unified model
terms are used, entered through an explicit design matrix in which the
CRP x Rx and TMB/BRAF x Rx product terms are built before fitting, so
each ridge interaction coefficient is the single treatment-by-biomarker
contrast that the Firth model estimates (issue \#153). Only the
biomarker main effects (CRP and TMB/BRAF) are penalized; Age, sex, Rx
and the two interaction terms carry a penalty factor of zero. The
penalty lambda is chosen by 10-fold cross-validation of the
partial-likelihood deviance (`glmnet::cv.glmnet`, `lambda.min`), with
fold assignment seeded for reproducibility.

| Interaction   | Firth (95% CI)   | Ridge (point est.) |
|:--------------|:-----------------|:-------------------|
| CRP x Rx      | 0.65 (0.19-2.60) | 0.40               |
| TMB/BRAF x Rx | 0.95 (0.32-2.82) | 0.86               |

Overall Survival: Firth vs Ridge Regression

*Note:* Ridge: glmnet L2 penalty (alpha = 0) on the CRP and TMB/BRAF
main effects only; interaction product terms built before fitting;
lambda by 10-fold CV (lambda.min = 86.857). glmnet gives no standard
errors. A substantial shift vs Firth = estimate sensitive to
regularization.

| Interaction   | Firth (95% CI)   | Ridge (point est.) |
|:--------------|:-----------------|:-------------------|
| CRP x Rx      | 0.42 (0.11-1.68) | 0.33               |
| TMB/BRAF x Rx | 0.66 (0.21-2.05) | 0.51               |

Progression-Free Survival: Firth vs Ridge Regression

*Note:* Ridge: glmnet L2 penalty (alpha = 0) on the CRP and TMB/BRAF
main effects only; interaction product terms built before fitting;
lambda by 10-fold CV (lambda.min = 109.761). glmnet gives no standard
errors.

**Interpretation:** Firth provides the best unbiased point estimate
given the sample size and data structure. Ridge applies L2 shrinkage to
the biomarker main effects (HR toward 1.0) with the penalty strength
chosen by cross-validation; the interaction terms are unpenalized, so
any movement in them reflects redistribution of the shrunken main-effect
signal. Comparing the two methods reveals estimate stability: if Firth
and Ridge agree closely, the estimate is robust; if the interaction
moves substantially under ridge, the effect is sensitive to how the
biomarker main effects are estimated in this sample.



# Exploratory: TLR as a Predictive Biomarker (Responder Analysis)

This section presents an **exploratory responder analysis** that
includes TLR (Tumor Lesion Reduction) as an effect modifier in a Cox
model with the same structural form as the primary analysis. Because TLR
is measured post-randomization (at the first on-treatment CT scan), this
is not a treatment-selection analysis: TLR status is itself influenced
by treatment, and the resulting TLR x Rx interaction is informative only
as a hypothesis-generating, responder-stratified contrast. It cannot be
used to guide the decision to add immunotherapy, which is made before
TLR is observed. The primary causal-DAG-informed analysis remains the
inferential anchor for predictive biomarker claims.

The model adjusts for the pre-immunotherapy biomarkers **CRP** and
**TMB/BRAF** as prognostic main effects (their treatment interactions
belong to the primary analysis and are not re-estimated here). In the
ridge sensitivity model, CRP and TMB/BRAF remain **unpenalized**
(well-established prognostic adjustment covariates), while the **TLR
main effect and TLR x Rx interaction are penalized**. These are the
post-randomization terms whose stability we want to stress-test against
the small effective sample size and the partially tautological TLR-PFS
association.

**Formula:**
`Surv(time, event) ~ Age + sex + Rx + CRP + TMB/BRAF + TLR + TLR:Rx`

## Descriptive Context

### TLR Prevalence by Treatment Arm

| Treatment Arm            |   N | TLR-positive |    % |
|:-------------------------|----:|-------------:|-----:|
| Control (FLOX)           |  29 |           22 | 75.9 |
| Experimental (FLOX/nivo) |  36 |           19 | 52.8 |

TLR Prevalence by Treatment Arm (Descriptive)

*Note:* TLR = Tumour Lesion Reduction \>= 10% at first CT scan
(post-treatment).

### Kaplan-Meier Curves Stratified by TLR

<img src="clinical_effectiveness_files/figure-commonmark/km-tlr-1.png"
style="width:100.0%" data-fig-align="center" />

**Note:** These curves are purely descriptive and stratify by a
post-treatment response variable. The observed difference in survival by
TLR status reflects treatment response, not a pre-treatment patient
characteristic.



## Standard Cox vs Firth

The standard Cox model is the unpenalized maximum-likelihood reference.
Divergence between standard Cox and Firth estimates flags small-sample
or near-separation instability that Firth corrects.

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| TLR x Rx          | 2.43 (0.74-7.95)         | 2.47 (0.75-7.75)  |
| TLR (main effect) | 0.21 (0.08-0.54)         | 0.21 (0.08-0.55)  |
| Rx (Experimental) | 0.60 (0.23-1.55)         | 0.59 (0.24-1.55)  |
| CRP               | 0.52 (0.25-1.06)         | 0.53 (0.25-1.05)  |
| TMB/BRAF          | 0.86 (0.48-1.52)         | 0.86 (0.48-1.52)  |
| Age               | 1.01 (0.98-1.04)         | 1.01 (0.98-1.04)  |
| Sex               | 1.48 (0.81-2.68)         | 1.47 (0.82-2.67)  |

TLR Responder Model - Overall Survival: Standard Cox vs Firth

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| TLR x Rx          | 1.47 (0.45-4.81)         | 1.49 (0.45-4.67)  |
| TLR (main effect) | 0.24 (0.09-0.61)         | 0.24 (0.10-0.62)  |
| Rx (Experimental) | 0.75 (0.29-1.96)         | 0.73 (0.30-1.97)  |
| CRP               | 0.46 (0.23-0.92)         | 0.47 (0.23-0.92)  |
| TMB/BRAF          | 0.55 (0.30-1.00)         | 0.56 (0.30-0.99)  |
| Age               | 0.99 (0.96-1.02)         | 0.99 (0.96-1.02)  |
| Sex               | 1.11 (0.62-2.01)         | 1.11 (0.62-2.01)  |

TLR Responder Model - Progression-Free Survival: Standard Cox vs Firth

## Firth-Corrected Cox Model

### Overall Survival

| Term                       | HR (95% CI)      | PLRT p-value |
|:---------------------------|:-----------------|:-------------|
| RxExperimental arm:tlr_num | 2.47 (0.75-7.75) | 0.135        |
| Age                        | 1.01 (0.98-1.04) | 0.626        |
| sex1                       | 1.47 (0.82-2.67) | 0.201        |
| RxExperimental arm         | 0.59 (0.24-1.55) | 0.270        |
| crp_num                    | 0.53 (0.25-1.05) | 0.071        |
| tmb_braf_num               | 0.86 (0.48-1.52) | 0.609        |
| tlr_num                    | 0.21 (0.08-0.55) | 0.002        |

TLR Responder Model: Overall Survival - Full Coefficient Table

### Progression-Free Survival

| Term                       | HR (95% CI)      | PLRT p-value |
|:---------------------------|:-----------------|:-------------|
| RxExperimental arm:tlr_num | 1.49 (0.45-4.67) | 0.504        |
| Age                        | 0.99 (0.96-1.02) | 0.357        |
| sex1                       | 1.11 (0.62-2.01) | 0.728        |
| RxExperimental arm         | 0.73 (0.30-1.97) | 0.522        |
| crp_num                    | 0.47 (0.23-0.92) | 0.027        |
| tmb_braf_num               | 0.56 (0.30-0.99) | 0.047        |
| tlr_num                    | 0.24 (0.10-0.62) | 0.005        |

TLR Responder Model: Progression-Free Survival - Full Coefficient Table



## Ridge Sensitivity Analysis (TLR Terms Penalized)

Ridge regression is applied here with a targeted penalty: only the TLR
main effect and the TLR x Rx interaction are L2-shrunk. CRP, TMB/BRAF,
Age, sex, and Rx are estimated without penalty, since they enter the
model as adjustment covariates whose effects are not the inferential
target.

| Term              | Firth (95% CI)   | Ridge HR (SE)   |
|:------------------|:-----------------|:----------------|
| TLR x Rx          | 2.47 (0.75-7.75) | 2.12 (SE: 0.58) |
| TLR (main effect) | 0.21 (0.08-0.55) | 0.24 (SE: 0.47) |
| Age               | 1.01 (0.98-1.04) | 1.01 (SE: 0.02) |
| Sex               | 1.47 (0.82-2.67) | 1.49 (SE: 0.30) |
| Rx (Experimental) | 2.47 (0.75-7.75) | 0.66 (SE: 0.48) |
| CRP               | 0.53 (0.25-1.05) | 0.52 (SE: 0.36) |
| TMB/BRAF          | 0.86 (0.48-1.52) | 0.86 (SE: 0.29) |

TLR Responder Model - Overall Survival: Firth vs Ridge, all coefficients

*Note:* Ridge applies L2 shrinkage to the TLR main effect and the TLR x
Rx interaction only; Age, sex, Rx, CRP, and TMB/BRAF are unpenalized.
survival::ridge() does not produce confidence intervals; standard errors
are shown for reference.

| Term              | Firth (95% CI)   | Ridge HR (SE)   |
|:------------------|:-----------------|:----------------|
| TLR x Rx          | 1.49 (0.45-4.67) | 1.33 (SE: 0.58) |
| TLR (main effect) | 0.24 (0.10-0.62) | 0.26 (SE: 0.47) |
| Age               | 0.99 (0.96-1.02) | 0.99 (SE: 0.02) |
| Sex               | 1.11 (0.62-2.01) | 1.12 (SE: 0.30) |
| Rx (Experimental) | 1.49 (0.45-4.67) | 0.80 (SE: 0.48) |
| CRP               | 0.47 (0.23-0.92) | 0.46 (SE: 0.36) |
| TMB/BRAF          | 0.56 (0.30-0.99) | 0.55 (SE: 0.30) |

TLR Responder Model - Progression-Free Survival: Firth vs Ridge, all
coefficients

*Note:* Ridge applies L2 shrinkage to the TLR main effect and the TLR x
Rx interaction only; Age, sex, Rx, CRP, and TMB/BRAF are unpenalized.
survival::ridge() does not produce confidence intervals; standard errors
are shown for reference.

**Interpretation:** The Firth HR is the best unbiased point estimate of
the TLR x Rx interaction under penalized maximum likelihood. The ridge
HR applies L2 shrinkage targeted at the TLR terms; if the ridge estimate
is substantially closer to HR = 1.0 than the Firth estimate, the
responder-stratified interaction is sensitive to regularization,
consistent with a small-sample, partly tautological signal. Findings
here are exploratory and should not be used for treatment selection,
since TLR is not measurable before nivolumab starts. They may motivate
landmark or formal causal-mediation analyses in future work.



# Landmark Analysis: TLR at Week 9

The responder analysis above treats TLR status as if it were known at
baseline. It is not: TLR is measured at the first on-treatment CT scan,
so a patient can only be classified TLR-positive or TLR-negative *after*
surviving (and remaining progression-free) long enough to reach that
scan. Conditioning survival on a post-baseline response therefore
induces **guarantee-time (immortal-time) bias** — TLR-classified
patients are, by construction, a more favourable risk set than the full
randomized population.

A **landmark analysis** removes this bias by (i) fixing a landmark time,
(ii) restricting to patients still at risk after the landmark, and (iii)
measuring survival **from the landmark onward**. The primary landmark is
**week 9**, the protocol time of the second CT (the first on-treatment
response assessment): TLR compares this scan to the baseline CT, so a
patient’s TLR status does not exist until that scan. The actual first
on-treatment scans, however, were not all at week 9: in the trial
extract they fell between weeks 6.7 and 12.1 after inclusion (median
8.6). A fixed week-9 landmark therefore does not guarantee that every
retained patient had reached the scan, and it does not remove every
patient whose progression was recorded at the first scan. Both
departures are counted below, and two sensitivity landmarks are
reported: a **per-patient landmark at each patient’s own first-scan
date**, which by construction retains only patients progression-free
after their TLR was read, and a **fixed week-12 landmark**, after the
latest first scan (issue \#155). This section **complements** the
responder analysis above (which is retained for comparison); it does not
replace it.

The three landmark definitions live in `R/tlr_landmark.R` and are shared
with Figure 1 of the clinical effectiveness paper and the DAG
association tests, so all three use one cohort rule: a patient enters an
endpoint-specific cohort only if the endpoint time is strictly after the
landmark. Because the small trial loses patients with events or
censoring before the landmark, the cohort sizes and event counts are
reported explicitly, and all model fits are wrapped so a too-small
cohort degrades gracefully rather than aborting the render.

## Landmark Cohort and Attrition

| Endpoint | N in TLR-complete cohort | N at risk at week 9 | N excluded (event/censor \< week 9) | Events from week 9 |
|:---|---:|---:|---:|---:|
| Overall survival | 65 | 65 | 0 | 56 |
| Progression-free survival | 65 | 60 | 5 | 55 |

Week-9 landmark cohort: attrition and remaining events

*Note:* Patients with death, progression, or censoring before week 9 are
excluded from the corresponding endpoint. Survival is measured from week
9 onward.

## Descriptive Context (Landmark Cohort)

### Characteristics of the Week-9 Landmark Cohorts by TLR Status

The two tables below describe the patients who enter the landmark
analysis, overall and split by the TLR status that becomes defined at
the week-9 scan. The at-risk set is endpoint-specific: the
**overall-survival** cohort keeps everyone alive at week 9, whereas the
**progression-free-survival** cohort additionally drops patients who had
already progressed by week 9. The PFS cohort is therefore smaller.
Reporting both fixes the denominators for the OS and PFS models that
follow and makes the composition of the two TLR groups directly
comparable.

#### Overall-survival at-risk set (OSwk \> 9)

| Characteristic                    | At-risk cohort | TLR-positive | TLR-negative |
|:----------------------------------|---------------:|-------------:|-------------:|
| n                                 |             65 |           41 |           24 |
| Age, years, mean (SD)             |    64.0 (10.0) |   65.0 (9.9) |  62.2 (10.2) |
| Female sex                        |     30 (46.2%) |   19 (46.3%) |   11 (45.8%) |
| Control arm (FLOX)                |     29 (44.6%) |   22 (53.7%) |    7 (29.2%) |
| Experimental arm (FLOX/nivolumab) |     36 (55.4%) |   19 (46.3%) |   17 (70.8%) |
| Deaths (OS events)                |     56 (86.2%) |   32 (78.0%) |  24 (100.0%) |
| PFS events (progression or death) |     60 (92.3%) |   36 (87.8%) |  24 (100.0%) |
| CRP-positive                      |     22 (33.8%) |   18 (43.9%) |    4 (16.7%) |
| TMB/BRAF-positive                 |     29 (44.6%) |   18 (43.9%) |   11 (45.8%) |

Characteristics of the week-9 OS landmark cohort (alive at week 9), by
TLR status

Because no patient dies or is censored before week 9, the OS at-risk set
coincides with the full TLR-complete cohort (n = 65; the three early
deaths without a TLR value are already outside this cohort).

#### Progression-free-survival at-risk set (PFSwk \> 9)

| Characteristic                    | At-risk cohort | TLR-positive | TLR-negative |
|:----------------------------------|---------------:|-------------:|-------------:|
| n                                 |             60 |           41 |           19 |
| Age, years, mean (SD)             |     64.4 (9.9) |   65.0 (9.9) |  63.0 (10.2) |
| Female sex                        |     28 (46.7%) |   19 (46.3%) |    9 (47.4%) |
| Control arm (FLOX)                |     28 (46.7%) |   22 (53.7%) |    6 (31.6%) |
| Experimental arm (FLOX/nivolumab) |     32 (53.3%) |   19 (46.3%) |   13 (68.4%) |
| Deaths (OS events)                |     51 (85.0%) |   32 (78.0%) |  19 (100.0%) |
| PFS events (progression or death) |     55 (91.7%) |   36 (87.8%) |  19 (100.0%) |
| CRP-positive                      |     22 (36.7%) |   18 (43.9%) |    4 (21.1%) |
| TMB/BRAF-positive                 |     26 (43.3%) |   18 (43.9%) |    8 (42.1%) |

Characteristics of the week-9 PFS landmark cohort (alive and
progression-free at week 9), by TLR status

The PFS at-risk set drops the 5 patient(s) who had already progressed,
died, or were censored by week 9, leaving n = 60.

*Definitions:* percentages are column percentages; the “Deaths” and “PFS
events” rows count events over the whole follow-up, not only after the
landmark. CRP-positive = CRP \< 5 mg/L; TMB/BRAF-positive = TMB \>= 9
mut/Mb or BRAF mutation; TLR-positive = tumour lesion reduction \>= 10%
at the first on-treatment CT.

**Note:** In both cohorts the two TLR groups are well matched on age,
sex, and TMB/BRAF status, but differ sharply on treatment arm and CRP:
the TLR-positive group is enriched for control-arm patients and for
CRP-positivity, while every TLR-negative patient died during follow-up.
These imbalances are the descriptive counterpart of the interaction
estimates below and reinforce that TLR status is realised after
randomization rather than being a baseline characteristic.

### TLR Prevalence by Treatment Arm at Week 9

| Treatment Arm            |   N | TLR-positive |    % |
|:-------------------------|----:|-------------:|-----:|
| Control (FLOX)           |  29 |           22 | 75.9 |
| Experimental (FLOX/nivo) |  36 |           19 | 52.8 |

TLR Prevalence by Treatment Arm, Week-9 Landmark Cohort (OS at-risk set)

*Note:* Restricted to patients alive at week 9. TLR = Tumour Lesion
Reduction \>= 10% at first CT scan.

### Landmark Kaplan-Meier Curves Stratified by TLR

<img
src="clinical_effectiveness_files/figure-commonmark/km-lm-tlr-1.png"
style="width:100.0%" data-fig-align="center" />

**Note:** Follow-up time is measured from the week-9 landmark. All
patients shown were alive (OS) or alive and progression-free (PFS) at
the landmark, so the curves are free of the guarantee-time bias
affecting the responder-analysis curves above.



## Standard Cox vs Firth (Landmark)

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| TLR x Rx          | 2.43 (0.74-7.95)         | 2.47 (0.75-7.75)  |
| TLR (main effect) | 0.21 (0.08-0.54)         | 0.21 (0.08-0.55)  |
| Rx (Experimental) | 0.60 (0.23-1.55)         | 0.59 (0.24-1.55)  |
| CRP               | 0.52 (0.25-1.06)         | 0.53 (0.25-1.05)  |
| TMB/BRAF          | 0.86 (0.48-1.52)         | 0.86 (0.48-1.52)  |
| Age               | 1.01 (0.98-1.04)         | 1.01 (0.98-1.04)  |
| Sex               | 1.48 (0.81-2.68)         | 1.47 (0.82-2.67)  |

Landmark (Week 9) TLR Model - Overall Survival: Standard Cox vs Firth

| Term              | Standard Cox HR (95% CI) | Firth HR (95% CI) |
|:------------------|:-------------------------|:------------------|
| TLR x Rx          | 1.86 (0.52-6.64)         | 1.88 (0.52-6.38)  |
| TLR (main effect) | 0.22 (0.08-0.61)         | 0.22 (0.09-0.63)  |
| Rx (Experimental) | 0.55 (0.19-1.61)         | 0.54 (0.20-1.63)  |
| CRP               | 0.53 (0.26-1.10)         | 0.54 (0.26-1.09)  |
| TMB/BRAF          | 0.44 (0.23-0.84)         | 0.45 (0.23-0.83)  |
| Age               | 0.99 (0.96-1.03)         | 0.99 (0.96-1.03)  |
| Sex               | 1.02 (0.55-1.91)         | 1.02 (0.55-1.91)  |

Landmark (Week 9) TLR Model - Progression-Free Survival: Standard Cox vs
Firth

## Firth-Corrected Cox Model (Landmark)

### Overall Survival

| Term                       | HR (95% CI)      | PLRT p-value |
|:---------------------------|:-----------------|:-------------|
| RxExperimental arm:tlr_num | 2.47 (0.75-7.75) | 0.135        |
| Age                        | 1.01 (0.98-1.04) | 0.626        |
| sex1                       | 1.47 (0.82-2.67) | 0.201        |
| RxExperimental arm         | 0.59 (0.24-1.55) | 0.270        |
| crp_num                    | 0.53 (0.25-1.05) | 0.071        |
| tmb_braf_num               | 0.86 (0.48-1.52) | 0.609        |
| tlr_num                    | 0.21 (0.08-0.55) | 0.002        |

Landmark TLR Model: Overall Survival - Full Coefficient Table

### Progression-Free Survival

| Term                       | HR (95% CI)      | PLRT p-value |
|:---------------------------|:-----------------|:-------------|
| RxExperimental arm:tlr_num | 1.88 (0.52-6.38) | 0.327        |
| Age                        | 0.99 (0.96-1.03) | 0.703        |
| sex1                       | 1.02 (0.55-1.91) | 0.946        |
| RxExperimental arm         | 0.54 (0.20-1.63) | 0.264        |
| crp_num                    | 0.54 (0.26-1.09) | 0.085        |
| tmb_braf_num               | 0.45 (0.23-0.83) | 0.011        |
| tlr_num                    | 0.22 (0.09-0.63) | 0.006        |

Landmark TLR Model: Progression-Free Survival - Full Coefficient Table



## Ridge Sensitivity Analysis (Landmark, TLR Terms Penalized)

Ridge regression is applied to the landmark cohorts with the same
targeted penalty as the responder analysis: only the TLR main effect and
the TLR x Rx interaction are L2-shrunk, while Age, sex, Rx, CRP, and
TMB/BRAF remain unpenalized.

| Term              | Firth (95% CI)   | Ridge HR (SE)   |
|:------------------|:-----------------|:----------------|
| TLR x Rx          | 2.47 (0.75-7.75) | 2.12 (SE: 0.58) |
| TLR (main effect) | 0.21 (0.08-0.55) | 0.24 (SE: 0.47) |
| Age               | 1.01 (0.98-1.04) | 1.01 (SE: 0.02) |
| Sex               | 1.47 (0.82-2.67) | 1.49 (SE: 0.30) |
| Rx (Experimental) | 2.47 (0.75-7.75) | 0.66 (SE: 0.48) |
| CRP               | 0.53 (0.25-1.05) | 0.52 (SE: 0.36) |
| TMB/BRAF          | 0.86 (0.48-1.52) | 0.86 (SE: 0.29) |

Landmark TLR Model - Overall Survival: Firth vs Ridge, all coefficients

*Note:* Week-9 landmark cohort. Ridge applies L2 shrinkage to the TLR
main effect and the TLR x Rx interaction only; Age, sex, Rx, CRP, and
TMB/BRAF are unpenalized.

| Term              | Firth (95% CI)   | Ridge HR (SE)   |
|:------------------|:-----------------|:----------------|
| TLR x Rx          | 1.88 (0.52-6.38) | 1.63 (SE: 0.62) |
| TLR (main effect) | 0.22 (0.09-0.63) | 0.25 (SE: 0.50) |
| Age               | 0.99 (0.96-1.03) | 0.99 (SE: 0.02) |
| Sex               | 1.02 (0.55-1.91) | 1.03 (SE: 0.32) |
| Rx (Experimental) | 1.88 (0.52-6.38) | 0.61 (SE: 0.54) |
| CRP               | 0.54 (0.26-1.09) | 0.53 (SE: 0.37) |
| TMB/BRAF          | 0.45 (0.23-0.83) | 0.45 (SE: 0.33) |

Landmark TLR Model - Progression-Free Survival: Firth vs Ridge, all
coefficients

*Note:* Week-9 landmark cohort. Ridge applies L2 shrinkage to the TLR
main effect and the TLR x Rx interaction only; Age, sex, Rx, CRP, and
TMB/BRAF are unpenalized.

**Interpretation:** The two endpoints behave very differently under the
landmark, and the contrast is itself informative.

*Overall survival* is clean: 0 patient(s) are excluded (no deaths occur
before week 9), so the landmark OS cohort is the full TLR-complete
population and the landmark Firth TLR x Rx HR (2.47 (0.75-7.75)) is
**identical** to the responder estimate (2.47 (0.75-7.75)) — a Cox model
is invariant to a common shift of the time origin when no one leaves the
risk set. The OS responder signal is therefore not an early-death
guarantee-time artefact.

*Progression-free survival* exposes a deeper problem. 5 patients are
excluded, of whom 5 progressed, 5 of them **TLR-negative**. This is not
a coincidence: progression and TLR are read from the *same* first
on-treatment scan, so a patient whose tumour grows at that scan is
simultaneously classified as a progression and as TLR-negative. In total
8 patients progressed at their first scan, all TLR-negative. The fixed
week-9 cut removes only the 5 whose scan fell before week 9 and **keeps
3** whose scan fell after it, with landmark times of 0.3 to 1.0 weeks;
it also keeps 19 patients whose TLR was read after week 9, so for them
the landmark covariate is measured after the landmark. Even so, the
landmark shifts the Firth TLR x Rx HR from 1.49 (0.45-4.67) (responder)
to 1.88 (0.52-6.38) (week-9 landmark). The responder PFS association is
thus **partly tautological** — TLR-negativity and early progression are
the same measurement — which the landmark makes explicit by excluding
the coupled events. The sensitivity landmarks below remove all 8
first-scan progressors.

### Landmark Sensitivity: Scan-Date and Week-12 Landmarks

| Landmark | N (OS) | N (PFS) | PFS events | First-scan progressors retained | TLR read after landmark | OS TLR x Rx HR (95% CI) | PLRT p | PFS TLR x Rx HR (95% CI) | PLRT p |
|:---|---:|---:|---:|---:|---:|:---|:---|:---|:---|
| Fixed week 9 (primary) | 65 | 60 | 55 | 3 | 19 | 2.47 (0.75-7.75) | 0.135 | 1.88 (0.52-6.38) | 0.327 |
| Per-patient first on-treatment CT date | 65 | 57 | 52 | 0 | 0 | 2.42 (0.73-7.58) | 0.144 | 1.88 (0.48-6.85) | 0.350 |
| Fixed week 12 | 65 | 57 | 52 | 0 | 1 | 2.47 (0.75-7.75) | 0.135 | 1.97 (0.51-7.16) | 0.318 |

TLR x Rx interaction under three landmark definitions (Firth Cox, full
adjustment set)

*Note:* Cohorts are patients with the endpoint time strictly after the
landmark; time is measured from the landmark. First on-treatment scans
fell between weeks 6.7 and 12.1 after inclusion. ‘First-scan progressors
retained’ counts patients whose progression was recorded at the first
scan (TLR-negative by construction) but who remain in the PFS cohort;
‘TLR read after landmark’ counts retained PFS patients whose first scan
came after the landmark. The per-patient landmark sets both to zero by
construction.

The per-patient scan-date landmark and the week-12 landmark both exclude
all 8 first-scan progressors (PFS cohort n = 57 and 57 versus 60 at week
9) and give a PFS TLR x Rx HR of 1.88 (0.48-6.85) (PLRT p 0.350) and
1.97 (0.51-7.16) (p 0.318), against 1.88 (0.52-6.38) (p 0.327) at week
9. No patient dies before week 12, so the OS estimate is unchanged under
every landmark. The direction of the PFS estimate is therefore not an
artefact of the three first-scan progressors that the week-9 cut
retains, but the confidence intervals are wide under every definition.

Neither analysis resolves the most fundamental issue: TLR lies on the
causal path of treatment (Rx -\> TLR -\> outcome), so the TLR x Rx
contrast is not a baseline patient-selection estimate under any time
origin; disentangling it would require formal causal mediation, not a
landmark. All estimates remain exploratory and imprecise (week-9 PLRT p
0.135 OS, 0.327 PFS) and must not guide treatment selection.



# Discussion

## Primary interaction results

Neither biomarker-treatment interaction reached conventional
significance in the unified model (Firth HR with 95% CI in parentheses).
For **overall survival**, CRP × Rx yielded HR 0.65 (0.19-2.60) (PLRT p =
0.520) and TMB/BRAF × Rx HR 0.95 (0.32-2.82) (p = 0.919) — point
estimates close to null with very wide confidence intervals. For
**progression-free survival**, the CRP × Rx interaction was
directionally favourable but imprecise (HR 0.42 (0.11-1.68), p = 0.213):
CRP-positive patients in the experimental arm had a lower estimated
hazard of progression or death than CRP-positive controls, though the CI
spans more than an order of magnitude and includes 1. The TMB/BRAF × Rx
PFS interaction was HR 0.66 (0.21-2.05) (p = 0.469). These results are
consistent with the trial being underpowered to detect treatment-effect
heterogeneity (n = 68, 63 PFS events (progression or death), 59 deaths):
even a moderate subgroup effect (HR ~0.5) would require far larger
samples for reliable estimation.

The PFS pattern is worth noting as a hypothesis-generating finding. The
CRP direction (HR 0.42) is consistent with the DAG association test
showing CRP → PFS (unadjusted Firth HR 0.41, p = 0.001; see the DAG
associations report) and the biological hypothesis that low CRP
(reflecting lower systemic inflammation) identifies patients more likely
to respond to immunotherapy. However, the OS CRP interaction estimate
(HR 0.65) is closer to null, suggesting that any PFS benefit does not
clearly translate to an OS benefit in this dataset.

## Proportional hazards

The PH assumption was well supported in overall survival (Schoenfeld
global p = 0.254). In progression-free survival the global test
indicated some departure from proportionality (p = 0.119), and this is
carried by the biomarker main effects rather than the interaction terms:
CRP main effect p = 0.130, TMB/BRAF main effect p = 0.084, versus CRP ×
Rx p = 0.272 and TMB/BRAF × Rx p = 0.203. This suggests that the
prognostic hazard ratios for CRP and TMB/BRAF are not constant over time
in PFS, which may reflect different kinetics of early vs late events
(now including deaths) across biomarker-defined subgroups. The PFS
main-effect estimates for CRP and TMB/BRAF should therefore be read as
time-averaged effects; a time-varying Cox model or landmark analysis
would be more appropriate for characterizing them. The interaction
estimates themselves are not flagged by this diagnostic.

## Ridge regression sensitivity

The ridge model penalizes the CRP and TMB/BRAF main effects with a
cross-validated lambda (OS lambda.min = 86.857, PFS lambda.min =
109.761) and estimates the interaction terms as pre-built product terms,
so the ridge and Firth interaction coefficients are the same
treatment-by-biomarker contrast. For OS, the CRP × Rx estimate changes
substantially under this penalization (Firth HR 0.65 → Ridge HR 0.40),
and the TMB/BRAF × Rx estimate changes little (0.95 → 0.86). For PFS,
the CRP × Rx estimate changes little (0.42 → 0.33) and the TMB/BRAF × Rx
estimate changes little (0.66 → 0.51). A shift of at least 1.5-fold in
the HR is labelled “substantial”. Because only the main effects are
penalized, any shift in an interaction term reflects redistribution of
the shrunken main-effect signal rather than independent support for a
particular effect size: an interaction that moves substantially is
sensitive to how the biomarker main effects are estimated and should not
be over-interpreted, while one that changes little is not driven by the
main-effect estimates.

## TLR responder analysis (exploratory)

The exploratory responder model —
`Surv ~ Age + sex + Rx + CRP + TMB/BRAF + TLR + TLR:Rx` — estimates the
differential treatment effect between TLR-positive and TLR-negative
subgroups while adjusting for the pre-immunotherapy biomarkers as
prognostic main effects. The Firth-estimated TLR × Rx HR was 2.47
(0.75-7.75) for overall survival and 1.49 (0.45-4.67) for
progression-free survival. The OS estimate sits **above 1.0**, which
directionally implies the (TLR-positive vs TLR-negative) hazard contrast
is *worse* on the experimental arm than on the control arm — the
opposite of a “TLR-positive predicts immunotherapy benefit” pattern. PFS
is closer to null. Ridge shrinkage of the TLR terms gives 2.12 (SE:
0.58) (OS) and 1.33 (SE: 0.58) (PFS): ridge moves the OS estimate
modestly toward null but preserves the directional pattern, and barely
changes the PFS estimate. PLRT p-values (0.135 OS, 0.504 PFS) do not
reach conventional significance, consistent with the wide
profile-likelihood CIs in this small sample.

The directional pattern is consistent with the descriptive imbalance:
TLR-positive prevalence is 75.9% in the control arm and 52.8% in the
experimental arm. If TLR cleanly captured immunotherapy response one
would expect higher TLR-positivity in the experimental arm, not lower.
The observed prevalence pattern, together with the OS direction (HR
\> 1) of the TLR × Rx interaction, is more consistent with TLR acting as
a **prognostic** marker (identifying chemotherapy-responsive disease,
which is more prevalent in the control arm by chance or by timing) than
as a predictive marker for immunotherapy benefit. Possible explanations
for the imbalance include small-sample randomization imbalance, a timing
artefact (TLR measured at a fixed CT that may occur during different
treatment phases across arms), or genuine heterogeneity in the patient
mix. Because TLR is itself influenced by treatment, the TLR × Rx
interaction is a responder-stratified contrast — useful for hypothesis
generation, but not a causal predictive-biomarker estimate.

## TLR landmark analysis (exploratory)

A **week-9 landmark analysis** conditions on survival to the protocol
time of the second (first on-treatment) CT — the scan at which TLR is
defined — and measures survival from that point, to test whether the
responder signal is a guarantee-time artefact. The two endpoints diverge
instructively. For **OS**, 0 patients are excluded (no early deaths), so
the landmark estimate is identical to the responder estimate (2.47
(0.75-7.75)): the OS signal is not driven by early-death immortal time.
For **PFS**, 5 patients are excluded, with 5 of the 5 excluded
progressors being TLR-negative, because progression and TLR are
ascertained at the same scan; removing these coupled events moves the
Firth TLR × Rx HR from 1.49 (0.45-4.67) to 1.88 (0.52-6.38). The PFS
responder association is therefore partly tautological. Because the
actual first scans fell between weeks 6.7 and 12.1, the fixed week-9 cut
retains 3 of the 8 first-scan progressors and 19 patients whose TLR was
read after the landmark; a per-patient landmark at the scan date and a
fixed week-12 landmark, which exclude all first-scan progressors, give
PFS TLR × Rx HRs of 1.88 (0.48-6.85) and 1.97 (0.51-7.16), so the
direction of the estimate does not depend on the landmark choice. Under
any time origin TLR remains on the causal path of treatment, so a
definitive predictive-biomarker estimate would require causal mediation
rather than a landmark. Both analyses remain underpowered and
exploratory, and neither supports baseline treatment selection on TLR.



# Summary

**Study:** METIMMOX-1 biomarker subgroup analysis. N = 68 complete
cases, 59 OS events, 63 PFS events (progression or death).

**Methods:** Two Firth-corrected Cox analyses are presented. **Primary
(DAG-informed):**
`Surv ~ Age + sex + Rx + CRP + TMB/BRAF + CRP:Rx + TMB/BRAF:Rx` —
pre-immunotherapy biomarkers (week-4 CRP before the first nivolumab
dose; baseline TMB/BRAF) as effect modifiers; TLR omitted because it is
measured on treatment. **Exploratory responder analysis:**
`Surv ~ Age + sex + Rx + CRP + TMB/BRAF + TLR + TLR:Rx` — adds TLR and
its treatment interaction, with CRP and TMB/BRAF retained as prognostic
main-effect adjustments. Profile likelihood CIs, PLRT for interaction
testing. Ridge sensitivity analysis penalizes the CRP/TMB main effects
in the primary model and the TLR terms (main effect + TLR:Rx) in the
exploratory model.

**Proportional hazards:** No violations in OS (global p = 0.254). In PFS
the global test gives p = 0.119, carried by the biomarker main effects
(CRP p = 0.130; TMB/BRAF p = 0.084) rather than the interaction terms
(CRP × Rx p = 0.272; TMB/BRAF × Rx p = 0.203); PFS biomarker main-effect
estimates should be read as time-averaged.

**Interaction estimates:**

| Interaction   | OS HR (95% CI)   | PLRT p | PFS HR (95% CI)  | PLRT p |
|---------------|------------------|--------|------------------|--------|
| CRP × Rx      | 0.65 (0.19-2.60) | 0.520  | 0.42 (0.11-1.68) | 0.213  |
| TMB/BRAF × Rx | 0.95 (0.32-2.82) | 0.919  | 0.66 (0.21-2.05) | 0.469  |

No interaction is statistically significant. CRP × Rx PFS (HR 0.42) is
the strongest directional signal; under ridge regularization of the
biomarker main effects (glmnet, CV lambda) it changes little (Ridge HR
0.33). The TMB/BRAF × Rx OS interaction changes little (Ridge HR 0.86 vs
Firth 0.95). The PFS PH departure concerns the biomarker main effects,
not the interaction terms.

**TLR responder analysis (exploratory):** Firth TLR × Rx HR 2.47
(0.75-7.75) for OS and 1.49 (0.45-4.67) for PFS. Ridge-penalized TLR
terms: 2.12 (SE: 0.58) (OS), 1.33 (SE: 0.58) (PFS). PLRT p-values (0.135
OS, 0.504 PFS) do not reach conventional significance in this small
sample. TLR-positive prevalence: control 75.9%, experimental 52.8% —
markedly higher TLR-positivity in the control arm. The OS direction (HR
\> 1) of TLR × Rx, combined with the prevalence pattern, is more
consistent with TLR acting as a prognostic marker for chemo-responsive
disease than as a predictive marker for immunotherapy benefit. Because
TLR is post-randomization, this is a responder-stratified,
hypothesis-generating contrast rather than a causal predictive-biomarker
estimate.

**TLR week-9 landmark analysis (exploratory):** A landmark analysis
conditioning on survival to the second (first on-treatment) CT in week 9
— the scan at which TLR is defined — was added. For OS, 0 patients are
excluded (no early deaths), so the landmark HR 2.47 (0.75-7.75) is
identical to the responder value: the OS signal is not an early-death
immortal-time artefact. For PFS, 5 patients are excluded, 5 of the 5
excluded progressors being TLR-negative because progression and TLR are
read from the same scan; this shifts the HR from 1.49 (0.45-4.67)
(responder) to 1.88 (0.52-6.38) (landmark), exposing the PFS association
as partly tautological (PLRT p 0.135 OS, 0.327 PFS). Because first scans
fell between weeks 6.7 and 12.1, the week-9 cut retains 3 first-scan
progressors and 19 patients scanned after the landmark; sensitivity
landmarks at each patient’s scan date (PFS HR 1.88 (0.48-6.85)) and at
week 12 (1.97 (0.51-7.16)) exclude all first-scan progressors and leave
the direction unchanged. Both analyses are retained for comparison; both
are underpowered and exploratory, and the residual confounding from TLR
being treatment-influenced needs causal mediation, not a landmark.
Neither supports baseline treatment selection on TLR.

------------------------------------------------------------------------

**Report completed on:** 2026-09-22  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 3.7
