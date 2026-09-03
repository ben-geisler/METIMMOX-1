# DAG Association Tests
Ben Geisler
2026-09-03

- [Introduction](#introduction)
- [Methods](#methods)
  - [Data and Variables](#data-and-variables)
  - [Statistical Tests](#statistical-tests)
  - [Multiple Testing](#multiple-testing)
- [Direct Edge Tests](#direct-edge-tests)
  - [Covariate edges (Age/Sex -\>
    TMB_BRAF)](#covariate-edges-agesex---tmb_braf)
  - [Latent edges (U -\> \*)](#latent-edges-u---)
  - [Biomarker edges (TMB_BRAF, CRP -\> PFS, OS,
    TLR)](#biomarker-edges-tmb_braf-crp---pfs-os-tlr)
  - [Treatment edges (T -\> PFS, OS,
    TLR)](#treatment-edges-t---pfs-os-tlr)
  - [Interaction edges (TxTMB -\> PFS,
    OS)](#interaction-edges-txtmb---pfs-os)
  - [Mediator edge (TLR -\> PFS)](#mediator-edge-tlr---pfs)
  - [Outcome chain (PFS -\> OS)](#outcome-chain-pfs---os)
  - [Edge summary table](#edge-summary-table)
- [Conditional Independence Tests](#conditional-independence-tests)
  - [Randomisation balance (marginal CIs involving
    T)](#randomisation-balance-marginal-cis-involving-t)
  - [Structural CIs (conditional)](#structural-cis-conditional)
  - [CI summary table](#ci-summary-table)
- [Consolidated Summary Table](#consolidated-summary-table)
- [Discussion](#discussion)
- [Conclusion](#conclusion)

# Introduction

This report provides the empirical companion to `dag.qmd`. The DAG
encodes the assumed causal structure for METIMMOX-1; this document asks
whether the observable parts of that structure are broadly compatible
with the data available in the repository.

The emphasis is deliberately modest. These tests do not prove causal
direction, and failure to reject a conditional independence does not
verify the DAG. Instead, the report checks whether the observable arrows
imply associations in the expected places and whether the DAG-implied
conditional independencies look plausible in the observed sample.

# Methods

## Data and Variables

Only scripts `02_setup_and_global_variables.R` and
`03_biomarker_strategies.R` were sourced. No survival-model fitting
scripts were run.

The available repository extract contains 74 patients after scripts 02
and 03. The main complete-case dataset used for non-TLR survival and
interaction models contains 68 patients. TLR in the DAG corresponds to
`tlr` in the data, and the TLR subset with observed `tlr` contains 68
patients.

Patients without follow-up CT assessment have `tlr = NA` and were
excluded from TLR analyses. Those analyses therefore rely on a
missing-at-random assumption conditional on the observed covariates
included in each model.

| Dataset | N | Deaths | PFS_events |
|:---|:---|:---|:---|
| Source dataset after scripts 02 and 03 | 74 | 65 | 69 |
| Main complete-case dataset for survival and interaction models | 68 | 59 | 63 |
| TLR subset with observed tlr | 68 | 59 | 63 |

Analysis datasets used in the DAG association report

For the interaction node, `T` was coded as `T_num = 1` for the
experimental arm and `0` for control, and `TxTMB` was defined as
`T_num * tmb_braf`.

## Statistical Tests

Each observable edge was tested using the prespecified method from the
plan:

- Wilcoxon rank-sum for continuous-versus-binary comparisons
- Fisher’s exact test for 2 x 2 tables
- Mantel-Haenszel tests for binary associations conditional on
  `tmb_braf`
- Firth-corrected Cox regression (`coxphf`) for time-to-event endpoints
- Firth-corrected logistic regression (`logistf`) for TLR models
- A conditional permutation test (`coin::independence_test`) for
  `Age _||_ TxTMB | {TMB_BRAF}`

The two arrows into `TxTMB` (`T -> TxTMB` and `TMB_BRAF -> TxTMB`) are
definitional because `TxTMB` is a constructed interaction node. The four
edges from the latent node `U` are not empirically testable and are
noted only in prose.

## Multiple Testing

No multiplicity correction was applied in this version of the report.
All statistical summaries therefore use raw p-values only, and the
findings should be interpreted as exploratory rather than confirmatory.

# Direct Edge Tests

## Covariate edges (Age/Sex -\> TMB_BRAF)

Age and sex are the only baseline covariates with direct arrows into
`TMB_BRAF`. The age association was tested with a Wilcoxon rank-sum test
and the sex association with Fisher’s exact test.

## Latent edges (U -\> \*)

The DAG includes four latent edges: `U -> TMB_BRAF`, `U -> CRP`,
`U -> PFS`, and `U -> OS`. Because `U` is unobserved, these edges are
not testable in the trial data and are therefore not assigned p-values.

## Biomarker edges (TMB_BRAF, CRP -\> PFS, OS, TLR)

Direct biomarker effects on `PFS` and `OS` were tested with univariable
Firth Cox models in `data_main`. Direct effects on `TLR` were tested in
the `tlr` subset using Firth logistic regression.

## Treatment edges (T -\> PFS, OS, TLR)

Treatment effects were evaluated with `T_num` as a binary predictor,
preserving the randomised-arm coding from the trial extract.

## Interaction edges (TxTMB -\> PFS, OS)

The `TxTMB` node was tested using the interaction coefficient from
`Surv(...) ~ T_num + tmb_braf + T_num:tmb_braf` for both `PFS` and `OS`.
These are the same interaction-scale questions emphasized in
`clinical_effectiveness.qmd`, but used here as DAG edge checks rather
than biomarker-model comparisons.

## Mediator edge (TLR -\> PFS)

The `TLR -> PFS` edge was estimated in patients with observed `tlr`.
This is a subset analysis because TLR is undefined when no follow-up
scan is available.

## Outcome chain (PFS -\> OS)

The `PFS -> OS` edge was represented by a Firth Cox model for `OS` with
`PFSwk` entered as a continuous predictor.

## Edge summary table

| Edge | Test | N | Effect | p-value | Conclusion |
|:---|:---|:---|:---|:---|:---|
| Age -\> TMB_BRAF | Wilcoxon rank-sum | 69 | Median Age: TMB/BRAF=0=61.0, TMB/BRAF=1=65.0 | 0.450 | No association detected (p \>= 0.05) |
| Age -\> OS | Firth Cox | 68 | 1.00 (95% CI 0.98 to 1.03) | 0.861 | No association detected (p \>= 0.05) |
| Sex -\> TMB_BRAF | Fisher’s exact | 69 | 0.87 (95% CI 0.30 to 2.48) | 0.812 | No association detected (p \>= 0.05) |
| TMB_BRAF -\> PFS | Firth Cox | 68 | 0.49 (95% CI 0.28 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| TMB_BRAF -\> OS | Firth Cox | 68 | 0.68 (95% CI 0.40 to 1.13) | 0.139 | No association detected (p \>= 0.05) |
| TMB_BRAF -\> TLR | Firth logistic | 65 | 0.92 (95% CI 0.34 to 2.52) | 0.877 | No association detected (p \>= 0.05) |
| CRP -\> PFS | Firth Cox | 68 | 0.41 (95% CI 0.23 to 0.71) | 0.001 | Association detected (p \< 0.05) |
| CRP -\> OS | Firth Cox | 68 | 0.49 (95% CI 0.27 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| CRP -\> TLR | Firth logistic | 68 | 3.48 (95% CI 1.14 to 12.59) | 0.028 | Association detected (p \< 0.05) |
| T -\> PFS | Firth Cox | 68 | 0.80 (95% CI 0.49 to 1.33) | 0.392 | No association detected (p \>= 0.05) |
| T -\> OS | Firth Cox | 68 | 0.96 (95% CI 0.58 to 1.61) | 0.873 | No association detected (p \>= 0.05) |
| T -\> TLR | Firth logistic | 68 | 0.33 (95% CI 0.11 to 0.90) | 0.030 | Association detected (p \< 0.05) |
| TLR -\> PFS | Firth Cox (TLR subset) | 68 | 0.31 (95% CI 0.19 to 0.53) | \<0.001 | Association detected (p \< 0.05) |
| PFS -\> OS | Firth Cox | 68 | 0.98 (95% CI 0.98 to 0.99) | \<0.001 | Association detected (p \< 0.05) |
| TxTMB -\> PFS | Firth Cox interaction | 68 | 0.36 (95% CI 0.13 to 1.00) | 0.051 | No association detected (p \>= 0.05) |
| TxTMB -\> OS | Firth Cox interaction | 68 | 0.63 (95% CI 0.22 to 1.79) | 0.385 | No association detected (p \>= 0.05) |

Direct edge tests

# Conditional Independence Tests

## Randomisation balance (marginal CIs involving T)

The four CI statements involving `T` double as randomisation checks. If
any of them were significant, that would suggest baseline imbalance in
this realized sample rather than a structural failure of the DAG.

| Item              | p-value | Conclusion                            |
|:------------------|:--------|:--------------------------------------|
| Age *\|\|* T      | 0.172   | Compatible with randomization balance |
| CRP *\|\|* T      | 0.023   | Possible randomization imbalance      |
| Sex *\|\|* T      | 0.254   | Compatible with randomization balance |
| T *\|\|* TMB_BRAF | 1.000   | Compatible with randomization balance |

Randomization balance checks implied by the DAG

One randomization-balance check, `CRP _||_ T`, was significant on the
raw scale (`p = 0.023`). The other three randomization checks remained
non-significant.

## Structural CIs (conditional)

The DAG implies 13 observable conditional independencies after excluding
the latent node `U`. These were tested using the prespecified
conditional models or stratified tests.

## CI summary table

| CI statement | Test | N | Effect | p-value | Randomization check | Conclusion |
|:---|:---|:---|:---|:---|:---|:---|
| Age *\|\|* CRP | Wilcoxon rank-sum | 71 | Median Age: CRP=0=65.0, CRP=1=64.5 | 0.961 | No | Compatible with DAG-implied CI |
| Age *\|\|* Sex | Wilcoxon rank-sum | 74 | Median Age: Sex=0=65.5, Sex=1=64.5 | 0.961 | No | Compatible with DAG-implied CI |
| Age *\|\|* T | Wilcoxon rank-sum (randomization check) | 74 | Median Age: T=0=65.5, T=1=60.5 | 0.172 | Yes | Compatible with randomization balance |
| Age *\|\|* TLR \| {CRP, TMB_BRAF} | Firth logistic | 65 | 1.03 (95% CI 0.98 to 1.08) | 0.292 | No | Compatible with DAG-implied CI |
| Age *\|\|* TxTMB \| {TMB_BRAF} | Conditional permutation test | 69 | Standardized Z = 0.77 | 0.457 | No | Compatible with DAG-implied CI |
| CRP *\|\|* Sex | Fisher’s exact | 71 | 0.74 (95% CI 0.25 to 2.23) | 0.619 | No | Compatible with DAG-implied CI |
| CRP *\|\|* T | Fisher’s exact (randomization check) | 71 | 3.51 (95% CI 1.12 to 12.10) | 0.023 | Yes | Possible randomization imbalance |
| CRP *\|\|* TxTMB \| {TMB_BRAF} | Mantel-Haenszel | 68 | 21.67 (95% CI 2.23 to 210.11) | 0.007 | No | CI contradicted by data |
| Sex *\|\|* T | Fisher’s exact (randomization check) | 74 | 0.58 (95% CI 0.21 to 1.59) | 0.254 | Yes | Compatible with randomization balance |
| Sex *\|\|* TxTMB \| {TMB_BRAF} | Mantel-Haenszel | 69 | 0.52 (95% CI 0.12 to 2.17) | 0.592 | No | Compatible with DAG-implied CI |
| Sex *\|\|* TLR \| {CRP, TMB_BRAF} | Firth logistic | 65 | 1.13 (95% CI 0.40 to 3.21) | 0.812 | No | Compatible with DAG-implied CI |
| TLR *\|\|* TxTMB \| {T, TMB_BRAF} | Firth logistic | 65 | 2.31 (95% CI 0.29 to 19.27) | 0.426 | No | Compatible with DAG-implied CI |
| T *\|\|* TMB_BRAF | Fisher’s exact (randomization check) | 69 | 0.96 (95% CI 0.33 to 2.76) | 1.000 | Yes | Compatible with randomization balance |

Conditional-independence tests

# Consolidated Summary Table

| Domain | Item | Effect | p-value | Conclusion |
|:---|:---|:---|:---|:---|
| Direct edge | Age -\> TMB_BRAF | Median Age: TMB/BRAF=0=61.0, TMB/BRAF=1=65.0 | 0.450 | No association detected (p \>= 0.05) |
| Direct edge | Age -\> OS | 1.00 (95% CI 0.98 to 1.03) | 0.861 | No association detected (p \>= 0.05) |
| Direct edge | Sex -\> TMB_BRAF | 0.87 (95% CI 0.30 to 2.48) | 0.812 | No association detected (p \>= 0.05) |
| Direct edge | TMB_BRAF -\> PFS | 0.49 (95% CI 0.28 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| Direct edge | TMB_BRAF -\> OS | 0.68 (95% CI 0.40 to 1.13) | 0.139 | No association detected (p \>= 0.05) |
| Direct edge | TMB_BRAF -\> TLR | 0.92 (95% CI 0.34 to 2.52) | 0.877 | No association detected (p \>= 0.05) |
| Direct edge | CRP -\> PFS | 0.41 (95% CI 0.23 to 0.71) | 0.001 | Association detected (p \< 0.05) |
| Direct edge | CRP -\> OS | 0.49 (95% CI 0.27 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| Direct edge | CRP -\> TLR | 3.48 (95% CI 1.14 to 12.59) | 0.028 | Association detected (p \< 0.05) |
| Direct edge | T -\> PFS | 0.80 (95% CI 0.49 to 1.33) | 0.392 | No association detected (p \>= 0.05) |
| Direct edge | T -\> OS | 0.96 (95% CI 0.58 to 1.61) | 0.873 | No association detected (p \>= 0.05) |
| Direct edge | T -\> TLR | 0.33 (95% CI 0.11 to 0.90) | 0.030 | Association detected (p \< 0.05) |
| Direct edge | TLR -\> PFS | 0.31 (95% CI 0.19 to 0.53) | \<0.001 | Association detected (p \< 0.05) |
| Direct edge | PFS -\> OS | 0.98 (95% CI 0.98 to 0.99) | \<0.001 | Association detected (p \< 0.05) |
| Direct edge | TxTMB -\> PFS | 0.36 (95% CI 0.13 to 1.00) | 0.051 | No association detected (p \>= 0.05) |
| Direct edge | TxTMB -\> OS | 0.63 (95% CI 0.22 to 1.79) | 0.385 | No association detected (p \>= 0.05) |
| Conditional independence | Age *\|\|* CRP | Median Age: CRP=0=65.0, CRP=1=64.5 | 0.961 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* Sex | Median Age: Sex=0=65.5, Sex=1=64.5 | 0.961 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* T | Median Age: T=0=65.5, T=1=60.5 | 0.172 | Compatible with randomization balance |
| Conditional independence | Age *\|\|* TLR \| {CRP, TMB_BRAF} | 1.03 (95% CI 0.98 to 1.08) | 0.292 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* TxTMB \| {TMB_BRAF} | Standardized Z = 0.77 | 0.457 | Compatible with DAG-implied CI |
| Conditional independence | CRP *\|\|* Sex | 0.74 (95% CI 0.25 to 2.23) | 0.619 | Compatible with DAG-implied CI |
| Conditional independence | CRP *\|\|* T | 3.51 (95% CI 1.12 to 12.10) | 0.023 | Possible randomization imbalance |
| Conditional independence | CRP *\|\|* TxTMB \| {TMB_BRAF} | 21.67 (95% CI 2.23 to 210.11) | 0.007 | CI contradicted by data |
| Conditional independence | Sex *\|\|* T | 0.58 (95% CI 0.21 to 1.59) | 0.254 | Compatible with randomization balance |
| Conditional independence | Sex *\|\|* TxTMB \| {TMB_BRAF} | 0.52 (95% CI 0.12 to 2.17) | 0.592 | Compatible with DAG-implied CI |
| Conditional independence | Sex *\|\|* TLR \| {CRP, TMB_BRAF} | 1.13 (95% CI 0.40 to 3.21) | 0.812 | Compatible with DAG-implied CI |
| Conditional independence | TLR *\|\|* TxTMB \| {T, TMB_BRAF} | 2.31 (95% CI 0.29 to 19.27) | 0.426 | Compatible with DAG-implied CI |
| Conditional independence | T *\|\|* TMB_BRAF | 0.96 (95% CI 0.33 to 2.76) | 1.000 | Compatible with randomization balance |

Consolidated summary of direct edge and conditional-independence tests

# Discussion

The direct-edge analysis found 7 of 16 tested edges with p \< 0.05.
Direct edges with p \< 0.05: TMB_BRAF -\> PFS; CRP -\> PFS; CRP -\> OS;
CRP -\> TLR; T -\> TLR; TLR -\> PFS; PFS -\> OS. The strongest supported
pattern sits along the prognostic and disease-course portion of the DAG:
low CRP was associated with lower hazards for both PFS and OS, TMB_BRAF
was associated with lower PFS hazard, TLR was strongly associated with
subsequent PFS, and longer PFS was tightly linked to longer OS. In other
words, the data support the idea that baseline biomarkers and
intermediate disease response are more informative than age or sex for
the downstream survival relationships represented here.

Treatment-related edges were still mixed. T -\> PFS, T -\> OS, and both
TxTMB interaction terms remained non-significant, so the report still
does not provide strong evidence for treatment effect modification
through the TMB/BRAF interaction node. However, the raw-p analysis now
shows p \< 0.05 for both CRP -\> TLR (p = 0.028) and T -\> TLR (p =
0.030), which suggests that tumour lesion reduction may be one route
through which baseline inflammation and treatment assignment connect to
later outcomes.

Among the 13 conditional-independence tests, 11 were compatible with the
DAG at p \>= 0.05 and 2 were contradicted by the data. Conditional
independencies contradicted by the data at p \< 0.05: CRP *\|\|* T; CRP
*\|\|* TxTMB \| {TMB_BRAF}. In practice, both violated statements
involve CRP: CRP *\|\|* T had p = 0.023, suggesting possible baseline
imbalance by treatment arm, and CRP *\|\|* TxTMB \| {TMB_BRAF} had p =
0.007, suggesting residual dependence around the
CRP-treatment-interaction part of the graph. That makes CRP the main
place where the current DAG looks least secure in this dataset.

Interpretation should remain cautious. The latent node U is unobserved,
TLR analyses depend on non-missing follow-up imaging, the available
repository extract is modest in size, and several tests are necessarily
low-power. These checks therefore speak to empirical compatibility, not
causal proof.

# Conclusion

This report adds a data-facing complement to the causal DAG by testing
its observable edges and implied conditional independencies in the
METIMMOX-1 repository extract. Overall, the empirical signal is stronger
for prognostic and disease-course relationships than for
treatment-effect-modification relationships: CRP, TMB_BRAF, TLR, and PFS
show the clearest links to later outcomes, while the direct treatment
and TxTMB interaction edges remain uncertain in this small sample.

Taken together, the observed data support the core survival pathway
(`TLR -> PFS -> OS`) and several biomarker-survival links, but they also
surface two CRP-related tensions with the assumed graph. The most
defensible takeaway is therefore that the DAG captures much of the
prognostic structure in the data, while the CRP-related balance and
interaction structure may need refinement or stronger justification in
future versions.

------------------------------------------------------------------------

**Report completed on:** 2026-09-03  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.0
