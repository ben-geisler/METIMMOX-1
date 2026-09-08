# DAG Association Tests
Ben Geisler
2026-09-07

- [Introduction](#introduction)
- [Methods](#methods)
  - [Data and Variables](#data-and-variables)
  - [Graph-derived test lists](#graph-derived-test-lists)
  - [Statistical Tests](#statistical-tests)
  - [Multiple Testing](#multiple-testing)
- [Direct Edge Tests](#direct-edge-tests)
  - [Covariate edges (Age/Sex -\>
    TMB_BRAF)](#covariate-edges-agesex---tmb_braf)
  - [Untested edges (latent and
    definitional)](#untested-edges-latent-and-definitional)
  - [Biomarker edges (TMB_BRAF, CRP -\> PFS, OS,
    TLR)](#biomarker-edges-tmb_braf-crp---pfs-os-tlr)
  - [Treatment edges (T -\> PFS, OS,
    TLR)](#treatment-edges-t---pfs-os-tlr)
  - [Interaction edges (TxTMB, TxCRP -\> PFS,
    OS)](#interaction-edges-txtmb-txcrp---pfs-os)
  - [Mediator edge (TLR -\> PFS)](#mediator-edge-tlr---pfs)
  - [Outcome chain (PFS -\> OS)](#outcome-chain-pfs---os)
  - [Edge summary table](#edge-summary-table)
- [Conditional Independence Tests](#conditional-independence-tests)
  - [Arm balance (marginal CIs involving
    T)](#arm-balance-marginal-cis-involving-t)
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
and 03. The main complete-case dataset used for survival and interaction
models contains 68 patients. TLR in the DAG corresponds to `tlr` in the
data, and the TLR subset with observed `tlr` contains 68 patients. The
TLR landmark PFS cohort (patients progression-free after the week 9
landmark, see below) contains 63 patients.

Patients without follow-up CT assessment have `tlr = NA` and were
excluded from TLR analyses. Those analyses therefore rely on a
missing-at-random assumption conditional on the observed covariates
included in each model.

| Dataset | N | Deaths | PFS_events |
|:---|:---|:---|:---|
| Source dataset after scripts 02 and 03 | 74 | 65 | 69 |
| Main complete-case dataset for survival and interaction models | 68 | 59 | 63 |
| TLR subset with observed tlr | 68 | 59 | 63 |
| TLR landmark PFS cohort (progression-free after week 9) | 63 | 54 | 58 |

Analysis datasets used in the DAG association report

For the interaction nodes, `T` was coded as `T_num = 1` for the
experimental arm and `0` for control, `TxTMB` was defined as
`T_num * tmb_braf` and `TxCRP` as `T_num * crp`.

## Graph-derived test lists

The graph tested here is the canonical `dag` object from
`dag_helpers.R`, the same object rendered in `dag.qmd`. Both the edge
list and the list of implied conditional independencies are derived from
that object at render time (`dagitty::edges()` and
`dagitty::impliedConditionalIndependencies()`), and the report stops if
the DAG implies a statement that has no prespecified test or a test
targets a statement the DAG no longer implies. An earlier version typed
both lists by hand from a graph without the `TxCRP` node and therefore
tested 13 of the 19 implied independencies and neither `TxCRP` edge.

The DAG has 26 edges. 4 originate from the latent node `U` and are not
testable. 4 point into the constructed interaction nodes (`T -> TxTMB`,
`TMB_BRAF -> TxTMB`, `T -> TxCRP`, `CRP -> TxCRP`) and hold by
definition. The remaining 18 edges are tested. The DAG implies 19
conditional independencies; 2 of them relate the two interaction nodes
given `T` and hold by construction (given `T = 0` both nodes are zero,
and given `T = 1` and the conditioned biomarker one of them is
constant), leaving 17 testable statements, of which 4 are marginal
statements involving `T` that double as arm-balance checks.

## Statistical Tests

Each testable edge or independence was tested using a prespecified
method:

- Wilcoxon rank-sum for continuous-versus-binary comparisons
- Fisher’s exact test for 2 x 2 tables
- Mantel-Haenszel tests for binary associations conditional on one
  binary stratifier
- Firth-corrected Cox regression (`coxphf`) for time-to-event endpoints
- Firth-corrected logistic regression (`logistf`) for TLR models
- A conditional permutation test (`coin::independence_test`) for
  `Age _||_ TxTMB | {TMB_BRAF}`

Two edge tests need a time origin other than randomisation (issue
\#155):

- **`PFS -> OS`**: progression is entered as a time-dependent indicator
  in a counting-process Firth Cox model for death (`survival::tmerge` on
  the progression-only trial event, `ProgressionExit` / `TTPwk`).
  Entering PFS time as a fixed covariate, as an earlier version did, is
  tautological because `PFSwk <= OSwk` by construction and PFS time is
  unknown at baseline (immortal-time bias).
- **`TLR -> PFS`**: TLR is read at the first on-treatment CT, so a test
  from randomisation conditions on surviving progression-free to that
  scan (guarantee-time bias). The edge is therefore tested on the week 9
  landmark PFS cohort from `tlr_landmark.R`, the same cohort and time
  origin used by `clinical_effectiveness.qmd` and Figure 1 of the
  clinical effectiveness paper: patients progression-free after the
  landmark, with time measured from the landmark. Of the 8 patients
  whose progression is recorded at the first scan itself (all
  TLR-negative because progression and TLR are read from the same scan),
  3 remain in the week-9 cohort because their scan fell after week 9,
  and 19 retained patients had their TLR read after the landmark. The
  clinical effectiveness report repeats the landmark analysis at each
  patient’s own scan date and at week 12 as sensitivity analyses.

## Multiple Testing

No multiplicity correction was applied in this version of the report.
All statistical summaries therefore use raw p-values only, and the
findings should be interpreted as exploratory rather than confirmatory.

# Direct Edge Tests

## Covariate edges (Age/Sex -\> TMB_BRAF)

Age and sex are the only baseline covariates with direct arrows into
`TMB_BRAF`. The age association was tested with a Wilcoxon rank-sum test
and the sex association with Fisher’s exact test.

## Untested edges (latent and definitional)

| Edge               | Reason                                           |
|:-------------------|:-------------------------------------------------|
| U -\> CRP          | Latent parent (U is unobserved)                  |
| U -\> OS           | Latent parent (U is unobserved)                  |
| U -\> PFS          | Latent parent (U is unobserved)                  |
| U -\> TMB_BRAF     | Latent parent (U is unobserved)                  |
| CRP -\> TxCRP      | Definitional (interaction node is T x biomarker) |
| T -\> TxCRP        | Definitional (interaction node is T x biomarker) |
| T -\> TxTMB        | Definitional (interaction node is T x biomarker) |
| TMB_BRAF -\> TxTMB | Definitional (interaction node is T x biomarker) |

DAG edges without an empirical test

## Biomarker edges (TMB_BRAF, CRP -\> PFS, OS, TLR)

Direct biomarker effects on `PFS` and `OS` were tested with univariable
Firth Cox models in the main dataset. Direct effects on `TLR` were
tested in the `tlr` subset using Firth logistic regression.

## Treatment edges (T -\> PFS, OS, TLR)

Treatment effects were evaluated with `T_num` as a binary predictor,
preserving the randomised-arm coding from the trial extract.

## Interaction edges (TxTMB, TxCRP -\> PFS, OS)

Each interaction node was tested using the interaction coefficient from
`Surv(...) ~ T_num + biomarker + T_num:biomarker` for both `PFS` and
`OS`. These are the same interaction-scale questions emphasised in
`clinical_effectiveness.qmd`, but used here as DAG edge checks rather
than biomarker-model comparisons. The `TxCRP` edges are the ones that
drive the CRP-guided strategy in the economic model.

## Mediator edge (TLR -\> PFS)

The `TLR -> PFS` edge was estimated in the week 9 landmark PFS cohort (n
= 63) with time measured from the landmark, as described in Methods.

## Outcome chain (PFS -\> OS)

The `PFS -> OS` edge was represented by a Firth Cox model for death with
a time-dependent progression indicator (n = 68 patients, 116
counting-process intervals). The hazard ratio compares the death hazard
after a recorded progression with the hazard before it.

## Edge summary table

| Edge | Test | N | Effect | p-value | Conclusion |
|:---|:---|:---|:---|:---|:---|
| Age -\> OS | Firth Cox | 68 | 1.00 (95% CI 0.98 to 1.03) | 0.861 | No association detected (p \>= 0.05) |
| Age -\> TMB_BRAF | Wilcoxon rank-sum | 69 | Median Age: TMB/BRAF=0=61.0, TMB/BRAF=1=65.0 | 0.450 | No association detected (p \>= 0.05) |
| CRP -\> OS | Firth Cox | 68 | 0.49 (95% CI 0.27 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| CRP -\> PFS | Firth Cox | 68 | 0.41 (95% CI 0.23 to 0.71) | 0.001 | Association detected (p \< 0.05) |
| CRP -\> TLR | Firth logistic | 68 | 3.48 (95% CI 1.14 to 12.59) | 0.028 | Association detected (p \< 0.05) |
| PFS -\> OS | Firth Cox, time-dependent progression | 68 | 3.75 (95% CI 2.09 to 7.12) | \<0.001 | Association detected (p \< 0.05) |
| Sex -\> TMB_BRAF | Fisher’s exact | 69 | 0.87 (95% CI 0.30 to 2.48) | 0.812 | No association detected (p \>= 0.05) |
| T -\> OS | Firth Cox | 68 | 0.96 (95% CI 0.58 to 1.61) | 0.873 | No association detected (p \>= 0.05) |
| T -\> PFS | Firth Cox | 68 | 0.80 (95% CI 0.49 to 1.33) | 0.392 | No association detected (p \>= 0.05) |
| T -\> TLR | Firth logistic | 68 | 0.33 (95% CI 0.11 to 0.90) | 0.030 | Association detected (p \< 0.05) |
| TLR -\> PFS | Firth Cox (week 9 landmark PFS cohort) | 63 | 0.36 (95% CI 0.21 to 0.65) | \<0.001 | Association detected (p \< 0.05) |
| TMB_BRAF -\> OS | Firth Cox | 68 | 0.68 (95% CI 0.40 to 1.13) | 0.139 | No association detected (p \>= 0.05) |
| TMB_BRAF -\> PFS | Firth Cox | 68 | 0.49 (95% CI 0.28 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| TMB_BRAF -\> TLR | Firth logistic | 65 | 0.92 (95% CI 0.34 to 2.52) | 0.877 | No association detected (p \>= 0.05) |
| TxCRP -\> OS | Firth Cox interaction | 68 | 0.63 (95% CI 0.19 to 2.41) | 0.480 | No association detected (p \>= 0.05) |
| TxCRP -\> PFS | Firth Cox interaction | 68 | 0.26 (95% CI 0.08 to 0.91) | 0.036 | Association detected (p \< 0.05) |
| TxTMB -\> OS | Firth Cox interaction | 68 | 0.63 (95% CI 0.22 to 1.79) | 0.385 | No association detected (p \>= 0.05) |
| TxTMB -\> PFS | Firth Cox interaction | 68 | 0.36 (95% CI 0.13 to 1.00) | 0.051 | No association detected (p \>= 0.05) |

Direct edge tests

# Conditional Independence Tests

## Arm balance (marginal CIs involving T)

The 4 CI statements involving `T` without a conditioning set double as
arm-balance checks. If any of them were significant, that would suggest
imbalance between arms in this realised sample rather than a structural
failure of the DAG. Age, Sex and TMB/BRAF are baseline characteristics,
so their checks are randomisation checks in the strict sense. CRP is the
week-4 value (before the first nivolumab dose, after two FLOX cycles
common to both arms), so its check asks whether the arms were balanced
on CRP at the point where the immunotherapy decision is made.

| Item              | p-value | Conclusion                       |
|:------------------|:--------|:---------------------------------|
| Age *\|\|* T      | 0.172   | Compatible with arm balance      |
| CRP *\|\|* T      | 0.023   | Arm imbalance in realised sample |
| Sex *\|\|* T      | 0.254   | Compatible with arm balance      |
| T *\|\|* TMB_BRAF | 1.000   | Compatible with arm balance      |

Arm balance checks implied by the DAG

One arm-balance check, `CRP _||_ T`, was significant on the raw scale (p
= 0.023): week-4 CRP-positivity was 17/36 in the experimental arm and
7/35 in the control arm. The baseline (cycle 1 day 1) CRP, which is not
used in the analysis, was balanced (8/38 vs 7/35; Fisher p = 1), so
randomisation itself was not at fault. Because both arms receive
identical FLOX up to week 4 and nivolumab has not yet been given, the
week-4 difference cannot be a treatment effect either; it is read as
chance imbalance on a week-4 measurement in a small trial. The other
checks, all on baseline characteristics, remained non-significant.

## Structural CIs (conditional)

The DAG implies 19 observable conditional independencies after excluding
the latent node `U`. 2 hold by construction (see Methods) and are listed
without a test; the remaining 17 were tested using the prespecified
conditional models or stratified tests. Statements that condition on `T`
and a biomarker and involve the matching interaction node (for example
`TLR _||_ TxTMB | {T, TMB_BRAF}`) are tested parametrically as the
interaction coefficient of a Firth logistic model, because within strata
of `T` and the biomarker the interaction node is constant.

## CI summary table

| CI statement | Test | N | Effect | p-value | Arm balance check | Conclusion |
|:---|:---|:---|:---|:---|:---|:---|
| Age *\|\|* CRP | Wilcoxon rank-sum | 71 | Median Age: CRP=0=65.0, CRP=1=64.5 | 0.961 | No | Compatible with DAG-implied CI |
| Age *\|\|* Sex | Wilcoxon rank-sum | 74 | Median Age: Sex=0=65.5, Sex=1=64.5 | 0.961 | No | Compatible with DAG-implied CI |
| Age *\|\|* T | Wilcoxon rank-sum (randomization check) | 74 | Median Age: T=0=65.5, T=1=60.5 | 0.172 | Yes | Compatible with arm balance |
| Age *\|\|* TLR \| {CRP, TMB_BRAF} | Firth logistic | 65 | 1.03 (95% CI 0.98 to 1.08) | 0.292 | No | Compatible with DAG-implied CI |
| Age *\|\|* TxCRP | Wilcoxon rank-sum | 71 | Median Age: TxCRP=0=65.0, TxCRP=1=68.0 | 0.557 | No | Compatible with DAG-implied CI |
| Age *\|\|* TxTMB \| {TMB_BRAF} | Conditional permutation test | 69 | Standardized Z = 0.77 | 0.457 | No | Compatible with DAG-implied CI |
| CRP *\|\|* Sex | Fisher’s exact | 71 | 0.74 (95% CI 0.25 to 2.23) | 0.619 | No | Compatible with DAG-implied CI |
| CRP *\|\|* T | Fisher’s exact (arm balance at week 4) | 71 | 3.51 (95% CI 1.12 to 12.10) | 0.023 | Yes | Arm imbalance in realised sample |
| CRP *\|\|* TxTMB \| {TMB_BRAF} | Mantel-Haenszel | 68 | 21.67 (95% CI 2.23 to 210.11) | 0.007 | No | CI contradicted by data |
| Sex *\|\|* T | Fisher’s exact (randomization check) | 74 | 0.58 (95% CI 0.21 to 1.59) | 0.254 | Yes | Compatible with arm balance |
| Sex *\|\|* TLR \| {CRP, TMB_BRAF} | Firth logistic | 65 | 1.13 (95% CI 0.40 to 3.21) | 0.812 | No | Compatible with DAG-implied CI |
| Sex *\|\|* TxCRP | Fisher’s exact | 71 | 0.66 (95% CI 0.19 to 2.27) | 0.578 | No | Compatible with DAG-implied CI |
| Sex *\|\|* TxTMB \| {TMB_BRAF} | Mantel-Haenszel | 69 | 0.52 (95% CI 0.12 to 2.17) | 0.592 | No | Compatible with DAG-implied CI |
| T *\|\|* TMB_BRAF | Fisher’s exact (randomization check) | 69 | 0.96 (95% CI 0.33 to 2.76) | 1.000 | Yes | Compatible with arm balance |
| TLR *\|\|* TxCRP \| {CRP, T} | Firth logistic | 68 | 1.25 (95% CI 0.01 to 20.93) | 0.898 | No | Compatible with DAG-implied CI |
| TLR *\|\|* TxTMB \| {T, TMB_BRAF} | Firth logistic | 65 | 2.31 (95% CI 0.29 to 19.27) | 0.426 | No | Compatible with DAG-implied CI |
| TMB_BRAF *\|\|* TxCRP \| {CRP} | Mantel-Haenszel | 68 | 7.14 (95% CI 0.68 to 75.22) | 0.203 | No | Compatible with DAG-implied CI |
| TxCRP *\|\|* TxTMB \| {T, TMB_BRAF} | Definitional (constant within every stratum) | 74 | – | – | No | Holds by construction |
| TxCRP *\|\|* TxTMB \| {CRP, T} | Definitional (constant within every stratum) | 74 | – | – | No | Holds by construction |

Conditional-independence tests

# Consolidated Summary Table

| Domain | Item | Effect | p-value | Conclusion |
|:---|:---|:---|:---|:---|
| Direct edge | Age -\> OS | 1.00 (95% CI 0.98 to 1.03) | 0.861 | No association detected (p \>= 0.05) |
| Direct edge | Age -\> TMB_BRAF | Median Age: TMB/BRAF=0=61.0, TMB/BRAF=1=65.0 | 0.450 | No association detected (p \>= 0.05) |
| Direct edge | CRP -\> OS | 0.49 (95% CI 0.27 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| Direct edge | CRP -\> PFS | 0.41 (95% CI 0.23 to 0.71) | 0.001 | Association detected (p \< 0.05) |
| Direct edge | CRP -\> TLR | 3.48 (95% CI 1.14 to 12.59) | 0.028 | Association detected (p \< 0.05) |
| Direct edge | PFS -\> OS | 3.75 (95% CI 2.09 to 7.12) | \<0.001 | Association detected (p \< 0.05) |
| Direct edge | Sex -\> TMB_BRAF | 0.87 (95% CI 0.30 to 2.48) | 0.812 | No association detected (p \>= 0.05) |
| Direct edge | T -\> OS | 0.96 (95% CI 0.58 to 1.61) | 0.873 | No association detected (p \>= 0.05) |
| Direct edge | T -\> PFS | 0.80 (95% CI 0.49 to 1.33) | 0.392 | No association detected (p \>= 0.05) |
| Direct edge | T -\> TLR | 0.33 (95% CI 0.11 to 0.90) | 0.030 | Association detected (p \< 0.05) |
| Direct edge | TLR -\> PFS | 0.36 (95% CI 0.21 to 0.65) | \<0.001 | Association detected (p \< 0.05) |
| Direct edge | TMB_BRAF -\> OS | 0.68 (95% CI 0.40 to 1.13) | 0.139 | No association detected (p \>= 0.05) |
| Direct edge | TMB_BRAF -\> PFS | 0.49 (95% CI 0.28 to 0.84) | 0.009 | Association detected (p \< 0.05) |
| Direct edge | TMB_BRAF -\> TLR | 0.92 (95% CI 0.34 to 2.52) | 0.877 | No association detected (p \>= 0.05) |
| Direct edge | TxCRP -\> OS | 0.63 (95% CI 0.19 to 2.41) | 0.480 | No association detected (p \>= 0.05) |
| Direct edge | TxCRP -\> PFS | 0.26 (95% CI 0.08 to 0.91) | 0.036 | Association detected (p \< 0.05) |
| Direct edge | TxTMB -\> OS | 0.63 (95% CI 0.22 to 1.79) | 0.385 | No association detected (p \>= 0.05) |
| Direct edge | TxTMB -\> PFS | 0.36 (95% CI 0.13 to 1.00) | 0.051 | No association detected (p \>= 0.05) |
| Conditional independence | Age *\|\|* CRP | Median Age: CRP=0=65.0, CRP=1=64.5 | 0.961 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* Sex | Median Age: Sex=0=65.5, Sex=1=64.5 | 0.961 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* T | Median Age: T=0=65.5, T=1=60.5 | 0.172 | Compatible with arm balance |
| Conditional independence | Age *\|\|* TLR \| {CRP, TMB_BRAF} | 1.03 (95% CI 0.98 to 1.08) | 0.292 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* TxCRP | Median Age: TxCRP=0=65.0, TxCRP=1=68.0 | 0.557 | Compatible with DAG-implied CI |
| Conditional independence | Age *\|\|* TxTMB \| {TMB_BRAF} | Standardized Z = 0.77 | 0.457 | Compatible with DAG-implied CI |
| Conditional independence | CRP *\|\|* Sex | 0.74 (95% CI 0.25 to 2.23) | 0.619 | Compatible with DAG-implied CI |
| Conditional independence | CRP *\|\|* T | 3.51 (95% CI 1.12 to 12.10) | 0.023 | Arm imbalance in realised sample |
| Conditional independence | CRP *\|\|* TxTMB \| {TMB_BRAF} | 21.67 (95% CI 2.23 to 210.11) | 0.007 | CI contradicted by data |
| Conditional independence | Sex *\|\|* T | 0.58 (95% CI 0.21 to 1.59) | 0.254 | Compatible with arm balance |
| Conditional independence | Sex *\|\|* TLR \| {CRP, TMB_BRAF} | 1.13 (95% CI 0.40 to 3.21) | 0.812 | Compatible with DAG-implied CI |
| Conditional independence | Sex *\|\|* TxCRP | 0.66 (95% CI 0.19 to 2.27) | 0.578 | Compatible with DAG-implied CI |
| Conditional independence | Sex *\|\|* TxTMB \| {TMB_BRAF} | 0.52 (95% CI 0.12 to 2.17) | 0.592 | Compatible with DAG-implied CI |
| Conditional independence | T *\|\|* TMB_BRAF | 0.96 (95% CI 0.33 to 2.76) | 1.000 | Compatible with arm balance |
| Conditional independence | TLR *\|\|* TxCRP \| {CRP, T} | 1.25 (95% CI 0.01 to 20.93) | 0.898 | Compatible with DAG-implied CI |
| Conditional independence | TLR *\|\|* TxTMB \| {T, TMB_BRAF} | 2.31 (95% CI 0.29 to 19.27) | 0.426 | Compatible with DAG-implied CI |
| Conditional independence | TMB_BRAF *\|\|* TxCRP \| {CRP} | 7.14 (95% CI 0.68 to 75.22) | 0.203 | Compatible with DAG-implied CI |
| Conditional independence | TxCRP *\|\|* TxTMB \| {T, TMB_BRAF} | – | – | Holds by construction |
| Conditional independence | TxCRP *\|\|* TxTMB \| {CRP, T} | – | – | Holds by construction |

Consolidated summary of direct edge and conditional-independence tests

# Discussion

The direct-edge analysis found 8 of 18 tested edges with p \< 0.05.
Direct edges with p \< 0.05: CRP -\> OS; CRP -\> PFS; CRP -\> TLR; PFS
-\> OS; T -\> TLR; TLR -\> PFS; TMB_BRAF -\> PFS; TxCRP -\> PFS. Among
the prognostic biomarker edges, those with p \< 0.05 were: CRP -\> PFS;
CRP -\> OS; TMB_BRAF -\> PFS. Among the edges into TLR: CRP -\> TLR; T
-\> TLR. Among the covariate edges: none.

The two disease-course edges are no longer tautological. PFS -\> OS, now
estimated with a time-dependent progression indicator, gives HR 3.75
(95% CI 2.09 to 7.12) (p = \<0.001) for death after versus before a
recorded progression; this is an estimate of how much a progression
raises the subsequent death hazard, not a restatement of PFS \<= OS. TLR
-\> PFS on the week 9 landmark cohort gives HR 0.36 (95% CI 0.21 to
0.65) (p = \<0.001, n = 63). Because progression and TLR are read from
the same scan, TLR-negativity and first-scan progression are partly the
same measurement; the week-9 landmark removes 5 of the 8 first-scan
progressors and keeps 3, so the estimate should be read alongside the
scan-date and week-12 landmark sensitivity analyses in the clinical
effectiveness report.

Treatment-related edges: those with p \< 0.05 were none for the direct
treatment edges and TxCRP -\> PFS for the interaction edges. TxCRP -\>
PFS gives HR 0.26 (95% CI 0.08 to 0.91) (p = 0.036) and TxCRP -\> OS HR
0.63 (95% CI 0.19 to 2.41) (p = 0.480); TxTMB -\> PFS gives HR 0.36 (95%
CI 0.13 to 1.00) (p = 0.051) and TxTMB -\> OS HR 0.63 (95% CI 0.22 to
1.79) (p = 0.385). These are the same treatment-by-biomarker contrasts
as the interaction terms of the economic survival model, estimated here
without covariate adjustment, and the intervals are wide.

Among the 17 tested conditional-independence statements, 15 were
compatible with the DAG at p \>= 0.05 and 2 were contradicted by the
data. Conditional independencies contradicted by the data at p \< 0.05:
CRP *\|\|* T; CRP *\|\|* TxTMB \| {TMB_BRAF}. CRP *\|\|* T had p =
0.023, reflecting a chance arm imbalance in week-4 CRP (baseline CRP was
balanced, and no nivolumab has been given by week 4, so this is neither
a randomisation failure nor a treatment effect), and CRP *\|\|* TxTMB \|
{TMB_BRAF} had p = 0.007, which follows from the same arm imbalance
because TxTMB is a function of T. The new statements involving TxCRP
(Age *\|\|* TxCRP, Sex *\|\|* TxCRP, TMB_BRAF *\|\|* TxCRP \| {CRP}, TLR
*\|\|* TxCRP \| {CRP, T}) had p = 0.557, 0.578, 0.203 and 0.898. The
week-4 CRP arm imbalance therefore remains the main place where the
current DAG looks least secure in this dataset.

Interpretation should remain cautious. The latent node U is unobserved,
TLR analyses depend on non-missing follow-up imaging and on the landmark
definition, the available repository extract is modest in size, and
several tests are necessarily low-power. These checks therefore speak to
empirical compatibility, not causal proof.

# Conclusion

This report adds a data-facing complement to the causal DAG by testing,
from the same `dag` object that `dag.qmd` renders, every observable edge
and every implied conditional independence of the graph. The prognostic
biomarker edges and the two disease-course edges, now estimated with a
time-dependent progression indicator and a landmark cohort rather than
tautological baseline covariates, are where the empirical signal is
clearest. The direct treatment edges are uncertain; of the interaction
edges, `TxCRP -> PFS` is the only one below p = 0.05 in this small
sample, and it should be read together with the adjusted interaction
estimates in the clinical effectiveness report.

The data also surface 2 conditional-independence violations (CRP *\|\|*
T; CRP *\|\|* TxTMB \| {TMB_BRAF}) against the assumed graph, traceable
to the chance arm imbalance in week-4 CRP. The most defensible takeaway
is therefore that the DAG captures much of the prognostic structure in
the data, while the CRP-related balance and interaction structure may
need refinement or stronger justification in future versions.

------------------------------------------------------------------------

**Report completed on:** 2026-09-07  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.2
