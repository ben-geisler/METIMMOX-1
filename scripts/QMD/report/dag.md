# Causal DAG
Ben Geisler
2026-03-27

- [Introduction](#introduction)
- [Node Definitions](#node-definitions)
- [Causal Structure](#causal-structure)
- [Causal Analysis](#causal-analysis)
  - [Adjustment Sets](#adjustment-sets)
  - [Implied Conditional
    Independencies](#implied-conditional-independencies)

# Introduction

This report presents the directed acyclic graph (DAG) for the METIMMOX-1
cost-effectiveness analysis. The DAG encodes the assumed causal
structure between patient characteristics, biomarkers, treatment,
intermediate outcomes, and survival endpoints. It was specified in
DAGitty format and visualised using the **ggdag** package in R.

Three node types are distinguished:

- **Exposure** (T): the randomised treatment (alternating FLOX +
  nivolumab vs FLOX alone)
- **Outcomes** (PFS, OS): the clinical endpoints of interest
- **Latent** (U): unmeasured patient characteristics that act as common
  causes of biomarkers and/or survival

# Node Definitions

| Node | Display | Type | Description |
|:---|:---|:---|:---|
| Age | Age | Covariate | Patient age at baseline |
| Sex | Sex | Covariate | Patient sex |
| U | U | Latent | Unmeasured prognostic factors (latent) |
| TMB_BRAF | TMB/BRAF | Biomarker | Combined biomarker: TMB \>= 9 mut/MB or BRAF mutation |
| CRP | CRP | Biomarker | C-reactive protein (positive: CRP \< 5 mg/L) |
| T | T | Exposure | Treatment: alternating FLOX + nivolumab vs FLOX alone |
| TxTMB | T x TMB | Effect modifier | Interaction between treatment and TMB/BRAF status |
| ETS | ETS | Mediator | Early tumour shrinkage (radiological response) |
| PFS | PFS | Outcome | Progression-free survival |
| OS | OS | Outcome | Overall survival |

DAG node definitions

# Causal Structure

![Causal DAG for the METIMMOX-1 analysis. Node colour indicates role:
blue = exposure (T), red = outcomes (PFS, OS), green = latent variable
(U), grey = covariates and intermediate variables. Arrows denote assumed
causal effects. Node abbreviations: T = Treatment; TxTMB = Treatment x
TMB/BRAF interaction; ETS = Early Tumour Shrinkage; U = Unmeasured
confounders.](dag_files/figure-commonmark/dag-plot-1.png)

# Causal Analysis

Because T is a randomised exposure with no parents in the DAG, there are
no backdoor paths from T to PFS or OS. The sections below confirm this
using dagitty’s formal algorithms and list the implied conditional
independencies encoded by the DAG.

## Adjustment Sets

**T -\> PFS:** Minimal adjustment sets: {}

**T -\> OS:** Minimal adjustment sets: {}

## Implied Conditional Independencies

| Implied conditional independence (X \_\|\|\_ Y \| {Z}) |
|:-------------------------------------------------------|
| Age \_\|\|\_ CRP                                       |
| Age \_\|\|\_ ETS \| {CRP, TMB_BRAF}                    |
| Age \_\|\|\_ Sex                                       |
| Age \_\|\|\_ T                                         |
| Age \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| CRP \_\|\|\_ Sex                                       |
| CRP \_\|\|\_ T                                         |
| CRP \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| ETS \_\|\|\_ Sex \| {CRP, TMB_BRAF}                    |
| ETS \_\|\|\_ TxTMB \| {T, TMB_BRAF}                    |
| Sex \_\|\|\_ T                                         |
| Sex \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| T \_\|\|\_ TMB_BRAF                                    |

Conditional independencies implied by the causal DAG

------------------------------------------------------------------------

**Report completed on:** 2026-03-27  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 1.0
