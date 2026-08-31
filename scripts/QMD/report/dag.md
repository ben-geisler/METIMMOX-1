# Causal DAG
Ben Geisler
2026-08-31

- [Introduction](#introduction)
- [Node Definitions](#node-definitions)
- [Causal Structure](#causal-structure)
- [Causal Analysis](#causal-analysis)
  - [Adjustment Sets](#adjustment-sets)
  - [Implied Conditional
    Independencies](#implied-conditional-independencies)
- [Sensitivity Analysis: Excluding Age/Sex →
  TMB/BRAF](#sensitivity-analysis-excluding-agesex--tmbbraf)
  - [Sensitivity DAG](#sensitivity-dag)
  - [Sensitivity Causal Analysis](#sensitivity-causal-analysis)
    - [Adjustment Sets](#adjustment-sets-1)
    - [Implied Conditional
      Independencies](#implied-conditional-independencies-1)
- [Simplified DAG](#simplified-dag)

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
| TxCRP | T x CRP | Effect modifier | Interaction between treatment and CRP status |
| TLR | TLR | Mediator | Tumour lesion reduction (radiological response) |
| PFS | PFS | Outcome | Progression-free survival |
| OS | OS | Outcome | Overall survival |

DAG node definitions

# Causal Structure

<img src="dag_files/figure-commonmark/dag-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Causal DAG for the METIMMOX-1 analysis. Node colour indicates role: blue = exposure (T), red = outcomes (PFS, OS), green = latent variable (U), grey = covariates and intermediate variables. Solid arrows denote assumed causal effects; dashed arrows indicate the hypothesised interaction effects of T x TMB/BRAF and T x CRP on PFS and OS. Node abbreviations: T = Treatment; TxTMB = Treatment x TMB/BRAF interaction; TxCRP = Treatment x CRP interaction; TLR = Tumour Lesion Reduction; U = Unmeasured confounders." />

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
| Age \_\|\|\_ Sex                                       |
| Age \_\|\|\_ T                                         |
| Age \_\|\|\_ TLR \| {CRP, TMB_BRAF}                    |
| Age \_\|\|\_ TxCRP                                     |
| Age \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| CRP \_\|\|\_ Sex                                       |
| CRP \_\|\|\_ T                                         |
| CRP \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| Sex \_\|\|\_ T                                         |
| Sex \_\|\|\_ TLR \| {CRP, TMB_BRAF}                    |
| Sex \_\|\|\_ TxCRP                                     |
| Sex \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| T \_\|\|\_ TMB_BRAF                                    |
| TLR \_\|\|\_ TxCRP \| {CRP, T}                         |
| TLR \_\|\|\_ TxTMB \| {T, TMB_BRAF}                    |
| TMB_BRAF \_\|\|\_ TxCRP \| {CRP}                       |
| TxCRP \_\|\|\_ TxTMB \| {T, TMB_BRAF}                  |
| TxCRP \_\|\|\_ TxTMB \| {CRP, T}                       |

Conditional independencies implied by the causal DAG

# Sensitivity Analysis: Excluding Age/Sex → TMB/BRAF

This sensitivity DAG removes the `Age -> TMB_BRAF` and `Sex -> TMB_BRAF`
edges while keeping the remainder of the structure identical to the
primary specification. The primary DAG retains these arrows on
biological grounds; this sensitivity analysis shows the downstream
implications of omitting them.

## Sensitivity DAG

<img src="dag_files/figure-commonmark/dag-sens-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Sensitivity DAG for the METIMMOX-1 analysis, excluding the Age -&gt; TMB_BRAF and Sex -&gt; TMB_BRAF edges while retaining all other assumed causal relationships. Solid arrows denote assumed causal effects; dashed arrows indicate the hypothesised interaction effects of T x TMB/BRAF and T x CRP on PFS and OS." />

## Sensitivity Causal Analysis

The treatment node remains randomised in the sensitivity DAG, so the
formal adjustment recommendations for `T -> PFS` and `T -> OS` are
unchanged. The main difference is in the implied conditional
independencies, which now reflect the omission of direct Age/Sex effects
on TMB/BRAF status.

### Adjustment Sets

**T -\> PFS:** Minimal adjustment sets: {}

**T -\> OS:** Minimal adjustment sets: {}

### Implied Conditional Independencies

| Implied conditional independence (X \_\|\|\_ Y \| {Z}) |
|:-------------------------------------------------------|
| Age \_\|\|\_ CRP                                       |
| Age \_\|\|\_ PFS                                       |
| Age \_\|\|\_ Sex                                       |
| Age \_\|\|\_ T                                         |
| Age \_\|\|\_ TLR                                       |
| Age \_\|\|\_ TMB_BRAF                                  |
| Age \_\|\|\_ TxCRP                                     |
| Age \_\|\|\_ TxTMB                                     |
| CRP \_\|\|\_ Sex                                       |
| CRP \_\|\|\_ T                                         |
| CRP \_\|\|\_ TxTMB \| {TMB_BRAF}                       |
| OS \_\|\|\_ Sex                                        |
| PFS \_\|\|\_ Sex                                       |
| Sex \_\|\|\_ T                                         |
| Sex \_\|\|\_ TLR                                       |
| Sex \_\|\|\_ TMB_BRAF                                  |
| Sex \_\|\|\_ TxCRP                                     |
| Sex \_\|\|\_ TxTMB                                     |
| T \_\|\|\_ TMB_BRAF                                    |
| TLR \_\|\|\_ TxCRP \| {CRP, T}                         |
| TLR \_\|\|\_ TxTMB \| {T, TMB_BRAF}                    |
| TMB_BRAF \_\|\|\_ TxCRP \| {CRP}                       |
| TxCRP \_\|\|\_ TxTMB \| {T, TMB_BRAF}                  |
| TxCRP \_\|\|\_ TxTMB \| {CRP, T}                       |

Conditional independencies implied by the sensitivity DAG

# Simplified DAG

The following DAG presents the same causal assumptions in a condensed
form intended for a clinical audience. The interaction terms (T ×
TMB/BRAF and T × CRP) and the latent common-cause node (U) have been
omitted; their roles are described in the text. Dashed arrows from T to
PFS and OS indicate that the magnitude of the treatment effect on
survival outcomes is hypothesised to be modified by TMB/BRAF and CRP
status — the central question addressed by the biomarker analysis.

<img src="dag_files/figure-commonmark/dag-simple-plot-1.png"
style="width:100.0%" data-fig-align="center"
alt="Figure 3. Simplified causal DAG for METIMMOX-1. Solid arrows denote assumed causal effects. The dashed arrow indicates that the treatment effect on survival is hypothesised to be modified by TMB/BRAF and CRP status. Interaction and latent-variable nodes are shown in the full DAG. Abbreviations: CRP = C-reactive protein, OS = overall survival, PFS = progression-free survival, T = treatment, TMB = tumour mutational burden." />

------------------------------------------------------------------------

**Report completed on:** 2026-08-31  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 2.2
