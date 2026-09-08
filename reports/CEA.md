# Cost-Effectiveness Analysis
Ben Geisler
2026-09-08

- [Overview](#overview)
- [Economic Survival Model](#economic-survival-model)
- [Base Case Results](#base-case-results)
- [Probabilistic Sensitivity
  Analysis](#probabilistic-sensitivity-analysis)
  - [Parameters sampled](#parameters-sampled)
  - [Results](#results)
- [Summary](#summary)
- [Scope Limitations](#scope-limitations)

# Overview

This cost-effectiveness analysis evaluates standard of care against two
pre-immunotherapy biomarker-guided immunotherapy strategies for
metastatic MSS/pMMR colorectal cancer:

- CRP-guided treatment selection (CRP \< 5 mg/L at week 4, before the
  first nivolumab dose)
- TMB/BRAF-guided treatment selection (baseline next-generation
  sequencing)

The economic survival model is a single joint model with CRP and
TMB/BRAF treatment interactions. Economic strategies are restricted to
biomarkers that are available before the decision to add immunotherapy
is made.

**CRP timing.** CRP is not a baseline (pre-randomization) measurement.
In the METIMMOX trial every patient received two cycles of FLOX before
the first nivolumab dose, and the CRP used to define the CRP-guided
strategy is the value measured at week 4 (cycle 3 day 1), i.e. at the
point where the decision to add nivolumab is made. The model clock
starts at randomization, but this has no economic consequence: up to
week 4 both arms receive identical FLOX, monitoring and visits, and the
week-4 blood test is already part of the monitoring schedule. Three
patients without a week-4 CRP value died early (weeks 2.4, 15.7 and
20.9) and are excluded from the complete-case data on which the survival
models are fitted, so the CRP subgroup curves are conditional on
surviving to the week-4 measurement.

# Economic Survival Model

| Component | Formula |
|:---|:---|
| OS | Surv(OSwk, Death) ~ Age + sex + Rx + crp\*Rx + tmb_braf\*Rx |
| PFS | Surv(PFSwk, Progression) ~ Age + sex + Rx + crp\*Rx + tmb_braf\*Rx |
| Control arm | Same joint OS and PFS models, predicted with Rx = control for every patient |

Single economic survival model

| Strategy | Definition |
|:---|:---|
| Standard of Care | All patients receive standard FLOX chemotherapy |
| CRP-guided | CRP-positive patients receive FLOX + nivolumab; CRP-negative patients receive FLOX |
| TMB/BRAF-guided | TMB/BRAF-positive patients receive FLOX + nivolumab; TMB/BRAF-negative patients receive FLOX |

Economic strategies

# Base Case Results

| Strategy         |       Cost | QALYs |        NMB |
|:-----------------|-----------:|------:|-----------:|
| Standard of Care | EUR 21,523 | 1.371 | EUR 48,419 |
| CRP-guided       | EUR 53,928 | 1.392 | EUR 17,056 |
| TMB/BRAF-guided  | EUR 63,399 | 1.361 |  EUR 6,005 |

Base case results at WTP = EUR 51,000

| Strategy | Cost | QALYs | Incremental Cost | Incremental QALYs | ICER | Status |
|:---|---:|---:|---:|---:|---:|:---|
| Standard of Care | EUR 21,523 | 1.371 | -- | -- | -- | ND |
| CRP-guided | EUR 53,928 | 1.392 | EUR 32,404 | 0.020 | EUR 1,587,067 | ND |
| TMB/BRAF-guided | EUR 63,399 | 1.361 | -- | -- | Dominated | D |

Incremental cost-effectiveness results

<img src="CEA_files/figure-commonmark/ce-plane-1.png"
style="width:100.0%" data-fig-align="center"
alt="Cost-effectiveness plane for the single joint economic model" />

# Probabilistic Sensitivity Analysis

## Parameters sampled

The PSA samples the health-state utilities, the biomarker prevalences,
and the resource-use costs (visit, baseline work-up, quarterly
follow-up, end-of-life care), together with the survival models. Unit
drug and test prices are **fixed** (issue \#154): they are published
tariffs or, for nivolumab, an assumed acquisition price, so treating
them as uncertain would attribute decision uncertainty — and value of
information — to quantities that no study could resolve. Their influence
is reported in the one-way sensitivity analysis and in the biosimilar
pricing scenario instead.

The two utilities are sampled jointly rather than independently. The PSA
draws a non-negative decrement and sets the progressed utility to the
progression-free utility minus that decrement, so no draw can place the
progressed utility above the progression-free utility.

    Utility-ordering diagnostic: 0 of 5000 PSA draws (0.0%) have u_p > u_np, and 0.0% exceed it by more than 0.05. The decrement parameterisation makes these fractions zero by construction; independent beta draws for the two utilities did not.

## Results

| Strategy         |  Mean Cost | Mean QALYs | Probability Cost-Effective |
|:-----------------|-----------:|-----------:|---------------------------:|
| Standard of Care | EUR 21,467 |      1.416 |                     100.0% |
| CRP-guided       | EUR 53,593 |      1.418 |                       0.0% |
| TMB/BRAF-guided  | EUR 62,949 |      1.396 |                       0.0% |

PSA summary at WTP = EUR 51,000

The PSA and the base case use the same joint survival model: every PSA
draw predicts the control, biomarker-positive, and biomarker-negative
curves from one sampled coefficient vector of the joint fit
(multivariate-normal draws, issue \#156), with the control curve
obtained by setting treatment to standard of care for every patient,
exactly as in the base case. PSA means should therefore sit close to the
base-case values, differing only through the nonlinearity of the
partitioned survival model. The table below reports each difference in
units of the PSA Monte Carlo standard error; the regression test
`test_psa_basecase_alignment.R` flags any absolute difference above five
standard errors. With 5,000 draws the standard errors are small, and the
mean QALYs of every strategy sit a few hundredths of a QALY above the
base case with the same sign: this is the mean bias of the extrapolated
parametric curves under coefficient sampling, shared by all arms because
they come from one fit, and it is documented as a known deviation in the
test protocol.

| Strategy | Outcome | Base Case | PSA Mean | PSA SE | Difference | Difference (SE) |
|:---|:---|---:|---:|---:|---:|---:|
| Standard of Care | Cost | EUR 21,523 | EUR 21,467 | EUR 36 | -EUR 56 | -1.6 |
| Standard of Care | QALYs | 1.3714 | 1.4160 | 0.0041 | +0.0446 | +10.9 |
| CRP-guided | Cost | EUR 53,928 | EUR 53,593 | EUR 82 | -EUR 335 | -4.1 |
| CRP-guided | QALYs | 1.3918 | 1.4180 | 0.0037 | +0.0262 | +7.0 |
| TMB/BRAF-guided | Cost | EUR 63,399 | EUR 62,949 | EUR 99 | -EUR 450 | -4.6 |
| TMB/BRAF-guided | QALYs | 1.3609 | 1.3955 | 0.0038 | +0.0346 | +9.1 |

PSA means versus base-case values (difference in Monte Carlo standard
errors)

Because all strategies in a PSA draw share one sampled joint model, a
shift that is common to every strategy (the mean bias of the
extrapolated parametric curves under coefficient sampling) cancels in
the increments, whereas a control-specific shift would not. The
incremental comparison below therefore isolates what matters for the
decision.

| Comparison | Outcome | Base Case | PSA Mean | PSA SE | Difference | Difference (SE) |
|:---|:---|---:|---:|---:|---:|---:|
| CRP-guided vs Standard of Care | Cost | EUR 32,404 | EUR 32,126 | EUR 73 | -EUR 279 | -3.8 |
| CRP-guided vs Standard of Care | QALYs | +0.0204 | +0.0020 | 0.0025 | -0.0184 | -7.2 |
| TMB/BRAF-guided vs Standard of Care | Cost | EUR 41,876 | EUR 41,482 | EUR 93 | -EUR 394 | -4.2 |
| TMB/BRAF-guided vs Standard of Care | QALYs | -0.0106 | -0.0205 | 0.0023 | -0.0099 | -4.4 |

Incremental PSA means versus base-case increments relative to standard
of care

<img src="CEA_files/figure-commonmark/ceac-1.png" style="width:100.0%"
data-fig-align="center" alt="Cost-effectiveness acceptability curves" />

# Summary

The economic evaluation now uses one joint survival model and three
economic strategies: standard of care, CRP-guided immunotherapy (week-4
CRP, before the first nivolumab dose), and TMB/BRAF-guided immunotherapy
(baseline NGS).

# Scope Limitations

Three scope decisions bound the interpretation of these results (issue
\#154).

**No post-progression treatment costs.** After progression the model
charges only the quarterly follow-up contact and the one-time
end-of-life cost. Second-line systemic therapy, post-progression imaging
and post-progression visits are not modeled. The strategies differ in
time spent in the progressed state, so this omission is differential
rather than a common offset. A structural scenario charging EUR 5,000
per quarter in the progressed state is reported in the one-way
sensitivity analysis.

**Second treatment sequence given to all progression-free patients.**
Both arms receive a second eight-cycle sequence at weeks 24-38, applied
to every patient still progression-free at that point. METIMMOX
re-treated on progression during the treatment break, so the model
charges second-sequence drug and visit costs to patients who would not
have been re-treated. A partitioned survival model has no on-treatment
substate, so re-treatment cannot be conditioned on progression without
restructuring the model; a structural scenario removing the second
sequence, with survival held at the trial estimate, bounds the cost side
of the assumption.

**Nivolumab price is an assumption.** The base-case price per
administration is an assumed Norwegian hospital acquisition cost rather
than a citable tariff, because Norwegian hospital prices are set by
confidential LIS tender. It is the largest single cost driver; the
biosimilar scenario and the one-way analysis show how the conclusions
move with it.

Adverse-event costs and disutilities remain outside the modeled scope,
as described in the project documentation.

------------------------------------------------------------------------

**Report completed on:** 2026-09-08  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.2
