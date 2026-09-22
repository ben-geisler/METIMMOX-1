# Cost-Effectiveness Analysis
Ben Geisler
2026-09-21

- [Overview](#overview)
- [Economic Survival Model](#economic-survival-model)
- [Base Case Results](#base-case-results)
- [Probabilistic Sensitivity
  Analysis](#probabilistic-sensitivity-analysis)
  - [Parameters sampled](#parameters-sampled)
  - [Results](#results)
- [Summary](#summary)
- [Common Target Population](#common-target-population)
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

QALYs and ongoing progressed-state cost rates use trapezoidal
integration over the weekly grid: 521 points cover exactly 520 weeks.
Scheduled treatment, diagnostic tests and visits, baseline costs and
end-of-life events retain full charges at their modeled time points
(issue \#163). See the [input-parameter report](input_parameters.md) for
the integration convention.

| Component | Formula |
|:---|:---|
| OS | Surv(OSwk, Death) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx |
| PFS | Surv(PFSwk, Progression) ~ Age + sex + Rx + crp*Rx + tmb_braf*Rx |
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
| Standard of Care | EUR 21,523 | 1.364 | EUR 48,061 |
| CRP-guided       | EUR 53,948 | 1.385 | EUR 16,684 |
| TMB/BRAF-guided  | EUR 62,685 | 1.353 |  EUR 6,312 |

Base case results at WTP = EUR 51,000

| Strategy | Cost | QALYs | Incremental Cost | Incremental QALYs | ICER | Status |
|:---|---:|---:|---:|---:|---:|:---|
| Standard of Care | EUR 21,523 | 1.364 | – | – | – | ND |
| CRP-guided | EUR 53,948 | 1.385 | EUR 32,424 | 0.021 | EUR 1,578,614 | ND |
| TMB/BRAF-guided | EUR 62,685 | 1.353 | – | – | Dominated | D |

Incremental cost-effectiveness results

<img src="CEA_files/figure-commonmark/ce-plane-1.png"
style="width:100.0%" data-fig-align="center"
alt="Cost-effectiveness plane for the single joint economic model" />

# Probabilistic Sensitivity Analysis

## Parameters sampled

The PSA samples the health-state utilities, the joint biomarker cell
probabilities, and the resource-use costs (visit, baseline work-up,
quarterly follow-up, end-of-life care), together with the survival
models. Unit drug and test prices are **fixed** (issue \#154): they are
published tariffs or, for nivolumab, an assumed acquisition price, so
treating them as uncertain would attribute decision uncertainty — and
value of information — to quantities that no study could resolve. Their
influence is reported in the one-way sensitivity analysis and in the
biosimilar pricing scenario instead.

The two utilities are sampled jointly rather than independently. The PSA
draws a non-negative decrement and sets the progressed utility to the
progression-free utility minus that decrement, so no draw can place the
progressed utility above the progression-free utility.

    Utility-ordering diagnostic: 0 of 5000 PSA draws (0.0%) have u_p > u_np, and 0.0% exceed it by more than 0.05. The decrement parameterisation makes these fractions zero by construction; independent beta draws for the two utilities did not.

## Results

| Strategy         |  Mean Cost | Mean QALYs | Probability Cost-Effective |
|:-----------------|-----------:|-----------:|---------------------------:|
| Standard of Care | EUR 21,497 |      1.401 |                     100.0% |
| CRP-guided       | EUR 53,531 |      1.407 |                       0.0% |
| TMB/BRAF-guided  | EUR 62,113 |      1.383 |                       0.0% |

PSA summary at WTP = EUR 51,000

The PSA and the base case use the same survival formulas and target
population. Since issue \#159, the OS and PFS coefficient vectors are
drawn together from a multivariate normal centred at their
maximum-likelihood estimates. Paired patient bootstrap fits estimate
their cross-endpoint dependence; block whitening and recolouring
preserve each original fit’s marginal covariance. The bootstrap is used
only to estimate dependence, never to supply PSA coefficient draws. See
the [joint sampling specification](technical/joint_survival_sampling.md)
for the covariance construction and diagnostics.

The deterministic result evaluates survival at the fitted coefficients;
the PSA averages nonlinear survival predictions over parameter
uncertainty. Those are different quantities: even correctly centred
coefficient draws need not give outcome means equal to the deterministic
result. Joint sampling changes the dependence and the ordering
correction, but, with fixed marginal coefficient distributions, cannot
remove the raw-curve mean shift in expectation. No draws or outcomes are
recentered. With all coefficient covariance set to zero and economic
parameters and population weights fixed, `test_psa_zero_uncertainty.R`
reproduces every deterministic cost and QALY result. The table below
retains the original five-Monte-Carlo-standard-error diagnostic in
`test_psa_basecase_alignment.R`; any continuing numerical failure is
reported explicitly.

| Strategy | Outcome | Base Case | PSA Mean | PSA SE | Difference | Difference (SE) |
|:---|:---|---:|---:|---:|---:|---:|
| Standard of Care | Cost | EUR 21,523 | EUR 21,497 | EUR 36 | -EUR 27 | -0.7 |
| Standard of Care | QALYs | 1.3644 | 1.4012 | 0.0042 | +0.0368 | +8.8 |
| CRP-guided | Cost | EUR 53,948 | EUR 53,531 | EUR 88 | -EUR 417 | -4.7 |
| CRP-guided | QALYs | 1.3849 | 1.4075 | 0.0038 | +0.0225 | +5.9 |
| TMB/BRAF-guided | Cost | EUR 62,685 | EUR 62,113 | EUR 94 | -EUR 572 | -6.1 |
| TMB/BRAF-guided | QALYs | 1.3529 | 1.3832 | 0.0039 | +0.0303 | +7.8 |

PSA means versus base-case values (difference in Monte Carlo standard
errors)

The nonlinear mean shift can differ across strategies, so it need not
cancel in increments even though every strategy shares the same
coefficient draw. The incremental comparison below reports this
remaining difference directly. Decisions under uncertainty use expected
net monetary benefit from the PSA; the deterministic results describe
the fitted-parameter scenario.

| Comparison | Outcome | Base Case | PSA Mean | PSA SE | Difference | Difference (SE) |
|:---|:---|---:|---:|---:|---:|---:|
| CRP-guided vs Standard of Care | Cost | EUR 32,424 | EUR 32,034 | EUR 80 | -EUR 390 | -4.9 |
| CRP-guided vs Standard of Care | QALYs | +0.0205 | +0.0063 | 0.0027 | -0.0143 | -5.3 |
| TMB/BRAF-guided vs Standard of Care | Cost | EUR 41,161 | EUR 40,616 | EUR 86 | -EUR 545 | -6.4 |
| TMB/BRAF-guided vs Standard of Care | QALYs | -0.0115 | -0.0180 | 0.0025 | -0.0065 | -2.6 |

Incremental PSA means versus base-case increments relative to standard
of care

<img src="CEA_files/figure-commonmark/ceac-1.png" style="width:100.0%"
data-fig-align="center" alt="Cost-effectiveness acceptability curves" />

# Summary

The economic evaluation now uses one joint survival model and three
economic strategies: standard of care, CRP-guided immunotherapy (week-4
CRP, before the first nivolumab dose), and TMB/BRAF-guided immunotherapy
(baseline NGS).

# Common Target Population

Control and both guided strategies are evaluated in the same
complete-case trial population (n = 68). The base case uses its observed
joint CRP and TMB/BRAF distribution. The PSA draws all four joint cell
probabilities together using a Dirichlet distribution with the observed
cell counts; marginal prevalences are derived from that draw. The age
and sex distributions within each cell remain fixed. Every draw’s common
population weights apply to control and both guided strategies. This
represents uncertainty in the trial population, not transport to an
external population. The prevalence DSA changes one marginal while
retaining the other marginal and joint odds ratio, and reweights every
strategy consistently (issue \#166).

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

**Report completed on:** 2026-09-21  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.5
