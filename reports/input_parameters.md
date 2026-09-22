# Input Parameters
Ben Geisler
2026-09-21

- [Model Configuration](#model-configuration)
- [Biomarker Prevalence](#biomarker-prevalence)
- [Utilities](#utilities)
- [Costs](#costs)
  - [Post-progression treatment
    costs](#post-progression-treatment-costs)
  - [Second treatment sequence](#second-treatment-sequence)
- [Structural Parameters](#structural-parameters)

# Model Configuration

State utilities and ongoing progressed-state costs are integrated over
the intervals between weekly grid points using the trapezoidal rule
(half weight at weeks 0 and 520). Thus 521 points span exactly 520
weeks: with no mortality, unit utilities and no discounting, the model
returns 10 QALYs. Quarterly follow-up and post-progression cost rates
accrue as four times the quarterly rate per year of progressed-state
occupancy. Discount weights are applied at each grid point before
integration. Drug administrations, scheduled tests and visits, baseline
screening and end-of-life events retain their full costs at their
modeled time points; they are not half-cycle corrected (issue \#163).

Diagnostic prices are read directly from `c_test_CRP` and `c_test_NGS`
on every calculation, including direct model and enriched-population
calls. A legacy `c_test_biomarker` lookup is ignored; editing a scalar
price needs no separate synchronization step (issue \#173).

| Item | Value |
|:---|:---|
| Economic strategies | control, crp, tmb_braf |
| Economic biomarkers | crp, tmb_braf |
| OS formula | Surv(OSwk, Death) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |
| PFS formula | Surv(PFSwk, Progression) ~ Age + sex + Rx + crp \* Rx + tmb_braf \* Rx |
| Control arm | The joint OS and PFS models predicted with Rx = control for every patient (no separate control formula; issue \#151) |

Single economic model configuration

# Biomarker Prevalence

| Biomarker         | Prevalence | Source                |
|:------------------|-----------:|:----------------------|
| CRP-positive      |      33.8% | METIMMOX-1 trial data |
| TMB/BRAF-positive |      44.1% | METIMMOX-1 trial data |

Economic biomarker prevalences

Every economic strategy uses the same complete-case population (n = 68).
Its observed joint CRP and TMB/BRAF distribution defines the base case.
In the PSA, the four joint cell probabilities are drawn together from a
Dirichlet distribution with the observed cell counts, and the two
marginal prevalences are derived. This represents uncertainty in the
trial population’s joint biomarker proportions, with age and sex
distributions within each cell held fixed. The population weights change
control and both guided strategies together. The one-way prevalence
analysis rakes this same joint distribution to the changed marginal,
retaining the other marginal and the joint odds ratio (issue \#166).

# Utilities

The two health-state utilities are not sampled independently in the PSA.
Independent beta draws placed the progressed utility above the
progression-free utility in a substantial minority of draws (about 16%
under the CORRECT utilities and about 49% under the IPD utilities),
which widened the PSA with orderings that cannot occur. The PSA instead
draws a non-negative decrement `u_decrement` (gamma, CV 0.15, mean equal
to the base-case difference `u_np - u_p`) and sets
`u_p = u_np - u_decrement`, so `u_p <= u_np` holds in every draw and the
marginal mean of `u_p` is unchanged (issue \#154). The realised reversal
fraction is reported in the cost-effectiveness report.

| Parameter                | Symbol | Base Case | DSA Low | DSA High |
|:-------------------------|:-------|----------:|--------:|---------:|
| Progression-free utility | u_np   |     0.730 |   0.584 |    0.876 |
| Post-progression utility | u_p    |     0.590 |   0.472 |    0.708 |

Health-state utilities

# Costs

Unit drug and test costs are **fixed in the PSA** and varied
deterministically only (issue \#154). They are published tariffs,
laboratory prices, or — for nivolumab — an assumed acquisition price;
none of them is a quantity that further research could resolve, so
sampling them attributed decision uncertainty — and, under the earlier
nearest-neighbour value-of-information estimator, research value — to
prices that are known. Fixing them roughly halved total EVPI (from EUR
45.2 to EUR 22.1 per patient) without changing any base-case result.
Visit, baseline, follow-up and end-of-life costs remain probabilistic
because they bundle genuine resource-use uncertainty, and the
deterministic influence of the prices is undiminished: the nivolumab
price is still the second-largest one-way driver for both guided
strategies.

The nivolumab price of EUR 13,923 per administration is an **assumption,
not a citable tariff**: it represents an assumed Norwegian hospital
acquisition price (about 40% below list). Norwegian hospital prices are
set by confidential LIS tender and cannot be published; international
list prices are higher. Eight administrations give about EUR 111,384 per
treated patient, making this the largest single cost driver in the
model. It is varied by +/-20% in the one-way analysis and by a factor of
three in the biosimilar scenario.

| Parameter | Symbol | Base Case | DSA Low | DSA High | PSA |
|:---|:---|---:|---:|---:|:---|
| Nivolumab per administration | c_drug_nivo | EUR 13,923 | EUR 11,138 | EUR 16,708 | Fixed |
| FLOX per administration | c_drug_FLOX | EUR 427 | EUR 342 | EUR 512 | Fixed |
| CT scan | c_test_CT | EUR 386 | EUR 309 | EUR 463 | Fixed |
| Routine blood monitoring | c_test_blood | EUR 16 | EUR 13 | EUR 19 | Fixed |
| CRP biomarker test | c_test_CRP | EUR 16 | EUR 13 | EUR 19 | Fixed |
| NGS test | c_test_NGS | EUR 2,518 | EUR 2,014 | EUR 3,022 | Fixed |
| Visit | c_other_visit | EUR 33 | EUR 26 | EUR 40 | Sampled (gamma) |
| Baseline other cost | c_other_baseline | EUR 530 | EUR 424 | EUR 636 | Sampled (gamma) |
| Follow-up other cost | c_other_follow | EUR 33 | EUR 26 | EUR 40 | Sampled (gamma) |
| Post-progression treatment (per quarter) | c_other_pp | EUR 0 | EUR 0 | EUR 0 | Fixed |
| End-of-life cost | c_other_last | EUR 13,803 | EUR 11,042 | EUR 16,564 | Sampled (gamma) |

Cost inputs

*Notes:* Fixed parameters are not sampled in the PSA and therefore carry
no value-of-information; they are varied in the one-way analysis.
c_other_pp is zero in the base case: post-progression treatment is
outside the modeled scope (see below). Its deterministic range is
degenerate, so it is tested as a structural scenario rather than in the
one-way analysis.

## Post-progression treatment costs

The base case applies **no post-progression treatment cost**. Once a
patient has progressed, the model charges only the quarterly follow-up
contact (`c_other_follow`) and the one-time end-of-life cost
(`c_other_last`); second-line systemic therapy, post-progression imaging
and post-progression outpatient visits are not modeled. Because the
strategies differ in the time patients spend in the progressed state,
this omission is **differential** rather than a common offset, and its
direction depends on which strategy accrues more progressed time. The
parameter `c_other_pp` makes the omission explicit and is varied in a
structural scenario (EUR 5,000 per quarter in the progressed state)
reported in the one-way sensitivity analysis.

## Second treatment sequence

The treatment schedules give a second eight-cycle sequence at weeks
24-38 to **every patient still progression-free** at that point, in both
arms. In METIMMOX the second sequence was started on progression during
the treatment break, so the model charges second-sequence drug and visit
costs to progression-free patients who would not have been re-treated. A
partitioned survival model has no on-treatment substate, so re-treatment
cannot be triggered on progression without restructuring the model; the
assumption is instead bounded by a structural scenario that removes the
second sequence entirely while holding survival at the trial estimate,
giving a cost-side bound.

# Structural Parameters

| Parameter                    |                             Value |
|:-----------------------------|----------------------------------:|
| Cycle length                 | 1 week (1/52 year = 0.0192 years) |
| Time horizon                 |              520 weeks (10 years) |
| Cost discount rate           |                              4.0% |
| Effect discount rate         |                              4.0% |
| Willingness-to-pay threshold |                        EUR 51,000 |
| PSA sample size              |                             5,000 |

Structural model parameters

------------------------------------------------------------------------

**Report completed on:** 2026-09-21  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.4
