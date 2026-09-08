# Input Parameters
Ben Geisler
2026-09-08

- [Model Configuration](#model-configuration)
- [Biomarker Prevalence](#biomarker-prevalence)
- [Utilities](#utilities)
- [Costs](#costs)
  - [Post-progression treatment
    costs](#post-progression-treatment-costs)
  - [Second treatment sequence](#second-treatment-sequence)
- [Structural Parameters](#structural-parameters)

# Model Configuration

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
| TMB/BRAF-positive |      44.9% | METIMMOX-1 trial data |

Economic biomarker prevalences

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

<table>

<caption>

Cost inputs
</caption>

<thead>

<tr>

<th style="text-align:left;">

Parameter
</th>

<th style="text-align:left;">

Symbol
</th>

<th style="text-align:right;">

Base Case
</th>

<th style="text-align:right;">

DSA Low
</th>

<th style="text-align:right;">

DSA High
</th>

<th style="text-align:left;">

PSA
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Nivolumab per administration
</td>

<td style="text-align:left;">

c_drug_nivo
</td>

<td style="text-align:right;">

EUR 13,923
</td>

<td style="text-align:right;">

EUR 11,138
</td>

<td style="text-align:right;">

EUR 16,708
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

FLOX per administration
</td>

<td style="text-align:left;">

c_drug_FLOX
</td>

<td style="text-align:right;">

EUR 427
</td>

<td style="text-align:right;">

EUR 342
</td>

<td style="text-align:right;">

EUR 512
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

CT scan
</td>

<td style="text-align:left;">

c_test_CT
</td>

<td style="text-align:right;">

EUR 386
</td>

<td style="text-align:right;">

EUR 309
</td>

<td style="text-align:right;">

EUR 463
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

Routine blood monitoring
</td>

<td style="text-align:left;">

c_test_blood
</td>

<td style="text-align:right;">

EUR 16
</td>

<td style="text-align:right;">

EUR 13
</td>

<td style="text-align:right;">

EUR 19
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

CRP biomarker test
</td>

<td style="text-align:left;">

c_test_CRP
</td>

<td style="text-align:right;">

EUR 16
</td>

<td style="text-align:right;">

EUR 13
</td>

<td style="text-align:right;">

EUR 19
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

NGS test
</td>

<td style="text-align:left;">

c_test_NGS
</td>

<td style="text-align:right;">

EUR 2,518
</td>

<td style="text-align:right;">

EUR 2,014
</td>

<td style="text-align:right;">

EUR 3,022
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

Visit
</td>

<td style="text-align:left;">

c_other_visit
</td>

<td style="text-align:right;">

EUR 33
</td>

<td style="text-align:right;">

EUR 26
</td>

<td style="text-align:right;">

EUR 40
</td>

<td style="text-align:left;">

Sampled (gamma)
</td>

</tr>

<tr>

<td style="text-align:left;">

Baseline other cost
</td>

<td style="text-align:left;">

c_other_baseline
</td>

<td style="text-align:right;">

EUR 530
</td>

<td style="text-align:right;">

EUR 424
</td>

<td style="text-align:right;">

EUR 636
</td>

<td style="text-align:left;">

Sampled (gamma)
</td>

</tr>

<tr>

<td style="text-align:left;">

Follow-up other cost
</td>

<td style="text-align:left;">

c_other_follow
</td>

<td style="text-align:right;">

EUR 33
</td>

<td style="text-align:right;">

EUR 26
</td>

<td style="text-align:right;">

EUR 40
</td>

<td style="text-align:left;">

Sampled (gamma)
</td>

</tr>

<tr>

<td style="text-align:left;">

Post-progression treatment (per quarter)
</td>

<td style="text-align:left;">

c_other_pp
</td>

<td style="text-align:right;">

EUR 0
</td>

<td style="text-align:right;">

EUR 0
</td>

<td style="text-align:right;">

EUR 0
</td>

<td style="text-align:left;">

Fixed
</td>

</tr>

<tr>

<td style="text-align:left;">

End-of-life cost
</td>

<td style="text-align:left;">

c_other_last
</td>

<td style="text-align:right;">

EUR 13,803
</td>

<td style="text-align:right;">

EUR 11,042
</td>

<td style="text-align:right;">

EUR 16,564
</td>

<td style="text-align:left;">

Sampled (gamma)
</td>

</tr>

</tbody>

</table>

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

**Report completed on:** 2026-09-08  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.2
