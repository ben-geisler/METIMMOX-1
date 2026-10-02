# Input Parameters
Ben Geisler
2026-10-02

- [Model Configuration](#model-configuration)
- [Biomarker Prevalence](#biomarker-prevalence)
- [Estimand of the CRP-Guided
  Strategy](#estimand-of-the-crp-guided-strategy)
- [Utilities](#utilities)
- [Costs](#costs)
  - [Post-progression treatment
    costs](#post-progression-treatment-costs)
  - [Second treatment sequence](#second-treatment-sequence)
- [Structural Parameters](#structural-parameters)
- [Input-to-Output Traceability](#input-to-output-traceability)
  - [Validation warnings](#validation-warnings)

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

# Estimand of the CRP-Guided Strategy

The economic model estimates a **randomisation-time policy**. Every
strategy is assigned at randomisation (t = 0), and costs and QALYs are
counted from then. The CRP-guided strategy is assigned using the week-4
CRP value (cycle 3 day 1), which is available before the first nivolumab
dose at modeled week 4. The fitted treatment and biomarker-by-treatment
effects act from t = 0, although all strategies receive identical
treatment until that dose. The 3 patients without a week-4 CRP are
excluded from the complete-case cohort of 68 patients (issue \#175).

The technical report `reports/technical/crp_week4_estimand.qmd`
quantifies both choices deterministically. It restores the patients
without a week-4 CRP under 4 combinations of assigned CRP status and
assigned TMB/BRAF status (for the patients whose TMB/BRAF status is
unknown):

- **CRP-guided:** pairwise ICERs from EUR 468,730 to EUR 815,624 in 4
  scenarios (base case: EUR 864,492).
- **TMB/BRAF-guided:** pairwise ICERs from EUR 456,666 to EUR 2,739,306
  in 3 scenarios; dominated by standard of care in 1 (base case:
  dominated by standard of care).
- Standard of Care has the highest net monetary benefit at WTP = EUR
  51,000 in the base case and in every scenario.

Before the first nivolumab dose, the CRP-guided strategy accrues
0.000013 more QALYs than standard of care. This is 0.035% of its
full-horizon QALY increment.

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
fraction is reported in the cost-effectiveness report. The reported
results use the CORRECT utilities (`UTILITY_SOURCE = 1`). The
alternative IPD mode (`UTILITY_SOURCE = 0`, issue \#188) derives both
utilities in the pipeline from the METIMMOX EQ-5D-5L responses with the
value set chosen by `EQ5D_VALUE_SET` (Danish by default: 0.912 and
0.897), and replaces the assumed CV 0.15 of `u_np` and `u_decrement`
with patient-clustered bootstrap standard errors (Danish set: 0.013 and
0.016).

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

**Outpatient visits** (issue \#36) are costed from the Norwegian DRG
(ISF) system: an outpatient contact carries DRG weight 0.047, and the
2023 unit price for somatic care is NOK 49,484 (EUR 4,541.60), giving
NOK 2,325.75 = EUR 213.46 per visit. The weight does not differ between
visit types, so the baseline visit (`c_other_baseline`), visits at
administrations and surveillance scans (`c_other_visit`) and the
quarterly follow-up visit in the progressed state (`c_other_follow`) all
take this value; the three parameters are kept separate so that each can
be varied on its own in the one-way analysis. The earlier values (EUR
33, 530 and 33) had no documented source.

**Monitoring and visits** (issue \#164). Progression-free patients have
a CT scan at baseline and every 12 weeks, and a blood test at baseline
and every 4 weeks, up to the horizon. A visit is charged at baseline, at
every drug administration and at every CT scan; during treatment each
scan already coincides with an administration, so the CT rule adds a
surveillance visit every 12 weeks after the last administration (week
38). Blood tests between scans carry no separate visit. Progressed
patients accrue one follow-up visit and one CT scan per quarter,
integrated over progressed-state occupancy.

**Dispersion of sampled parameters** (issue \#56). The four resource-use
costs are gamma distributed with a coefficient of variation (CV) of
0.20, and `u_np` (beta) and `u_decrement` (gamma) with a CV of 0.15.
Both CVs are **assumptions**: no study reports the dispersion of these
costs, and the utility sources give no usable standard error, so
conventional moderate values were chosen. They are neither estimated
from data nor taken from the literature, and they are not varied in a
scenario. A CV of 0.20 gives a gamma 95% interval of about 0.65 to 1.43
times the mean; a CV of 0.15 gives a beta 95% interval for `u_np` = 0.73
of about 0.49 to 0.91.

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
| Outpatient visit | c_other_visit | EUR 213 | EUR 171 | EUR 256 | Sampled (gamma) |
| Baseline visit | c_other_baseline | EUR 213 | EUR 171 | EUR 256 | Sampled (gamma) |
| Follow-up visit, progressed (per quarter) | c_other_follow | EUR 213 | EUR 171 | EUR 256 | Sampled (gamma) |
| Post-progression treatment (per quarter) | c_other_pp | EUR 0 | EUR 0 | EUR 0 | Fixed |
| End-of-life cost | c_other_last | EUR 13,803 | EUR 11,042 | EUR 16,564 | Sampled (gamma) |

Cost inputs

*Notes:* Fixed parameters are not sampled in the PSA and therefore carry
no value-of-information; they are varied in the one-way analysis.
Sampled costs are gamma distributed with an assumed CV of 0.20 (no
source; see text). Visit costs: DRG weight 0.047 times the 2023 ISF unit
price (NOK 49,484). The CT scan price is also charged once per quarter
in the progressed state. c_other_pp is zero in the base case:
post-progression treatment is outside the modeled scope (see below). Its
deterministic range is degenerate, so it is tested as a structural
scenario rather than in the one-way analysis.

## Post-progression treatment costs

The base case applies **no post-progression treatment cost**. Once a
patient has progressed, the model charges only the quarterly follow-up
visit (`c_other_follow`), a quarterly CT scan (`c_test_CT`, issue \#164)
and the one-time end-of-life cost (`c_other_last`); second-line systemic
therapy and any further post-progression visits are not modeled. Because
the strategies differ in the time patients spend in the progressed
state, this omission is **differential** rather than a common offset,
and its direction depends on which strategy accrues more progressed
time. The parameter `c_other_pp` makes the omission explicit and is
varied in a structural scenario (EUR 5,000 per quarter in the progressed
state) reported in the one-way sensitivity analysis.

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

# Input-to-Output Traceability

The table maps every field of the base-case parameter list
`l_params_base` to the function and calculation step that consumes it
(issue \#165). `model_fun()` resolves the grid, the population weights
and the survival curves of each subgroup and calls
`partitioned_survival_states()` and `calculate_outcomes()`, which
converts state occupancy into discounted costs and QALYs. The table is
built from the field names of `l_params_base`, and rendering stops if a
field has no entry or an entry names a field that does not exist, so it
cannot fall out of step with the model inputs.

| Field | Function | Step |
|:---|:---|:---|
| cl | model_fun(), discount_weights(), calculate_outcomes() | Cycle length in years: the year of each grid point for discounting, the duration of each interval for QALYs (trapezoidal weights), and the conversion of quarterly cost rates (4 x cl quarters per interval). |
| time_horizon | model_fun() | Number of weekly intervals; fixes the time_horizon + 1 grid points against which schedules and survival curves are checked and on which PSA curves are predicted. |
| dr_costs | discount_weights() | Annual discount rate of the cost weights applied to every cost stream. |
| dr_effects | discount_weights() | Annual discount rate of the QALY weights. |
| u_np | calculate_outcomes() | Utility weight of progression-free occupancy. |
| u_p | calculate_outcomes() | Utility weight of progressed occupancy (derived as u_np - u_decrement, floored at 0, in every PSA draw). |
| c_drug_nivo | calculate_outcomes() | Drug cost per nivolumab administration (l_nivo), biomarker-positive subgroups only. |
| c_drug_FLOX | calculate_outcomes() | Drug cost per FLOX administration: l_FLOX_exp in biomarker-positive subgroups, l_FLOX_control otherwise. |
| c_test_CT | calculate_outcomes() | CT price: scheduled scans (l_CT) while progression-free, and one scan per quarter of progressed occupancy. |
| c_test_blood | calculate_outcomes() | Blood-test price on the l_blood schedule. |
| c_test_CRP | calculate_outcomes() | One-time diagnostic test at week 0 for every patient of the CRP-guided strategy (also the screening cost of the enriched-population analysis). |
| c_test_NGS | calculate_outcomes() | One-time diagnostic test at week 0 for every patient of the TMB/BRAF-guided strategy (also the screening cost of the enriched-population analysis). |
| c_other_visit | calculate_outcomes() | Outpatient visit price on the l_visit schedule. |
| c_other_baseline | calculate_outcomes() | One-time baseline visit at week 0. |
| c_other_follow | calculate_outcomes() | Follow-up visit rate per quarter of progressed occupancy. |
| c_other_pp | calculate_outcomes() | Post-progression treatment rate per quarter of progressed occupancy (0 in the base case). |
| c_other_last | calculate_outcomes() | End-of-life cost charged on each interval’s increase in dead occupancy. |
| l_nivo | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| l_FLOX_exp | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| l_FLOX_control | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| l_CT | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| l_blood | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| l_visit | calculate_outcomes() | Indicator of a scheduled event at each weekly grid point; multiplied by the unit price and by progression-free occupancy. |
| p_os | model_fun(), partitioned_survival_states() | OS curve of each subgroup (control, biomarker-positive, biomarker-negative): progressed = max(OS - PFS, 0), dead = 1 - OS. The \*\_weighted curves are for reporting only. |
| p_pfs | model_fun(), partitioned_survival_states() | PFS curve of each subgroup: progression-free occupancy. |
| u_decrement | apply_derived_psa_parameters() | Not read by model_fun(); the sampled PSA decrement from which u_p is derived. |
| p_crp | model_fun(), model_population_weights() | Share of the CRP-guided cohort in the biomarker-positive subgroup: weights its costs and QALYs against the negative subgroup. |
| p_tmb_braf | model_fun(), model_population_weights() | Share of the TMB/BRAF-guided cohort in the biomarker-positive subgroup: weights its costs and QALYs against the negative subgroup. |
| population_curves | model_fun(), standardize_population_curves() | Patient-level curves, re-averaged when the population weights change. |
| prediction_population | model_fun() | Complete-case patients whose predictions are averaged, deterministically and in every PSA draw. |
| population_weights | model_fun() | Weights behind the current curves; compared with the weights implied by the prevalence fields. |
| joint_counts | configure_parameter_distributions() | Not read by model_fun(); complete-case joint cell counts, the Dirichlet parameters of the PSA prevalence draws. |
| p_joint_00 | model_population_weights() | Joint biomarker cell masses of the common target population; a change re-standardises all curves, control included. |
| p_joint_01 | model_population_weights() | Joint biomarker cell masses of the common target population; a change re-standardises all curves, control included. |
| p_joint_10 | model_population_weights() | Joint biomarker cell masses of the common target population; a change re-standardises all curves, control included. |
| p_joint_11 | model_population_weights() | Joint biomarker cell masses of the common target population; a change re-standardises all curves, control included. |

Fields of l_params_base and the model step that consumes each

------------------------------------------------------------------------

**Report completed on:** 2026-10-02  
**Repository:** ben-geisler/METIMMOX-1  
**Report version:** 4.8

## Validation warnings

No warnings recorded during rendering.
