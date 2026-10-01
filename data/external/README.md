# Norwegian population mortality

`ssb_07902_2025.csv` contains the public 2025 period life-table probabilities
of death between exact ages x and x+1, by sex, from Statistics Norway table
[07902](https://www.ssb.no/en/statbank/table/07902). Retrieved 23 September 2026
from the full **Life tables** table on the official
[Deaths statistics page](https://www.ssb.no/en/befolkning/fodte-og-dode/statistikk/dode).
Table 07902 was updated 12 March 2026. The published probabilities per 1,000
were divided by 1,000 to obtain `qx`; ages 0 through 106 and both sexes are retained.
The source's three-decimal precision per 1,000 is preserved. These are public
population statistics, not trial data.

The extrapolation report uses piecewise constant annual hazards
`-log(1-qx)` at attained age, holding 2025 mortality rates fixed through follow-up.
It does not forecast mortality improvements or change economic model inputs.

# Published survival benchmarks (issue #162)

`survival_benchmarks.csv` holds the external survival estimates against which
`reports/technical/survival_external_validation.qmd` checks the extrapolated
standard-of-care curves. Every value is transcribed from the text or tables of a
published source (retrieved 1 October 2026 through PubMed/PMC); nothing was
digitised from a figure. The file was committed before any model comparison was
computed, so the benchmarks and the tolerance below are fixed in advance.

| Source | Use | DOI |
|---|---|---|
| Tveit et al. 2012, NORDIC-VII primary analysis, arm A (Nordic FLOX, n = 185) | median OS and PFS with reported 95% CIs; the CIs were taken from the tabulation of the same trial in Oncotarget 2016 Table 2 (doi:10.18632/oncotarget.8644) because the JCO full text was not accessible | 10.1200/JCO.2011.38.0915 |
| Guren et al. 2017, NORDIC-VII final survival analysis | 12 of 185 arm-A patients alive at the censoring date, each with at least 6.5 years of follow-up | 10.1038/bjc.2017.93 |
| Aasebo et al. 2019, Scandinavian population-based mCRC cohort | MSS patients after first-line chemotherapy (Table 2): combination chemotherapy and any chemotherapy | 10.1002/cam4.2205 |
| Central Norway population-based study, 2001-2015 (Acta Oncologica 2025) | primary palliative chemotherapy, synchronous and metachronous metastases (Table 2) | 10.2340/1651-226X.2025.42985 |
| Sorbye et al. 2013, Nordic cancer registries | Norwegian synchronous mCRC diagnosed 2006-2008 (abstract values; the full text was not accessible) | 10.1093/annonc/mdt197 |

The METIMMOX primary publication (Ree et al. 2024, doi:10.1038/s41416-024-02696-6)
reports the same trial and is therefore not an independent benchmark.

**Roles.** `comparable` rows are first-line oxaliplatin-based chemotherapy cohorts
close enough to the METIMMOX control arm for a two-sided check. `lower_bound` rows
are unselected or real-world populations (older patients, poorer performance
status, monotherapy, untreated patients in the registries), whose survival a
trial population receiving FLOX should not fall below; they give a one-sided check
only.

**Tolerance (fixed in advance).** The primary model curve is the economic
standard-of-care curve (all complete-case patients predicted with `Rx = control`).
For a `comparable` row the model value must lie within the benchmark's reported
95% confidence interval (inclusive). For a `lower_bound` row the model value must
be at least the benchmark estimate. Medians are compared in months (weeks x 12 / 52)
and survival proportions at `time_years`. A row the model cannot evaluate (for
example a median beyond the horizon) is "Not evaluable". A failed check is
reported and counted in the report; it does not stop the render. The rule is
implemented in `evaluate_benchmark()` in `R/survival_validation.R`.

**Known differences** recorded in the `note` column: some sources measure OS from
the diagnosis of metastatic disease rather than from randomisation (which lengthens
their OS), cohorts predate METIMMOX (2018-2023) by 10-20 years, and only Aasebo
2019 is restricted to MSS tumours.

**Headline criterion added after review (1 October 2026).** The point-estimate
rule above stays as pre-specified and is still reported. After the comparison was
run, it was judged too strict as the headline, because it ignores the uncertainty
of a model fitted to 68 patients. The report therefore leads with consistency
under uncertainty:

- a comparable benchmark is consistent if its published estimate lies inside the
  95% coefficient-draw interval of the parametric survival model;
- a lower bound is consistent if the upper limit of that interval reaches the bound.

This rule is implemented in `benchmark_consistent()` in `R/survival_validation.R`.
