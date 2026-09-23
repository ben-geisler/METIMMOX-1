# Issue #179: extrapolation plausibility

The new [technical report](../../reports/technical/psa_extrapolation_plausibility.md)
reads the current sampling, paired PSA and ordering-diagnostic caches. It sources
02-06 with the report's read-only sampling guard; no analysis cache is rebuilt.
The public 2025 Norwegian life-table extract and source details are in
[`data/external`](../../data/external/README.md).

Validation:

- `Rscript tests/test_psa_extrapolation_diagnostics.R`: pass. Constant hazards,
  sex mapping, fractional birthdays, terminal mortality, RMST-ranked contributions
  and stale-cache rejection without mutation are tested on synthetic inputs.
- `Rscript tests/test_report_contracts.R`: 563 checks pass; one expected skip
  for the empty zero-EVPI table.
- Both report formats reconstruct all 5,000 saved PSA rows using their
  `model_idx` and joint-population weights. QALYs agree with cached outcomes to
  8.88e-16; every raw subgroup crossing flag agrees with the ordering cache.
  The join includes simulation ID and curve as well as model index to avoid
  multiplying replacement-model rows with different population weights.
- The report includes central 50/80/95% ribbons, 50 evenly spaced spaghetti
  curves, deterministic and PSA means, the fixed matched-population reference,
  OS/PFS RMST distributions, 1/5/10% tail accounting and diagnostic trimmed
  means, plus crossing associations. The population comparison reports timing
  and magnitude as well as any-week and week-520 counts.
- PDF rendered successfully first, followed by GFM. PDF text, Markdown results,
  figure links and the survival graphic were inspected. All 11 existing sampling,
  PSA, EVPPI, scenario and ordering cache files retain their pre-task MD5 hashes.

The top 1% by OS RMST contributes 1.42%, 1.28% and 1.32% of QALYs for SoC,
CRP and TMB/BRAF. Removing those rows leaves PSA-minus-deterministic gaps of
0.03089, 0.01863 and 0.02585. No draw exceeds population survival at week 520,
although early exceedances are frequent; these must not be mistaken for
implausibly long tails. All exceedances lie in weeks 1-27, with a maximum
vertical excess of 0.2922 percentage points and no exceedance after year one.
Mean positive excess areas are 0.0000254, 0.0000260 and 0.0000227 life-years.
The known five-MCSE alignment failure is unchanged.
Trimming is a diagnostic, not a corrected economic result. A mortality hazard
floor remains outside the paper's scope.
