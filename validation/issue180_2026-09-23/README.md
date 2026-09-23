# Primary probabilistic reporting (#180)

The CEA report and publication Table 5 now lead with PSA means and 95%
percentile uncertainty intervals. Incremental intervals use paired differences
within each retained draw versus standard of care. ICERs divide mean incremental
cost by mean incremental QALYs; they never average per-draw ratios. The CSV
contains unrounded numeric results and pairwise status, while the display uses
dominance labels. Tables 4/S5 remain deterministic point evaluations anchoring
the one-way sensitivity analysis.

The historical five-MCSE outcome-mean alignment criterion was retired because
E[f(theta)] need not equal f(E[theta]). It was not relaxed. The report preserves
the level and incremental comparisons as methods diagnostics. The structural
joint-control test remains, and zero-uncertainty equality is the correctness
criterion. Joint OS/PFS sampling is unchanged.

Validation on 23 September 2026:

- `Rscript tests/test_psa_summary.R`: passed synthetic paired-percentile,
  ratio-of-means, reference, dominance, zero-effect and invalid-input checks.
- `Rscript tests/test_psa_basecase_alignment.R`: passed the structural
  joint-model prediction check with treatment set to control for every patient.
- `Rscript tests/test_psa_zero_uncertainty.R`: passed on trial data; all
  deterministic costs reproduced within EUR 1e-8 and QALYs within 1e-12.
- `Rscript tests/test_report_contracts.R`: 575 checks passed; one inapplicable
  empty-EVPPI regression-row block skipped. This includes independent arithmetic
  checks for every numeric Table 5 column against the existing PSA cache.
- Tables 4/S5 and Table 5 rendered through their generating vignettes and
  registered in the publication manifest. Tables 4/S5 CSV values are unchanged.
- CEA rendered to PDF followed by GFM; checked the Markdown results and PDF
  text. The report is version 4.6.

No sampling, PSA, EVPPI or scenario cache regeneration was required. No model
calculation or economic parameter changed. Standard of care remains optimal in
all 5,000 PSA draws at EUR 51,000/QALY. CRP mean incremental cost is EUR 32,034.23
and mean incremental QALYs 0.00628881 (ratio EUR 5,093,844.42/QALY); TMB/BRAF is
dominated by standard of care on mean outcomes.

Manuscript #160 remains a downstream writing task: `docs/manuscript/` contains
only a placeholder. Its headline values should come from the probabilistic
Table 5 export.
