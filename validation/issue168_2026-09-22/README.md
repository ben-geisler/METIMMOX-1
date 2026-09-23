# Issue #168: publication artifacts and provenance

Implemented on 22 September 2026; final validation completed on 23 September.

Table S7 and Figure 5 now retain the zero-EVPI base case. Scenario presence and
total EVPI come from each scenario's PSA object and WTP. Exactly zero EVPI
bounds each configured group at zero (SE 0); percentages at zero EVPI are
undefined. Missing or failed estimates at positive EVPI remain `NA` with a
reason. Table S7 carries `error` and `method` columns, and Figure 5 labels zeros
and missing estimates explicitly. The scenario report separately identifies
total EVPI and excludes it from the EVPPI plot.

Every publication vignette declares its outputs before setup, removes its
predecessors, verifies CSV serialization against the generating object and
records completed outputs in [the manifest](../../outputs/manifest.csv).
Empty renders remove obsolete files or write an explicit empty-state output.
The manifest records generating vignettes, loaded cache fingerprints, UTC render
times, output hashes, source/input digests and write/verification status.
Its scope is the 37 vignette-owned publication exports; the three `traces_*.png`
files are separate script-07 pipeline diagnostics.

All 26 vignettes rendered successfully. All 17 tracked CSVs were compared with
their generating objects and the original Git tree: 16 are value-identical;
Table S7 now has 42 rows, including seven base-case rows. Table 5 already matched
the current PSA before this fix: SoC EUR 21,497 / 1.40 QALYs / 100.0%, CRP
EUR 53,531 / 1.41 / 0.0%, and TMB/BRAF EUR 62,113 / 1.38 / 0.0%. These rounded
values now have an independent regression check against the PSA.
See [csv_audit.csv](csv_audit.csv).

## Snapshots and cache provenance

The requested commands were run unchanged:

```powershell
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' analysis/13_save_snapshot.R 168 baseline
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' analysis/13_save_snapshot.R 168 fixed
```

Both snapshot pairs precede the fix commit and carry HEAD `b8e1b02`. The baseline
was saved at 11:25:26 on 22 September before editing sources. The fixed snapshot
was saved at 12:04:45 after regenerating all 5,000 PSA draws and EVPPI under the
snapshot script's existing cache-before-HEAD rule. Deterministic results are
identical. Maximum absolute differences across all PSA rows are EUR 2.91e-11
and 8.88e-16 QALYs; EVPI remains zero. Sampling indices and PSA input fingerprints
are unchanged. See [snapshot_comparison.json](snapshot_comparison.json).

[recompile_scenarios.R](recompile_scenarios.R) verifies that only
`compile_evppi_results()` changed in the parsed scenario source, reconstructs
the original implementation identity from `b8e1b02`, and validates the old
scenario cache, its embedded PSA pairs and result hash. It then rebuilds the
changed compilation and population columns from those validated scenario
results, saves the new identity and validates the result. Estimation functions
and inputs were unchanged, so their existing results were reused; obsolete
results were not merely relabelled with a new fingerprint.

## Validation

| Check | Result |
|---|---|
| `tests/test_publication_artifacts.R` | Pass: exact zero, small positive and failed estimates; real Figure 5 plotting code with zero/NA/empty fixtures; CSV object agreement; stale CSV/PNG removal and manifest lifecycle |
| `tests/test_report_contracts.R` | 550 checks pass; one inapplicable regression-row block skipped because the base EVPPI cache is empty at EVPI 0 |
| `tests/test_cache_provenance.R` | Pass |
| `tests/test_psa_zero_uncertainty.R` | Pass |
| `tests/test_psa_basecase_alignment.R` | Structural check passes; historical numerical criterion still fails for four of six level comparisons, unchanged |
| `tests/compare_snapshots.R 168` | Pass; comparison outputs in ignored `outputs/temp/` |
| `check_snapshots.R` | Pass: deterministic and PSA invariance within stated tolerances |
| `audit_csvs.R` | All 17 tracked CSVs verified |

The unchanged five-SE failures are control QALYs (+8.77 SE), CRP QALYs
(+5.94 SE), TMB/BRAF costs (-6.11 SE) and TMB/BRAF QALYs (+7.80 SE).
No threshold was relaxed or coefficient draws recentered. The runtime emits
pre-existing locale and package-build-version warnings; no new test warnings
were observed.

The scenario and bug-fix impact reports were rendered to PDF first and GFM
second. Figure 5 was inspected visually after its final render. `AGENTS.md`
records the output contract and the updated acceptance-test results.

## Reproduce

Run from the repository root with the confidential data and valid caches:

```powershell
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' validation/issue168_2026-09-22/recompile_scenarios.R
& ./validation/issue168_2026-09-22/render_artifacts.ps1
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' validation/issue168_2026-09-22/audit_csvs.R
& 'C:\Program Files\R\R-4.3.2\bin\x64\Rscript.exe' validation/issue168_2026-09-22/check_snapshots.R
```

Render vignettes sequentially: the manifest has one writer.
