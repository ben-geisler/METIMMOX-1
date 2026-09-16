# Cache identity and reader isolation: issues #167 and #169

Implemented locally on 2026-09-16. Full production regeneration and report renders
are deferred by user decision until the other adversarial-review repairs are ready.
No model equations, economic assumptions or survival sampling distributions changed.

## Changes

- Sampling hashes cover joint-formula variables, endpoints and ID, preserving row
  order, factor levels and contrasts. Serialization version 2 and SHA-256 hash
  canonical values. Local caches retain canonical inputs, per-column hashes and
  runtime/library diagnostics; console diagnostics contain hashes, not patient values.
- Script 06 permits regeneration in the analysis pipeline. Report setup temporarily
  disables it, and knitting also disables it. Tests can redirect managed caches with
  `options(metimmox.cache_dir = ...)`. The index-validation test now uses synthetic
  fixtures and returns nonzero on assertion failure.
- PSA identities include calculation source and loaded function definitions,
  dependency versions, both current prediction populations and the existing economic
  inputs. Outcome and parameter files share generation metadata and content hashes;
  loaders verify exact ordered draw IDs and model indices.
- EVPPI identities independently cover the PSA generation/results, WTP, parameters,
  groups, estimator settings/code, seed and population scaling. The report loader
  recalculates EVPI from current PSA outcomes and WTP.
- Scenario identities cover the expected PSA inputs, scenario definitions and
  estimator/population inputs. Embedded PSA pairs and saved scenario result contents
  are checked. Reports and publication consumers use the validated loading paths.

## Verification

| Check | Result |
|---|---|
| `test_cache_provenance.R` | Pass: independent-session hashes; read-only sampling guard; interactive code edits; missing, shuffled and mixed PSA files; outcome corruption; duplicate draw IDs; EVPPI/scenario input changes; EVPI reconciliation; forced index-test failure returns nonzero |
| `test_sim_idx_validation.R` | Pass, all six cases; no trial data or analysis scripts sourced |
| `test_mvn_sampling.R` | Pass |
| `test_rng_reproducibility.R` | Pass |
| `test_sampling_rework.R` | Pass |
| `test_sampling_failure_fallback.R` | Pass |
| Optional 50-draw integration, scripts 06/10/11/12 | Pass: PSA, EVPPI and all six scenarios; repeated sampling/PSA loads reuse the cache; validated loaders accept the new artifacts |
| Report contracts on temporary generated caches | Pass, 284 checks; two skips: zero-EVPI row checks and the absent 5,000-draw sampling file (the integration uses 50 draws) |
| Missing sampling cache through actual `setup_report()` | Correctly stops without creating a cache |
| Production cache files before/after integration | Every `.rds`/`.RData` MD5 unchanged |
| Syntax | All 61 live R scripts and 44 Quarto documents parse |
| Production-cache report contracts | Expected failure on three legacy contracts: sampling diagnostics/schema, PSA pair binding, independent EVPPI identity; scenario cache absent |

The optional [integration.R](integration.R) reproduces the small pipeline run using
the local confidential trial dataset and a temporary cache directory. Run it with
Rscript from the repository root; it does not render or publish reports. The syntax
check emitted existing `xfun::attr()` deprecation warnings; the R session also emitted
locale/package-build warnings. Neither prevented the checks above.

## Historical data-hash investigation

The review-session hash `6cc5a549826e9837d4496533f4bb6243` was reproduced from the
current complete-case data using the previous hashing procedure. Sourcing script 05,
materialising/rebuilding all columns and serialising/unserialising the canonical
column list did not change it. The alternate pre-review hash `13f4cb24...` was not
reproduced, so the historical difference is not conclusively explained.

A separate representation defect is reproducible with synthetic values:
`rlang::hash(1:10)` gives `7f7765390f3eb9d1851af5858c6277d4`, whereas
`rlang::hash(c(1:10))` gives `d27012ad27a63c1cf259289712d92483` despite identical
integer values. The new serialization-based hash agrees for both representations.
This supports the hardening change but does not establish the cause of the historical
trial-data mismatch. The newly stored canonical inputs and runtime metadata make a
future discrepancy diagnosable.

## Deferred migration

Keep the existing production caches as historical results. After completing the
remaining model repairs, run setup scripts 02–05, then 06, 10, 11 and 12. Regenerate
the affected reports and publication vignettes afterwards, rendering PDF before GFM.
Do not copy new fingerprints onto legacy results to bypass validation. Until that
regeneration, strict readers and the live-cache contract test intentionally reject
the old caches. No GitHub comments, pushes or publications were made.
